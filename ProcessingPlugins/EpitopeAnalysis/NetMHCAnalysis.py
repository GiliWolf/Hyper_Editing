__author__ = 'Hillel'
# =====================imports=======================#
# region Builtin Imports
import argparse
import logging
import os
from _csv import reader
from collections import namedtuple
from sys import path

# endregion

# region Internal Imports
if __name__ == "__main__":
    path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from Outputs.Outputers.CSVOutputer import CSVOutputer
from ProcessingPlugins.EpitopeAnalysis.annovar_analysis import CHR, POSITION, P_ID_SEP
from ProcessingPlugins.EpitopeAnalysis.EpitopeAnalysisConsts import *
from Commons.general_functions import init_logging_dict
# endregion

# =====================classes=======================#
"""This namedtuple hold the data from one line of netMHC output file.
per_allele_affinity: the affinity of the peptide per HLA allele checked.
p_length: the length of the checked peptide.
chromosome: the chromosome of the mismacth site.
position: the position of the mismatch
avg_rank: the average rank of the peptide.
binders_num: the number of allele the peptide binds to."""
NetMHCOutputRec = namedtuple("NetMHCOutputRec", "per_allele_affinity p_length chromosome position avg_rank binders_num")
# =====================functions=====================#


def filter_net_mhc4_results(neo_net_mhc_tab, orig_net_mhc_tab, neopeptides_csv, output_path):
    """
    This function filters and processes the output from run of NetMHC on the mutated peptide and it's original, based
    on output from annovar_analysis.py
    :param neo_net_mhc_tab: the path to the NetMHC output of the mutated peptides.
    :param orig_net_mhc_tab: the path to the NetMHC output of the original peptides.
    :param neopeptides_csv:  the output summery file of the mutated peptides from epitope_analysis
    :param output_path: output path.
    :return: None
    """
    with open(neo_net_mhc_tab) as mhc:
        neo_mhc_data = [l for l in reader(mhc, delimiter="\t")]
    with open(orig_net_mhc_tab) as mhc:
        orig_mhc_data = [l for l in reader(mhc, delimiter="\t")]
    with open(neopeptides_csv) as neo_p:
        neo_p_data = [l for l in reader(neo_p)]

    alleles_line = neo_mhc_data[0]
    headers_line = neo_mhc_data[1]
    neo_p_alleles_i, neo_p_offset_i, neo_p_id_i, neo_p_affinities_i, neo_window_p_i, neo_avg_rank_i, neo_binders_num_i \
        = get_net_mhc_indexes(alleles_line, headers_line)

    alleles_line = orig_mhc_data[0]
    headers_line = orig_mhc_data[1]
    orig_p_alleles_i, orig_p_offset_i, orig_p_id_i, orig_p_affinities_i, orig_window_p_i, orig_avg_rank_i, \
    orig_binders_num_i = get_net_mhc_indexes(alleles_line, headers_line)

    change_positions = {}
    for line in neo_p_data[1:]:
        if line == "":
            continue
        change_positions.setdefault(line[CHR_I], {})[line[POSITION_I]] = line[CHANGE_P_I]

    neo_mhc_filtered_recs = parse_net_mhc_output(change_positions, neo_mhc_data, neo_p_affinities_i, neo_p_alleles_i,
                                                 neo_p_id_i,
                                                 neo_p_offset_i, neo_window_p_i, neo_avg_rank_i, neo_binders_num_i)
    orig_mhc_filtered_recs = parse_net_mhc_output(change_positions, orig_mhc_data, orig_p_affinities_i,
                                                  orig_p_alleles_i, orig_p_id_i, orig_p_offset_i, orig_window_p_i,
                                                  orig_avg_rank_i, orig_binders_num_i)

    out_recs, diff_recs, add_headers = get_net_mhc_output_recs(neo_mhc_filtered_recs, orig_mhc_filtered_recs)
    outputer = CSVOutputer()
    outputer.output([output_path], add_headers, out_recs)
    outputer.output([output_path + ".diff.csv"], add_headers, diff_recs)


def get_net_mhc_output_recs(neo_mhc_filtered_recs, orig_mhc_filtered_recs):
    """
    This function builds the output recs.
    :param neo_mhc_filtered_recs: The parsed records of the mutated peptides.
    :param orig_mhc_filtered_recs: The parsed records of the original peptides.
    :return: all records, records that were different between the originals and mutated ones, and the headers.
    """
    headers = []
    first_headers = []
    h_flags = {}
    out_recs = []
    diff_recs = []

    for chromosome in neo_mhc_filtered_recs:
        for pos in neo_mhc_filtered_recs[chromosome]:
            identical = True
            for p_length in neo_mhc_filtered_recs[chromosome][pos]:
                for offset in neo_mhc_filtered_recs[chromosome][pos][p_length]:
                    type_tag = NEO_PEP_TAG
                    avg_rank = neo_mhc_filtered_recs[chromosome][pos][p_length][offset].avg_rank
                    binders_num = neo_mhc_filtered_recs[chromosome][pos][p_length][offset].binders_num
                    rec = {CHR: chromosome, POSITION: pos, type_tag + RANKS_AVG_HEADER: avg_rank,
                           type_tag + BINDERS_NUM_HEADER: binders_num, "Offset": str(offset),
                           P_LENGTH: str(p_length)}
                    if CHR not in h_flags:
                        first_headers.extend([CHR, POSITION, P_LENGTH, OFFSET_FORMAT, type_tag + RANKS_AVG_HEADER,
                                              type_tag + BINDERS_NUM_HEADER])
                        h_flags[CHR] = True
                    neo_rec = build_net_mhc_affinity_rec(chromosome=chromosome, pos=pos, h_flags=h_flags,
                                                         headers=headers, mhc_filtered_recs=neo_mhc_filtered_recs,
                                                         offset=offset, p_length=p_length, type_tag=type_tag)
                    type_tag = ORIG_PEP_TAG
                    avg_rank = orig_mhc_filtered_recs[chromosome][pos][p_length][offset].avg_rank
                    binders_num = orig_mhc_filtered_recs[chromosome][pos][p_length][offset].binders_num
                    rec.update({type_tag + RANKS_AVG_HEADER: avg_rank, type_tag + BINDERS_NUM_HEADER: binders_num})
                    if type_tag + RANKS_AVG_HEADER not in h_flags:
                        first_headers.extend([type_tag + RANKS_AVG_HEADER, type_tag + BINDERS_NUM_HEADER])
                        h_flags[type_tag + RANKS_AVG_HEADER] = True
                    orig_rec = build_net_mhc_affinity_rec(chromosome=chromosome, pos=pos, h_flags=h_flags,
                                                          headers=headers, mhc_filtered_recs=orig_mhc_filtered_recs,
                                                          offset=offset, p_length=p_length, type_tag=type_tag)

                    if sorted(neo_rec.values()) != sorted(orig_rec.values()):
                        identical = False
                    rec.update(neo_rec)
                    rec.update(orig_rec)
                    out_recs.append(rec)
                    if not identical:
                        diff_recs.append(rec)

    return out_recs, diff_recs, first_headers + sorted(headers)


def build_net_mhc_affinity_rec(chromosome, pos, h_flags, headers, mhc_filtered_recs, offset, p_length, type_tag):
    """
    This function is a help functio to build a single output rec.
    :param chromosome: the chromosome
    :param pos: The position on the chr.
    :param h_flags: a dict of flag to know what header were already added.
    :param headers: the collection of headers to add the headers to.
    :param mhc_filtered_recs: The parsed records
    :param offset: current offset
    :param p_length: current peptide length.
    :param type_tag: the type of the record -  neo or original.
    :return: the updated record (by reference).
    """
    rec = {}
    for hla, aff_rec in mhc_filtered_recs[chromosome][pos][p_length][offset].per_allele_affinity.iteritems():
        core_h = CORE_FORMAT % (hla, type_tag)
        rank_h = RANK_FORMAT % (hla, type_tag)
        nm_aff_h = NM_AFFINITY_FORMAT % (hla, type_tag)
        aff_strength_h = BIND_STRENGTH_FORMAT % (hla, type_tag)
        rank = float(aff_rec[RANK_HEADER])
        if rank <= STRONG_BINDING_PEPTIDES_THRESHOLD:
            aff_strength = STROG_BIND_VAL
        elif rank <= WEAK_BINDING_PEPTIDES_THRESHOLD:
            aff_strength = INTER_BIND_VAL
        else:
            aff_strength = WEAK_BIND_VAL
        if core_h not in h_flags:
            headers.extend([core_h, rank_h, nm_aff_h, aff_strength_h])
            h_flags[core_h] = True
        rec.update({core_h: aff_rec[CORE_HEADER], rank_h: aff_rec[RANK_HEADER],
                    nm_aff_h: aff_rec[NM_AFIINITY_HEADER], aff_strength_h: aff_strength})
    return rec


def parse_net_mhc_output(change_positions, mhc_data, p_affinities_i, p_alleles_i, p_id_i,
                         p_offset_i, window_p_i, avg_rank_i, binders_num_i):
    """
    This function parses file of NetMHC output and filters out all non localized peptides
     (that don't contain the mutation).
    :param change_positions: the position of changes in all peptides as a dict of
        [chromosome]=>[position]=>[change index]
    :param mhc_data: The data of the file.
    :param p_affinities_i: the indexes of each alleles affinities results in each record..
    :param p_alleles_i: the index of all alleles in each record..
    :param p_id_i: the index of the peptide id in each record.
    :param p_offset_i: the index of the offset in each record.
    :param window_p_i: the index of the window in each record.
    :param avg_rank_i: the index of the avg_rank_i in each record.
    :param binders_num_i: the index of the binders_num_i in each record.
    :return: the parsed records.
    """
    mhc_filtered_recs = {}

    for line in mhc_data[2:]:
        if line == "":
            continue

        offset = int(line[p_offset_i])
        chromosome, position = line[p_id_i].split(P_ID_SEP)[:2]

        try:
            change_p = int(change_positions[chromosome][position])
        except KeyError:
            logging.exception(FILTER_NET_MHC_WINDOWA_ERROR)

        p_length = len(line[window_p_i])

        if offset > change_p or change_p >= offset + p_length:
            continue  # this peptide doesn't include the changed AA.

        per_allele_aff = {}

        for hla in p_alleles_i.itervalues():
            nm_aff = line[p_affinities_i[hla][NM_AFIINITY_HEADER]]
            rank = line[p_affinities_i[hla][RANK_HEADER]]
            core = line[p_affinities_i[hla][CORE_HEADER]]
            per_allele_aff[hla] = {NM_AFIINITY_HEADER: nm_aff, RANK_HEADER: rank, CORE_HEADER: core}

        avg_rank = line[avg_rank_i]
        binders_num = line[binders_num_i]
        mhc_filtered_recs.setdefault(chromosome, {}).setdefault(position, {}).setdefault(p_length, {})[offset] \
            = NetMHCOutputRec(
            per_allele_aff,
            p_length,
            chromosome,
            position, avg_rank, binders_num)

    return mhc_filtered_recs


def get_net_mhc_indexes(alleles_line, headers_line):
    """
    This function gets the indexes of all atts in each record.
    :param alleles_line: The line from the NetMCH output containing the alelles (currently the first).
    :param headers_line: The line from the NetMCH output containing the headers (currently the second).
    :return: The indexes of all attributes, also connected to the corresponding HLA allele.
    """
    alleles_i = {}
    for i, cell_vall in enumerate(alleles_line):
        if cell_vall:
            alleles_i[i] = cell_vall
    offset_i = headers_line.index(OFFSET_HEADER)
    window_p_i = headers_line.index(WINDOW_CONTENT_HEADER)
    p_id_i = headers_line.index(P_ID_HEADER)
    try:
        avg_rank_i = headers_line.index(RANKS_AVG_HEADER)
        binders_num_i = headers_line.index(BINDERS_NUM_HEADER)
    except ValueError:
        avg_rank_i = headers_line.index(NET_MHC_PAN_RANKS_AVG_HEADER)
        binders_num_i = headers_line.index(NET_MHC_PAN_BINDERS_NUM_HEADER)

    affinites_i = {}
    iterate_alleles_i = sorted(alleles_i)
    for i, hla_i in enumerate(iterate_alleles_i):
        if i < len(iterate_alleles_i) -1:
            allele_headers = headers_line[hla_i:iterate_alleles_i[i + 1]]
        else:
            allele_headers = headers_line[hla_i:]
        nm_aff_i = allele_headers.index(NM_AFIINITY_HEADER) + hla_i
        rank_i = allele_headers.index(RANK_HEADER) + hla_i
        try:
            core_i = allele_headers.index(CORE_HEADER) + hla_i
        except ValueError:
            core_i = allele_headers.index(CORE_HEADER.lower()) + hla_i

        affinites_i[alleles_i[hla_i]] = {NM_AFIINITY_HEADER: nm_aff_i, RANK_HEADER: rank_i, CORE_HEADER: core_i}

    return alleles_i, offset_i, p_id_i, affinites_i, window_p_i, avg_rank_i, binders_num_i


if __name__ == '__main__':
    desc = 'A Script for extracting the correlated output of both mutated peptides and originals from NetMHC output'


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='NetMhcAnalysis', description=desc, formatter_class=MyFormatter)
    parser.add_argument('-in', '--neo_net_mhc_tab', metavar="neo NetMHC tab file", dest='neo_net_mhc_tab', nargs='?',
                        required=True, help='The input tab file of the mutated peptides.')
    parser.add_argument('-io', '--orig_net_mhc_tab', metavar="original NetMHC tab file", dest='orig_net_mhc_tab', nargs='?',
                        required=True, help='The input tab file of the original peptides.')
    parser.add_argument('-c', '--neopeptides_csv', metavar="neopeptides csv", dest='neopeptides_csv',
                        nargs='?', required=True,
                        help='The path to the csv summery file of the neo peptides from annovar_analysis.py.')
    parser.add_argument('-o', '--output_path', metavar="output path", dest='output_path', nargs='?', required=True,
                        help='The path for the output.')

    args = parser.parse_args()

    init_logging_dict(os.path.join(os.path.dirname(args.output_path), LOG_FILE))
    filter_net_mhc4_results(neo_net_mhc_tab=args.neo_net_mhc_tab,
                            orig_net_mhc_tab=args.orig_net_mhc_tab,
                            neopeptides_csv=args.neopeptides_csv,
                            output_path=args.output_path)

