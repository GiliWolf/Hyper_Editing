__author__ = 'Hillel'
#  This is made ot fit the current results, further generalization is required.

# =====================imports=====================#
import argparse
import logging
import os
from _csv import reader
from sys import path

if __name__ == "__main__":
    path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from Outputs.Outputers.CSVOutputer import CSVOutputer
from Commons.general_functions import init_logging_dict
from Commons.consts import SENSE_STRAND, MHC_CLASS_I, GROUP, MHC_ALLELE, MHC_CLASS, MHC_LOCUS, SAMPLE, SAMPLE_PATH,\
    DATA_SOURCE
from Commons.data_structs import Site, MHCAllele
from Commons.help_functions import determine_sep

from DataObjects.DataRecords.PeptidesAffinities import PeptidesAffinities
from DataObjects.KnownProperties.MHCAlleles import MHCAlleles
from DataObjects.BaseClasses.Sample import Sample

from ProcessingPlugins.EpitopeAnalysis.NetMHCAnalysis import CHR, POSITION
from ProcessingPlugins.EpitopeAnalysis.EpitopeAnalysisConsts import *

# =====================constants=====================#
# for defining which affinity is stronger than the other
AFFINITIES_ENUM = [WEAK_BIND_VAL, INTER_BIND_VAL, STROG_BIND_VAL, ]

AFFINITY_ASSESSMENT_SOURCE_NAME_ORIG = "AffinityAssessmentOrig"
AFFINITY_ASSESSMENT_SOURCE_NAME_NEO = "AffinityAssessmentNeo"
AFFINITY_ASSESSMENT_MHC_SOURCE = "AffinityAssessmentMHCs"

HLA_RANK_KEY = "HLAs"
CORE_PEPTIDE_KEY = "CorePep"
AFF_KEY = "Affinity"
ORIG_AFF_KEY = "OrigAffinity"
BIND_I = 0
CORE_I = 1

HLAMINER_ALLELE_FORMAT = "%(hla_class)s*%(allele)s"

HLA_HEADER = "HLA_Allele"
ORIG_BIND_HEADER = "OrigBindStrength"
NEO_BIND_HEADER = "NeoBindStrength"
BIND_DIFF_HEADER = "BindingDelta_Orig-Neo"
BIND_DIFF_PC_HEADER = "BindingDelta_PC_Orig-Neo"
BIND_DIFF_AVG_HEADER = "BindingDelta_Orig-Neo_AVG"
ORIG_CORE_HEADER = "OrigCore"
NEO_CORE_HEADER = "NeoCore"
GROUP_SIZE_FORMAT = "NumOfSamplesWithAlleleIn-%s"

SUM_ORIG_BIND_HEADER_S = "-".join([ORIG_BIND_HEADER, STROG_BIND_VAL])
SUM_ORIG_BIND_HEADER_I = "-".join([ORIG_BIND_HEADER, INTER_BIND_VAL])
SUM_ORIG_BIND_HEADER_W = "-".join([ORIG_BIND_HEADER, WEAK_BIND_VAL])
SUM_NEO_BIND_HEADER_S = "-".join([NEO_BIND_HEADER, STROG_BIND_VAL])
SUM_NEO_BIND_HEADER_I = "-".join([NEO_BIND_HEADER, INTER_BIND_VAL])
SUM_NEO_BIND_HEADER_W = "-".join([NEO_BIND_HEADER, WEAK_BIND_VAL])

SUM_REC_SEP = ";"

SPEC_FILE_NAME = "BindStrengthPerAllelePerGroupSpec.BindingAtLeast%s.csv"
SUM_FILE_NAME = "BindStrengthPerAllelePerGroupSum.BindingAtLeast%s.csv"

#  ---- logging and errors------

LOG_FILE = 'affinityAssessment.log'
HLA_FILE_UNAVAILABLE = 'Cannot Open HLA File! (%s)'
NET_MHC_FILE_UNAVAILABLE = 'Cannot Open NetNHC Analysis File! (%s)'
MISSING_HLA_ALLELE_MSG = "HLA Allele %s Was Found In NetMHC Analysis But Not In Samples, Skipping!"


# =====================functions=====================#


def get_hla_allels_per_sample(hlas_per_sample_file, merge_sub_strs_dict={}, merge_ignores=[]):
    """
    This function associate each sample with its HLA alleles
    :param hlas_per_sample_file: the file containing the HLA alleles of each sample (from HLAMinerAnalyser)
    :return: a dictionary of HLA allele to groups to samples.
    :rtype: C{dict} of C{str} to C{dict} of C{str} to C{list} of C{str}
    """
    # TODO: convert to entity load...
    try:
        data = [line for line in reader(open(hlas_per_sample_file))]
    except IOError:
        logging.exception(HLA_FILE_UNAVAILABLE % hlas_per_sample_file)
        return None

    group_i = data[0].index(GROUP)
    allele_i = data[0].index(MHC_ALLELE)
    allele_class_i = data[0].index(MHC_CLASS)
    allele_locus_i = data[0].index(MHC_LOCUS)
    sample_i = data[0].index(SAMPLE)


    hla_allels = {}
    for line in data[1:]:
        allele = line[allele_i]
        allele_class = line[allele_class_i]
        allele_locus = line[allele_locus_i]


        mhc_allele = MHCAllele(allele_class, allele_locus, allele)

        hla_allels.setdefault(mhc_allele, {}).setdefault(line[group_i], []).append(line[sample_i])

        if line[group_i] in merge_ignores:
            continue
        for sub_str in merge_sub_strs_dict:
            if sub_str in line[group_i]:
                merged_group = merge_sub_strs_dict[sub_str]
                hla_allels.setdefault(mhc_allele, {}).setdefault(merged_group, []).append(line[sample_i])
    return hla_allels


def get_affinities(net_mhc_analysis_output):
    neo_affinities = PeptidesAffinities(AFFINITY_ASSESSMENT_SOURCE_NAME_NEO)
    orig_affinities = PeptidesAffinities(AFFINITY_ASSESSMENT_SOURCE_NAME_ORIG)
    try:
        data = open(net_mhc_analysis_output).readlines()
    except IOError:
        logging.exception(NET_MHC_FILE_UNAVAILABLE % net_mhc_analysis_output)
        return None

    sep = determine_sep(data[0])
    data = [line.split(sep) for line in data]

    hla_aff = {}

    chr_i, pos_i, hla_indexes, offset_i = get_indexes(data)

    for line in data[1:]:
        chrom = line[chr_i]
        pos = int(line[pos_i])
        offset = int(line[offset_i])
        for hla in hla_indexes:
            """:type hla MHCAllele"""
            hla_orig_aff_rank = line[hla_indexes[hla][HLA_RANK_KEY][ORIG_PEP_TAG]]
            hla_neo_aff_rank = line[hla_indexes[hla][HLA_RANK_KEY][NEO_PEP_TAG]]
            core_orig = line[hla_indexes[hla][CORE_PEPTIDE_KEY][ORIG_PEP_TAG]]
            core_neo = line[hla_indexes[hla][CORE_PEPTIDE_KEY][NEO_PEP_TAG]]
            orig_aff = float(line[hla_indexes[hla][AFF_KEY][ORIG_PEP_TAG]])
            neo_aff = float(line[hla_indexes[hla][AFF_KEY][NEO_PEP_TAG]])
            # orig_neo_aff_delta = float(line[hla_indexes[hla][AFF_KEY][ORIG_PEP_TAG]]) - \
            #                      float(line[hla_indexes[hla][AFF_KEY][NEO_PEP_TAG]])
            # if affinity_rank(hla_orig_aff_rank, hla_neo_aff_rank, minimal_aff):
            orig_affinities.add_peptide_aff(region=chrom,
                                           start=pos - 1,
                                           end=pos,
                                           strand=SENSE_STRAND, # TODO: strand...
                                           mhc_class=hla.mhc_class,
                                           locus_name=hla.locus,
                                           allele_name=hla.allele,
                                           offset=offset,
                                           ligand_length=len(core_orig),
                                           rank=0, # TODO: get the right input for this
                                           affinity_rank=hla_orig_aff_rank,
                                           affinity=orig_aff,
                                           binding_core=core_orig)
            neo_affinities.add_peptide_aff(region=chrom,
                                           start=pos - 1,
                                           end=pos,
                                           strand=SENSE_STRAND, # TODO: strand...
                                           mhc_class=hla.mhc_class,
                                           locus_name=hla.locus,
                                           allele_name=hla.allele,
                                           offset=offset,
                                           ligand_length=len(core_neo),
                                           rank=0, # TODO: get the right input for this
                                           affinity_rank=hla_neo_aff_rank,
                                           affinity=neo_aff,
                                           binding_core=core_neo)
            # hla_aff.setdefault(hla, {}).setdefault(chrom, {}).setdefault(pos, {}) \
            #         .setdefault(offset, {ORIG_PEP_TAG: [hla_orig_aff_rank, core_orig], NEO_PEP_TAG: [hla_neo_aff_rank, core_neo],
            #                              AFF_KEY: orig_neo_aff_delta, ORIG_AFF_KEY: line[hla_indexes[hla][AFF_KEY][ORIG_PEP_TAG]]})
    return hla_aff


def get_indexes(data):
    """
    This function derive the indexes of each value type in the data
    :param data: the data from NetMHCAnalysis output.
    :return:  the indexes for chromosome, position, peptide and every hla allele.
    """
    hla_bind_rank_format = BIND_STRENGTH_FORMAT.split('_')[0]
    core_format = CORE_FORMAT.split('_')[0]
    aff_format = NM_AFFINITY_FORMAT.split('_')[0]

    hla_indexes = {}
    chr_i = data[0].index(CHR)
    pos_i = data[0].index(POSITION)
    offset_i = data[0].index(OFFSET_FORMAT)

    for i, header in enumerate(data[0]):
        if hla_bind_rank_format in header:
            type_tag = ORIG_PEP_TAG if ORIG_PEP_TAG in header else NEO_PEP_TAG
            hla = transform_allele(header.replace(hla_bind_rank_format, '').replace(ORIG_PEP_TAG, '').replace(NEO_PEP_TAG, '').replace("_", ""))
            hla_indexes.setdefault(hla, {}).setdefault(HLA_RANK_KEY, {})[type_tag] = i
        elif core_format in header:
            type_tag = ORIG_PEP_TAG if ORIG_PEP_TAG in header else NEO_PEP_TAG
            hla = transform_allele(header.replace(core_format, '').replace(ORIG_PEP_TAG, '').replace(NEO_PEP_TAG, '').replace("_", ""))
            hla_indexes.setdefault(hla, {}).setdefault(CORE_PEPTIDE_KEY, {})[type_tag] = i
        elif aff_format in header:
            type_tag = ORIG_PEP_TAG if ORIG_PEP_TAG in header else NEO_PEP_TAG
            hla = transform_allele(header.replace(aff_format, '').replace(ORIG_PEP_TAG, '').replace(NEO_PEP_TAG, '').replace("_", ""))
            hla_indexes.setdefault(hla, {}).setdefault(AFF_KEY, {})[type_tag] = i

    return chr_i, pos_i, hla_indexes, offset_i


def transform_allele(hla_allele):
    """
    This function "translates" from HLA standard format to HLAMiner format.
    :param hla_allele: The hla allele
    :return: the HLAMiner format.
    """
    allele = hla_allele.split("-")[1]
    hla_locus = "HLA_" + allele[0]
    allele_class = MHC_CLASS_I # TODO: add detctabilty (no need if loading is done)
    # hla_locus = allele[0]

    return MHCAllele(allele_class, hla_locus, allele[0] + "*" + allele[1:])


def asses_affinity(hlas_per_sample_file, net_mhc_analysis_outputs, output_dir="."):
    hlas = get_hla_allels_per_sample(hlas_per_sample_file)
    # for minimal_affinity in [STROG_BIND_VAL, INTER_BIND_VAL, WEAK_BIND_VAL]:
    create_output(hlas,  net_mhc_analysis_outputs, output_dir)


def create_output(hlas,  net_mhc_analysis_outputs, output_dir):
    affinities = {}
    for net_mhc_analysis_output in net_mhc_analysis_outputs:
        get_affinities(net_mhc_analysis_output)
    orig_affs = PeptidesAffinities(AFFINITY_ASSESSMENT_SOURCE_NAME_ORIG)
    neo_affs = PeptidesAffinities(AFFINITY_ASSESSMENT_SOURCE_NAME_NEO)
    dump_path_orig = os.path.join(output_dir, "OrigAffs.csv")
    dump_path_neo = os.path.join(output_dir, "NeoAffs.csv")
    orig_affs.dump_me(dump_path_orig)
    PeptidesAffinities.remove_entity(orig_affs.record_id)
    neo_affs.dump_me(dump_path_neo)
    # PeptidesAffinities.remove_entity(orig_affs.record_id)

    # outputer = CSVOutputer()
    # spec_headers = [CHR, POSITION, OFFSET_FORMAT, HLA_HEADER, BIND_DIFF_HEADER, ORIG_BIND_HEADER, NEO_BIND_HEADER,
    #                 ORIG_CORE_HEADER, NEO_CORE_HEADER]
    # sum_headers = [CHR, POSITION, OFFSET_FORMAT, HLA_HEADER, BIND_DIFF_AVG_HEADER, BIND_DIFF_HEADER, BIND_DIFF_PC_HEADER,
    #                SUM_ORIG_BIND_HEADER_S, SUM_ORIG_BIND_HEADER_I, SUM_ORIG_BIND_HEADER_W, SUM_NEO_BIND_HEADER_S,
    #                SUM_NEO_BIND_HEADER_I, SUM_NEO_BIND_HEADER_W]
    # added_groups = {}
    # spec_recs = []
    # sum_recs = []
    # for hla, chroms in affinities.iteritems():
    #     for chrom, positions in chroms.iteritems():
    #         for pos, offsets in positions.iteritems():
    #             sum_offset = SUM_REC_SEP.join(sorted(offsets.keys()))
    #             sum_bind_orig = {STROG_BIND_VAL: 0, INTER_BIND_VAL: 0, WEAK_BIND_VAL: 0}
    #             sum_bind_neo = {STROG_BIND_VAL: 0, INTER_BIND_VAL: 0, WEAK_BIND_VAL: 0}
    #             deltas = []
    #             deltas_percentages = []
    #             groups_counter = {}
    #             for offset in sorted(offsets.keys()):
    #                 pep_aff = offsets[offset]
    #                 trans_hla = transform_allele(hla)
    #                 orig_bind = pep_aff[ORIG_PEP_TAG][BIND_I]
    #                 neo_bind = pep_aff[NEO_PEP_TAG][BIND_I]
    #                 orig_core = pep_aff[ORIG_PEP_TAG][CORE_I]
    #                 neo_core = pep_aff[NEO_PEP_TAG][CORE_I]
    #                 orig_neo_aff_delta = pep_aff[AFF_KEY]
    #                 delta_precentage = 100 * float(pep_aff[AFF_KEY]) /float(pep_aff[ORIG_AFF_KEY])
    #
    #                 sum_bind_neo[neo_bind] += 1
    #                 sum_bind_orig[orig_bind] += 1
    #
    #                 deltas.append(orig_neo_aff_delta)
    #                 deltas_percentages.append(delta_precentage)
    #
    #                 rec = {CHR: chrom, POSITION: pos, OFFSET_FORMAT: offset, BIND_DIFF_HEADER: str(orig_neo_aff_delta),
    #                        HLA_HEADER: trans_hla, ORIG_BIND_HEADER: orig_bind, NEO_BIND_HEADER: neo_bind,
    #                        ORIG_CORE_HEADER: orig_core, NEO_CORE_HEADER: neo_core}
    #                 # import pdb;pdb.set_trace()
    #                 if trans_hla not in hlas:
    #                     logging.warning(MISSING_HLA_ALLELE_MSG % trans_hla)
    #                     continue
    #
    #                 t_rec = rec.copy()
    #                 for group, samples in hlas[trans_hla].iteritems():
    #                     group_header = GROUP_SIZE_FORMAT % group
    #
    #                     if added_groups.get(group_header, True):
    #                         spec_headers.append(group_header)
    #                         sum_headers.append(group_header)
    #                         added_groups[group_header] = False
    #
    #                     t_rec[group_header] = str(len(samples))
    #
    #                     groups_counter.setdefault(group, 0)
    #                     groups_counter[group] = + len(samples)
    #
    #                 spec_recs.append(t_rec)
    #             sum_rec = {CHR: chrom, POSITION: pos, OFFSET_FORMAT: sum_offset, HLA_HEADER: trans_hla,
    #                        BIND_DIFF_AVG_HEADER: str(sum(deltas) / float(len(deltas))),
    #                        BIND_DIFF_HEADER: SUM_REC_SEP.join([str(d) for d in deltas]),
    #                        BIND_DIFF_PC_HEADER: SUM_REC_SEP.join([str(d) for d in deltas_percentages]),
    #                        SUM_ORIG_BIND_HEADER_S: str(sum_bind_orig[STROG_BIND_VAL]),
    #                        SUM_ORIG_BIND_HEADER_I: str(sum_bind_orig[INTER_BIND_VAL]),
    #                        SUM_ORIG_BIND_HEADER_W: str(sum_bind_orig[WEAK_BIND_VAL]),
    #                        SUM_NEO_BIND_HEADER_S: str(sum_bind_neo[STROG_BIND_VAL]),
    #                        SUM_NEO_BIND_HEADER_I: str(sum_bind_neo[INTER_BIND_VAL]),
    #                        SUM_NEO_BIND_HEADER_W: str(sum_bind_neo[WEAK_BIND_VAL])}
    #             for group in groups_counter:
    #                 sum_rec[GROUP_SIZE_FORMAT % group] = str(groups_counter[group])
    #             sum_recs.append(sum_rec)
    # outputer.output([os.path.join(output_dir, SPEC_FILE_NAME % minimal_affinity)], spec_headers, spec_recs)
    # outputer.output([os.path.join(output_dir, SUM_FILE_NAME %minimal_affinity)], sum_headers, sum_recs)


if __name__ == '__main__':
    desc = 'A script for assesing the affinty per sample from the netMHCAnalysis output'


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='AffinityAssesment', description=desc, formatter_class=MyFormatter)
    parser.add_argument('-i', '--net_mhc_analysis_outputs', metavar="NetMHCAnalysis outputs",
                        dest='net_mhc_analysis_outputs', nargs='*', required=True,
                        help='The input file(s). Files from NetMHCAnalysis')
    parser.add_argument('-hla', '--hlas_per_sample_file', metavar="hla allles per sample file",
                        dest='hlas_per_sample_file',
                        nargs='?', required=True, help='A file containing hla alleles per sample.')
    parser.add_argument('-o', '--output_dir', metavar="output dir", dest='output_dir', nargs='?', required=True,
                        help='The path for the output.')

    args = parser.parse_args()

    init_logging_dict(os.path.join(args.output_dir, LOG_FILE))
    asses_affinity(hlas_per_sample_file=args.hlas_per_sample_file,
                   net_mhc_analysis_outputs=args.net_mhc_analysis_outputs,
                   output_dir=args.output_dir)

