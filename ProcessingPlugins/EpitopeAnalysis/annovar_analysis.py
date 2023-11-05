__author__ = 'Hillel'
# =====================imports=====================#
# TODO: move most to converters...
import argparse
import logging
import os
import pickle
import re
from collections import namedtuple
from csv import reader
from sys import exit, path

if __name__ == "__main__":
    path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from Outputs.Outputers.CSVOutputer import CSVOutputer
from Outputs.Outputers.RawOutputer import RawOutputer
from Commons.general_functions import init_logging_dict
from Commons.consts import MismatchesAndRefsEnum
# =====================consts=====================#
NA = "NA"
EXON_ID_SEP = "."

# indexes for readability and import simplicity
GENE_ID_I = 0
GENE_NAME_I = 1
GENE_MRNA_EXON_I = 2

# ---- extract_protein_data_from_ncbi_gpffs constants ----
REC_REG = re.compile(
    r"LOCUS *(NP_\d+) *(\d+) aa.*?VERSION *(NP_[\d\.]+).*?DBSOURCE    REFSEQ: accession (NM_[\d\.]+).*?ORIGIN(.*?)//",
    re.DOTALL)
SEQ_SUB_REG = re.compile(r"[^a-zA-Z]")

PROTEIN_ID_I = 0
NUM_OF_AA_I = 1
PROTEIN_VERSION_ID_I = 2
GENE_ID = 3
PROTEIN_SEQ_I = 4

# ---- parse_annovar_change_str constants ----
ANNOVAR_CHANGE_STR_SEP = ":"
CHANGED_VERSION_SEP1 = ";"
CHANGED_VERSION_SEP2 = ","
EXON = "exon"

P_CHANGE_RE = re.compile(r"([A-z])(\d+)([A-z])")

# ---- get_annovar_changed_peptides constants ----
CHR = "Chr"
POSITION = "Position"
CHANGE_POSITION = "ChangePosition"
ORIG_AA = "OriginalAA"
NEO_AA = "NeoAA"
PEPTIDE = "Peptide"
GENE_NAME = "GeneName"
CHANGE_STR = "AAChange.refgene"
CHANGE_TYPE = "ChangeType"
P_ID_SEP = "-"
NEO_PEPTIDES_HEADERS = [CHR, POSITION, GENE_NAME, PEPTIDE, CHANGE_POSITION, ORIG_AA, NEO_AA]
EXONIC_SITES_NOT_NONSYN_HEADERS = [CHR, POSITION, GENE_NAME, CHANGE_TYPE, CHANGE_STR]
EXONIC_SITES_NOT_NONSYN_FILE = os.path.join(r"%s","ExonicSitesNotNonSyn.csv")
NEO_PEPTIDES_FILE_FORMAT = os.path.join(r"%(mismatch)s","WindowSize%(window_size)s","NeoPeptides.%(suffix)s")
NEO_PEPTIDES_ORG_FILE_FORMAT = os.path.join(r"%(mismatch)s","WindowSize%(window_size)s","OrigPeptides.%(suffix)s")
NON_SYNONYMOUS = "nonsynonymous SNV"

FAILED_P = "ProcessingFailure"

#  ---- get_exonic_recs_from_annovar_csv constants ----
SITE_FUNCTION_HEADER = "Func.refgene"
EXONIC_FUNC_HEADER = "ExonicFunc.refgene"
CHANGE_STR_HEADER = "AAChange.refgene"
GENE_NAME_HEADER = "Gene.refgene"
CHR_HEADER = "Chr"
POS_HEADER = "Start"
REF_HEADER = "Ref"
ALT_HEADER = "Alt"

EXONIC = "exonic"

# ---- logging msgs and consts ----
GPPF_PROCESSING_STARTED = "Beginning Loading Protein Data From NCBI GPFFs In Folder %s"
GPPF_PROCESSING_DONE = "Folder %s Was Successfully Processed. Pickle Saved At %s"
GPPF_FILE_TRY = "Trying To Load GPFF File %s!"
GPPF_FILE_FAILED = "Failed To Load GPFF File %s. Exiting!"
GPPF_PICKLE_FAIL = "Failed Pickling Data!"

P_PICKLE_LOAD_FAILED = "Could Not Unpickle Proteins Data Pickle From %s, Exiting!"
P_PICKLE_FILE_ERROR = "Could Not Open Proteins Data Pickle, Exiting!"

ANNO_CSV_PARSE_FAILED = "Failed Parsing Annovar CSV Outputer (At %s), Exiting!"

CHANGED_PEPTIDE_RET_FAIL = "Failed To Retrieve Peptide For Change String %s"
PROTEIN_PROCESSING_ERROR = "Failed To Implant Change In Using %s, Position Is Off Protein!"

NCBI_DIR_WRONG = "The Folder Provided To Load NCBI From Doesn't Exist! (%s)"

LOG_FILE = "neo_epitopes_analysis.log"
# =====================classes=======================#

# This named tuple hold the data from the ncbi gpffs (all are strings).
# protein_id: The ncbi protein acces ID.
# num_of_aa: The number of amino acids in the protein.
# protein_version_id: The protein ncbi access ID including version suffix.
# gene_id: RefSeq gene name.
# gene_version: RefSeq gene name including version suffix.
# protein_seq: The amino acids sequence of the protein.
ProteinData = namedtuple("ProteinData", "protein_id num_of_aa protein_version_id gene_id gene_version protein_seq")

# This named tuple hold the data from the annovar csv output.
# chr: The site chromosome.
# position: The position in the chromosome.
# exonic_function: The function of change (e.g. synonymous).
# change_str: The change string from annovar . e.g. HIST2H2BE:NM_003528:exon1:c.T19G:p.S7A
AnnovarExonRec = namedtuple("AnnovarExonRec", "chr position gene_name exonic_function change_str")


# =====================functions=====================#


def extract_protein_data_from_ncbi_gpffs(gpffs_folder, pickle_output_path):
    """
    This function process the gpff files from (on 19/05/16) ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot to
     extract from them proteins data and saves the data as pickle file.
    :param gpffs_folder: The root folder containing the gpffs.
    :param pickle_output_path: The path to save the pickle at.
    :return: True if succeeded,False if not.
    """
    prot = {}
    logging.info(GPPF_PROCESSING_STARTED % gpffs_folder)
    for root, dirs, files in os.walk(gpffs_folder):
        for f in files:
            logging.debug(GPPF_FILE_TRY % f)
            try:
                data = open(os.path.join(root, f)).read()
                recs = REC_REG.findall(data)
                for rec in recs:
                    gene_id = get_id_from_version_id(rec[GENE_ID].split(".")[0])
                    prot[gene_id] = ProteinData(rec[PROTEIN_ID_I], rec[NUM_OF_AA_I], rec[PROTEIN_VERSION_ID_I],
                                                gene_id, rec[GENE_ID], SEQ_SUB_REG.sub("", rec[PROTEIN_SEQ_I]))
            except:
                logging.exception(GPPF_FILE_FAILED % f)
                return False
    try:
        with open(pickle_output_path, 'wb') as p:
            pickle.dump(prot, p)
    except:
        logging.exception(GPPF_PICKLE_FAIL)
        return False
    logging.info(GPPF_PROCESSING_DONE % (gpffs_folder, pickle_output_path))

    return True


def get_id_from_version_id(exon_id):
    """
    This function extracts the access ID from the exon ID.
    :param exon_id: The exon ID to process.
    :return: The ID.
    """
    return exon_id.split(EXON_ID_SEP)[0]


def implant_amino_acid_change(gene_id, position, change, protein_data):
    """
    This function change an amino acid in a protein sequence.
    :param gene_id: The RefSeq gene name of the protein.
    :param position: Position in the protein, as present in annovar output (first aa index is 1).
    :param change: What to change the aa to.
    :param protein_data: A dictionary holding the gene names to their proteins' data.
    :type protein_data: C{dict} of C{str} to L{ProteinData}
    :return: The changed protein. If errors where encountered - None.
    """
    if gene_id not in protein_data:
        return None

    if int(protein_data[gene_id].num_of_aa) < position:
        return None

    seq = protein_data[gene_id].protein_seq

    seq = seq[:position - 1] + change + seq[position:]

    return seq.upper()


def parse_annovar_change_str(annovar_change_str):
    """
    This function parses the annovar change string (under AAChange.refgene) and returns the data from it.
    :param annovar_change_str: The change string from annovar . e.g. HIST2H2BE:NM_003528:exon1:c.T19G:p.S7A
    :return: The gene name C{str}, RefSeq name C{str}, exon number C{int}, protein change position C{int} and
    protein change C{str}.
    """
    gene_name, ref_seq_name, exon_number, _, raw_change = annovar_change_str.split(ANNOVAR_CHANGE_STR_SEP)
    exon_number = int(exon_number.replace(EXON, ""))
    orig, position, to = P_CHANGE_RE.findall(raw_change)[0]

    return gene_name, ref_seq_name, exon_number, int(position), to


def get_peptide_changed_seq(annovar_change_str, protein_data, window_size=9):
    """
    This function generates the changed peptide from the change string (e.g. HIST2H2BE:NM_003528:exon1:c.T19G:p.S7A).
    :param annovar_change_str:  The change string from annovar . e.g. HIST2H2BE:NM_003528:exon1:c.T19G:p.S7A
    :param protein_data: A dictionary holding the gene names to their proteins' data.
    :type protein_data: C{dict} of C{str} to L{ProteinData}
    :param window_size: The number of amino acid from each size of the change to retrieve.
    :return: The peptide with and without the change of length window_size*2 +1 and the index of the changed AA.
     If fails - None.
    """
    try:
        gene_name, ref_seq_name, exon_number, position, change = parse_annovar_change_str(annovar_change_str)
    except:
        logging.exception(CHANGED_PEPTIDE_RET_FAIL % annovar_change_str)
        return None, None

    whole_seq = implant_amino_acid_change(gene_id=ref_seq_name, position=position, change=change,
                                          protein_data=protein_data)

    if whole_seq:
        # no need to check again for sanity for the original seq
        orig_seq = protein_data[ref_seq_name].protein_seq

        if position - window_size - 1 >= 0:
            win_start_i = position - window_size - 1
            change_pos = window_size
        else:
            win_start_i = 0
            change_pos = position
        win_end_i = position + window_size if position + window_size <= len(whole_seq) else len(whole_seq)

        return whole_seq[win_start_i:win_end_i], orig_seq[win_start_i:win_end_i].upper(), change_pos
    else:
        logging.warning(PROTEIN_PROCESSING_ERROR % annovar_change_str)
        return None, None, None


def get_annovar_changed_peptides(annovar_csv_path, protein_data_pickle, output_dir, window_size=9, split_fasta=0,
                                 not_stranded=False):
    try:
        exonic_non_syn_recs, exonic_other_recs = get_exonic_recs_from_annovar_csv(annovar_csv_path, not_stranded)
    except (ValueError, IndexError):
        logging.exception(ANNO_CSV_PARSE_FAILED % annovar_csv_path)
        return
    try:
        with open(protein_data_pickle, 'rb') as p_pkl:
            protein_data_pickle = pickle.load(p_pkl)
    except pickle.UnpicklingError:
        logging.exception(P_PICKLE_LOAD_FAILED % p_pkl)
    except IOError:
        logging.exception(P_PICKLE_FILE_ERROR)

    outputer = CSVOutputer()
    raw_outputer = RawOutputer(pprint_flag=False)
    for mismatch in exonic_non_syn_recs:
        orig_recs = []
        recs = []
        for exon_rec in exonic_non_syn_recs[mismatch]:
            if not isinstance(exon_rec, AnnovarExonRec):
                continue
            neo_p, orig_p, change_position = get_peptide_changed_seq(annovar_change_str=exon_rec.change_str,
                                                                     protein_data=protein_data_pickle,
                                                                     window_size=window_size)

            if None is neo_p:
                neo_p = FAILED_P
                orig_p = FAILED_P
                orig_aa = FAILED_P
                new_aa = FAILED_P
            else:
                orig_aa = orig_p[change_position]
                new_aa = neo_p[change_position]

            recs.append({CHR: exon_rec.chr, POSITION: exon_rec.position, GENE_NAME: exon_rec.gene_name, PEPTIDE: neo_p,
                         CHANGE_POSITION: str(change_position), ORIG_AA: orig_aa, NEO_AA: new_aa})
            orig_recs.append({CHR: exon_rec.chr, POSITION: exon_rec.position, GENE_NAME: exon_rec.gene_name,
                              PEPTIDE: orig_p, CHANGE_POSITION: str(change_position)})
        out_file = os.path.join(output_dir, NEO_PEPTIDES_FILE_FORMAT % dict(mismatch=mismatch,
                                                                            window_size=str(window_size), suffix="csv"))
        outputer.output([out_file], NEO_PEPTIDES_HEADERS, recs)

        out_file = os.path.join(output_dir, NEO_PEPTIDES_ORG_FILE_FORMAT % dict(mismatch=mismatch,
                                                                                window_size=str(window_size),
                                                                                suffix="csv"))
        outputer.output([out_file], NEO_PEPTIDES_HEADERS, orig_recs)

        fasta = []
        orig_fasta = []
        for rec in recs:
            if rec[PEPTIDE] == FAILED_P:
                continue
            fasta.append(">" + P_ID_SEP.join([rec[CHR], rec[POSITION], rec[GENE_NAME]]) + "\n" + rec[PEPTIDE])
        for rec in orig_recs:
            if rec[PEPTIDE] == FAILED_P:
                continue
            orig_fasta.append(">" + P_ID_SEP.join([rec[CHR], rec[POSITION], rec[GENE_NAME]]) + "\n" + rec[PEPTIDE])

        if split_fasta == 0:
            out_file = os.path.join(output_dir, NEO_PEPTIDES_FILE_FORMAT % dict(mismatch=mismatch,
                                                                                window_size=str(window_size),
                                                                                suffix="fasta"))
            raw_outputer.output([out_file], "\n".join(fasta))
            out_file = os.path.join(output_dir, NEO_PEPTIDES_ORG_FILE_FORMAT % dict(mismatch=mismatch,
                                                                                    window_size=str(window_size),
                                                                                    suffix="fasta"))
            raw_outputer.output([out_file], "\n".join(orig_fasta))
        else:
            for i in xrange(1 + len(fasta) / split_fasta):
                out_file = os.path.join(output_dir, NEO_PEPTIDES_FILE_FORMAT % dict(mismatch=mismatch,
                                                                                    window_size=str(window_size),
                                                                                    suffix="[%s].fasta" % i))
                raw_outputer.output([out_file], "\n".join(fasta[split_fasta * i:split_fasta * (i + 1)]))
            for i in xrange(1 + len(orig_fasta) / split_fasta):
                out_file = os.path.join(output_dir, NEO_PEPTIDES_ORG_FILE_FORMAT % dict(mismatch=mismatch,
                                                                                        window_size=str(window_size),
                                                                                        suffix="[%s].fasta" % i))
                raw_outputer.output([out_file], "\n".join(orig_fasta[split_fasta * i:split_fasta * (i + 1)]))

    for mismatch in exonic_other_recs:
        recs = []
        for exon_rec in exonic_other_recs[mismatch]:
            if not isinstance(exon_rec, AnnovarExonRec):
                continue

            recs.append({CHR: exon_rec.chr, POSITION: exon_rec.position, GENE_NAME: exon_rec.gene_name,
                         CHANGE_TYPE: exon_rec.exonic_function, CHANGE_STR: exon_rec.change_str})

        out_file = os.path.join(output_dir, EXONIC_SITES_NOT_NONSYN_FILE % mismatch)
        outputer.output([out_file], EXONIC_SITES_NOT_NONSYN_HEADERS, recs)


def get_exonic_recs_from_annovar_csv(annovar_csv_path, not_stranded):
    """
    This function parses the annovar csv output and returns only exonic records.
    :param annovar_csv_path: The path to the output.
    :param not_stranded: If True will treat the data as not stranded uniting complementary mismatches.
    :return: A list of the records with only the relevant data.
    :rtype: C{list} of L{AnnovarExonRec}
    """
    exonic_non_syn_recs = {}
    exonic_other_recs = {}

    with open(annovar_csv_path) as annovar_csv:
        csv_r = reader(annovar_csv)
        is_headers = True
        for rec in csv_r:
            if is_headers:
                is_headers = False
                site_function_i = rec.index(SITE_FUNCTION_HEADER)
                exonic_func_i = rec.index(EXONIC_FUNC_HEADER)
                change_str_i = rec.index(CHANGE_STR_HEADER)
                chr_i = rec.index(CHR_HEADER)
                pos_i = rec.index(POS_HEADER)
                gene_name_i = rec.index(GENE_NAME_HEADER)
                ref_i = rec.index(REF_HEADER)
                alt_i = rec.index(ALT_HEADER)
                continue

            if rec[site_function_i] != EXONIC:
                continue
            # The first record is the newer (so far)
            if CHANGED_VERSION_SEP1 in rec[change_str_i]:
                gene_name = rec[gene_name_i].split(CHANGED_VERSION_SEP1)[0]
            else:
                gene_name = rec[gene_name_i].split(CHANGED_VERSION_SEP2)[0]

            if CHANGED_VERSION_SEP1 in rec[change_str_i]:
                change_str = rec[change_str_i].split(CHANGED_VERSION_SEP1)[0]
            else:
                change_str = rec[change_str_i].split(CHANGED_VERSION_SEP2)[0]

            if CHANGED_VERSION_SEP1 in rec[exonic_func_i]:
                exonic_func = rec[exonic_func_i].split(CHANGED_VERSION_SEP1)[0]
            else:
                exonic_func = rec[exonic_func_i].split(CHANGED_VERSION_SEP2)[0]

            chromosome = "chr" + rec[chr_i] if "chr" not in rec[chr_i] else rec[chr_i]
            pos = rec[pos_i]

            mm = MismatchesAndRefsEnum.MISMATCH_TYPE_PARSE.get(rec[ref_i] + rec[alt_i], "NA")
            if not_stranded:
                mm = MismatchesAndRefsEnum.HE_NOT_STRANDED_TRANSFORM_DICT[mm]

            if exonic_func != NON_SYNONYMOUS:
                exonic_other_recs.setdefault(mm, []).append(
                    AnnovarExonRec(chromosome, pos, gene_name, exonic_func, change_str))
            else:
                exonic_non_syn_recs.setdefault(mm, []).append(
                    AnnovarExonRec(chromosome, pos, gene_name, exonic_func, change_str))

    return exonic_non_syn_recs, exonic_other_recs


if __name__ == '__main__':
    desc = 'A Script for extracting the change peptides from annovar changes prediction for downstream analysis.'


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='AnnovarAnalysis', description=desc, formatter_class=MyFormatter)
    parser.add_argument('-i', '--input_csv', metavar="input csv path", dest='input_csv', nargs='?', required=True,
                        help='The input csv file (The annovar output).')
    parser.add_argument('-p', '--protein_data_pickle', metavar="protein data pickle path", dest='protein_data_pickle',
                        nargs='?', required=True,
                        help='The path to the pickle containing the proteins data (From the NCBI GPFFs).')
    parser.add_argument('-o', '--output_dir', metavar="output path", dest='output_dir', nargs='?', required=True,
                        help='The dir for the output.')
    parser.add_argument('-w', '--window_size', metavar="peptide window size", dest='window_size', nargs='?',
                        required=False, default=9, type=int,
                        help='The number of amino acids to get from each side of the change.')
    parser.add_argument('-l', '--load_pickle', dest='load_pickle', required=False, default=False, action="store_true",
                        help='If set, will try to load the proteins data base from the path provided'
                             ' with <ncbi_db_dir> and will save it at <protein_data_pickle> ')
    parser.add_argument('-n', '--ncbi_db_dir', metavar="ncbi db dir", dest='ncbi_db_dir', nargs='?', required=False,
                        help='The root folder containing the gpffs from NCBI.')
    parser.add_argument('-s', '--split_fasta', metavar="split peptide fasta size", dest='split_fasta', nargs='?',
                        required=False, type=int, default=0,
                        help='If given will split the fasta file to number of recs in each file..')
    parser.add_argument('--not_stranded', dest='not_stranded', required=False, default=False, action="store_true",
                        help='If set, will try treat the data as not stranded, uniting complemetary mismatches')

    args = parser.parse_args()

    if args.load_pickle:
        assert os.path.isdir(args.ncbi_db_dir), NCBI_DIR_WRONG % args.ncbi_db_dir

        if not extract_protein_data_from_ncbi_gpffs(args.ncbi_db_dir, args.protein_data_pickle):
            exit(1)
    init_logging_dict(os.path.join(args.output_dir, LOG_FILE))
    get_annovar_changed_peptides(annovar_csv_path=args.input_csv, protein_data_pickle=args.protein_data_pickle,
                                 output_dir=args.output_dir, window_size=args.window_size, split_fasta=args.split_fasta,
                                 not_stranded=args.not_stranded)
