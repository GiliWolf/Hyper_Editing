__author__ = 'Hillel'
"""
This module contains the consts for the wrapper of Hagit Porath's HE tool.
"""
# =====================imports===================#
import re

# =====================constants===================#

# region input filename regexes
# ---- regexes for file names ---- TODO:change
FILENAME_FASTQ_RE = re.compile(r"([\w_\-]+?)_?([12]?)\.fastq$")
FILENAME_UNMAPPED_RE = re.compile(r"([\w_\-]+?)_?unmapped.out.mate([12]?)", re.IGNORECASE)  # TODO: flex regex
SAMPLE_NAME_MERGED_CLUSTERS_RE = re.compile(r"/([^/]*?)\.UE\.bed_files/")

ANALYSIS_SAMPLE_RE = re.compile(r"Analysis.*?of (/.*?)/([\w\-]*?)\n(.*?)(?=(?:Analysis|$))", re.DOTALL)

NAME_TYPE_TO_RE = {"fastq": FILENAME_FASTQ_RE, "unmapped": FILENAME_UNMAPPED_RE}
# endregion

# region Outputer Files Consts
# region Outputer Files Names and Formats
INPUT_OUTPUT_FILENAME = "%(date)s_HE_inputs_outputs_%(run_type)s_%(param_hash)s.tab"
DEFAULT_SITES_OUTPUT_FOLDER = "./HE_sites"
SITES_OUTPUT_FILE = "%s_HE_sites.bed"
SITES_REDITOOLS_OUTPUT_FILE = "%s_HE_sites.ForRedITools.bed"
SITES_TO_SAMPLE_OUTPUT_FILE = "%s_HE_sites.tab"
SPEC_FILE_NAME = "HESummerySpecification_%(data_type)s_%(stranded)s_%(run_type)s_%(params)s.csv"
SUMMERY_FILENAME = "HESummery_%s.csv"
LOG_FILE = "HE_run_log.log"
# endregion

# region Outputer Headers
SAMPLE_HEADER_FORMAT = "%(data_type)s-%(run_type)s-%(sample)s"
TOTAL_SITES_HEADER_FORAMT = "%s_TotalSites"
UNQ_SITES_HEADER_FORAMT = "%s_UniqueSites"
CLUSTERS_HEADER_FORAMT = "%s_ClustersCount"
MERGED_CLUSTERS_HEADER_FORAMT = "%s_MergedClustersCount"
MMF_NORM_TOTAL_SITES_HEADER_FORAMT = "%s_TotalSitesNormPerMillionMappedFrags"
MMF_NORM_UNQ_SITES_HEADER_FORAMT = "%s_UniqueSitesNormPerMillionMappedFrags"
MMF_NORM_CLUSTERS_HEADER_FORAMT = "%s_ClustersCountNormPerMillionMappedFrags"
MMF_NORM_MERGED_CLUSTERS_HEADER_FORAMT = "%s_MergedClustersCountNormPerMillionMappedFrags"
SAMPLE = "Sample"
PATH = "Path"
SAMPLES = "Samples"
PATHS = "Paths"
TOTAL_FRAGMENTS = "TotalFragments"
CHROMOSOME = "Chromosome"
START = "Start"
END = "End"
MISMATCH = "Mismatch"
STRAND = "Strand"
EXPRESSING_SAMPLES = "ExpressingSamples"
NEIGHBORS = "Neighbors"
NUM_OF_READS = "NumOfReads"
HE_UNMAPPED_FRAGMENTS = "HE_BWA_UnmappedFragments"
HE_MAPPED_FRAGMENTS = "HE_BWA_MappedFragments"
ALIGNER_UNMAPPED_FRAGMENTS = "STAR_UnmappedFragments"
ALIGNER_MAPPED_FRAGMENTS = "STAR_MappedFragments"
OUTPUT_HEADERS = (TOTAL_SITES_HEADER_FORAMT, UNQ_SITES_HEADER_FORAMT, CLUSTERS_HEADER_FORAMT,
                  MERGED_CLUSTERS_HEADER_FORAMT, MMF_NORM_TOTAL_SITES_HEADER_FORAMT, MMF_NORM_UNQ_SITES_HEADER_FORAMT,
                  MMF_NORM_CLUSTERS_HEADER_FORAMT, MMF_NORM_MERGED_CLUSTERS_HEADER_FORAMT)
# endregion

#: this is the delimiter between names of samples and paths in the output files.
SAMPLES_AND_PATHS_SEP = ";"
# endregion

# These consts are used to delete intermediate files of the Hyper Editing Tool. Change them if the output changes.
# Note: Sorting (instead of deleting) requires "samtools" to be in the $PATH
# region Bash Foramts
MERGE_BED_BASH_FORMAT = "bedtools intersect -wo -s -b %(raw_bed)s -a %(merged_bed)s |" \
                        " awk 'OFS=\"\t\"{print $1,$2,$3,$10,$11,$12}' > %(out_file)s"

SORT_BAMS_FORMAT = "%(samtools_path)s sort -l 9 -m 500M -@ 30 -o %(bam)s %(bam)s"

DEL_BAMS_FOLDERS_FORMAT = "rm -r %(output_dir)s/unMap %(output_dir)s/TransRun"

DEL_FASTQS_FORMAT = "find %(output_dir)s/unMap %(output_dir)s/TransRun/*/* -name '*.fastq'|xargs -P 5 -I %% rm %%"
# endregion

# Note: The following section is dependant on the Hyper Editing tool's input format - change accordingly.
# region Hyper Editing Run Params And Help Strings

# region Hyper Editing Tool Default Run Params
#: (generated automatically, don't change) a file containing for each sample a line of the format:
#:          <sample path><\t><output_name_prefix>
INPUT_OUTPUT_FILE = "generated automatically"

#: dirname + prefix of the bwa indexed genome so that there are 5 files found there:
#:        <dirname>/<prefix>.amb, <dirname>/<prefix>.ann, <dirname>/<prefix>.bwt, <dirname>/<prefix>.pac, <dirname>/<prefix>.sa
GENOME_BWA_INDEX_PREFIX = r"/private/common/Data/Genomes/Human/hg19.bwa_index/hg19"

#: same as GENOME_BWA_INDEX_PREFIX, but for the transformed genome created by the designated script of the Hyper Editing.
GENOME_TRANSFORMED_BWA_INDEX_PREFIX = r"/private/common/Data/Genomes/Human/hg19-transformed.bwa_index/hg19"

#: full path to the genome's fasta.
GENOME_FASTA = "/private/common/Data/Genomes/Human/hg19/hg19.fa"

#: base quality for PHRED score(33 or 64).
PHRED_SCORE_OFFSET = "33"

#: (generated automatically, don't change) run as single end (SE) or paired (PE).
RUN_PAIRED = "0"

#: max gap allowed (in bp) between mismatches on mates when running as PE.
MAX_GAP_SIZE = "500000"

#: (generated automatically, don't change) the output path.
OUTPUT_DIR = "generated automatically"

#: A flag. If not set the Hyper Editing tool will not run BWA mapping if analyseMM file exists.
RUN_BWA_MAPPING = "1"

#: A flag of values (0, 1, 2):
#:       0 => just detect UE without analysis.
#:       1=> analyse each src file (or PE) results.
#:       2=> analyse each src file (or PE) results and combine and analyse all the files together.
ANALYSIS_MODE = "2"

#: minimal percentage of hyper edited sites (i.e. mismatches of A->G) from a read's length to count it as hyper edited.
HE_SITES_MIN_PC = "0.05"

#: minimal ratio A->G to all other mismatches in a read to count it as hyper edited.
HE_TO_MM_RATIO = "0.6"

#: minimal quality of a mismatch site to count it as editing
SEQ_MIN_Q = "30"

#: maximal ratio of a single letter from the entire read (if passed the read is ignored).
MAX_SAME_LETTER_RATIO = "0.6"

#: minimal fraction (in percentage) of a cluster from the length of a read, to be considered as a hyper edited cluster.
MIN_PC_CLUSTER_LEN = "0.1"

#: The last position (in percentages of read length) of the *first* detected editing event in a read,
#  to be considered as a hyper edited cluster.
LAST_CLUSTER_START_INDEX = "0.8"

#: The first position (in percentages of read length) of the *last* detected editing event in a read,
#  to be considered as a hyper edited cluster.
FIRST_CLUSTER_END_INDEX = "0.2"

#: The full path of the main of the Hyper Editing tool (currently: 'run_hyper_editing.sh')
HE_SCRIPT_PATH = r"/home/alu/hagitpt/scripts/hyper_editing_scripts/run_hyper_editing.sh"
# endregion

# The following consts section (RUN_HE_FORMAT_DICT Key Names) contains the key name in RUN_HE_FORMAT_DICT for
# usage in the main script. *** Change them only according to RUN_HE_FORMAT_DICT! ***
# region RUN_HE_FORMAT_DICT Key Names
INPUT_OUTPUT_FILE_KEY = 'input_output_file'
GENOME_BWA_INDEX_PREFIX_KEY = 'genome_bwa_index_prefix'
GENOME_TRANSFORMED_BWA_INDEX_PREFIX_KEY = 'genome_transformed_bwa_index_prefix'
GENOME_FASTA_KEY = 'genome_fasta'
PHRED_SCORE_OFFSET_KEY = 'phred_score_offset'
RUN_PAIRED_KEY = 'run_paired'
MAX_GAP_SIZE_KEY = 'max_gap_size'
OUTPUT_DIR_KEY = 'output_dir'
RUN_BWA_MAPPING_KEY = 'run_bwa_mapping'
ANALYSIS_MODE_KEY = 'analysis_mode'
HE_SITES_MIN_PC_KEY = 'he_sites_min_pc'
HE_TO_MM_RATIO_KEY = 'he_to_mm_ratio'
SEQ_MIN_Q_KEY = 'seq_min_q'
MAX_SAME_LETTER_RATIO_KEY = 'max_same_letter_ratio'
MIN_PC_CLUSTER_LEN_KEY = 'min_pc_cluster_len'
LAST_CLUSTER_START_INDEX_KEY = 'last_cluster_start_index'
FIRST_CLUSTER_END_INDEX_KEY = 'first_cluster_end_index'
HE_SCRIPT_PATH_KEY = 'HE_script_path'

#: All the detection args (this args is only for the string format). It is also dependat on Hyper Editing tool's output
# (the format of the otuput dirs containing the UE BED files).
UE_DETECT_ARGS = (HE_SITES_MIN_PC_KEY, HE_TO_MM_RATIO_KEY, SEQ_MIN_Q_KEY, MAX_SAME_LETTER_RATIO_KEY,
                  MIN_PC_CLUSTER_LEN_KEY, LAST_CLUSTER_START_INDEX_KEY, FIRST_CLUSTER_END_INDEX_KEY)

# endregion

#: A dict containing the params for the run of the Hyper Editing tool, initialized with the default param values.
RUN_HE_FORMAT_DICT = {INPUT_OUTPUT_FILE_KEY: INPUT_OUTPUT_FILE,
                      GENOME_BWA_INDEX_PREFIX_KEY: GENOME_BWA_INDEX_PREFIX,
                      GENOME_TRANSFORMED_BWA_INDEX_PREFIX_KEY: GENOME_TRANSFORMED_BWA_INDEX_PREFIX,
                      GENOME_FASTA_KEY: GENOME_FASTA,
                      PHRED_SCORE_OFFSET_KEY: PHRED_SCORE_OFFSET,
                      RUN_PAIRED_KEY: RUN_PAIRED,
                      MAX_GAP_SIZE_KEY: MAX_GAP_SIZE,
                      OUTPUT_DIR_KEY: OUTPUT_DIR,
                      RUN_BWA_MAPPING_KEY: RUN_BWA_MAPPING,
                      ANALYSIS_MODE_KEY: ANALYSIS_MODE,
                      HE_SITES_MIN_PC_KEY: HE_SITES_MIN_PC,
                      HE_TO_MM_RATIO_KEY: HE_TO_MM_RATIO,
                      SEQ_MIN_Q_KEY: SEQ_MIN_Q,
                      MAX_SAME_LETTER_RATIO_KEY: MAX_SAME_LETTER_RATIO,
                      MIN_PC_CLUSTER_LEN_KEY: MIN_PC_CLUSTER_LEN,
                      LAST_CLUSTER_START_INDEX_KEY: LAST_CLUSTER_START_INDEX,
                      FIRST_CLUSTER_END_INDEX_KEY: FIRST_CLUSTER_END_INDEX,
                      HE_SCRIPT_PATH_KEY: HE_SCRIPT_PATH}

#: a dict containing the documentation of the RUN_HE_FORMAT_DICT variables for user help display. *** If changing make
#  sure keys are identical to RUN_HE_FORMAT_DICT ***
BASH_HELP_DICT = {'input_output_file': "(generated automatically, don't change) a file containing for each sample a "
                                       "line of the format: <sample path><\t>"
                                       "<output_name_prefix>",
                  'genome_bwa_index_prefix': "dirname + prefix of the bwa indexed genome so that there are 5 files "
                                             "found there: <dirname>/<prefix>.amb, <dirname>/<prefix>.ann, <dirname>"
                                             "/<prefix>.bwt, <dirname>/<prefix>.pac, <dirname>/<prefix>.sa",
                  'genome_transformed_bwa_index_prefix': "same as genome_bwa_index_prefix, but for the transformed"
                                                         " genome created by the designated script of the Hyper Editing"
                                                         "tool.",
                  'genome_fasta': "full path to the genome's fasta.",
                  'phred_score_offset': "base quality for PHRED score(33 or 64).",
                  'run_paired': "(generated automatically, don't change) run as single end (SE) or paired (PE).",
                  'max_gap_size': "max gap allowed (in bp) between mismatches on mates when running as PE.",
                  'output_dir': "(generated automatically, don't change) the output path.",
                  'run_bwa_mapping': "A flag. If not set the Hyper Editing tool will not run BWA mapping if"
                                     " analyseMM file exists.",
                  'analysis_mode': '''A flag of values (0, 1, 2): 0 => just detect UE without analysis;
                   1=> analyse each src file (or PE) results; 2=> analyse each src file (or PE) results and combine and
                    analyse all the files together.''',
                  'he_sites_min_pc': "minimal percentage of hyper edited sites (i.e. mismatches of A->G) from a read's "
                                     "length to count it as hyper edited.",
                  'he_to_mm_ratio': "minimal ratio A->G to all other mismatches in a read to count it as hyper edited.",
                  'seq_min_q': "minimal quality of a mismatch site to count it as editing",
                  'max_same_letter_ratio': "maximal ratio of a single letter from the entire read (if passed the read "
                                           "is ignored).",
                  'min_pc_cluster_len': "minimal fraction (in percentage) of a cluster from the length of a read, to be"
                                        " considered as a hyper edited cluster.",
                  'last_cluster_start_index': "The last position (in percentages of read length) of the FIRSTdetected"
                                              " editing event in a read, to be considered as a hyper edited cluster.",
                  'first_cluster_end_index': "The first position (in percentages of read length) of the LAST detected"
                                             " editing event in a read, to be considered as a hyper edited cluster.",
                  'HE_script_path': "The full path of the main of the Hyper Editing tool (currently: "
                                    "'run_hyper_editing.sh')"
                  }

BASH_HELP_STR = "\t\r\n".join(["%(att)s: %(h_str)s\r\n\t\tdefault: %(d_val)s" %
                               dict(att=key, h_str=BASH_HELP_DICT[key], d_val=value)
                               for key, value in RUN_HE_FORMAT_DICT.iteritems()])

#: this is an ordered list to create the bash and help string automatically with the key names
RUN_HE_BASH_ARGS = (HE_SCRIPT_PATH_KEY, GENOME_BWA_INDEX_PREFIX_KEY, GENOME_TRANSFORMED_BWA_INDEX_PREFIX_KEY,
                    GENOME_FASTA_KEY, PHRED_SCORE_OFFSET_KEY, RUN_PAIRED_KEY, MAX_GAP_SIZE_KEY, OUTPUT_DIR_KEY,
                    RUN_BWA_MAPPING_KEY, ANALYSIS_MODE_KEY, HE_SITES_MIN_PC_KEY, HE_TO_MM_RATIO_KEY, SEQ_MIN_Q_KEY,
                    MAX_SAME_LETTER_RATIO_KEY, MIN_PC_CLUSTER_LEN_KEY, LAST_CLUSTER_START_INDEX_KEY,
                    FIRST_CLUSTER_END_INDEX_KEY, INPUT_OUTPUT_FILE_KEY)
# the foramt for the bash running the Hyper Editing tool
RUN_HE_BASH_FORMAT = "%(" + ")s %(".join(RUN_HE_BASH_ARGS) + ")s"

# values for the flags used in the script
HE_TRUE_VAL = "1"
HE_FALSE_VAL = "0"

# endregion

# These consts are used to extract data from the Hyper Editing Tool. Change them if the output's format changes.
# region Hyper Editing Tool-Related Consts
# ---- stats summery consts ---
STATS_FILE_DATA_START = 4
STATS_FILE_DATA_SEP = '\t'
AG_HE_PC_COL = "%UE_of_tot"
PRE_G_PC_COL = "G-up"
PAST_G_PC_COL = "G-down"
SAMPLE_PATH_COL = "src_path"
SAMPLE_NAME_COL = "sample"
MAPPED_READS_COL = "Map_reads"
UNMAPPED_READS_COL = "unMap_reads"
SENSE_STRAND = "+"
HE_EMPTY_VAL = "-"
MM_EDIT_TYPE = "Edit_type"
# ---- BED files consts ---
SITES_BED_CHR_I = 0
SITES_BED_START_I = 1
SITES_BED_END_I = 2
SITES_BED_STRAND_I = -1
OUTPUT_STAR_ALL_READS = 0
OUTPUT_STAR_UNMAPPED_READS = 1
OUTPUT_STAR_MAPPED_READS = 2
# ---- BED files formats consts ---
HE_OUTFILES_PARAMS_FORMAT = "%(paired)s_%(" + ")s_%(".join(UE_DETECT_ARGS) + ")s"

BED_CLUSTERS_MERGED_FILE_FORMAT = r"%(output_dir)s/UEdetect.%(params)s/%(sample)s.UE.bed_files/UEmergeSd20.%(mismatch)s.bed"
BED_CLUSTERS_FILE_FORMAT = r"%(output_dir)s/UEdetect.%(params)s/%(sample)s.UE.bed_files/%(mismatch)s.bed"
UNQ_SITES_BED_FILE_FORMAT = r"%(output_dir)s/UEdetect.%(params)s/%(sample)s.ES.bed_files/ESuniqS.%(mismatch)s.bed"
SITES_BED_FILE_FORMAT = r"%(output_dir)s/UEdetect.%(params)s/%(sample)s.ES.bed_files/%(mismatch)s.bed"
HE_STATS_FILE_FORMAT = r"%(output_dir)s/statistic/analyse.%(params)s.A2G"
ALL_SAMPLES = "all"
TRANS_FOLDER_FORMAT = r"%(output_dir)s/TransRun"
UNMAP_FOLDER_FORMAT = r"%(output_dir)s/unMap"
# endregion

# region Logging Messages
RUN_PARAMS_DEBUG = "Starting Run with params: input_dir = %(input_dir)s, output_dir = %(output_dir)s, files_suffix = " \
                   "%(files_suffix)s, output_files_dir = %(output_files_dir)s, sites_folder = %(sites_folder)s, " \
                   "run_single_end = %(run_single_end)s, run_paired_end = %(run_paired_end)s, repeat_on_noise = %(" \
                   "repeat_on_noise)s, STAR_output = %(STAR_output)s, STAR_recursion_depth = %(" \
                   "STAR_recursion_depth)s, stranded = %(stranded)s, calc_output_only = %(calc_output_only)s, " \
                   "keep_bams = %(keep_bams)s, max_pre_g_pc = %(max_pre_g_pc)s, min_past_g_pc = %(min_past_g_pc)s, " \
                   "strict_check = %(strict_check)s, folder_names_only = %(folder_names_only)s, run_bash_args = %(" \
                   "run_bash_args)s"
RUN_HE_MSG = "Running Hyper Editing Tool With The Following: \r\n\tBash Line: '%s'\r\n\t(" + " ".join(RUN_HE_BASH_ARGS) \
             + ")\r\n\t Not As Defaults Params: %s "
RUN_RECALC_HE_MSG = "Recalculating Summery Outputs For Hyper Editing Run With The Following Bash Line: '%s'"
NOISY_RERUN_MSG = "%s Of The Samples Were Found Too Noisy (Lacking Motif %s Low A->G Mismatches Ratio) Rerunning With " \
                  "Stricter Params"
CONCAT_STATS_WRN = "Encountered A Non Numerical Line In the Hyper Editing Stats File (%s), Using Later Outputer! " \
                   " (This usually means there's a concatenation from a previous run)"
SAMPLE_FAILED_NOISE_CHECK_WRN = "Sample %(sample)s Was Found Noisy! (Checking Motif %(strict_operator)s Mismatch Ratio)" \
                                "\r\n\tMotif Thresholds:\r\n\t\tG-Up Threshold: %(g_up_threshold)s," \
                                " Found %(g_up_found)s (%(g_up_res)s)\r\n\t\tG-Down Threshold: %(g_dwn_threshold)s," \
                                " Found %(g_dwn_found)s (%(g_dwn_res)s)" \
                                "\r\n\tMismatch Ratio Threshold: %(mm_r)s, Found %(mm_r_found)s (%(mm_r_res)s)"
COND_SUCCESS = "Passed"
COND_FAILED = "Failed"
BAD_PATH_WARN = "Path %s Is Problematic, Can't Add It To Analysis!"
HE_START = "Started Running HE."
SUMMERY_CREATION_MSG = "Creating Summery File: %s"
BED_AND_SPEC_CREATION_MSG = "Creating BED and Sites Specification Files."
MISSING_HE_SE_MATE_OUTFILE_DBG = "No 2nd Mate Outputer File (%s - %s) For Sample %s, Will Skip It!"
MISSING_HE_OUTFILE_WRN = "Missing Outputer File (%s - %s) For Sample %s!"
MISSING_HE_OUT_SKIPPING_WRN = "Missing Outputer For Sample %s!"
CLEANUP_MSG = "Removing Temp Fastqs and %s Intermediate BAM Files"
CLEANUP_DEL_FASTQS_MSG = "Removing Temp Fastqs: %s"
CLEANUP_BAMS_MSG = "%s Intermediate BAMs: %s"
# endregion

STRANDED = "Stranded"
NOT_STRANDED = "NotStranded"
DATA_TYPE_CLUSTERS = "Clusters"
DATA_TYPE_MERGED_CLUSTERS = "MergedClusters"
DATA_TYPE_SITES = "Sites"
DATA_TYPE_UNQ_SITES = "UniqueSites"

MATE1_SUFF = "-1"
MATE2_SUFF = "-2"

RUN_TYPE_SE = "SE"
RUN_TYPE_PE = "PE"

#: this is the mismatch summed in the output files.
# TODO: consider adding all others
MISMATCH_FOR_SUMMERY = "A2G"
