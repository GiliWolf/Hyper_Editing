__author__ = 'Hillel'
"""
This module contains the consts for the wrapper of Hagit Porath's HE tool.
"""

from os.path import dirname, join

# TODO: add path locator
SCRIPT_INSTALL_DIR = dirname(__file__)

# region Hyper Editing Scripts
# "/home/alu/hagitpt/scripts/hyper_editing_scripts/stranded_sample.sh"
HE_SETUP_SCRIPT = "/home/alu/fulther/Scripts/bash/HyperEditing/setup_file_management.sh"
ANALYSE_UE_SCRIPT = "/home/alu/hagitpt/scripts/hyper_editing_scripts/Analyse_UE_basic_4pl.sh"
STAT_SUMMARY_SCRIPT = "/home/alu/hagitpt/scripts/hyper_editing_scripts/stat_summary_basic.pl"
COMBINE_ALL_SCRIPT = "/home/alu/hagitpt/scripts/hyper_editing_scripts/combine_all.sh"
UNMAPPED_SCRIPT = "/home/alu/hagitpt/scripts/hyper_editing_scripts/pre_unmapped_4pl.sh"
TRANS_RUN_SCRIPT = "/home/alu/hagitpt/scripts/hyper_editing_scripts/TransRun_4pl.sh"
ANALYSE_MM_SCRIPT = "/home/alu/hagitpt/scripts/hyper_editing_scripts/analyse_mm_4pl.sh"
DETECT_UE_SCRIPT = "/home/alu/hagitpt/scripts/hyper_editing_scripts/detect_ue_4pl.pl"
COMBINE_PE_SCRIPT = "/home/alu/hagitpt/scripts/hyper_editing_scripts/combine_PE.sh"
FULL_ANALYSIS_ALL_SCRIPT = "/home/alu/fulther/Scripts/bash/HyperEditing/combined_analysis_all_files.sh"
# endregion

# region Hyper Editing Configs
CONFIGS_DIR = "%s/Configs" % SCRIPT_INSTALL_DIR
BWA_ALIGNMENTS_CONF = join(CONFIGS_DIR, "runHE_bwaOnly.conf")
RUN_DETECTION_CONF = join(CONFIGS_DIR, "runHE_analysisOnly.conf")
RUN_ANALYZE_MM2BAM_CONF = join(CONFIGS_DIR, "runHE_analyzeMMtoBAM.conf")
RUN_ANALYZE_MM2BAM_MERGE_CONF = join(CONFIGS_DIR, "runHE_analyzeMMtoBAM_merge.conf")
# endregion

# region Hyper Editing Params And Help Strings
NO_2_MATE_INPUT_STR = "No2Mate"

ANALYSIS_MODE_FULL = "FullOutput"
ANALYSIS_MODE_NO_SUM = "NoTotal"
ANALYSIS_MODE_JUST_DETECT = "DetectOnly"

#: A flag of values (0, 1, 2):
#:       0 => just detect UE without analysis.
#:       1=> analyse each src file (or PE) results.
#:       2=> analyse each src file (or PE) results and combine and analyse all the files together.
ANALYSIS_MODE_DICT = {
    ANALYSIS_MODE_FULL: 2,
    ANALYSIS_MODE_NO_SUM: 1,
    ANALYSIS_MODE_JUST_DETECT: 0
}

# values for the flags used in the script
HE_TRUE_VAL = "1"
HE_FALSE_VAL = "0"

HE_BOOLEAN_TRANS = {
    True: HE_TRUE_VAL,
    False: HE_FALSE_VAL
}

HE_SE_VAL = "SE"
HE_PE_VAL = "PE"

HE_SITES_VAL = "ES"
HE_CLUSTERS_VAL = "UE"

ANALYSIS_MODE_HELP_STR = "%s - just detect UE without analysis, %s - output analysis for each sample, and %s - " \
                         "output combined summary too." % (
                             ANALYSIS_MODE_JUST_DETECT, ANALYSIS_MODE_NO_SUM, ANALYSIS_MODE_FULL)

HE_DETECT_ARGS_FORMAT = "_".join([r"%(he_sites_min_pc)s", r"%(he_to_mm_ratio)s", r"%(seq_min_q)s",
                                  r"%(max_same_letter_ratio)s", r"%(min_pc_cluster_len)s",
                                  r"%(last_cluster_start_index)s", r"%(first_cluster_end_index)s"])

HE_SETUP_ARGS_FORMAT = " ".join([r"%(he_output_dir)s", r"%(run_paired_end)s", r"%(force_run_bwa_mapping)s",
                                 r"%(analysis_mode)s", STAT_SUMMARY_SCRIPT, r"%(he_detect_args)s"])

UNCMPRS_TMP_DIR_FORMAT = r"%s/uncompress_tmp_%s"
# endregion

# region Hyper Editing Folders Formats
HE_DETECT_DIR_FORMAT = r"%(he_output_dir)s/UEdetect.%(detect_args)s"
HE_UNMAP_DIR_FORMAT = r"%(he_output_dir)s/unMap"
HE_TRANS_RUN_DIR_FORMAT = r"%(he_output_dir)s/TransRun"
HE_ANALYZE_RUN_DIR_FORMAT = r"%(he_output_dir)s/AnalyseMM"
HE_STATISTICS_DIR_FORMAT = r"%(he_output_dir)s/statistic/analyse.%(detect_args)s"
HE_STATISTICS_GENERAL_FILE_FORMAT = r"%(he_output_dir)s/statistic/general"

HE_ANALYZE_MM_FILE_SUFF = "analyseMM"
# this is what's nested inside HE_DETECT_DIR_FORMAT
HE_DETECT_DIR_PER_SAMPLE_BEDS_FORMAT = "%(sample_name)s.%(data_type)s.bed_files"
# endregion

# endregion


# region Outputer Files Consts
# region Outputer Files Names and Formats
INPUT_OUTPUT_FILENAME = r"%(date)s_HE_inputs_outputs_%(run_type)s_%(param_hash)s.tab"
DEFAULT_SITES_OUTPUT_FOLDER = r"./HE_sites"
SITES_OUTPUT_FILE = "%s_HE_sites.bed"
SITES_REDITOOLS_OUTPUT_FILE = "%s_HE_sites.ForRedITools.bed"
SITES_TO_SAMPLE_OUTPUT_FILE = "%s_HE_sites.tab"
SPEC_FILE_NAME = "HESummerySpecification_%(data_type)s_%(stranded)s_%(run_type)s_%(params)s.csv"
SUMMERY_FILENAME = "HESummery_%s.csv"
LOG_FILE = "RunHyperEditing.log"
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

# region Logging
RUN_PARAMS_DEBUG_MSG = "Starting Run with Params: %s"
RUN_SETUP_MSG = "Running Setup For %s With Params: %s"
CREATE_UNCMPRS_DIR_MSG = "Creating Temp Uncompress Directory at %s"
UNCMPRS_DIR_EXISTS_WARN_MSG = "Temp Uncompress Directory at %s Already Exists, Make Sure You don't Have Conflicting Runs!"
RUN_ALIGN_PIPELINE_MSG = "\n\nRunning Alignments Using Logging Directory %s"
RUN_HE_DETECT_PIPELINE_MSG = "\n\nRunning Hyper Editing Detection Using Logging Directory %s"
RUN_HE_FULL_ANALYSIS_MSG = "\n\nRunning Hyper Editing Full Analysis Using Command %s"
# endregion


# region System Commands
HE2BAM_CMD = "%(python3)s %(HE2BAM_prog)s %(HE_dir)s %(params)s %(PE)s %(file_name)s %(sam_output_prefix)s %(header_sample_bam)s %(samtools)s"
# endregion
