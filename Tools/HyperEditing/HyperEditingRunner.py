__author__ = 'Hillel & Roni'
"""
This module wraps the run of Hagit's HE tool.
"""
# region Imports
# region Builtin Imports
from collections import namedtuple, OrderedDict
from csv import reader
from datetime import datetime
from operator import and_, or_
from pprint import pformat

import argparse
import glob
import inspect
import logging
import multiprocessing
import os
import subprocess

# endregion

# region Imports From Internal Modules
if __name__ == "__main__":
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from Outputs.Outputers.CSVOutputer import CSVOutputer
from Outputs.Outputers.RawOutputer import RawOutputer

from DataConverters.HyperEditingBEDConverter import HyperEditingBEDConverter, HyperEditingBEDRec, DICT_KEY_SEP

from Session.BlackBoxesWrapperSctipts.STAR242_functios import get_STAR_total_reads, get_STAR_unmapped, \
    get_STAR_uniquely_aligned_reads

from Commons.data_structs import SortingHelpFormatter
from Commons.general_functions import init_logging_dict, convert_args_to_dict, remove_files, \
    add_get_paths_function_to_argparser
from Commons.help_functions import reverse_strand
from Commons.consts import *

from Session import PipelineManger

from Tools.HyperEditing.HyperEditingConsts import *
import path_locator


# endregion

# endregion

def run_hyper_editing_setup(analysis_mode, force_run_bwa_mapping, he_detect_args, he_output_dir, run_mode_str,
                            dont_run):
    # type: (int, bool, str, str, str, bool) -> tuple
    he_detect_args = "_".join([run_mode_str, he_detect_args])

    # region Setup
    he_setup_args = HE_SETUP_ARGS_FORMAT % dict(
        he_output_dir=he_output_dir,
        run_paired_end=HE_FALSE_VAL,
        force_run_bwa_mapping=HE_BOOLEAN_TRANS[force_run_bwa_mapping],
        analysis_mode=analysis_mode,
        he_detect_args=he_detect_args
    )
    uncompress_temp_dir = UNCMPRS_TMP_DIR_FORMAT % (he_output_dir, he_detect_args)

    if dont_run:
        return he_detect_args, uncompress_temp_dir

    logging.info(RUN_SETUP_MSG % (run_mode_str, he_setup_args))
    os.system(HE_SETUP_SCRIPT + " " + he_setup_args)

    if not os.path.exists(uncompress_temp_dir):
        logging.info(CREATE_UNCMPRS_DIR_MSG % uncompress_temp_dir)
        os.makedirs(uncompress_temp_dir)
    else:
        logging.warn(UNCMPRS_DIR_EXISTS_WARN_MSG % uncompress_temp_dir)
    # endregion

    return he_detect_args, uncompress_temp_dir


def run_hyper_editing_bwa_alignments(bwa_genome_index_pre, exclude_paths, exclude_paths_operator, files_suffix,
                                     follow_links, force_run_bwa_mapping, genome_fasta, he_output_dir, include_paths,
                                     include_paths_operator, is_pe, logs_dir, mate1_suff, mate2_suff, mates_prefix,
                                     phred_score_offset, recursion_depth, root_dir, transformed_bwa_genome_index_pre,
                                     uncompress_cmd, uncompress_temp_dir, decrease_parallel_run):
    run_bwa_logs_dir = os.path.join(logs_dir, "Run_BWA_Alignments")
    logging.info(RUN_ALIGN_PIPELINE_MSG % run_bwa_logs_dir)

    truncs = get_pipeline_truncs(mate1_suff, mates_prefix)
    conf_override_args = dict(genome_bwa_ind=bwa_genome_index_pre,
                              genome_trans_bwa_ind=transformed_bwa_genome_index_pre,
                              genome_fasta=genome_fasta,
                              dir_pre=he_output_dir,
                              bwa_run=HE_BOOLEAN_TRANS[force_run_bwa_mapping],
                              Q=str(phred_score_offset),
                              PE=HE_BOOLEAN_TRANS[is_pe],
                              uncompress_fastq_cmd=uncompress_cmd,
                              suffix=files_suffix,
                              mateprefix=mates_prefix,
                              mate1=mate1_suff,
                              mate2=mate2_suff,
                              unmapped_script=UNMAPPED_SCRIPT,
                              TransRun_script=TRANS_RUN_SCRIPT,
                              analyse_mm_script=ANALYSE_MM_SCRIPT,
                              limited_alignment_ps=str(decrease_parallel_run))
    try:
        PipelineManger.main(root_dir=root_dir,
                            output_dir=uncompress_temp_dir,
                            include_paths=include_paths,
                            include_paths_operator=include_paths_operator,
                            exclude_paths=exclude_paths,
                            exclude_paths_operator=exclude_paths_operator,
                            recursion_depth=recursion_depth,
                            follow_links=follow_links,
                            config_file=BWA_ALIGNMENTS_CONF,
                            log_path=run_bwa_logs_dir,
                            files_suffix=files_suffix,
                            truncs=truncs,
                            trunc_after_join=True,
                            config_args=conf_override_args)
        logging.info("Done Running Alignments! (Check for failed samples)")
    except Exception, e:
        logging.exception("Running Alignments Pipeline Has Failed! The Program Will Now Exit!")
        sys.exit(11)


def get_pipeline_truncs(mate1_suff, mates_prefix):
    num_of_periods = (mate1_suff + mates_prefix).count(".")
    return [num_of_periods, len(mates_prefix + mate1_suff)]


def run_hyper_editing_detection(analysis_mode, exclude_paths, exclude_paths_operator, files_suffix, follow_links,
                                he_output_dir, include_paths, include_paths_operator, is_pe, logs_dir, mate1_suff,
                                mate2_suff, mates_prefix, max_gap_size, recursion_depth, root_dir, detect_args,
                                uncompress_temp_dir, uncompress_cmd):
    run_detection_logs_dir = os.path.join(logs_dir, "Run_HE_Detection_" + detect_args)
    logging.info(RUN_HE_DETECT_PIPELINE_MSG % run_detection_logs_dir)
    truncs = get_pipeline_truncs(mate1_suff, mates_prefix)
    strip_val = HE_PE_VAL if is_pe else HE_SE_VAL
    conf_override_args = dict(dir_pre=he_output_dir,
                              PE=HE_BOOLEAN_TRANS[is_pe],
                              uncompress_fastq_cmd=uncompress_cmd,
                              suffix=files_suffix,
                              mateprefix=mates_prefix,
                              mate1=mate1_suff,
                              mate2=mate2_suff,
                              detectUE_script=DETECT_UE_SCRIPT,
                              AnalyseUE_script=ANALYSE_UE_SCRIPT,
                              stat_summary_script=STAT_SUMMARY_SCRIPT,
                              combine_PE_script=COMBINE_PE_SCRIPT,
                              GAP=str(max_gap_size),
                              analyse=str(analysis_mode),
                              ue_detect_args=detect_args.lstrip(strip_val + "_").replace("_", " "),
                              args=detect_args,
                              )
    try:
        PipelineManger.main(root_dir=root_dir,
                            output_dir=uncompress_temp_dir,
                            include_paths=include_paths,
                            include_paths_operator=include_paths_operator,
                            exclude_paths=exclude_paths,
                            exclude_paths_operator=exclude_paths_operator,
                            recursion_depth=recursion_depth,
                            follow_links=follow_links,
                            config_file=RUN_DETECTION_CONF,
                            log_path=run_detection_logs_dir,
                            files_suffix=files_suffix,
                            truncs=truncs,
                            trunc_after_join=True,
                            config_args=conf_override_args)
        logging.info("Done Running Alignments! (Check for failed samples)")
    except Exception, e:
        logging.exception("Running Alignments Pipeline Has Failed! The Program Will Now Exit!")
        sys.exit(11)

    return run_detection_logs_dir


def run_hyper_editing_full_analysis(he_output_dir, is_pe, detect_args, run_detection_logs_dir):
    analysis_cmd = " ".join([FULL_ANALYSIS_ALL_SCRIPT,
                             HE_BOOLEAN_TRANS[is_pe],
                             COMBINE_ALL_SCRIPT,
                             ANALYSE_UE_SCRIPT,
                             STAT_SUMMARY_SCRIPT,
                             os.path.join(run_detection_logs_dir, "flags"),
                             HE_DETECT_DIR_FORMAT % dict(he_output_dir=he_output_dir,
                                                         detect_args=detect_args),
                             HE_STATISTICS_DIR_FORMAT % dict(
                                 he_output_dir=he_output_dir,
                                 detect_args=detect_args),
                             HE_STATISTICS_GENERAL_FILE_FORMAT % dict(he_output_dir=
                                                                      he_output_dir)])
    logging.info(RUN_HE_FULL_ANALYSIS_MSG % analysis_cmd)
    os.system(analysis_cmd)


def create_he_bam(star_dir, bam_files_suffix, he_output_dir, he_detect_args, is_pe, bam_out_dir, follow_links,
                  python3_path, logs_dir, samtools_path, decrease_parallel, include_paths, include_paths_operator,
                  exclude_paths, exclude_paths_operator):
    script_dir = path_locator.module_path()
    unmapped_dir = HE_UNMAP_DIR_FORMAT % dict(he_output_dir=he_output_dir)
    header_sample_bam = None

    # region Choose a BAM for the headers
    if os.path.exists(unmapped_dir):
        logging.info("Trying to use BAM from unmapped folder (Unmapped dir %s is used)" % unmapped_dir)
        for root, dirs, files in os.walk(unmapped_dir):
            if header_sample_bam:
                break
            for f in files:
                if f.endswith("mem.bam"):
                    header_sample_bam = os.path.join(root, f)
                    break
        if header_sample_bam is None:
            logging.warning("Your unmapped dir doesn't contain any BAMs!")

    elif star_dir:
        logging.info("Trying to use BAM from STAR folder (STAR dir %s is used)" % star_dir)
        for root, dirs, files in os.walk(star_dir):
            if header_sample_bam:
                break
            for f in files:
                if f.endswith(bam_files_suffix):
                    header_sample_bam = os.path.join(root, f)
                    break
        if header_sample_bam is None:
            logging.warning("Your STAR dir doesn't contain any BAMs!")

    if header_sample_bam is None:
        logging.info("Using header from empty hg38 BAM file!")
        header_sample_bam = os.path.join(script_dir, "GeneralHumanHg38Header.bam")
    # endregion

    conf_override_args = dict(dir_pre=he_output_dir,
                              PE="PE" if is_pe else "SE",
                              he2bam_prog=os.path.join(script_dir, "HE2BAM.py"),
                              params=he_detect_args,
                              HE_dir=he_output_dir,
                              header_sample_bam=header_sample_bam,
                              python3=python3_path,
                              samtools=samtools_path,
                              limited_alignment_ps=str(decrease_parallel))
    try:
        PipelineManger.main(root_dir=HE_ANALYZE_RUN_DIR_FORMAT % dict(he_output_dir=he_output_dir),
                            output_dir=bam_out_dir,
                            include_paths=include_paths,
                            include_paths_operator=include_paths_operator,
                            exclude_paths=exclude_paths,
                            exclude_paths_operator=exclude_paths_operator,
                            follow_links=follow_links,
                            config_file=RUN_ANALYZE_MM2BAM_CONF,
                            log_path=os.path.join(logs_dir, "HE2BAM_Conversion"),
                            files_suffix=HE_ANALYZE_MM_FILE_SUFF,
                            config_args=conf_override_args,
                            truncs=[0, 2] if is_pe else [0, 0],
                            no_extension=True)
        logging.info("Done Running HE files to BAM conversion! (Check for failed samples)")
    except:
        logging.exception("Running HE files to BAM conversion Pipeline Has Failed! The Program Will Now Exit!")
        sys.exit(11)


def merge_he_bam(star_dir, bam_files_suffix, he_detect_args, bam_out_dir, follow_links,
                 samtools_path, logs_dir, decrease_parallel, remove_recurring_reads, include_paths,
                 include_paths_operator, exclude_paths, exclude_paths_operator):
    script_dir = path_locator.module_path()

    conf_override_args = dict(extract_unmapped_reads_prog=os.path.join(script_dir, "extract_read_from_bam.py"),
                              remove_recurring_reads=str(remove_recurring_reads),
                              samtools=samtools_path,
                              aligned_bam_suffix=bam_files_suffix,
                              HE_BAM_dir=bam_out_dir,
                              extra_formatting=he_detect_args,
                              limited_alignment_ps=str(decrease_parallel))
    try:

        PipelineManger.main(root_dir=star_dir,
                            output_dir=star_dir,
                            include_paths=include_paths,
                            include_paths_operator=include_paths_operator,
                            exclude_paths=exclude_paths,
                            exclude_paths_operator=exclude_paths_operator,
                            follow_links=follow_links,
                            config_file=RUN_ANALYZE_MM2BAM_MERGE_CONF,
                            log_path=os.path.join(logs_dir, "HE2BAM_Merge_" + he_detect_args),
                            files_suffix=bam_files_suffix,
                            truncs=[0, len(bam_files_suffix.split(".")[0])],
                            config_args=conf_override_args,
                            no_extension=True)
        logging.info("Done Running HE files to BAM merge! (Check for failed samples)")
    except:
        logging.exception("Running HE files to BAM merge Pipeline Has Failed! The Program Will Now Exit!")
        sys.exit(11)


def main(  # Inputs options
        root_dir, include_paths, include_paths_operator, exclude_paths, exclude_paths_operator, recursion_depth,
        follow_links, files_suffix, uncompress_cmd,
        # mates options
        mate1_suff, mate2_suff, mates_prefix,
        # Aligner (STAR) params
        star_dir, star_recursion_depth, bam_files_suffix,

        # Output Options
        he_output_dir, sum_out_dir, sites_folder, bam_out_dir,

        # Hyper-Editing Run Parameters:
        # Genome Options
        bwa_genome_index_pre, transformed_bwa_genome_index_pre, genome_fasta,
        # Run Mode
        run_single_end, run_paired_end, stranded,
        # Hyper-Editing Run Parameters
        phred_score_offset, max_gap_size, force_run_bwa_mapping, analysis_mode, he_sites_min_pc, he_to_mm_ratio,
        seq_min_q, max_same_letter_ratio, min_pc_cluster_len, last_cluster_start_index, first_cluster_end_index,

        # Rerun and Clean Setting
        repeat_on_noise, max_noisy_samples_ratio, a2g_cleanness_threshold, max_upstream_g_pc,
        min_downstream_g_pc, strict_check,

        # Miscellaneous
        calc_output_only, folder_names_only, keep_unmapped_bams, keep_trans_run, decrease_parallel_run, python3_path,
        samtools_path, dont_run, remove_recurring_reads):
    # type: (str, list, operator, list, operator, int, bool, str, str, str, str, str, str, int, str, str, str, str, str, str, str, str, bool, bool, bool, int, int, bool, int, float, float, int, str, float, float, float, bool, float, float, float, float, bool, bool, bool, bool, bool, int, str, str, bool) -> None
    assert len(mate1_suff) == len(
        mate2_suff) or mate2_suff == NO_2_MATE_INPUT_STR, "Mates Suffixes Must Be of the Same Length (For using " \
                                                          "SE data only use %s for mate 2)!" % NO_2_MATE_INPUT_STR

    logs_dir = os.path.join(he_output_dir, "Logs")
    init_logging_dict(os.path.join(logs_dir, datetime.today().isoformat() + "_" + LOG_FILE))
    args, _, _, values = inspect.getargvalues(inspect.currentframe())


    logging.debug(RUN_PARAMS_DEBUG_MSG % pformat(dict([(i, values[i]) for i in args])))

    he_detect_args = HE_DETECT_ARGS_FORMAT % dict(
        he_sites_min_pc=he_sites_min_pc,
        he_to_mm_ratio=he_to_mm_ratio,
        seq_min_q=seq_min_q,
        max_same_letter_ratio=max_same_letter_ratio,
        min_pc_cluster_len=min_pc_cluster_len,
        last_cluster_start_index=last_cluster_start_index,
        first_cluster_end_index=first_cluster_end_index)

    if bam_out_dir is None:
        bam_out_dir = os.path.join(he_output_dir, "HyperEditingBAM")

    # call setup (create dirs)
    ran_se_alignments = False
    if run_single_end:
        se_detect_args, se_uncompress_temp_dir = run_hyper_editing_setup(analysis_mode, force_run_bwa_mapping,
                                                                         he_detect_args, he_output_dir, HE_SE_VAL,
                                                                         dont_run)
        is_pe = False
        if not dont_run:
            run_hyper_editing_bwa_alignments(bwa_genome_index_pre, exclude_paths, exclude_paths_operator, files_suffix,
                                             follow_links, force_run_bwa_mapping, genome_fasta, he_output_dir,
                                             include_paths, include_paths_operator, is_pe, logs_dir, mate1_suff,
                                             mate2_suff,
                                             mates_prefix, phred_score_offset, recursion_depth, root_dir,
                                             transformed_bwa_genome_index_pre, uncompress_cmd, se_uncompress_temp_dir,
                                             decrease_parallel_run)
            ran_se_alignments = True

            se_run_detection_logs_dir = run_hyper_editing_detection(analysis_mode, exclude_paths,
                                                                    exclude_paths_operator,
                                                                    files_suffix, follow_links, he_output_dir,
                                                                    include_paths,
                                                                    include_paths_operator, is_pe, logs_dir, mate1_suff,
                                                                    mate2_suff, mates_prefix, max_gap_size,
                                                                    recursion_depth,
                                                                    root_dir, se_detect_args, se_uncompress_temp_dir,
                                                                    uncompress_cmd)
            remove_files([se_uncompress_temp_dir])

            if analysis_mode == ANALYSIS_MODE_DICT[ANALYSIS_MODE_FULL]:
                run_hyper_editing_full_analysis(he_output_dir, is_pe, se_detect_args, se_run_detection_logs_dir)

        if bam_out_dir:
            bam_out_dir_full = os.path.join(bam_out_dir, "SE_" + he_detect_args)
            create_he_bam(star_dir, bam_files_suffix, he_output_dir, he_detect_args, is_pe, bam_out_dir_full,
                          follow_links, python3_path, logs_dir, samtools_path, decrease_parallel_run,
                          include_paths, include_paths_operator, exclude_paths, exclude_paths_operator)

            if star_dir:
                merge_he_bam(star_dir, bam_files_suffix, "SE_" + he_detect_args, bam_out_dir_full, follow_links, samtools_path,
                             logs_dir, decrease_parallel_run, remove_recurring_reads, include_paths,
                             include_paths_operator, exclude_paths, exclude_paths_operator)

    if run_paired_end:
        pe_detect_args, pe_uncompress_temp_dir = run_hyper_editing_setup(analysis_mode, force_run_bwa_mapping,
                                                                         he_detect_args, he_output_dir, HE_PE_VAL,
                                                                         dont_run)
        is_pe = True
        if not dont_run:
            if ran_se_alignments:
                logging.info("Alignments were run for single-end, no re-run is needed for paired-end detection.")
                pe_force_run_bwa_mapping = False
            else:
                pe_force_run_bwa_mapping = force_run_bwa_mapping

            run_hyper_editing_bwa_alignments(bwa_genome_index_pre, exclude_paths, exclude_paths_operator, files_suffix,
                                             follow_links, pe_force_run_bwa_mapping, genome_fasta, he_output_dir,
                                             include_paths, include_paths_operator, is_pe, logs_dir, mate1_suff,
                                             mate2_suff,
                                             mates_prefix, phred_score_offset, recursion_depth, root_dir,
                                             transformed_bwa_genome_index_pre, uncompress_cmd, pe_uncompress_temp_dir,
                                             decrease_parallel_run)

            pe_run_detection_logs_dir = run_hyper_editing_detection(analysis_mode, exclude_paths,
                                                                    exclude_paths_operator,
                                                                    files_suffix, follow_links, he_output_dir,
                                                                    include_paths,
                                                                    include_paths_operator, is_pe, logs_dir, mate1_suff,
                                                                    mate2_suff, mates_prefix, max_gap_size,
                                                                    recursion_depth,
                                                                    root_dir, pe_detect_args, pe_uncompress_temp_dir,
                                                                    uncompress_cmd)
            remove_files([pe_uncompress_temp_dir])

            if analysis_mode == ANALYSIS_MODE_DICT[ANALYSIS_MODE_FULL]:
                run_hyper_editing_full_analysis(he_output_dir, is_pe, pe_detect_args, pe_run_detection_logs_dir)

        if bam_out_dir:
            bam_out_dir_full = os.path.join(bam_out_dir, "PE_" + he_detect_args)
            create_he_bam(star_dir, bam_files_suffix, he_output_dir, he_detect_args, is_pe, bam_out_dir_full,
                          follow_links, python3_path, logs_dir, samtools_path, decrease_parallel_run, include_paths,
                          include_paths_operator, exclude_paths, exclude_paths_operator)

            if star_dir:
                merge_he_bam(star_dir, bam_files_suffix, "PE_" + he_detect_args, bam_out_dir_full, follow_links, samtools_path,
                             logs_dir, decrease_parallel_run, remove_recurring_reads, include_paths,
                             include_paths_operator, exclude_paths, exclude_paths_operator)

    if not dont_run:
        to_del = []
        if not keep_unmapped_bams:
            to_del.append(HE_UNMAP_DIR_FORMAT % dict(he_output_dir=he_output_dir))
        if not keep_trans_run:
            to_del.append(HE_TRANS_RUN_DIR_FORMAT % dict(he_output_dir=he_output_dir))

        remove_files(to_del)


if __name__ == '__main__':
    desc = "A wrapper script running the HyperEditing analysis"

    parser = argparse.ArgumentParser(prog='RunHE', description=desc, formatter_class=SortingHelpFormatter)
    # region CMD Line - Parser Options
    # region Input Options
    inputs_g = parser.add_argument_group(title="Input Files Options")
    add_get_paths_function_to_argparser(parser=inputs_g)
    inputs_g.add_argument('-f', '--files_suffix', metavar="files_suffix", dest='files_suffix', nargs='?',
                          required=True, help="The suffix of the files to run on.")
    inputs_g.add_argument('-u', '--uncompress_cmd', metavar="uncompress command", dest='uncompress_cmd', nargs='?',
                          required=True, help="The command to use to extract the FASTQs (e.g. bzcat)")
    inputs_mates_g = parser.add_argument_group(
        title=r"Mates Recognition Options (These should match the format <sample name><mates prefix><mate1\mate2>")
    inputs_mates_g.add_argument('-m1', '--mate1_suffix', metavar="mate 1 suffix", dest='mate1_suff', nargs='?',
                                required=False,
                                help='The suffix of mate 1', default="1")
    inputs_mates_g.add_argument('-m2', '--mate2_suffix', metavar="mate 2 suffix", dest='mate2_suff', nargs='?',
                                required=False,
                                help='The suffix of mate 2. For SE data only use "%s".' % NO_2_MATE_INPUT_STR,
                                default="2")
    inputs_mates_g.add_argument('-mp', '--mates_prefix', metavar="mates prefix", dest='mates_prefix', nargs='?',
                                required=False,
                                help='The prefix of the mates')
    inputs_STAR_g = parser.add_argument_group(
        title=r"STAR directory options (for data normalization, if omitted BWA stats will be used)")
    inputs_STAR_g.add_argument('-so', '--STAR_dir', metavar="STAR alignments directory", dest='STAR_dir', nargs='?',
                               required=False,
                               help='Output path of your STAR run (where the BAM files are)',
                               default=None)
    inputs_STAR_g.add_argument('-srd', '--STAR_recursion_depth', metavar="STAR recursion depth",
                               dest='STAR_recursion_depth',
                               required=False, type=int, default=5,
                               help='The depth to enter recursively in STAR folder to look for the logs (>=1)')
    inputs_STAR_g.add_argument('-b', '--bam_files_suffix', metavar="bam files suffix", dest='bam_files_suffix',
                               nargs='?', required=False, default="Aligned.sortedByCoord.out.bam",
                               help="A suffix of the BAM files to run on (e.g. Sorted.By.Coord.bam). Should be "
                                    "the *full* suffix")
    # endregion

    # region Output Options
    outputs_g = parser.add_argument_group(title="Outputs Options (all output directories are generated automatically)")
    outputs_g.add_argument('-o', '--he_output_dir', metavar="hyper-editing output path", dest='he_output_dir',
                           nargs='?', required=True, help='Output directory for the hyper-editing tool main output.')
    outputs_g.add_argument('-os', '--summary_output_dir', metavar="summary output dir", dest='sum_out_dir', nargs='?',
                           required=False,
                           help='Summary output directory, if wanted (required for summary files output)',
                           default=None)
    outputs_g.add_argument('-ob', '--he_bam_output_dir', metavar="bam output dir", dest='bam_out_dir', nargs='?',
                           required=False,
                           help='HE BAMs output directory, if wanted (else will output to STAR folder or HE folder)',
                           default=None)
    outputs_g.add_argument('-sf', '--sites_folder', metavar="processed sites folder", dest='sites_folder', nargs='?',
                           required=False,  # TODO: upadate docs
                           help='Processed sites output directory, if wanted (required for sites summary'
                                ' files output)', default=None)
    # endregion

    # region Hyper-Editing Run Parameters

    # region Genome Options
    he_run_genomes_g = parser.add_argument_group(title="Hyper-Editing Genomes Options")
    he_run_genomes_g.add_argument('-gi', '--bwa_genome_index', dest='bwa_genome_index_pre', required=False,
                                  help='dirname + prefix of the bwa indexed genome so that there are 5 files found there:'
                                       ' <dirname>/<prefix>.amb, <dirname>/<prefix>.ann, <dirname>/<prefix>.bwt,'
                                       ' <dirname>/<prefix>.pac, <dirname>/<prefix>.sa',
                                  default="/private/common/Data/Genomes/Human/hg38.bwa_index/hg38")
    he_run_genomes_g.add_argument('-gti', '--transformed_bwa_genome_index', dest='transformed_bwa_genome_index_pre',
                                  required=False,
                                  help='dirname + prefix of the bwa indexed genome so that there are 5 files found there:'
                                       ' <dirname>/<prefix>.amb, <dirname>/<prefix>.ann, <dirname>/<prefix>.bwt,'
                                       ' <dirname>/<prefix>.pac, <dirname>/<prefix>.sa',
                                  default=r"/private/common/Data/Genomes/Human/hg38-transformed.bwa_index/hg38")
    he_run_genomes_g.add_argument('-gf', '--genome_fasta', dest='genome_fasta', required=False,
                                  help="full path to the genome's fasta",
                                  default=r"/private/common/Data/Genomes/Human/hg38/hg38.fa")
    # endregion

    # region Run Mode
    he_run_data_mode_g = parser.add_argument_group(title="Hyper-Editing Run Mode Options")
    he_run_data_mode_g.add_argument('-se', '--run_single_end', dest='run_single_end', required=False,
                                    action="store_true",
                                    help='If set, will run single end analysis (can be used with or without PE)',
                                    default=False)
    he_run_data_mode_g.add_argument('-pe', '--run_paired_end', dest='run_paired_end', required=False,
                                    action="store_true",
                                    help='If set, will run paired end analysis (can be used with or without SE)',
                                    default=False)
    he_run_data_mode_g.add_argument('--stranded', dest='stranded', required=False, action="store_true", default=False,
                                    help='If set, will treat data as stranded data, otherwise will transform all to +'
                                         ' strand. *CURRENTLY INACTIVE*')
    # endregion

    # region Hyper-Editing Run Parameters
    he_run_params_g = parser.add_argument_group(title="Hyper-Editing Run Parameters")
    he_run_params_g.add_argument('--phred_score_offset', dest='phred_score_offset', required=False,
                                 help="base quality for PHRED score(33 or 64).", default=33, type=int)
    he_run_params_g.add_argument('--max_gap_size', dest='max_gap_size', required=False,
                                 help="max gap allowed (in bp) between mismatches on mates when running as PE.",
                                 default=500000, type=int)
    he_run_params_g.add_argument('--force_run_bwa_mapping', dest='force_run_bwa_mapping', required=False,
                                 action="store_true",
                                 help="A flag. If set the Hyper Editing tool will run BWA mapping even if analyseMM file"
                                      " exists.")
    he_run_params_g.add_argument('--analysis_mode', dest='analysis_mode', required=False, choices=ANALYSIS_MODE_DICT,
                                 default=ANALYSIS_MODE_FULL, help=ANALYSIS_MODE_HELP_STR)
    he_run_params_g.add_argument('--he_sites_min_pc', dest='he_sites_min_pc', required=False,
                                 help="minimal percentage of hyper edited sites (i.e. mismatches of A->G) from a read's"
                                      " length to count it as hyper edited.", default=0.05, type=float)
    he_run_params_g.add_argument('--he_to_mm_ratio', dest='he_to_mm_ratio', required=False,
                                 help="minimal ratio A->G to all other mismatches in a read to count it as hyper edited.",
                                 default=0.6, type=float)
    he_run_params_g.add_argument('--seq_min_q', dest='seq_min_q', required=False,
                                 help="minimal quality of a mismatch site to count it as editing", default=30, type=int)
    he_run_params_g.add_argument('--max_same_letter_ratio', dest='max_same_letter_ratio', required=False,
                                 help="maximal ratio of a single letter from the entire read (if passed the read is"
                                      " ignored).", default=0.6, type=float)
    he_run_params_g.add_argument('--min_pc_cluster_len', dest='min_pc_cluster_len', required=False,
                                 help="minimal fraction (in percentage) of a cluster from the length of a read, to be "
                                      "considered as a hyper edited cluster.", default=0.1, type=float)
    he_run_params_g.add_argument('--last_cluster_start_index', dest='last_cluster_start_index', required=False,
                                 help="The last position (in percentages of read length) of the FIRST detected editing "
                                      "event in a read, to be considered as a hyper edited cluster.", default=0.8,
                                 type=float)
    he_run_params_g.add_argument('--first_cluster_end_index', dest='first_cluster_end_index', required=False,
                                 help="The first position (in percentages of read length) of the LAST detected editing"
                                      " event in a read, to be considered as a hyper edited cluster.", default=0.2,
                                 type=float)
    # endregion

    # endregion

    # region Rerun and Clean Setting
    wrapper_run_g = parser.add_argument_group(title="Rerun Options")

    wrapper_run_g.add_argument('--repeat_on_noise', dest='repeat_on_noise', required=False, action="store_true",
                               help='If` set, will re-run the script if <max_noise_ratio> %% of samples are below'
                                    ' <noise_level> ', default=False)
    wrapper_run_g.add_argument('-snr', '--noisy_ratio', metavar="noisy samples ratio", dest='max_noisy_samples_ratio',
                               required=False,
                               type=float,
                               default=0.1,
                               help='The ratio of "noisy" samples threshold, if passed will repeat run as PE')
    wrapper_run_g.add_argument('-cl', '--a2g_cleanness_threshold', metavar="A2G cleanness threshold",
                               dest='a2g_cleanness_threshold',
                               required=False,
                               type=float, default=0.80,
                               help='The minimal ratio of A to G mismatches (events) to others to define a sample as "clean"')
    wrapper_run_g.add_argument('-ug', '--upstream_g', metavar="max upstream g", dest='max_upstream_g_pc',
                               required=False,
                               type=float, default=0.20, help='The minimal ratio of  G upstream to the mismatch '
                                                              '(for motif check) to define a sample as "clean"')
    wrapper_run_g.add_argument('-pg', '--downstream_g', metavar="min downstream g", dest='min_downstream_g_pc',
                               required=False,
                               type=float, default=0.25, help='The minimal ratio of  G upstream to the mismatch '
                                                              '(for motif check) to define a sample as "clean"')
    wrapper_run_g.add_argument('--strict_check', dest='strict_check', required=False, action='store_true',
                               help=' If set, will require *both* motif and A2g mismatch ratio to pass, else - either.')
    # endregion

    # ---- Miscellaneous ----
    miscellaneous_options = parser.add_argument_group(title="Miscellaneous")

    miscellaneous_options.add_argument('-c', '--calc_output_only', dest='calc_output_only', required=False,
                                       action='store_true',
                                       help='If set, will only re-caluculte the summery output file.'
                                            ' *note* folders are according to given bash args with the run.')

    miscellaneous_options.add_argument('--keep_unmapped_BAMs', dest='keep_unmapped_bams', required=False,
                                       action='store_true',
                                       help='If set, will not delete the unMap directory')

    miscellaneous_options.add_argument('--keep_trans_run', dest='keep_trans_run', required=False,
                                       action='store_true',
                                       help='If set, will not delete the transformed runs BAMs')
    miscellaneous_options.add_argument('--no_recurring_reads', dest='no_recurring_reads', required=False,
                                       action='store_false',
                                       help='If set, will not remove the reads overlapping between HE and STAR BAM files (for decreasing computation time for HE run on unmapped)')
    miscellaneous_options.add_argument('-fn', '--folder_names_only', dest='folder_names_only', required=False,
                                       action='store_true',
                                       help='If set, will only get summery stats from every and each sample folder instead of the'
                                            ' "all" dir. Use when correlation between reads names and samples doesn\'t exist')
    miscellaneous_options.add_argument('--decrease_parallelization', dest='decrease_parallel', required=False,
                                       help='Set number of parallel runs (<=5)', default=5, type=int)
    miscellaneous_options.add_argument('--dont_run', dest='dont_run', required=False,
                                       action='store_true',
                                       help='If set, will not run analysis only conversions (and summaries)')
    miscellaneous_options.add_argument('--python3_path', dest='python3_path', required=False, default="python3.6",
                                       help='python 3 path')
    miscellaneous_options.add_argument('--samtools_path', dest='samtools_path', required=False, default="samtools-1.9",
                                       help='samtools path')
    # endregion

    options = parser.parse_args()
    assert options.decrease_parallel <= 5, "No more than 5 in parallel!!"
    main(
        # Inputs options
        root_dir=options.root_dir,
        include_paths=options.include_prefixes,
        include_paths_operator=options.include_operator,
        exclude_paths=options.exclude_prefixes,
        exclude_paths_operator=options.exclude_operator,
        recursion_depth=options.recursion_depth,
        follow_links=options.follow_links,
        files_suffix=options.files_suffix,
        uncompress_cmd=options.uncompress_cmd,
        # mates options
        mate1_suff=options.mate1_suff,
        mate2_suff=options.mate2_suff,
        mates_prefix=options.mates_prefix,
        # Aligner (STAR) params
        star_dir=options.STAR_dir,
        star_recursion_depth=options.STAR_recursion_depth,
        bam_files_suffix=options.bam_files_suffix,

        # Output Options
        he_output_dir=options.he_output_dir,
        sum_out_dir=options.sum_out_dir,
        sites_folder=options.sites_folder,
        bam_out_dir=options.bam_out_dir,

        # Hyper-Editing Run Parameters:
        # Genome Options
        bwa_genome_index_pre=options.bwa_genome_index_pre,
        transformed_bwa_genome_index_pre=options.transformed_bwa_genome_index_pre,
        genome_fasta=options.genome_fasta,
        # Run Mode
        run_single_end=options.run_single_end,
        run_paired_end=options.run_paired_end,
        stranded=options.stranded,
        # Hyper-Editing Run Parameters
        phred_score_offset=options.phred_score_offset,
        max_gap_size=options.max_gap_size,
        force_run_bwa_mapping=options.force_run_bwa_mapping,
        analysis_mode=ANALYSIS_MODE_DICT[options.analysis_mode],
        he_sites_min_pc=options.he_sites_min_pc,
        he_to_mm_ratio=options.he_to_mm_ratio,
        seq_min_q=options.seq_min_q,
        max_same_letter_ratio=options.max_same_letter_ratio,
        min_pc_cluster_len=options.min_pc_cluster_len,
        last_cluster_start_index=options.last_cluster_start_index,
        first_cluster_end_index=options.first_cluster_end_index,

        # Rerun and Clean Setting
        repeat_on_noise=options.repeat_on_noise,
        max_noisy_samples_ratio=options.max_noisy_samples_ratio,
        a2g_cleanness_threshold=options.a2g_cleanness_threshold,
        max_upstream_g_pc=options.max_upstream_g_pc,
        min_downstream_g_pc=options.min_downstream_g_pc,
        strict_check=options.strict_check,

        # Miscellaneous
        calc_output_only=options.calc_output_only,
        folder_names_only=options.folder_names_only,
        keep_unmapped_bams=options.keep_unmapped_bams,
        keep_trans_run=options.keep_trans_run,
        decrease_parallel_run=options.decrease_parallel,
        python3_path=options.python3_path,
        samtools_path=options.samtools_path,
        dont_run=options.dont_run,
        remove_recurring_reads=options.no_recurring_reads
    )
