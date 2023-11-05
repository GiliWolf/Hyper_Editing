__author__ = 'Hillel'
"""
This module wraps the run of Hagit's HE tool.
"""
# TODO: change paths to fit also windows...
# =====================imports=====================#
# region Python Builtin Imports
import argparse
import glob
import logging
import os
import multiprocessing
import subprocess
from collections import namedtuple, OrderedDict
from csv import reader
from datetime import date
from operator import and_, or_

# endregion

# region Imports From Internal Modules
if __name__ == "__main__":
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from Outputs.Outputers.CSVOutputer import CSVOutputer
from Outputs.Outputers.RawOutputer import RawOutputer

from DataConverters.HyperEditingBEDConverter import HyperEditingBEDConverter, HyperEditingBEDRec, DICT_KEY_SEP

from STAR242_functios import get_STAR_total_reads, get_STAR_unmapped, get_STAR_uniquely_aligned_reads

from Commons.general_functions import init_logging_dict, convert_args_to_dict, remove_files
from Commons.help_functions import reverse_strand
from Commons.consts import *

from Session.BlackBoxesWrapperSctipts.HEWrapperConsts import *

# endregion

# =====================classes===================#
'''
This namedtuple holds the data about a single mismatch type of HE
'''
HERunSitesData = namedtuple("HERunSitesData", "sites unq_sites clusters merged_clusters")

'''
This namedtuple holds the mapping data for a sample form the HE run
'''
HERunSampleMapData = namedtuple("HERunSampleMapData", "mapped unmapped")


class HyperEditingConvertingException(Exception):
    pass


# =====================functions===================#


# region Help Functions
def get_stricter_mm_pc(curr_pc):
    """
    This function returns a stricter percentage of he_to_mm_ratio (see help string for meaning) for biuldin stricter
    re-run in case of noise.
    :param curr_pc: The previous he_to_mm_ratio.
    :return: The stricter value
    """
    return '0.8' if float(curr_pc) < 0.8 else str(float(curr_pc) * 1.2)


def merge_stats_from_samples(stats_dict, new_stats_dict):
    """
    This function merges the dicts sent to it (in depth merging).
    :param stats_dict: The "old" dict.
    :param new_stats_dict: The dict to merge it with.
    :return: The merged dicts
    """
    for chromosome in new_stats_dict:
        old_chr = stats_dict.setdefault(chromosome, {})
        for site in new_stats_dict[chromosome]:
            old_chr.setdefault(site, []).extend(new_stats_dict[chromosome][site])


def get_io_filename(bash_args, o_dir):
    """
    This function joins the bash arguments to a legal path and formats the path format with it.
    :param bash_args: The current argument for the HE run.
    :param o_dir: The HE output dir.
    :return: The inputs outputs file's path.
    """
    run_type = RUN_TYPE_SE if bash_args[RUN_PAIRED_KEY] == HE_FALSE_VAL else RUN_TYPE_PE
    bash_p_hash = hash(str(bash_args))
    io_filename = os.path.join(o_dir, INPUT_OUTPUT_FILENAME % dict(date=date.today().isoformat(),
                                                                   run_type=run_type, param_hash=bash_p_hash))
    return io_filename


def get_output_path(ifile, se_run, file_re=FILENAME_UNMAPPED_RE):
    """
    This function parse the file name extracting from it the name of the sample and the mate number if exists.
    *Note* If running SE with output names of the format <sample name>-<mate number>, when running PE the BAM files
    won't be created again, thus reducing run time.
    :param ifile:  The file's name.
    :param se_run: If set, will create output name with mate's artificial  number as a preperation for re-run.
    :param file_re:  The regex to parse the file name with.
    :return: The sample name and basic output name for the sample.
    """
    reg_res = file_re.findall(ifile)
    if reg_res:
        reg_res = reg_res[0]
    else:
        return None, None

    if reg_res[0]:  # the sample's name
        name = reg_res[0]
    else:  # if didn't match anything the file wasn't a sample.
        return None, None
    if reg_res[1]:  # the sample's mate number
        if se_run:  # attach the mate's number as a preparation for re-run of PE.
            name += "-" + reg_res[1]
    else:
        if se_run:
            name += MATE1_SUFF
        else:  # if mate number is not present it means that the data is not PE and there for cannot be ran this way.
            return None, None

    return name, reg_res[0]


def calc_noise_and_get_mapping_pc(a2g_stats, a2g_mismatch_min_ratio, max_pre_g_pc=0.20, min_past_g_pc=0.25,
                                  strict_check=False, files_regex=FILENAME_UNMAPPED_RE):
    """
    This function calculates the ratio of "noisy" samples to the whole. To save run time, this function also extracts
     the mapping data.
    :param a2g_stats: The handle to the A2G stats file.
    :param a2g_mismatch_min_ratio: The minimal ratio of A to G of all mismatches to count as "clean"
    :param max_pre_g_pc: The maximal percentage of G before sites (for motif check).
    :param min_past_g_pc: The minimal percentage of G after sites (for motif check).
    :param strict_check: If set, will require *both* motif and A2g mismatch ratio to pass, else - either.
    :param re.Pattern files_regex: The regex to use for extracting the data from sample path.
    :return: The noise ratio(the ratio of "noisy" samples to the whole).
    """
    logical_op = and_ if strict_check else or_
    lines = [l for l in reader(a2g_stats, delimiter=STATS_FILE_DATA_SEP)]
    stats_file_data_start = get_he_stats_file_start(lines)
    lines = lines[stats_file_data_start:]
    ag_he_pc_i, mapped_c_i, past_g_i, pre_g_i, sample_path_i, unmapped_c_i, sample_name_i = get_a2g_stats_indexes(
        lines[0])
    all_count = 0
    noise_count = 0

    samples_paths_to_mapping_data = {}

    for l in lines[1:]:
        # Santiy check for the line
        if not len(l) > ag_he_pc_i and not len(l) > pre_g_i and not len(l) > past_g_i:
            continue
        if not (l[pre_g_i] and l[past_g_i]):
            continue

        if l[sample_path_i] == HE_EMPTY_VAL:  # This is the summery line.
            continue
        # get mapping data
        sample_name = l[sample_name_i]
        sample_path = l[sample_path_i]
        sample_name_from_path = get_output_path(ifile=sample_path,
                                                se_run=MATE1_SUFF in sample_name or MATE2_SUFF in sample_name,
                                                file_re=files_regex)[0]
        if sample_name != sample_name_from_path:
            sample_path = sample_path.replace(get_output_path(sample_path, False, files_regex)[0],
                                              sample_name.replace(MATE1_SUFF, "").replace(MATE2_SUFF, ""))
        samples_paths_to_mapping_data[sample_path] = HERunSampleMapData(l[mapped_c_i], l[unmapped_c_i])

        # noise calc
        all_count += 1
        try:
            ag_he_pc = l[ag_he_pc_i] if l[ag_he_pc_i] else 0.0
            mm_cond = float(ag_he_pc) > a2g_mismatch_min_ratio
            g_up_cond = float(l[pre_g_i]) < max_pre_g_pc
            g_down_cond = float(l[past_g_i]) > min_past_g_pc
            motif_cond = g_up_cond and g_down_cond
        except ValueError:
            logging.warn(CONCAT_STATS_WRN % a2g_stats.name)
            try:  # Try to retrieve the column indexes in case they changed
                ag_he_pc_i, mapped_c_i, past_g_i, pre_g_i, sample_path_i, unmapped_c_i, sample_name_i = get_a2g_stats_indexes(
                    l)
            except ValueError:
                pass
            all_count = 0
            noise_count = 0
            samples_paths_to_mapping_data = dict()
            continue

        if not logical_op(mm_cond, motif_cond):
            # TODO: send a list of noisy samples
            logging.warn(SAMPLE_FAILED_NOISE_CHECK_WRN % dict(sample=l[sample_path_i],
                                                              strict_operator="And" if strict_check else "Or",
                                                              g_up_threshold=max_pre_g_pc, g_up_found=l[pre_g_i],
                                                              g_up_res=COND_SUCCESS if g_up_cond else COND_FAILED,
                                                              g_dwn_threshold=min_past_g_pc, g_dwn_found=l[past_g_i],
                                                              g_dwn_res=COND_SUCCESS if g_down_cond else COND_FAILED,
                                                              mm_r=a2g_mismatch_min_ratio, mm_r_found=l[ag_he_pc_i],
                                                              mm_r_res=COND_SUCCESS if g_up_cond else COND_FAILED, ))
            noise_count += 1
    try:
        noise_ratio = noise_count / float(all_count)
    except ZeroDivisionError:
        noise_ratio = 1

    return noise_ratio, samples_paths_to_mapping_data


def get_he_stats_file_start(lines):
    for i, line in enumerate(lines):
        try:
            line.index(SAMPLE_PATH_COL)
            return i
        except ValueError:
            continue
    return -1


def get_a2g_stats_indexes(line):
    sample_path_i = line.index(SAMPLE_PATH_COL)
    sample_name_i = line.index(SAMPLE_NAME_COL)
    mapped_c_i = line.index(MAPPED_READS_COL)
    unmapped_c_i = line.index(UNMAPPED_READS_COL)
    ag_he_pc_i = line.index(AG_HE_PC_COL)
    pre_g_i = line.index(PRE_G_PC_COL)
    past_g_i = line.index(PAST_G_PC_COL)
    return ag_he_pc_i, mapped_c_i, past_g_i, pre_g_i, sample_path_i, unmapped_c_i, sample_name_i


# endregion

# region Run Hyper Editing
def get_he_bash_config(inputs_outputs_file, output_dir, **run_he_params):
    """
    This function formats the bash string the runs the HE tool.
    :param inputs_outputs_file: The file containing the inputs outputs names for the Hyper Editing tool(see help)
    :param output_dir: The outputdir for the Hyper Editing tool(see help)
    :param run_he_params: different than default values to run HE with them (as depicted in the help).
    :return: The formatted bash string, and params representation string.
    """
    run_dict = dict(RUN_HE_FORMAT_DICT)
    run_dict.update(**run_he_params)
    run_dict.update({INPUT_OUTPUT_FILE_KEY: inputs_outputs_file, OUTPUT_DIR_KEY: output_dir})

    bash_str = RUN_HE_BASH_FORMAT % run_dict
    # NOTE: the following line is dependant on the name of the run_paired option in the RUN_HE_FORMAT_DICT
    run_dict.update(paired=RUN_TYPE_SE if run_dict[RUN_PAIRED_KEY] == HE_FALSE_VAL else RUN_TYPE_PE)
    params_str = HE_OUTFILES_PARAMS_FORMAT % run_dict

    return bash_str, params_str


def run_he(bash_args, max_noise_ratio, noise_level, output_dir, se_input_output, pe_input_output, repeat_on_noise,
           run_paired_end, run_single_end, calc_output_only, max_pre_g_pc, min_past_g_pc, strict_check, files_regex):
    # type: (dict, float, float, str, str, str, bool, bool, bool, bool, float, float, bool, re.Pattern) -> object
    """
    This function wraps the run of the HE tool Hagit wrote.
    :param dict[str, str] bash_args: Any arguments to set in the bash.
    :param float max_noise_ratio: The threshold percentage of "noisy" samples.
    :param float noise_level: The minimal ratio of A to G mismatches to others to define a sample as "clean", otherwise "noisy"
    :param str output_dir:  Outputer dir for the hyper editing output.
    :param str se_input_output: The string for 'inputs outputs file' of the  single end HE run.
    :param str pe_input_output: The string for 'inputs outputs file' of the  paired end HE run.
    :param bool repeat_on_noise: a flag, if set will check for noisy samples ratio (see next params), and if they pass the threshold will re-run with stricter params.
    :param bool run_paired_end: a flag, if set will run paired end version.
    :param bool run_single_end: a flag, if set will run single end version.
    :param bool calc_output_only: a flag, if set will no run the script but rather recalculate the output files.
    :param float max_pre_g_pc: The maximal percentage of G before sites (for motif check).
    :param float min_past_g_pc: The minimal percentage of G after sites (for motif check).
    :param bool strict_check: If set, will require *both* motif and A2g mismatch ratio to pass, else - either.
    :param re.Pattern files_regex: The regex to use to extract sample data from path.
    :return: The string representing the SE run params and the string representing the PE run params, and the mapping
    data.
    """
    se_params = {}
    pe_params = {}
    pe_samples_paths_to_mapping_data = se_samples_paths_to_mapping_data = None
    orig_he_mm_ratio = None
    bash_args = dict(**bash_args)
    o_file_path = os.path.dirname(output_dir)
    o_dir = os.path.basename(output_dir)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    ran_bwa = False  # a flag to know if alignment was done.

    if run_single_end and se_input_output:
        bash_args[RUN_PAIRED_KEY] = HE_FALSE_VAL
        io_filename = get_io_filename(bash_args, o_dir)

        if not calc_output_only:
            with open(os.path.join(o_file_path, io_filename), 'wb') as o:
                o.write(se_input_output)

        # get the bash string to run.
        run_str, se_params = get_he_bash_config(os.path.join(o_file_path, io_filename), output_dir, **bash_args)
        if not calc_output_only:
            logging.debug(RUN_HE_MSG % (run_str, str(bash_args)))
            subprocess.call(run_str, shell=True)
        else:
            logging.debug(RUN_RECALC_HE_MSG % str(bash_args))

        ran_bwa = True  # "remember" that bwa was ran.

        if repeat_on_noise:  # check that the run resulted in "clean" enough results.
            stats_file = HE_STATS_FILE_FORMAT % dict(params=se_params, output_dir=output_dir)
            with open(stats_file) as a2g_stats:
                noise_ratio, se_samples_paths_to_mapping_data = calc_noise_and_get_mapping_pc(a2g_stats, noise_level,
                                                                                              max_pre_g_pc,
                                                                                              min_past_g_pc,
                                                                                              strict_check,
                                                                                              files_regex)
                if noise_ratio > max_noise_ratio:  # if failed - re run with stricter params.
                    logging.debug(NOISY_RERUN_MSG % (noise_ratio, "Or" if strict_check else "And"))
                    bash_args[RUN_BWA_MAPPING_KEY] = HE_FALSE_VAL
                    # store the original he_to_mm_ratio for following PE run.
                    orig_he_mm_ratio = bash_args.get(HE_TO_MM_RATIO_KEY, HE_TO_MM_RATIO)
                    # get the stricter params
                    bash_args[HE_TO_MM_RATIO_KEY] = get_stricter_mm_pc(orig_he_mm_ratio)
                    io_filename = get_io_filename(bash_args, o_dir)
                    if not calc_output_only:
                        with open(os.path.join(o_file_path, io_filename), 'wb') as o:
                            o.write(se_input_output)
                    run_str, se_params = get_he_bash_config(os.path.join(o_file_path, io_filename),
                                                            output_dir, **bash_args)
                    if not calc_output_only:
                        subprocess.call(run_str, shell=True)

    if run_paired_end and pe_input_output:
        bash_args[RUN_PAIRED_KEY] = HE_TRUE_VAL
        bash_args[RUN_BWA_MAPPING_KEY] = HE_FALSE_VAL if ran_bwa else HE_TRUE_VAL
        if ran_bwa:
            if orig_he_mm_ratio:  # if SE was reran using stricter params - set this to the original value
                bash_args[HE_TO_MM_RATIO_KEY] = orig_he_mm_ratio

        io_filename = get_io_filename(bash_args, o_dir)

        if not calc_output_only:
            with open(os.path.join(o_file_path, io_filename), 'wb') as o:
                o.write(pe_input_output)

        run_str, pe_params = get_he_bash_config(os.path.join(o_file_path, io_filename), output_dir, **bash_args)

        if not calc_output_only:
            logging.debug(RUN_HE_MSG % (run_str, str(bash_args)))
            subprocess.call(run_str, shell=True)
        else:
            logging.debug(RUN_RECALC_HE_MSG % str(bash_args))

        if repeat_on_noise:
            stats_file = HE_STATS_FILE_FORMAT % dict(params=pe_params, output_dir=output_dir)
            with open(stats_file) as a2g_stats:
                noise_ratio, pe_samples_paths_to_mapping_data = calc_noise_and_get_mapping_pc(a2g_stats, noise_level,
                                                                                              max_pre_g_pc,
                                                                                              min_past_g_pc,
                                                                                              strict_check,
                                                                                              files_regex)

                if noise_ratio > max_noise_ratio:
                    logging.debug(NOISY_RERUN_MSG % (noise_ratio, "Or" if strict_check else "And"))
                    # get the stricter params
                    bash_args[HE_TO_MM_RATIO_KEY] = get_stricter_mm_pc(
                        bash_args.get(HE_TO_MM_RATIO_KEY, HE_TO_MM_RATIO))
                    bash_args[RUN_BWA_MAPPING_KEY] = HE_FALSE_VAL
                    io_filename = get_io_filename(bash_args, o_dir)
                    with open(os.path.join(o_file_path, io_filename), 'wb') as o:
                        if not calc_output_only:
                            o.write(pe_input_output)
                    run_str, pe_params = get_he_bash_config(os.path.join(o_file_path, io_filename),
                                                            output_dir, **bash_args)
                    if not calc_output_only:
                        subprocess.call(run_str, shell=True)

    return pe_params, pe_samples_paths_to_mapping_data, se_params, se_samples_paths_to_mapping_data


# endregion

# region Files Parsing and Processing


def get_overall_sites_data(output_dir, params, sample_names=None,
                           possible_mismatches=MismatchesAndRefsEnum.HE_NOT_STRANDED_TRANSFORM_DICT.values()):
    """
    This function parses the data from the HE output BED files.
    :param output_dir: The output HE dir.
    :param params: The params string representation.
    :param sample_names: If provided will extract the data from each sample dir instead of the main one.
    :param possible_mismatches: The possible mismatches to check for summery files.
    :rtype L{HERunSitesData}
    :return: HERunSitesData containing the data for all sites, unique sites, clusters and merged cluster.
    """
    sites = {}
    unq_sites = {}
    clusters = {}
    merged_clusters = {}
    run_from_all_folder = True if None is sample_names else False
    converter = HyperEditingBEDConverter(use_folder_names=(not run_from_all_folder))
    if run_from_all_folder:
        sample_names = [ALL_SAMPLES, ]  # TODO: what happens in PE?

    if RUN_TYPE_SE in params:
        run_from_all_folder = False  # TODO: If fastq mate can be derived from records enable this again.
    samples = []
    if not run_from_all_folder and RUN_TYPE_SE in params:
        for sample in sample_names:
            samples.append(sample + MATE1_SUFF)
            samples.append(sample + MATE2_SUFF)
    else:
        samples = sample_names[:]

    for sample in samples:
        for mismatch in possible_mismatches:
            sites_path = SITES_BED_FILE_FORMAT % {"output_dir": output_dir, "mismatch": mismatch, "params": params,
                                                  "sample": sample}
            if os.path.exists(sites_path):
                with open(sites_path, 'rb') as sites_file:
                    if run_from_all_folder:
                        sites[mismatch] = converter.convert(sites_file)
                    else:
                        merge_stats_from_samples(sites.setdefault(mismatch, {}),
                                                 converter.convert(sites_file, sample_name=sample))
                        # In both unique sites BED and merged clusters BED lack the samples supporting each record.
                        #  Thus here we BED intersect the unique sites BED and  merged clusters BED with their non
                        #  unique\unmerged corresponding files or use the not joined files.
            else:
                if RUN_TYPE_SE in params and MATE2_SUFF in sample:
                    logging.debug(MISSING_HE_SE_MATE_OUTFILE_DBG % ("Sites", sites_path, sample))
                    continue
                logging.warn(MISSING_HE_OUTFILE_WRN % ("Sites", sites_path, sample))
            unq_sites_path = UNQ_SITES_BED_FILE_FORMAT % {"output_dir": output_dir, "mismatch": mismatch,
                                                          "params": params,
                                                          "sample": sample}
            if os.path.exists(unq_sites_path):
                merged_path = unq_sites_path + "merged.bed"
                subprocess.call(MERGE_BED_BASH_FORMAT % dict(raw_bed=sites_path, merged_bed=unq_sites_path,
                                                             out_file=merged_path), shell=True)
                with open(merged_path, 'rb') as unq_sites_file:
                    if run_from_all_folder:
                        unq_sites[mismatch] = converter.convert(unq_sites_file)
                    else:
                        merge_stats_from_samples(unq_sites.setdefault(mismatch, {}),
                                                 converter.convert(unq_sites_file, sample_name=sample))
                os.remove(merged_path)
            else:
                logging.warn(MISSING_HE_OUTFILE_WRN % ("Unique Sites", unq_sites_path, sample))

            clusters_path = BED_CLUSTERS_FILE_FORMAT % {"output_dir": output_dir, "mismatch": mismatch,
                                                        "params": params,
                                                        "sample": sample}
            if os.path.exists(clusters_path):
                with open(clusters_path, 'rb') as clusters_file:
                    if run_from_all_folder:
                        clusters[mismatch] = converter.convert(clusters_file, no_neighbors=True)
                    else:
                        merge_stats_from_samples(clusters.setdefault(mismatch, {}),
                                                 converter.convert(clusters_file, sample_name=sample,
                                                                   no_neighbors=True))
            else:
                logging.warn(MISSING_HE_OUTFILE_WRN % ("Clusters", clusters_path, sample))
            if run_from_all_folder:
                merged_c_path_format = BED_CLUSTERS_MERGED_FILE_FORMAT % {"output_dir": output_dir,
                                                                          "mismatch": mismatch,
                                                                          "params": params, "sample": "*"}
                merged_c_paths = glob.glob(merged_c_path_format)
                for merged_c_path in merged_c_paths:
                    sample_mc_name = SAMPLE_NAME_MERGED_CLUSTERS_RE.findall(merged_c_path)[0]
                    if sample_mc_name == ALL_SAMPLES:
                        continue
                    sample_mc_name = re.sub(r"-[12]$", "", sample_mc_name)
                    with open(merged_c_path, 'rb') as merged_clusters_file:
                        merge_stats_from_samples(merged_clusters.setdefault(mismatch, {}),
                                                 converter.convert(merged_clusters_file, sample_name=sample_mc_name,
                                                                   no_neighbors=True, unique=True))
            else:
                merged_c_path = BED_CLUSTERS_MERGED_FILE_FORMAT % {"output_dir": output_dir, "mismatch": mismatch,
                                                                   "params": params, "sample": sample}
                if os.path.exists(merged_c_path):
                    with open(merged_c_path, 'rb') as merged_clusters_file:
                        merge_stats_from_samples(merged_clusters.setdefault(mismatch, {}),
                                                 converter.convert(merged_clusters_file, sample_name=sample,
                                                                   no_neighbors=True, unique=True))
                else:
                    logging.warn(MISSING_HE_OUTFILE_WRN % ("Merged Clusters", merged_c_path, sample))

    return HERunSitesData(sites, unq_sites, clusters, merged_clusters)


def get_bed_recs(samples_to_paths, he_data, sites_recs, sites_to_samples_recs, run_type):
    """
    This function calc and creates the records for overall sites BED file creation and sites to sample specification.
    :param samples_to_paths: A dict of samples to their corresponding paths.
    :param he_data: The data from the single end run for all sites, unique sites, clusters and merged cluster.
    :type he_data: L{HERunSitesData}
    :param sites_recs: Outputer variable containing the sites records.
    :param sites_to_samples_recs: Outputer variable containing the sites samples specification records.
    :param run_type: Either for RUN_TYPE_SE or RUN_TYPE_PE
    :return: None
    """
    for mismatch in he_data.unq_sites:
        for chromosome in he_data.unq_sites[mismatch]:
            for site in he_data.unq_sites[mismatch][chromosome]:
                overall_reads = 0
                samples = ""
                samples_paths = ""
                start, end, strand = site.split(DICT_KEY_SEP)
                he_mm = mismatch if strand == SENSE_STRAND else reverse_strand(seq=mismatch,
                                                                               complementary=False,
                                                                               ignore_non_bases=True)
                for he_rec in he_data.unq_sites[mismatch][chromosome][site]:
                    if not isinstance(he_rec, HyperEditingBEDRec):
                        continue
                    overall_reads += he_rec.reads_count
                    if run_type == RUN_TYPE_PE:
                        samples_paths += ";".join(samples_to_paths[he_rec.sample]) + SAMPLES_AND_PATHS_SEP
                    else:
                        if he_rec.sample.endswith(MATE1_SUFF):
                            samples_paths += samples_to_paths[he_rec.sample[:-2]][0] + SAMPLES_AND_PATHS_SEP
                        else:
                            samples_paths += samples_to_paths[he_rec.sample[:-2]][1] + SAMPLES_AND_PATHS_SEP

                    samples += he_rec.sample + SAMPLES_AND_PATHS_SEP

                sites_recs.append({CHROMOSOME: he_rec.chr, START: he_rec.start, END: he_rec.end, MISMATCH: he_mm,
                                   NUM_OF_READS: str(overall_reads), STRAND: he_rec.strand})
                sites_to_samples_recs.append({CHROMOSOME: he_rec.chr, START: he_rec.start, END: he_rec.end,
                                              MISMATCH: he_mm, SAMPLES: samples, PATHS: samples_paths})
    sites_recs.sort(key=lambda x: [x[CHROMOSOME], int(x[END]), x[STRAND]])


def convert_bed_for_annovar(filename, output_path):
    """
    This function converts the sites file to an annovar compatible format.
    :param filename: The sites BED's path.
    :param output_path:  The path for output.
    :return: None
    """
    anno_lines = []

    with open(filename) as i:
        lines = [line for line in reader(i, delimiter="\t")]
        for line in lines[1:]:
            rec = list()
            rec.append(line[0].replace("chr", ""))
            rec.append(line[2])
            rec.append(line[2])

            fro, to = line[3].split("2")

            rec.append(fro)
            rec.append(to)

            anno_lines.append("\t".join(rec))

    outputer = RawOutputer(pprint_flag=False)
    outputer.output([output_path], "\n".join(anno_lines))


# endregion


# region Cleanup
def sort_bams(bam_path, samtools_path, semaphore):
    """
    This is target for paralleling sorting.
    :param multiprocessing.Semaphore semaphore:
    :param str bam_path: The bam to sort
    :param str samtools_path: The invoke cmd of samtools
    :return: None
    """
    sort_cmd = SORT_BAMS_FORMAT % dict(samtools_path=samtools_path, bam=bam_path)
    subprocess.Popen(sort_cmd, shell=True)
    semaphore.release()


def remove_or_sort_he_files(output_dir, keep_bams):
    """
    This function remove the unnecessary files after the run to save on disk space.
    :param output_dir: The dir the script was run into (--output_dir, -o option)
    :param keep_bams: If set, will not erase the bam files need for re-run, instead sorting them only.
    :return: None
    """
    log_msg = CLEANUP_MSG % "Sorting" if keep_bams else "Removing"
    logging.debug(log_msg)
    trans_folder = TRANS_FOLDER_FORMAT % dict(output_dir=output_dir)
    unmap_folder = UNMAP_FOLDER_FORMAT % dict(output_dir=output_dir)
    del_fastq = glob.glob(os.path.join(unmap_folder, "*.fastq"))
    logging.debug(CLEANUP_DEL_FASTQS_MSG % del_fastq)
    remove_files(del_fastq)

    bam_files = [l.strip() for l in os.popen("find %s -name '*.bam'" % output_dir).read().split("\n")]

    if keep_bams:
        ps = list()
        logging.debug(CLEANUP_BAMS_MSG % ("Sorting", SORT_BAMS_FORMAT))
        semaphore = multiprocessing.Semaphore(5)
        for bam in bam_files:
            ps.append(multiprocessing.Process(target=sort_bams, args=("samtools18", bam, semaphore)))

        for p in ps:
            semaphore.acquire()
            p.start()

        for p in multiprocessing.active_children():
            p.join(3600)
    else:

        logging.debug(CLEANUP_BAMS_MSG % ("Removing", ";".join([trans_folder, unmap_folder])))
        remove_files(bam_files)


# endregion

# region Outputer Creation
def create_sites_bed(se_he_data, pe_he_data, sites_folder, samples_to_paths):
    """
    This function creates the overall sites BED file, the annovar input and the sties to samples specification file.
    :param se_he_data: The data from the single end run for all sites, unique sites, clusters and merged cluster.
    :type se_he_data: L{HERunSitesData}
    :param pe_he_data: The data from the paired end run for all sites, unique sites, clusters and merged cluster.
    :type pe_he_data: L{HERunSitesData}
    :param sites_folder: The folder to output the file to.
    :param samples_to_paths: A dict of samples to their corresponding paths.
    :return: None
    """
    outputer = CSVOutputer(delim="\t")
    reditools_outputer = CSVOutputer(delim="\t", write_headers=False)

    site_headers = [CHROMOSOME, START, END, MISMATCH, NUM_OF_READS, STRAND]
    site_reditools_headers = [CHROMOSOME, END, STRAND]
    se_sites_recs = []
    pe_sites_recs = []
    se_sites_to_samples_recs = []
    pe_sites_to_samples_recs = []

    if se_he_data:
        se_sites_file_name = os.path.join(sites_folder, SITES_OUTPUT_FILE % RUN_TYPE_SE)
        se_reditools_file_name = os.path.join(sites_folder, SITES_REDITOOLS_OUTPUT_FILE % RUN_TYPE_SE)
        get_bed_recs(samples_to_paths, se_he_data, se_sites_recs, se_sites_to_samples_recs, RUN_TYPE_SE)
        outputer.output([se_sites_file_name], site_headers, se_sites_recs)
        reditools_outputer.output([se_reditools_file_name], site_reditools_headers, se_sites_recs)

    if pe_he_data:
        pe_sites_file_name = os.path.join(sites_folder, SITES_OUTPUT_FILE % RUN_TYPE_PE)
        pe_reditools_file_name = os.path.join(sites_folder, SITES_REDITOOLS_OUTPUT_FILE % RUN_TYPE_PE)
        get_bed_recs(samples_to_paths, pe_he_data, pe_sites_recs, pe_sites_to_samples_recs, RUN_TYPE_PE)
        outputer.output([pe_sites_file_name], site_headers, pe_sites_recs)
        reditools_outputer.output([pe_reditools_file_name], site_reditools_headers, pe_sites_recs)

    sites_to_sample_headers = [CHROMOSOME, START, END, MISMATCH, SAMPLES, PATHS]
    sites_pe_file_name = os.path.join(sites_folder, SITES_TO_SAMPLE_OUTPUT_FILE % "PE_SamplesSpec")
    sites_se_file_name = os.path.join(sites_folder, SITES_TO_SAMPLE_OUTPUT_FILE % "SE_SamplesSpec")

    outputer.output([sites_se_file_name], sites_to_sample_headers, se_sites_to_samples_recs)
    outputer.output([sites_pe_file_name], sites_to_sample_headers, pe_sites_to_samples_recs)

    anno_sites_se_file_name = os.path.join(sites_folder, SITES_TO_SAMPLE_OUTPUT_FILE % "SE_SamplesSpecForAnnovar")
    anno_sites_pe_file_name = os.path.join(sites_folder, SITES_TO_SAMPLE_OUTPUT_FILE % "PE_SamplesSpecForAnnovar")

    convert_bed_for_annovar(sites_se_file_name, anno_sites_se_file_name)
    convert_bed_for_annovar(sites_pe_file_name, anno_sites_pe_file_name)


def create_summery_output(se_he_data, se_params, se_samples_paths_to_mapping_data, pe_he_data,
                          pe_samples_paths_to_mapping_data, pe_params, samples_to_paths, STAR_output,
                          STAR_recursion_depth, output_path, stranded):
    """
    This function creates the summery outputs.
    :param se_params:
    :param pe_params:
    :param se_he_data: The data from the single end run for all sites, unique sites, clusters and merged cluster.
    :type se_he_data: L{HERunSitesData}
    :param se_samples_paths_to_mapping_data: Dictionary of sample paths to thier mapping data from the HE tool for the SE run.
    :type se_samples_paths_to_mapping_data: C{dict} of C{str} to L{HERunSampleMapData}
    :param pe_he_data: The data from the paired end run for all sites, unique sites, clusters and merged cluster.
    :type pe_he_data: L{HERunSitesData}
    :param pe_samples_paths_to_mapping_data: Dictionary of sample paths to thier mapping data from the HE tool for the SE run.
    :type pe_samples_paths_to_mapping_data: C{dict} of C{str} to L{HERunSampleMapData}
    :param samples_to_paths: a dictionary containing sample names to their paths.
    :param STAR_output: The output dir of the STAR run to check for STAR alignment logs (for mapping stats).
    :param STAR_recursion_depth: Maximal recursion depth to check for STAR logs inside <STAR_output>.
    :param output_path: The path to output to.
    :param stranded: True if stranded.
    :return: None
    """
    assert (se_he_data is None or isinstance(se_he_data, HERunSitesData)) and \
           (isinstance(pe_he_data, HERunSitesData) or pe_he_data is None)

    outputer = CSVOutputer(delim=",")
    samples_to_recs = {}
    stranded_f = STRANDED if stranded else NOT_STRANDED
    summery_spec_full_path_format = os.path.join(output_path, SPEC_FILE_NAME)

    reads = {}

    # TODO: add other aligners (maybe through DB)
    if STAR_output:  # get mapped, unmapped, and total fragments for each sample from STAR logs.
        for sample in samples_to_paths:
            for i in xrange(STAR_recursion_depth):
                star_o_p = STAR_output + (os.sep + "*" + os.sep) * i + os.sep + sample
                if len(glob.glob(star_o_p)) == 1:
                    try:
                        reads[sample] = (
                            get_STAR_total_reads(glob.glob(star_o_p)[0])[0],
                            get_STAR_unmapped(glob.glob(star_o_p)[0])[0],
                            get_STAR_uniquely_aligned_reads(glob.glob(star_o_p)[0])[0])
                    except TypeError:
                        reads[sample] = (ERROR_VAL, ERROR_VAL, ERROR_VAL)
                    break
    if se_he_data:
        added_headers = {}
        sum_spec_recs_sites = {}
        sum_spec_recs_clusters = {}
        sum_spec_headers = []
        sum_spec_headers_clusters = []
        samples_to_recs[RUN_TYPE_SE] = {}
        run_type = RUN_TYPE_SE
        calc_spec_output_recs(added_headers, reads, run_type, samples_to_paths, samples_to_recs, se_he_data,
                              sum_spec_headers, sum_spec_headers_clusters, sum_spec_recs_clusters, sum_spec_recs_sites)
        logging.debug(SUMMERY_CREATION_MSG % summery_spec_full_path_format % dict(data_type=DATA_TYPE_UNQ_SITES,
                                                                                  stranded=stranded_f,
                                                                                  run_type=run_type, params=se_params))
        outputer.output([summery_spec_full_path_format % dict(data_type=DATA_TYPE_UNQ_SITES, stranded=stranded_f,
                                                              run_type=run_type, params=se_params)],
                        sum_spec_headers, sum_spec_recs_sites.values())

        logging.debug(SUMMERY_CREATION_MSG % summery_spec_full_path_format % dict(data_type=DATA_TYPE_MERGED_CLUSTERS,
                                                                                  stranded=stranded_f,
                                                                                  run_type=run_type, params=se_params))
        outputer.output([summery_spec_full_path_format % dict(data_type=DATA_TYPE_MERGED_CLUSTERS, stranded=stranded_f,
                                                              run_type=run_type, params=se_params)],
                        sum_spec_headers_clusters, sum_spec_recs_clusters.values())

    if pe_he_data:
        added_headers = {}
        sum_spec_recs_sites = {}
        sum_spec_recs_clusters = {}
        sum_spec_headers = []
        sum_spec_headers_clusters = []
        samples_to_recs[RUN_TYPE_PE] = {}
        run_type = RUN_TYPE_PE
        calc_spec_output_recs(added_headers, reads, run_type, samples_to_paths, samples_to_recs, pe_he_data,
                              sum_spec_headers, sum_spec_headers_clusters, sum_spec_recs_clusters, sum_spec_recs_sites)

        logging.debug(SUMMERY_CREATION_MSG % summery_spec_full_path_format % dict(data_type=DATA_TYPE_UNQ_SITES,
                                                                                  stranded=stranded_f,
                                                                                  run_type=run_type, params=pe_params))
        outputer.output([summery_spec_full_path_format % dict(data_type=DATA_TYPE_UNQ_SITES, stranded=stranded_f,
                                                              run_type=run_type, params=pe_params)],
                        sum_spec_headers, sum_spec_recs_sites.values())

        logging.debug(SUMMERY_CREATION_MSG % summery_spec_full_path_format % dict(data_type=DATA_TYPE_MERGED_CLUSTERS,
                                                                                  stranded=stranded_f,
                                                                                  run_type=run_type, params=pe_params))
        outputer.output([summery_spec_full_path_format % dict(data_type=DATA_TYPE_MERGED_CLUSTERS, stranded=stranded_f,
                                                              run_type=run_type, params=pe_params)],
                        sum_spec_headers_clusters, sum_spec_recs_clusters.values())

    # TODO: aligner data from DB
    output_headers = [SAMPLE, PATH, TOTAL_FRAGMENTS, HE_MAPPED_FRAGMENTS, HE_UNMAPPED_FRAGMENTS]
    if reads != {}:
        output_headers.extend([ALIGNER_UNMAPPED_FRAGMENTS, ALIGNER_MAPPED_FRAGMENTS, ])
    se_output_recs = {}
    pe_output_recs = {}
    if se_he_data:
        calc_summery_records(output_headers, reads, samples_to_paths, samples_to_recs, se_he_data, RUN_TYPE_SE,
                             se_samples_paths_to_mapping_data, se_output_recs)
    if pe_he_data:
        calc_summery_records(output_headers, reads, samples_to_paths, samples_to_recs, pe_he_data, RUN_TYPE_PE,
                             pe_samples_paths_to_mapping_data, pe_output_recs)
    logging.debug(SUMMERY_CREATION_MSG % os.path.join(output_path, SUMMERY_FILENAME) % RUN_TYPE_SE)
    outputer.output([os.path.join(output_path, SUMMERY_FILENAME) % RUN_TYPE_SE], output_headers,
                    se_output_recs.values())
    logging.debug(SUMMERY_CREATION_MSG % os.path.join(output_path, SUMMERY_FILENAME) % RUN_TYPE_PE)
    outputer.output([os.path.join(output_path, SUMMERY_FILENAME) % RUN_TYPE_PE], output_headers,
                    pe_output_recs.values())


def calc_summery_records(output_headers, reads, samples_to_paths, samples_to_recs, he_data, run_type,
                         samples_paths_to_mapping_data, output_recs):
    """
    This function builds the output records for the summery file.
    :param output_headers: The base headers list for the output.
    :param reads: A dict of sample to [mapped, unmapped, and total fragments] from STAR logs.
    :param samples_to_paths: A dictionary of sample name to its paths.
    :param samples_to_recs: A dict to populate with records of data for each sample to later save run time.
    :param he_data: HE data from the single end run for all sites, unique sites, clusters and merged cluster.
    :param run_type: RUN_TYPE_SE or RUN_TYPE_PE.
    :param samples_paths_to_mapping_data: dict containing paths to mapping data.
    :param output_recs: Outputer variable of the recs.
    :return: None
    :type output_recs: C{dict}
    """
    if None is he_data:
        return
    if samples_to_recs[run_type] == {}:
        return
    output_headers.extend([h % run_type for h in OUTPUT_HEADERS])
    output_samples_to_paths = {}
    for sample, paths in samples_to_paths.iteritems():
        if run_type == RUN_TYPE_SE:
            output_samples_to_paths[sample + MATE1_SUFF] = paths[0]
            if len(paths) > 1:
                output_samples_to_paths[sample + MATE2_SUFF] = paths[1]
        else:
            output_samples_to_paths[sample] = paths[0]

    for sample, path in output_samples_to_paths.iteritems():
        output_recs[sample] = {SAMPLE: sample, PATH: path}
        try:
            he_mapped_reads = samples_paths_to_mapping_data[path].mapped
            he_unmapped_reads = samples_paths_to_mapping_data[path].unmapped
            he_total_reads = he_mapped_reads + he_unmapped_reads
            output_recs[sample].update({TOTAL_FRAGMENTS: str(he_total_reads),
                                        HE_MAPPED_FRAGMENTS: str(he_mapped_reads),
                                        HE_UNMAPPED_FRAGMENTS: str(he_unmapped_reads)})
            num_of_sites = sum(
                [s.reads_count for s in samples_to_recs[run_type][MISMATCH_FOR_SUMMERY][sample][DATA_TYPE_SITES]])
            num_of_clusters = sum(
                [s.reads_count for s in samples_to_recs[run_type][MISMATCH_FOR_SUMMERY][sample][DATA_TYPE_CLUSTERS]])
            num_of_unq_sites = len(samples_to_recs[run_type][MISMATCH_FOR_SUMMERY][sample][DATA_TYPE_UNQ_SITES])
            num_of_merged_clusters = len(
                samples_to_recs[run_type][MISMATCH_FOR_SUMMERY][sample][DATA_TYPE_MERGED_CLUSTERS])
            if reads != {}:
                sample_name = sample.replace(MATE1_SUFF, "").replace(MATE2_SUFF, "")
                output_recs[sample][TOTAL_FRAGMENTS] = reads[sample_name][OUTPUT_STAR_ALL_READS]
                output_recs[sample][ALIGNER_UNMAPPED_FRAGMENTS] = reads[sample_name][OUTPUT_STAR_UNMAPPED_READS]
                output_recs[sample][ALIGNER_MAPPED_FRAGMENTS] = reads[sample_name][OUTPUT_STAR_MAPPED_READS]
                mapped_reads = float(reads[sample_name][OUTPUT_STAR_MAPPED_READS]) \
                    if reads[sample_name][OUTPUT_STAR_MAPPED_READS] != ERROR_VAL else -1
            else:
                mapped_reads = float(he_mapped_reads)

            norm_num_of_sites = 1000000 * num_of_sites / mapped_reads
            norm_num_of_clusters = 1000000 * num_of_clusters / mapped_reads
            norm_num_of_unq_sites = 1000000 * num_of_unq_sites / mapped_reads
            norm_num_of_merged_clusters = 1000000 * num_of_merged_clusters / mapped_reads
        except KeyError:
            logging.warn(MISSING_HE_OUT_SKIPPING_WRN % sample)
            num_of_sites = num_of_clusters = num_of_unq_sites = num_of_merged_clusters = ERROR_VAL
            norm_num_of_sites = norm_num_of_clusters = norm_num_of_unq_sites = norm_num_of_merged_clusters = ERROR_VAL
            if reads != {}:
                norm_num_of_sites = norm_num_of_clusters = norm_num_of_unq_sites = norm_num_of_merged_clusters = ERROR_VAL
        output_recs[sample].update({
            MMF_NORM_TOTAL_SITES_HEADER_FORAMT % run_type: str(norm_num_of_sites),
            MMF_NORM_UNQ_SITES_HEADER_FORAMT % run_type: str(norm_num_of_unq_sites),
            MMF_NORM_CLUSTERS_HEADER_FORAMT % run_type: str(norm_num_of_clusters),
            MMF_NORM_MERGED_CLUSTERS_HEADER_FORAMT % run_type: str(norm_num_of_merged_clusters)})
        output_recs[sample].update({TOTAL_SITES_HEADER_FORAMT % run_type: str(num_of_sites),
                                    UNQ_SITES_HEADER_FORAMT % run_type: str(num_of_unq_sites),
                                    CLUSTERS_HEADER_FORAMT % run_type: str(num_of_clusters),
                                    MERGED_CLUSTERS_HEADER_FORAMT % run_type: str(num_of_merged_clusters)})


def calc_spec_output_recs(added_headers, reads, run_type, samples_to_paths, samples_to_recs, he_data,
                          sum_spec_headers, sum_spec_headers_clusters, sum_spec_recs_clusters, sum_spec_recs_sites):
    """
    This function creates the records for the summery output.
    :param added_headers: A flag dict marking all headers added to a single output file so far.
    :param reads: a dict of sample to [mapped, unmapped, and total fragments] from STAR logs.
    :param run_type: RUN_TYPE_SE or RUN_TYPE_PE
    :param samples_to_paths: A dictionary of sample name to its paths.
    :param samples_to_recs: A dict to populate with records of data for each sample to later save run time.
    :param he_data: HE data from the single end run for all sites, unique sites, clusters and merged cluster.
    :param sum_spec_headers: Basic headers for the sites summery specification.
    :param sum_spec_headers_clusters: Basic headers for the sites summery specification.
    :param sum_spec_recs_clusters: Outputer variable for the cluster output recs.
    :param sum_spec_recs_sites: Outputer variable for the sites output recs.
    :return: None
    """
    spec_headers = [CHROMOSOME, START, END, MISMATCH, STRAND, EXPRESSING_SAMPLES, NEIGHBORS]
    sum_spec_headers_clusters.extend(spec_headers)
    sum_spec_headers.extend(spec_headers)

    data_type = DATA_TYPE_UNQ_SITES
    for mismatch in he_data.unq_sites:
        for chromosome in he_data.unq_sites[mismatch]:
            for site in he_data.unq_sites[mismatch][chromosome]:
                key = chromosome + mismatch + site  # create a key to join releveant data ont he same site.
                start, end, strand = site.split(DICT_KEY_SEP)
                sum_spec_recs_sites[key] = {CHROMOSOME: chromosome, START: start, END: end,
                                            MISMATCH: mismatch, STRAND: strand,
                                            EXPRESSING_SAMPLES: str(
                                                len(he_data.unq_sites[mismatch][chromosome][site]))}
                for he_rec in he_data.unq_sites[mismatch][chromosome][site]:
                    if not isinstance(he_rec, HyperEditingBEDRec):
                        continue

                    s_header = SAMPLE_HEADER_FORMAT % dict(run_type=run_type, data_type=data_type,
                                                           sample=he_rec.sample)
                    if added_headers.get(s_header, True):
                        sum_spec_headers.append(s_header)
                        added_headers[s_header] = False

                    sum_spec_recs_sites[key][s_header] = str(he_rec.reads_count)
                    if NEIGHBORS not in sum_spec_recs_sites[key]:
                        sum_spec_recs_sites[key][NEIGHBORS] = he_rec.neighbors

                    # Add the record to the dictionary for later.
                    samples_to_recs[run_type].setdefault(mismatch, {})
                    samples_to_recs[run_type][mismatch].setdefault(he_rec.sample, {})
                    samples_to_recs[run_type][mismatch][he_rec.sample].setdefault(DATA_TYPE_UNQ_SITES, []).append(
                        he_rec)

    # Add the rest of records to the dictionary for later.
    try:
        for mismatch in he_data.sites:
            for chromosome in he_data.sites[mismatch]:
                for site in he_data.sites[mismatch][chromosome]:
                    for rec in he_data.sites[mismatch][chromosome][site]:
                        samples_to_recs[run_type][mismatch][rec.sample].setdefault(DATA_TYPE_SITES, []).append(rec)
        for mismatch in he_data.clusters:
            for chromosome in he_data.clusters[mismatch]:
                for site in he_data.clusters[mismatch][chromosome]:
                    for rec in he_data.clusters[mismatch][chromosome][site]:
                        samples_to_recs[run_type][mismatch][rec.sample].setdefault(DATA_TYPE_CLUSTERS, []).append(rec)
        for mismatch in he_data.merged_clusters:
            for chromosome in he_data.merged_clusters[mismatch]:
                for site in he_data.merged_clusters[mismatch][chromosome]:
                    for rec in he_data.merged_clusters[mismatch][chromosome][site]:
                        samples_to_recs[run_type][mismatch][rec.sample].setdefault(DATA_TYPE_MERGED_CLUSTERS, []). \
                            append(rec)
    except KeyError:
        try:
            sample_name = rec.sample
        except (UnboundLocalError, NameError):
            sample_name = None
        if sample_name == "Unknown":  # TODO: make this const in converter...
            raise HyperEditingConvertingException(
                "Could Not Convert Outputer Files Correctly, Probably Read Names In The"
                " Fastq Files Don't Contain Sample Name. Try Using The -fn Option")
        else:
            raise HyperEditingConvertingException("Could Not Convert Outputer Files Correctly, Unknown Problem!")

    data_type = DATA_TYPE_MERGED_CLUSTERS
    for mismatch in he_data.merged_clusters:
        for chromosome in he_data.merged_clusters[mismatch]:
            for site in he_data.merged_clusters[mismatch][chromosome]:
                key = chromosome + mismatch + site
                start, end, strand = site.split(DICT_KEY_SEP)
                sum_spec_recs_clusters[key] = {CHROMOSOME: chromosome, START: start, END: end,
                                               MISMATCH: mismatch, STRAND: strand,
                                               EXPRESSING_SAMPLES: str(
                                                   len(he_data.merged_clusters[mismatch][chromosome][site]))}
                for he_rec in he_data.merged_clusters[mismatch][chromosome][site]:
                    if not isinstance(he_rec, HyperEditingBEDRec):
                        continue
                    s_header = SAMPLE_HEADER_FORMAT % dict(run_type=run_type, data_type=data_type,
                                                           sample=he_rec.sample)
                    if added_headers.get(s_header, True):
                        sum_spec_headers_clusters.append(s_header)
                        added_headers[s_header] = False

                    sum_spec_recs_clusters[key][s_header] = str(he_rec.reads_count)

    # Add records describing the path and mapped fragments for each sample.
    paths = {CHROMOSOME: "Paths"}
    frags = {CHROMOSOME: "Fragments"}
    for sample in samples_to_paths:
        s_header = SAMPLE_HEADER_FORMAT % dict(run_type=run_type, data_type=DATA_TYPE_UNQ_SITES,
                                               sample=sample)
        paths[s_header] = samples_to_paths[sample][0]
    if reads != {}:
        for sample in samples_to_paths:
            s_header = SAMPLE_HEADER_FORMAT % dict(run_type=run_type, data_type=DATA_TYPE_UNQ_SITES,
                                                   sample=sample)
            try:
                frags[s_header] = reads[sample][OUTPUT_STAR_MAPPED_READS]
            except KeyError:
                frags[s_header] = ERROR_VAL
    for sample in samples_to_paths:
        s_header = SAMPLE_HEADER_FORMAT % dict(run_type=run_type, data_type=DATA_TYPE_MERGED_CLUSTERS,
                                               sample=sample)
        paths[s_header] = samples_to_paths[sample][0]
    if reads != {}:
        for sample in samples_to_paths:
            s_header = SAMPLE_HEADER_FORMAT % dict(run_type=run_type, data_type=DATA_TYPE_MERGED_CLUSTERS,
                                                   sample=sample)
            try:
                frags[s_header] = reads[sample][OUTPUT_STAR_MAPPED_READS]
            except KeyError:
                frags[s_header] = ERROR_VAL
    sum_spec_recs_clusters["frags"] = frags
    sum_spec_recs_clusters["paths"] = paths
    sum_spec_recs_sites["frags"] = frags
    sum_spec_recs_sites["paths"] = paths


# endregion


#: Main
def run_he_script(input_dir, output_dir, files_suffix, output_files_dir=None, sites_folder=DEFAULT_SITES_OUTPUT_FOLDER,
                  run_single_end=True, run_paired_end=True, repeat_on_noise=True, max_noise_ratio=0.05,
                  noise_level=0.85, STAR_output=None, STAR_recursion_depth=5, calc_output_only=False,
                  stranded=False, keep_bams=False, max_pre_g_pc=0.20, min_past_g_pc=0.25, strict_check=False,
                  folder_names_only=False, **bash_args):
    """
    This is the main of the wrapper.
    :param max_pre_g_pc: The maximal percentage of G before sites (for motif check).
    :param min_past_g_pc: The minimal percentage of G after sites (for motif check).
    :param strict_check: If set, will require *both* motif and A2g mismatch ratio to pass, else - either.
    :param input_dir: The input dir to run on (where the fastqs to read are).
    :param output_dir:  Outputer dir for the hyper editing output.
    :param files_suffix: The suffix of the fastq files for the script.
    :param output_files_dir: Outputer summery file path, if wanted
    :param sites_folder: The path of the folder to output to all found sites ( in BED format).
    :param run_single_end: a flag, if set will run single end version.
    :param run_paired_end: a flag, if set will run paired end version.
    :param repeat_on_noise: a flag, if set will check for noisy samples ratio (see next params), and if they pass the
     threshold will re-run with stricter params.
    :param max_noise_ratio: The threshold percentage of "noisy" samples.
    :param noise_level: The minimal ratio of A to G mismatches to others to define a sample as "clean" otherwise "noisy"
    :param STAR_output: The output dir of the STAR run to check for STAR alignment logs (for mapping stats).
    :param STAR_recursion_depth: Maximal recursion depth to check for STAR logs inside <STAR_output>.
    :param calc_output_only: a flag, if set will no run the script but rather recalculate the output files.
    :param stranded: If set, will treat data as stranded.
    :param keep_bams: If set, will not erase the HE BAM files, but only sort them, enabling a quick re-run.
    :param folder_names_only: If set will extract data from every sample dir instead of the "all" dir.
    :param bash_args: Any arguments to set in the bash.
    :return: None
    """

    logging_args_dict = dict(input_dir=str(input_dir), output_dir=str(output_dir),
                             files_suffix=str(files_suffix), output_files_dir=str(output_files_dir),
                             sites_folder=str(sites_folder), run_single_end=str(run_single_end),
                             run_paired_end=str(run_paired_end), repeat_on_noise=str(repeat_on_noise),
                             STAR_output=str(STAR_output), STAR_recursion_depth=str(STAR_recursion_depth),
                             stranded=str(stranded), calc_output_only=str(calc_output_only),
                             keep_bams=str(keep_bams), max_pre_g_pc=str(max_pre_g_pc),
                             min_past_g_pc=str(min_past_g_pc), strict_check=str(strict_check),
                             folder_names_only=str(folder_names_only), run_bash_args=str(bash_args))
    logging.debug(RUN_PARAMS_DEBUG % logging_args_dict)

    se_input_output = []
    pe_input_output = []

    samples_to_paths = {}
    #  get the corresponding file names regex (for either fastqs or STAR unmapped files)
    # TODO: this somehow should be more generic
    file_re = NAME_TYPE_TO_RE[files_suffix.lower()]

    for root, dirs, files in os.walk(input_dir, followlinks=True):
        for ifile in files:
            if files_suffix.lower() in ifile.lower():
                full_path = os.path.join(root, ifile)
                # get the output path and sample name for SE run.
                out_path, s_name = get_output_path(ifile, True, file_re)
                if None is out_path:
                    logging.warn(BAD_PATH_WARN % ifile)
                    continue
                samples_to_paths.setdefault(s_name, []).append(full_path)
                se_input_output.append(full_path + "\t" + out_path)
    se_input_output.sort()
    se_input_output = "\n".join(se_input_output) + "\n"

    if run_paired_end:
        for sample, paths in samples_to_paths.iteritems():
            if len(paths) != 2:
                continue
            for path in paths:
                out_path, s_name = get_output_path(path, False, file_re)
                if None is out_path:
                    logging.warn(BAD_PATH_WARN % path)
                    continue
                pe_input_output.append(path + "\t" + out_path)
        pe_input_output.sort()
        pe_input_output = "\n".join(pe_input_output) + "\n"

    # run the HE tool
    pe_params, pe_samples_paths_to_mapping_data, se_params, se_samples_paths_to_mapping_data = \
        run_he(bash_args=bash_args, max_noise_ratio=max_noise_ratio, noise_level=noise_level, output_dir=output_dir,
               se_input_output=se_input_output, pe_input_output=pe_input_output, repeat_on_noise=repeat_on_noise,
               run_paired_end=run_paired_end, run_single_end=run_single_end, calc_output_only=calc_output_only,
               max_pre_g_pc=max_pre_g_pc, min_past_g_pc=min_past_g_pc, strict_check=strict_check, files_regex=file_re)

    se_he_data = pe_he_data = None
    sample_names = samples_to_paths.keys() if folder_names_only else None
    possible_mismatches = MismatchesAndRefsEnum.MISMATCH_TYPE_PARSE.values() if stranded else \
        MismatchesAndRefsEnum.HE_NOT_STRANDED_TRANSFORM_DICT.values()
    if se_params:
        se_he_data = get_overall_sites_data(output_dir, se_params, samples_to_paths.keys(),
                                            possible_mismatches)  # See in function todo.

    if pe_params:
        pe_he_data = get_overall_sites_data(output_dir, pe_params, sample_names, possible_mismatches)

    if output_files_dir:
        create_summery_output(se_he_data=se_he_data,
                              se_samples_paths_to_mapping_data=se_samples_paths_to_mapping_data,
                              se_params=se_params,
                              pe_he_data=pe_he_data,
                              pe_samples_paths_to_mapping_data=pe_samples_paths_to_mapping_data,
                              pe_params=pe_params,
                              samples_to_paths=samples_to_paths,
                              STAR_output=STAR_output,
                              STAR_recursion_depth=STAR_recursion_depth,
                              output_path=output_files_dir,
                              stranded=stranded
                              )
    logging.debug(BED_AND_SPEC_CREATION_MSG)
    create_sites_bed(se_he_data, pe_he_data, sites_folder, samples_to_paths)

    if not calc_output_only:
        remove_or_sort_he_files(output_dir, keep_bams)


def process_he_analysis_table(table_str, from_line=0):
    # type: (str, int) -> list
    data = list()
    lines = table_str.split("\n")[from_line:]

    indexes = {i: header for i, header in enumerate(lines[0].split(STATS_FILE_DATA_SEP)) if header}
    for line in lines[1:]:
        rec = OrderedDict()
        line_recs = line.split(STATS_FILE_DATA_SEP)
        for i in indexes:
            rec[indexes[i]] = line_recs[i]

        data.append(rec)

    return data


def get_data_from_analyse_results(analyse_filename, output_path):
    out_recs = list()
    with open(analyse_filename) as analyse_fh:
        data = analyse_fh.read()

        for sample_record in ANALYSIS_SAMPLE_RE.findall(data):
            res = OrderedDict()
            res[SAMPLE] = sample_record[1]
            res[SAMPLE_PATH] = os.path.join(sample_record[0], sample_record[1])
            try:
                he_data, strand_bias_data, motif_data, _ = sample_record[2].split("\n\n\n")
            except ValueError:
                he_data, strand_bias_data, motif_data = sample_record[2].split("\n\n\n")

            he_data = process_he_analysis_table(he_data.strip("\n"))
            combined_data_dict = OrderedDict()
            for rec in he_data:
                mm_type = rec.pop(MM_EDIT_TYPE)
                combined_data_dict[mm_type] = rec.copy()

            strand_bias_data = process_he_analysis_table(strand_bias_data.strip("\n"), from_line=1)
            for rec in strand_bias_data:
                mm_type = rec.pop(MM_EDIT_TYPE)
                for key in sorted(rec):
                    combined_data_dict[mm_type][key] = rec[key]

            motif_data = process_he_analysis_table(motif_data.strip("\n"), from_line=1)
            for rec in motif_data:
                mm_type = rec.pop(MM_EDIT_TYPE)
                for key in sorted(rec):
                    combined_data_dict[mm_type][key] = rec[key]

            for mm_type in combined_data_dict:
                full_rec = res.copy()
                full_rec[MISMATCH] = mm_type
                full_rec.update(combined_data_dict[mm_type])
                out_recs.append(full_rec)

    outputer = CSVOutputer()
    outputer.output(paths=[output_path, ],
                    headers=full_rec.keys(),
                    records=out_recs)


if __name__ == '__main__':
    desc = "A wrapper script running the HyperEditing analysis"


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='RunHE', description=desc, formatter_class=MyFormatter)

    # region CMD Line - Parser Options
    # ---- io options ----
    parser.add_argument('-i', '--input_dir', metavar="input path", dest='input_dir', nargs='?', required=True,
                        help='The input dir to run on (where the fastqs for the Hyper Editing tool are).')
    parser.add_argument('-o', '--output_dir', metavar="output path", dest='output_dir', nargs='?', required=True,
                        help='Outputer dir for the output')
    parser.add_argument('-f', '--files_suffix_re', metavar="files suffix regex", dest='files_suffix', nargs='?',
                        required=True, help='A regex for the fastq files for the script.',
                        choices=NAME_TYPE_TO_RE.keys())
    parser.add_argument('-of', '--output_files_dir', metavar="output files dir", dest='output_files_dir', nargs='?',
                        required=False,
                        help='Outputer summery dir path, if wanted', default=None)
    parser.add_argument('-sf', '--sites_folder', metavar="overall sites folder", dest='sites_folder', nargs='?',
                        required=False,
                        help='The path of the folder to output to all found sites (output in BED format).',
                        default=DEFAULT_SITES_OUTPUT_FOLDER)
    parser.add_argument('-so', '--STAR_output', metavar="STAR output", dest='STAR_output', nargs='?', required=False,
                        help='Outputer path of STAR run, to normalize reads', default=None)
    parser.add_argument('-srd', '--STAR_recursion_depth', metavar="STAR recursion depth", dest='STAR_recursion_depth',
                        required=False, type=int, default=5,
                        help='The depth to enter recursively in STAR output folder to look for the logs (>=1)')
    # ---- run options ----
    parser.add_argument('-se', '--run_single_end', dest='run_single_end', required=False, action="store_true",
                        help='If set, will run single end analysis (can be used with or without PE)', default=False)
    parser.add_argument('-pe', '--run_paired_end', dest='run_paired_end', required=False, action="store_true",
                        help='If set, will run paired end analysis (can be used with or without SE)', default=False)
    parser.add_argument('--stranded', dest='stranded', required=False, action="store_true", default=False,
                        help='If set, will treat data as stranded data, otherwise will transform all to + strand.')
    parser.add_argument('--keep_bams', dest='keep_bams', required=False, action="store_true", default=False,
                        help='If set, will not delete alignment and transformed BAMs but will instead only sort them.')
    parser.add_argument('-r', '--repeat_on_noise', dest='repeat_on_noise', required=False, action="store_true",
                        help='If` set, will re-run the script if <max_noise_ratio> %% of samples are below'
                             ' <noise_level> ', default=False)
    parser.add_argument('-snr', '--noise_ratio', metavar="noise_ratio", dest='max_noise_ratio', required=False,
                        type=float,
                        default=0.1, help='The ratio of "noisey" samples threshold, if passed will repeat run as PE')
    parser.add_argument('-cl', '--a2g_cleanness_level', metavar="A2G cleanness level", dest='noise_level',
                        required=False,
                        type=float, default=0.80,
                        help='The minimal ratio of A to G mismatches to others to define a sample as "clean"')
    parser.add_argument('-ug', '--upstream_g', metavar="max upstream g", dest='max_pre_g_pc', required=False,
                        type=float, default=0.20, help='The minimal ratio of  G upstream to the mismatch '
                                                       '(for motif check) to define a sample as "clean"')
    parser.add_argument('-pg', '--downstream_g', metavar="min downstream g", dest='min_past_g_pc', required=False,
                        type=float, default=0.25, help='The minimal ratio of  G upstream to the mismatch '
                                                       '(for motif check) to define a sample as "clean"')
    parser.add_argument('--strict_check', dest='strict_check', required=False, action='store_true',
                        help=' If set, will require *both* motif and A2g mismatch ratio to pass, else - either.')

    # ---- miscellaneous ----
    parser.add_argument('-b', '--bash_args', metavar="script bash args", dest='bash_args', required=False, default='',
                        nargs='*',
                        help='named args (in the format of <var>=\'<val>\',<var>=\'<val>\') to pass to the bash running '
                             'the Hyper Editing tool. Possible args:\r\n' +
                             BASH_HELP_STR)

    parser.add_argument('-c', '--calc_output_only', dest='calc_output_only', required=False, action='store_true',
                        help='If set, will only re-caluculte the summery output file.'
                             ' *note* folders are according to given bash args with the run.')

    parser.add_argument('-fn', '--folder_names_only', dest='folder_names_only', required=False, action='store_true',
                        help='If set, will only get summery stats from every and each sample folder instead of the'
                             ' "all" dir. Use when correlation between reads names and samples doesn\'t exist')
    # endregion

    args = parser.parse_args()
    b_args = convert_args_to_dict(args.bash_args)

    init_logging_dict(os.path.join(args.output_dir, date.today().isoformat() + "_" + LOG_FILE))

    run_he_script(input_dir=args.input_dir,
                  output_dir=args.output_dir,
                  files_suffix=args.files_suffix,
                  output_files_dir=args.output_files_dir,
                  sites_folder=args.sites_folder,
                  run_single_end=args.run_single_end,
                  run_paired_end=args.run_paired_end,
                  repeat_on_noise=args.repeat_on_noise,
                  STAR_output=args.STAR_output,
                  STAR_recursion_depth=args.STAR_recursion_depth,
                  stranded=args.stranded,
                  calc_output_only=args.calc_output_only,
                  keep_bams=args.keep_bams,
                  max_pre_g_pc=args.max_pre_g_pc,
                  min_past_g_pc=args.min_past_g_pc,
                  strict_check=args.strict_check,
                  folder_names_only=args.folder_names_only,
                  **b_args)
