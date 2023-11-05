__author__ = 'Hillel'
"""
This module runs editing index.
"""
# =====================imports=====================#
# region Builtin Import
import argparse
import logging
import multiprocessing
import operator
import os
import random
import sys
import threading
import cPickle
import gzip
import re

from collections import OrderedDict
from ConfigParser import ConfigParser
from datetime import datetime
from pprint import pformat

# endregion

# region Internal Imports
if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

import path_locator
from Tools.EditingIndex.EditingIndexConsts import *
from Commons.general_functions import init_logging_dict, convert_args_to_dict, add_get_paths_function_to_argparser, \
    get_paths, add_groups_file_to_argparser, load_groups_and_samples_from_file, get_groups_and_sample_names_dict, \
    remove_files, outof_temp_env
from Commons.consts import ALL, UNKNOWN_STRAND, SENSE_STRAND, ANTISENSE_STRAND, SAMPLE, SAMPLE_PATH, GROUP, NA, \
    POSSIBLE_STRANDS, MismatchesAndRefsEnum, GENOMIC_REGION, SITE_START, SITE_END
from Commons.data_structs import Site, RefSeqPosEnum, SortingHelpFormatter
from Tools.EditingIndex.DataConverters.GeneRegularExpressionCoverter import load_gene_expression_bed
from DataObjects.BaseClasses.Sample import Sample

from Tools.EditingIndex.DataConverters.CountPileupCoverter import parse_count_pileup, get_empty_region_counts

from Tools.EditingIndex.EIPipeline.EIPipelineManger import run_pipeline

from DataObjects.DataRecords.GeneRegularExpressionData import GeneRegularExpressionData
from Tools.EditingIndex.DataConverters.RefSeqGenesCoverter import load_refseq_bed
from Tools.EditingIndex.DataConverters.SNPsDBCoverter import load_snps_file

from Outputs.Outputers.CSVOutputer import CSVOutputer

# endregion

# =====================constants===================#
EDITING_INDEX_SOURCE = "EditingIndexScript"
PICKLE_FILE_FORMAT = "%(resource_file_name)s.%(genome_name)s.%(regions_name)s.Preprocessed.pkl.gz"
COUNTS_KEY = "MismatchesAndCoverageCounts"
STRANDS_KEY = "Strands"

# region Logging Messages
PATH_UNAVAILABLE_ERR = "Could Not Access %s at %s!"
FASTQ_DIR_SUFF_ERR = "Fastqs Dir and Fastq Suffix Must Both Be Absent Or Provided!"
SPECIAL_COUNT_FAIL_ERR = "Failed At Generating Special Count %s for Sample %s"
# endregion

# region Operational and Output Options
STRAND_DECIDING_BY_UNFILTERED_MM_SITES = "MMSitesThenRefSeq"
STRAND_DECIDING_BY_REFSEQ_AND_UNFILTERED_MM_SITES = "RefSeqThenMMSites"
STRAND_DECIDING_BY_REFSEQ_AND_RANDOM = "RefSeqThenRandom"
STRAND_DECIDING_BY_RANDOM = "Randomly"
STRAND_DECIDING_ALWAYS_BOTH = "AlwaysAverageOfBothStrands"
STRAND_DECIDING_STRANDED_DATA = "StrandedData"

ALL_MM = "AllMismatches"
UNSTRANDED_MM = "UnstrandedMismatches"
POSSIBLE_MM = {ALL_MM: MismatchesAndRefsEnum.ALL_MISMATCHES,
               UNSTRANDED_MM: MismatchesAndRefsEnum.UNSTRANDED_MISMATCHES
               }

PER_REGION_STRAND_COUNTS_RE = re.compile(r"NumOf\w{3}MMSitesAvg")


def mismatch_sites_supported_strand(counts, mm_type):
    sense_mm_sites_count = counts[NUM_OF_MM_SITES_FORMAT % mm_type]
    antisense_mm_sites_count = counts[NUM_OF_MM_SITES_FORMAT % MismatchesAndRefsEnum.COMPLEMENTARY_MM[mm_type]]
    has_mm = sense_mm_sites_count + antisense_mm_sites_count > 0

    if has_mm:
        strand = SENSE_STRAND if sense_mm_sites_count >= antisense_mm_sites_count else ANTISENSE_STRAND
    else:
        strand = UNKNOWN_STRAND

    return strand


def strand_deciding_by_mm_sites_first(counts, refseq_strand, mm_type):
    strand = mismatch_sites_supported_strand(counts, mm_type)

    if strand == UNKNOWN_STRAND:
        strand = refseq_strand

    return strand


def strand_deciding_by_refseq_first(counts, refseq_strand, mm_type):
    if refseq_strand != UNKNOWN_STRAND:
        return refseq_strand

    return mismatch_sites_supported_strand(counts, mm_type)


def strand_deciding_by_refseq_then_random(counts, refseq_strand, mm_type):
    if refseq_strand != UNKNOWN_STRAND:
        return refseq_strand

    return strand_deciding_by_random(counts, refseq_strand, mm_type)


def strand_deciding_by_random(counts, refseq_strand, mm_type):
    random_strand_dict = {0: SENSE_STRAND, 1: ANTISENSE_STRAND, 2: UNKNOWN_STRAND}

    return random_strand_dict[random.randint(0, 2)]


def strand_deciding_always_both_strands(counts, refseq_strand, mm_type):
    return UNKNOWN_STRAND


MAIN_STRAND_DECIDING_METHOD = STRAND_DECIDING_BY_REFSEQ_AND_UNFILTERED_MM_SITES

STRAND_DECIDING_OPTIONS = {
    # STRAND_DECIDING_BY_RANDOM: strand_deciding_by_random,
    # STRAND_DECIDING_ALWAYS_BOTH: strand_deciding_always_both_strands,
    STRAND_DECIDING_BY_REFSEQ_AND_UNFILTERED_MM_SITES: strand_deciding_by_refseq_first,
    # STRAND_DECIDING_BY_REFSEQ_AND_RANDOM: strand_deciding_by_refseq_then_random,
    # STRAND_DECIDING_BY_UNFILTERED_MM_SITES: strand_deciding_by_mm_sites_first
}

PER_REGION_OUTPUT_FILE_NAME = "StrandDerivingCountsPerRegion.csv"
EDITING_INDEX_OUTPUT_FILE_NAME = "EditingIndex.csv"

SENSE_POSITION = "SenseGenomicPosition"
ANTISENSE_POSITION = "AntisenseGenomicPosition"
SENSE_NAME = "SenseGeneCommonName"
SENSE_REFSEQ_ID = "SenseGeneRefSeqID"
SENSE_RFPKMS = "SenseGeneExpressionGeometricAvg"
ANTISENSE_NAME = "AntisenseGeneCommonName"
ANTISENSE_REFSEQ_ID = "AntisenseGeneRefSeqID"
ANTISENSE_RFPKMS = "AntisenseGeneExpressionGeometricAvg"


# endregion


# =====================functions===================#
def init_empty_index_dict():
    """
        This function inits an empty index dict.
        :return dict: an initialized index dict.
        """
    empty_index_counts = dict()
    for ref_pos in MismatchesAndRefsEnum.ALL_REFS:
        empty_index_counts[COVERAGE_AT_N_POSITIONS_FORMAT % ref_pos] = 0
        empty_index_counts[NUM_OF_N_SITES_COVERED_FORMAT % ref_pos] = 0

    for mm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
        empty_index_counts[NUM_OF_MM_FORMAT % mm_type] = 0
        empty_index_counts[NUM_OF_MM_SITES_FORMAT % mm_type] = 0
        empty_index_counts[SNPS_NUM_OF_MM_FORMAT % mm_type] = 0
        empty_index_counts[SNPS_NUM_OF_MM_SITES_FORMAT % mm_type] = 0

    indexed_mms = MismatchesAndRefsEnum.ALL_MISMATCHES  # if stranded else MismatchesAndRefsEnum.UNSTRANDED_MISMATCHES

    for mm_type in indexed_mms:
        empty_index_counts[INDEXED_EDITED_FORMAT % mm_type] = 0
        empty_index_counts[INDEXED_CANONICAL_FORMAT % mm_type] = 0
        empty_index_counts[NUM_OF_INDEXED_MM_SITES_FORMAT % mm_type] = 0
        empty_index_counts[NUM_OF_INDEXED_OVERALL_SITES_FORMAT % mm_type] = 0
        empty_index_counts[NUM_OF_REGIONS % mm_type] = 0
        for strand in POSSIBLE_STRANDS:
            empty_index_counts[STRANDED_INDEXED_EDITED_FORMAT % (mm_type, strand)] = 0
            empty_index_counts[STRANDED_INDEXED_CANONICAL_FORMAT % (mm_type, strand)] = 0
            empty_index_counts[STRANDED_NUM_OF_INDEXED_MM_SITES_FORMAT % (mm_type, strand)] = 0
            empty_index_counts[STRANDED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT % (mm_type, strand)] = 0
            empty_index_counts[STRANDED_NUM_OF_REGIONS % (mm_type, strand)] = 0
        for region in RefSeqPosEnum.ALL_ANNOTATIONS:
            empty_index_counts[REGIONED_INDEXED_EDITED_FORMAT % (mm_type, region)] = 0
            empty_index_counts[REGIONED_INDEXED_CANONICAL_FORMAT % (mm_type, region)] = 0
            empty_index_counts[REGIONED_NUM_OF_INDEXED_MM_SITES_FORMAT % (mm_type, region)] = 0
            empty_index_counts[REGIONED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT % (mm_type, region)] = 0
            empty_index_counts[REGIONED_NUM_OF_REGIONS % (mm_type, region)] = 0

    for base in MismatchesAndRefsEnum.ALL_REFS:
        empty_index_counts[NUM_OF_CANONICAL_FORMAT % base] = 0

    return empty_index_counts


def derive_strand_from_refseq(genes_exp, refseqs, region):
    sense_genomic_pos, genomic_region_sense = refseqs[region][SENSE_STRAND]
    antisense_genomic_pos, genomic_region_antisense = refseqs[region][ANTISENSE_STRAND]
    if sense_genomic_pos == antisense_genomic_pos:
        if sense_genomic_pos == RefSeqPosEnum.INTERGENIC:
            derived_strand = UNKNOWN_STRAND
        else:
            sense_rfpkms = genes_exp.get_gavg_rfkpm(region=genomic_region_sense.site.region,
                                                    start=genomic_region_sense.site.start,
                                                    end=genomic_region_sense.site.end,
                                                    strand=genomic_region_sense.site.strand,
                                                    connected_region=genomic_region_sense.site)
            antisense_rfpkms = genes_exp.get_gavg_rfkpm(region=genomic_region_antisense.site.region,
                                                        start=genomic_region_antisense.site.start,
                                                        end=genomic_region_antisense.site.end,
                                                        strand=genomic_region_antisense.site.strand,
                                                        connected_region=genomic_region_antisense.site)
            if sense_rfpkms > antisense_rfpkms:
                derived_strand = SENSE_STRAND
            else:
                if sense_rfpkms > antisense_rfpkms:
                    derived_strand = ANTISENSE_STRAND
                else:
                    derived_strand = UNKNOWN_STRAND
    else:
        stronger_annotation = RefSeqPosEnum.stronger_annotation(sense_genomic_pos, antisense_genomic_pos)
        if stronger_annotation == sense_genomic_pos:
            derived_strand = SENSE_STRAND
        else:
            derived_strand = ANTISENSE_STRAND
    return derived_strand


def derive_strand_from_refseq_p(genes_exp, refseqs, region, sema_decide, refseqs_strands):
    try:
        sense_genomic_pos, genomic_region_sense = refseqs[region][SENSE_STRAND]
        antisense_genomic_pos, genomic_region_antisense = refseqs[region][ANTISENSE_STRAND]
        if sense_genomic_pos == antisense_genomic_pos:
            if sense_genomic_pos == RefSeqPosEnum.INTERGENIC:
                derived_strand = UNKNOWN_STRAND
            else:
                sense_rfpkms = genes_exp.get_gavg_rfkpm(region=genomic_region_sense.site.region,
                                                        start=genomic_region_sense.site.start,
                                                        end=genomic_region_sense.site.end,
                                                        strand=genomic_region_sense.site.strand,
                                                        connected_region=genomic_region_sense.site)
                antisense_rfpkms = genes_exp.get_gavg_rfkpm(region=genomic_region_antisense.site.region,
                                                            start=genomic_region_antisense.site.start,
                                                            end=genomic_region_antisense.site.end,
                                                            strand=genomic_region_antisense.site.strand,
                                                            connected_region=genomic_region_antisense.site)
                if sense_rfpkms > antisense_rfpkms:
                    derived_strand = SENSE_STRAND
                else:
                    if sense_rfpkms > antisense_rfpkms:
                        derived_strand = ANTISENSE_STRAND
                    else:
                        derived_strand = UNKNOWN_STRAND
        else:
            stronger_annotation = RefSeqPosEnum.stronger_annotation(sense_genomic_pos, antisense_genomic_pos)
            if stronger_annotation == sense_genomic_pos:
                derived_strand = SENSE_STRAND
            else:
                derived_strand = ANTISENSE_STRAND

        refseqs_strands[region] = derived_strand
    except:
        logging.exception("Failed to Decide Strand for %s" % region)
    finally:
        sema_decide.release()


def rm_cmpileups(confs, pipeline_output_dir, is_stranded):
    # type: (list, str, bool) -> None
    """
    This function removes all cmpileups and the genome index
    :param list[dict[str,str]] confs: the loaded configs of all samples
    :param str pipeline_output_dir:
    :param bool is_stranded:
    :return: None
    """
    logging.info("Trying to Delete cmpileups")
    to_del_files = list()
    if is_stranded:
        for conf in confs:
            to_del_files.append(conf[REGIONS_PILEUP_COUNT_STRAND_1_OPTION])
            to_del_files.append(conf[REGIONS_PILEUP_COUNT_STRAND_2_OPTION])
    else:
        for conf in confs:
            to_del_files.append(conf[REGIONS_PILEUP_COUNT_OPTION])

    # not cloud opt: >> to_del_files.append(conf[GENOME_INDEX_PATH_OPTION])
    remove_files(set(to_del_files))
    os.system("find %s -type d -empty -delete" % pipeline_output_dir)


def get_index(refseq_strands, refseq_annotations, samples_configs, summery_dir, pipeline_output_dir, snps,
              get_regions_data, per_sample_output, delete_cmpileups, verbose, max_processes_sample,
              max_processes_strand_decision, is_stranded, possible_mismatches):
    # type: (dict, dict, dict, str, str, dict, bool, bool, bool, bool, int, int, bool, list) -> None
    """
    This is the main functional method, that calculates the index.
    :param dict[Site, str] refseq_annotations: The annotations for each region according to the refseq.
    :param dict[Site, str] refseq_strands: The most likely strand of the regions according to the annotations (by
    region).
    :param dict[str, str] samples_configs: A dictionary containing each pof the samples config instances.
    :param str summery_dir: The path of summery output dir (will be created)
    :param str pipeline_output_dir: The path of pipeline (where pileup etc are created) output dir (will be created)
    :param dict snps: A dictionary with the SNPs per mismatch type per position.
    :param bool get_regions_data: A flag. If set will also output the counts per region (huge output, use
    carefully)
    :param bool per_sample_output: A flag. If set will output the counts per region also per sample.
    :param bool delete_cmpileups: A flag. If set will delete the cmpileup after processing it.
    :param bool verbose: A flag. If set, will output a verbose output for the index, including all the counts
    :param int max_processes_sample: The maximal number of samples to process in parallel.
    :param int max_processes_strand_decision:  The maximal strand decisions per sample to process in parallel.
    :param bool is_stranded: A flag. If set, will treat the data as stranded.
    :param list possible_mismatches: what mismatch types to output.
    :return: None
    """

    empty_index_dict = init_empty_index_dict()
    empty_counter = get_empty_region_counts()
    confs = list()
    manager = multiprocessing.Manager()
    ps_pool = multiprocessing.Pool(processes=max_processes_sample)
    samples_groups_dict = get_groups_and_sample_names_dict(by_sample=True)
    samples = Sample.get_all_records(of_types=(Sample,))
    num_of_samples = len(samples)

    if get_regions_data:
        regions_headers = OrderedDict()
        for h in [GENOMIC_REGION, SITE_START, SITE_END, ANTISENSE_POSITION, ANTISENSE_NAME,
                  ANTISENSE_REFSEQ_ID, SENSE_POSITION, SENSE_NAME, SENSE_REFSEQ_ID]:
            regions_headers[h] = True

        get_regions_data_dict = manager.dict()
        threads = list()
        sema_populate = threading.Semaphore(5)
        for region, annotations in refseq_annotations.iteritems():
            tr = threading.Thread(target=populate_get_region_data_d,
                                  args=(annotations, empty_counter, get_regions_data_dict, manager,
                                        region, regions_headers, sema_populate),
                                  name="EditingIndexPopGetRegionDSubprocessRefSeq")
            tr.daemon = True

            threads.append(tr)
        # decide all strands (paralleling)
        for thr in threads:
            sema_populate.acquire()
            thr.start()
        for thr in multiprocessing.active_children():
            thr.join(3600)  # if waited for more than an hour something is wrong....
        logging.debug("Done Initializing Regions Counter For Get Regions Data")

        logging.debug("Added headers:" + ",\n".join(regions_headers.keys()))
    else:
        get_regions_data_dict = None

    for rec_id, sample in samples.iteritems():
        """type sample Sample"""
        sample_name = sample.sample_name
        # get the config dict to retrieve the cmpileup
        try:
            config_defaults = samples_configs[rec_id].defaults(formatted=True)  # get the sample's specific paths
            confs.append(config_defaults)
        except KeyError:
            logging.info("No Config Was Available For %s, Cannot Retrieve cmpileup!" % sample_name)
            continue  # No data for this sample or this is a group

        all_methods_sample_index_dict = {method: empty_index_dict.copy() for method in STRAND_DECIDING_OPTIONS}
        ps_pool.apply_async(func=load_and_process_sample,
                            args=(
                                all_methods_sample_index_dict,
                                config_defaults,
                                get_regions_data_dict,
                                refseq_annotations,
                                refseq_strands,
                                samples_groups_dict,
                                snps,
                                per_sample_output,
                                max_processes_strand_decision,
                                num_of_samples,
                                pipeline_output_dir,
                                summery_dir,
                                sample,
                                verbose,
                                is_stranded,
                                possible_mismatches))

    ps_pool.close()
    ps_pool.join()
    if get_regions_data:
        logging.info("Creating Region Counts Output For All Samples")
        recs = [OrderedDict([(key, str(val)) for key, val in get_regions_data_dict[region].iteritems()]) for region in
                get_regions_data_dict.keys()]
        out_format = os.path.join(summery_dir, PER_REGION_OUTPUT_FILE_NAME)
        csv_outputer = CSVOutputer()
        logging.debug("Writing Data for Headers: %s", ", ".join(recs[0].keys()))
        csv_outputer.output([out_format], regions_headers.keys(), recs)

    if delete_cmpileups:
        rm_cmpileups(confs, pipeline_output_dir, is_stranded)


def populate_get_region_data_d(annotations, empty_counter, get_regions_data_dict, manager, region,
                               regions_headers, sema_populate):
    antisense_position, antisense_refseq = annotations[ANTISENSE_STRAND]
    sense_position, sense_refseq = annotations[SENSE_STRAND]
    reg_dict = manager.dict({
        GENOMIC_REGION: region.region,
        SITE_START: region.start,
        SITE_END: region.end,
        ANTISENSE_POSITION: antisense_position,
        ANTISENSE_NAME: antisense_refseq.common_name,
        ANTISENSE_REFSEQ_ID: antisense_refseq.refseq_name,
        SENSE_POSITION: sense_position,
        SENSE_NAME: sense_refseq.common_name,
        SENSE_REFSEQ_ID: sense_refseq.refseq_name,
    })
    for count_type in sorted(empty_counter.iterkeys()):
        avg_header = PER_REGIONS_OUTPUT_COUNTS_FORMAT % dict(count_type=count_type, operation="Average")
        tot_header = PER_REGIONS_OUTPUT_COUNTS_FORMAT % dict(count_type=count_type, operation="Total")
        reg_dict[avg_header] = 0
        reg_dict[tot_header] = 0
        regions_headers[avg_header] = True
        regions_headers[tot_header] = True

    get_regions_data_dict[region] = reg_dict
    sema_populate.release()


def generate_per_region_output(get_regions_data_dict, num_of_samples, per_sample_output,
                               per_sample_output_dir, sample, regions_d, by_groups_dict):
    # type: (dict, int, bool, str, Sample, dict, dict) -> None
    """
    This function generates the per region and per sample output.
    :param get_regions_data_dict: dictionary with data on each region and for storage of counts
    :param num_of_samples: The total number of samples for calculation.
    :param per_sample_output: a flag, if set will output counts per sample.
    :param per_sample_output_dir: The output dir for per-sample output.
    :param sample: The sample processed.
    :param regions_d: The processed data of the sample (per region).
    :param by_groups_dict: The samples group dict.
    :return: None
    """
    per_sample_regions_recs = list()
    per_sample_dict = OrderedDict()
    csv_outputer = CSVOutputer(append=True, override=False)

    """:type sample Sample"""
    assert isinstance(sample, Sample)
    group = NA
    for group_p in by_groups_dict.get(sample.sample_name).keys():
        if group_p[1] <= 1:
            group = group_p[0]
            break

    sample_name = sample.sample_name

    logging.info("Generating Per Region Data for %s" % sample.sample_name)
    for region, region_data_d in regions_d.iteritems():
        if per_sample_output:
            per_sample_dict[GROUP] = group
            per_sample_dict[SAMPLE] = sample_name
            per_sample_dict[SAMPLE_PATH] = sample.sample_path
            per_sample_dict[GENOMIC_REGION] = region.region
            per_sample_dict[SITE_START] = str(region.start)
            per_sample_dict[SITE_END] = str(region.end)
            per_sample_dict[ANTISENSE_POSITION] = str(get_regions_data_dict[region][ANTISENSE_POSITION])
            per_sample_dict[ANTISENSE_NAME] = get_regions_data_dict[region][ANTISENSE_NAME]
            per_sample_dict[ANTISENSE_REFSEQ_ID] = get_regions_data_dict[region][ANTISENSE_REFSEQ_ID]
            per_sample_dict[SENSE_POSITION] = str(get_regions_data_dict[region][SENSE_POSITION])
            per_sample_dict[SENSE_NAME] = get_regions_data_dict[region][SENSE_NAME]
            per_sample_dict[SENSE_REFSEQ_ID] = get_regions_data_dict[region][SENSE_REFSEQ_ID]

        for count_type, count in region_data_d[COUNTS_KEY].iteritems():
            avg_header = PER_REGIONS_OUTPUT_COUNTS_FORMAT % dict(count_type=count_type, operation="Average")
            tot_header = PER_REGIONS_OUTPUT_COUNTS_FORMAT % dict(count_type=count_type, operation="Total")
            get_regions_data_dict[region].setdefault(avg_header, 0)
            get_regions_data_dict[region].setdefault(tot_header, 0)
            get_regions_data_dict[region][avg_header] += count / float(num_of_samples)
            get_regions_data_dict[region][tot_header] += count
            if per_sample_output:
                per_sample_dict[count_type] = str(count)

        strands_dict = region_data_d[STRANDS_KEY][MAIN_STRAND_DECIDING_METHOD]

        if per_sample_output:
            per_sample_dict[STRAND_DECIDING_METHOD] = MAIN_STRAND_DECIDING_METHOD
            per_sample_dict.update(strands_dict)

            for mm_type, strand in strands_dict.iteritems():
                if strand == SENSE_STRAND:
                    ref_base = mm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

                    stranded_edited = region_data_d[COUNTS_KEY][NUM_OF_MM_FORMAT % mm_type]
                    stranded_edited_sites = region_data_d[COUNTS_KEY][NUM_OF_MM_SITES_FORMAT % mm_type]
                    stranded_overall_sites = region_data_d[COUNTS_KEY][NUM_OF_N_SITES_COVERED_FORMAT % ref_base]
                    stranded_canonical = region_data_d[COUNTS_KEY][NUM_OF_CANONICAL_FORMAT % ref_base]

                elif strand == ANTISENSE_STRAND:
                    stranded_mm_type = MismatchesAndRefsEnum.COMPLEMENTARY_MM[mm_type]
                    ref_base = stranded_mm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

                    stranded_edited = region_data_d[COUNTS_KEY][NUM_OF_MM_FORMAT % stranded_mm_type]
                    stranded_edited_sites = region_data_d[COUNTS_KEY][NUM_OF_MM_SITES_FORMAT % stranded_mm_type]
                    stranded_overall_sites = region_data_d[COUNTS_KEY][
                        NUM_OF_N_SITES_COVERED_FORMAT % ref_base]
                    stranded_canonical = region_data_d[COUNTS_KEY][NUM_OF_CANONICAL_FORMAT % ref_base]
                else:
                    ref_base = mm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

                    stranded_edited = region_data_d[COUNTS_KEY][NUM_OF_MM_FORMAT % mm_type] * 0.5
                    stranded_edited_sites = region_data_d[COUNTS_KEY][NUM_OF_MM_SITES_FORMAT % mm_type] * 0.5
                    stranded_overall_sites = region_data_d[COUNTS_KEY][
                                                 NUM_OF_N_SITES_COVERED_FORMAT % ref_base] * 0.5
                    stranded_canonical = region_data_d[COUNTS_KEY][NUM_OF_CANONICAL_FORMAT % ref_base] * 0.5

                    stranded_mm_type = MismatchesAndRefsEnum.COMPLEMENTARY_MM[mm_type]
                    ref_base = stranded_mm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

                    stranded_edited += region_data_d[COUNTS_KEY][NUM_OF_MM_FORMAT % stranded_mm_type] * 0.5
                    stranded_edited_sites += region_data_d[COUNTS_KEY][
                                                 NUM_OF_MM_SITES_FORMAT % stranded_mm_type] * 0.5
                    stranded_overall_sites += region_data_d[COUNTS_KEY][
                                                  NUM_OF_N_SITES_COVERED_FORMAT % ref_base] * 0.5
                    stranded_canonical += region_data_d[COUNTS_KEY][NUM_OF_CANONICAL_FORMAT % ref_base] * 0.5

                per_sample_dict[INDEXED_EDITED_FORMAT % mm_type] = str(stranded_edited)

                per_sample_dict[INDEXED_CANONICAL_FORMAT % mm_type] = str(stranded_canonical)

                per_sample_dict[NUM_OF_INDEXED_MM_SITES_FORMAT % mm_type] = str(stranded_edited_sites)
                per_sample_dict[NUM_OF_INDEXED_OVERALL_SITES_FORMAT % mm_type] = str(stranded_overall_sites)

                try:
                    per_sample_dict[EDITING_INDEX_FORMAT % mm_type] = str(
                        100 * float(stranded_edited) / (stranded_edited + stranded_canonical))
                except ZeroDivisionError:
                    per_sample_dict[EDITING_INDEX_FORMAT % mm_type] = "0.0"

            per_sample_regions_recs.append(per_sample_dict.copy())

    logging.info("Done Summing Per Sample Regions Data for %s" % sample.sample_name)

    if per_sample_output:
        logging.info("Creating Per Sample Region Counts Output For %s" % sample.sample_name)
        out_format = os.path.join(per_sample_output_dir, PER_REGION_OUTPUT_FILE_NAME)
        csv_outputer.output([out_format, ], per_sample_regions_recs[0].keys(),
                            per_sample_regions_recs)


def generate_per_region_output_stranded(get_regions_data_dict, num_of_samples, per_sample_output, per_sample_output_dir,
                                        sample, refseq_annotations, regions_counts, by_groups_dict, split_by_strand):
    # type: (dict, int, bool, str, Sample, dict, dict, dict, bool) -> None
    """
    This methods  generates the per region and per sample output for stranded data
    :param get_regions_data_dict: dictionary with data on each region and for storage of counts
    :param num_of_samples: The total number of samples for calculation.
    :param per_sample_output: a flag, if set will output counts per sample.
    :param per_sample_output_dir: The output dir for per-sample output.
    :param sample: The sample processed.
    :param regions_counts: The processed data of the sample (per region).
    :param by_groups_dict: The samples group dict.
    :param split_by_strand: If set, will divide the counts according to the strand of origin
    :return: None
    """
    region_counts_trans_key = get_complement_region_counts()

    csv_outputer = CSVOutputer(append=True, override=False)

    """:type sample Sample"""
    assert isinstance(sample, Sample)
    group = NA
    for group_p in by_groups_dict.get(sample.sample_name).keys():
        if group_p[1] <= 1:
            group = group_p[0]
            break

    sample_name = sample.sample_name

    logging.info("Summing Per Sample Regions Data for %s" % sample.sample_name)
    per_sample_all_regions_dict = OrderedDict()
    per_sample_region_data_dict = OrderedDict()
    per_sample_regions_output_recs = dict()

    empty_index_dict = init_empty_index_dict()
    sample_index_dict = dict()
    for region, region_data_d in regions_counts[SENSE_STRAND].iteritems():
        sample_index_dict[region] = empty_index_dict.copy()
        per_sample_region_data_dict = OrderedDict()
        if per_sample_output:
            init_per_region_region_data_d(get_regions_data_dict, group, per_sample_region_data_dict, region, sample,
                                          sample_name)
            # if split_by_strand:
            # per_sample_region_data_dict[STRAND_OF_ORIGIN] = SENSE_STRAND

        per_region_count_dict = get_regions_data_dict[region]  # to overcome python2.7 multiprocessing nested proxy bug
        for count_type, count in region_data_d.iteritems():
            avg_header = PER_REGIONS_OUTPUT_COUNTS_FORMAT % dict(count_type=count_type, operation="Average")
            tot_header = PER_REGIONS_OUTPUT_COUNTS_FORMAT % dict(count_type=count_type, operation="Total")
            per_region_count_dict[avg_header] += count / float(num_of_samples)
            per_region_count_dict[tot_header] += count

            if per_sample_output:
                per_sample_region_data_dict[count_type] = count

        get_regions_data_dict[region] = per_region_count_dict  # to overcome python2.7 multiprocessing nested proxy bug

        annotations = refseq_annotations[region]
        antisense_position, _ = annotations[ANTISENSE_STRAND]
        sense_position, _ = annotations[SENSE_STRAND]

        if per_sample_output:
            per_sample_regions_output_recs[region] = per_sample_region_data_dict.copy()
            for count_type, count_val in region_data_d.iteritems():
                sample_index_dict[region][count_type] = count_val

            strand = SENSE_STRAND
            genomic_position = sense_position
            index_counts_stranded_data(genomic_position, region_data_d, sample_index_dict[region], strand)

    for region, region_data_d in regions_counts[ANTISENSE_STRAND].iteritems():
        if per_sample_output and region not in per_sample_regions_output_recs:
            sample_index_dict[region] = empty_index_dict.copy()
            per_sample_region_data_dict = OrderedDict()
            init_per_region_region_data_d(get_regions_data_dict, group, per_sample_region_data_dict, region, sample,
                                          sample_name)

            # per_sample_region_data_dict = per_sample_all_regions_dict.copy()

            per_sample_regions_output_recs[region] = per_sample_region_data_dict.copy()
            if split_by_strand:
                per_sample_region_data_dict[STRAND_OF_ORIGIN] = SENSE_STRAND
        per_region_count_dict = get_regions_data_dict[region]  # to overcome python2.7 multiprocessing nested proxy bug
        for count_type, count in region_data_d.iteritems():
            sense_count_type = region_counts_trans_key[count_type]
            avg_header = PER_REGIONS_OUTPUT_COUNTS_FORMAT % dict(count_type=sense_count_type, operation="Average")
            tot_header = PER_REGIONS_OUTPUT_COUNTS_FORMAT % dict(count_type=sense_count_type, operation="Total")
            per_region_count_dict[avg_header] += count / float(num_of_samples)
            per_region_count_dict[tot_header] += count
            if per_sample_output:
                if region not in regions_counts[SENSE_STRAND]:
                    per_sample_regions_output_recs[region].update(get_empty_region_counts())
                per_sample_regions_output_recs[region][sense_count_type] += count
        get_regions_data_dict[region] = per_region_count_dict
        if per_sample_output:
            for ref_pos in MismatchesAndRefsEnum.ALL_REFS:
                pos_ref = MismatchesAndRefsEnum.COMPLEMENTARY_REF[ref_pos]
                sample_index_dict[region][COVERAGE_AT_N_POSITIONS_FORMAT % pos_ref] += region_data_d[
                    COVERAGE_AT_N_POSITIONS_FORMAT % pos_ref]
                sample_index_dict[region][NUM_OF_N_SITES_COVERED_FORMAT % pos_ref] += region_data_d[
                    NUM_OF_N_SITES_COVERED_FORMAT % pos_ref]
                sample_index_dict[region][NUM_OF_CANONICAL_FORMAT % pos_ref] += region_data_d[
                    NUM_OF_CANONICAL_FORMAT % pos_ref]
            for mm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
                pos_mm_type = MismatchesAndRefsEnum.COMPLEMENTARY_MM[mm_type]
                sample_index_dict[region][NUM_OF_MM_FORMAT % pos_mm_type] += region_data_d[
                    NUM_OF_MM_FORMAT % pos_mm_type]
                sample_index_dict[region][NUM_OF_MM_SITES_FORMAT % pos_mm_type] += region_data_d[
                    NUM_OF_MM_SITES_FORMAT % pos_mm_type]
                sample_index_dict[region][SNPS_NUM_OF_MM_FORMAT % pos_mm_type] += region_data_d[
                    SNPS_NUM_OF_MM_FORMAT % pos_mm_type]
                sample_index_dict[region][SNPS_NUM_OF_MM_SITES_FORMAT % pos_mm_type] += region_data_d[
                    SNPS_NUM_OF_MM_SITES_FORMAT % pos_mm_type]

            annotations = refseq_annotations[region]
            antisense_position, _ = annotations[ANTISENSE_STRAND]

            strand = ANTISENSE_STRAND
            genomic_position = antisense_position
            index_counts_stranded_data(genomic_position, region_data_d, sample_index_dict[region], strand)

    for region in sample_index_dict:
        for mm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
            try:
                per_sample_regions_output_recs[region][EDITING_INDEX_FORMAT % mm_type] = str(
                    100 * float(sample_index_dict[region][
                                    INDEXED_EDITED_FORMAT % mm_type]) / (
                            sample_index_dict[region][
                                INDEXED_EDITED_FORMAT % mm_type] +
                            sample_index_dict[region][
                                INDEXED_CANONICAL_FORMAT % mm_type]))
            except ZeroDivisionError:
                per_sample_regions_output_recs[region][EDITING_INDEX_FORMAT % mm_type] = "0.0"
    logging.info("Done Summing Per Sample Regions Data for %s" % sample.sample_name)

    if per_sample_output:
        per_sample_regions_output_recs = [OrderedDict([(key, str(val)) for key, val in ps_rec.iteritems()])
                                          for ps_rec in per_sample_regions_output_recs.values()]
        logging.info("Creating Per Sample Region Counts Output For %s" % sample.sample_name)
        out_format = os.path.join(per_sample_output_dir, PER_REGION_OUTPUT_FILE_NAME)
        csv_outputer.output([out_format, ], per_sample_regions_output_recs[0].keys(),
                            per_sample_regions_output_recs)


def init_per_region_region_data_d(get_regions_data_dict, group, per_sample_region_data_dict, region, sample,
                                  sample_name):
    per_sample_region_data_dict[GROUP] = group
    per_sample_region_data_dict[SAMPLE] = sample_name
    per_sample_region_data_dict[SAMPLE_PATH] = sample.sample_path
    per_sample_region_data_dict[STRAND_DECIDING_METHOD] = STRAND_DECIDING_STRANDED_DATA
    per_sample_region_data_dict[GENOMIC_REGION] = region.region
    per_sample_region_data_dict[SITE_START] = str(region.start)
    per_sample_region_data_dict[SITE_END] = str(region.end)
    per_sample_region_data_dict[ANTISENSE_POSITION] = str(get_regions_data_dict[region][ANTISENSE_POSITION])
    per_sample_region_data_dict[ANTISENSE_NAME] = get_regions_data_dict[region][ANTISENSE_NAME]
    per_sample_region_data_dict[ANTISENSE_REFSEQ_ID] = get_regions_data_dict[region][ANTISENSE_REFSEQ_ID]
    per_sample_region_data_dict[SENSE_POSITION] = str(get_regions_data_dict[region][SENSE_POSITION])
    per_sample_region_data_dict[SENSE_NAME] = get_regions_data_dict[region][SENSE_NAME]
    per_sample_region_data_dict[SENSE_REFSEQ_ID] = get_regions_data_dict[region][SENSE_REFSEQ_ID]


def output_sample_indexes(summery_dir, sample, refseq_annotations, regions_d, by_groups_dict,
                          all_methods_sample_index_dict, verbose, possible_mismatches):
    # type: (str, Sample, dict, dict, dict, dict, bool, list) -> None
    """
    :param summery_dir: The dir to create summery outputs in.
    :param get_regions_data_dict: dictionary with data on each region and for storage of counts
    :param num_of_samples: The total number of samples for calculation.
    :param per_sample_output: a flag, if set will output counts per sample.
    :param per_sample_output_dir: The output dir for per-sample output.
    :param sample: The sample processed.
    :param refseq_annotations: A dictionary with the regions annotations.
    :param regions_d: The processed data of the sample (per region).
    :param by_groups_dict: The samples group dict.
    :param all_methods_sample_index_dict: a counter dict, sent as param to save run time
    :param bool verbose: A flag. If set, will output a verbose output for the index, including all the counts
    :param list possible_mismatches: what mismatch types to output.
    :return: None
    """
    overall_index_recs = list()
    per_sample_regions_dict = OrderedDict()
    csv_outputer = CSVOutputer(append=True, override=False)

    """:type sample Sample"""
    assert isinstance(sample, Sample)
    group = NA
    for group_p in by_groups_dict.get(sample.sample_name).keys():
        if group_p[1] <= 1:
            group = group_p[0]
            break

    sample_name = sample.sample_name
    per_sample_regions_dict[GROUP] = group
    per_sample_regions_dict[SAMPLE] = sample_name
    per_sample_regions_dict[SAMPLE_PATH] = sample.sample_path

    logging.info("Summing Data for %s" % sample.sample_name)
    for region, region_data_d in regions_d.iteritems():
        for strand_deciding_method, strands_dict in region_data_d[STRANDS_KEY].iteritems():
            if strand_deciding_method != MAIN_STRAND_DECIDING_METHOD and not verbose:
                continue

            sample_index_dict = all_methods_sample_index_dict[strand_deciding_method]

            for count_type, count_val in region_data_d[COUNTS_KEY].iteritems():
                sample_index_dict[count_type] += count_val

            annotations = refseq_annotations[region]
            antisense_position, _ = annotations[ANTISENSE_STRAND]
            sense_position, _ = annotations[SENSE_STRAND]

            for mm_type, strand in strands_dict.iteritems():
                if strand == SENSE_STRAND:
                    genomic_position = sense_position
                    ref_base = mm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

                    stranded_edited = region_data_d[COUNTS_KEY][NUM_OF_MM_FORMAT % mm_type]
                    stranded_edited_sites = region_data_d[COUNTS_KEY][NUM_OF_MM_SITES_FORMAT % mm_type]
                    stranded_overall_sites = region_data_d[COUNTS_KEY][NUM_OF_N_SITES_COVERED_FORMAT % ref_base]
                    stranded_canonical = region_data_d[COUNTS_KEY][NUM_OF_CANONICAL_FORMAT % ref_base]
                elif strand == ANTISENSE_STRAND:
                    genomic_position = antisense_position
                    stranded_mm_type = MismatchesAndRefsEnum.COMPLEMENTARY_MM[mm_type]
                    ref_base = stranded_mm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

                    stranded_edited = region_data_d[COUNTS_KEY][NUM_OF_MM_FORMAT % stranded_mm_type]
                    stranded_edited_sites = region_data_d[COUNTS_KEY][NUM_OF_MM_SITES_FORMAT % stranded_mm_type]
                    stranded_overall_sites = region_data_d[COUNTS_KEY][
                        NUM_OF_N_SITES_COVERED_FORMAT % ref_base]
                    stranded_canonical = region_data_d[COUNTS_KEY][NUM_OF_CANONICAL_FORMAT % ref_base]
                else:
                    genomic_position = RefSeqPosEnum.INTERGENIC
                    ref_base = mm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

                    stranded_edited = region_data_d[COUNTS_KEY][NUM_OF_MM_FORMAT % mm_type] * 0.5
                    stranded_edited_sites = region_data_d[COUNTS_KEY][NUM_OF_MM_SITES_FORMAT % mm_type] * 0.5
                    stranded_overall_sites = region_data_d[COUNTS_KEY][
                                                 NUM_OF_N_SITES_COVERED_FORMAT % ref_base] * 0.5
                    stranded_canonical = region_data_d[COUNTS_KEY][NUM_OF_CANONICAL_FORMAT % ref_base] * 0.5

                    stranded_mm_type = MismatchesAndRefsEnum.COMPLEMENTARY_MM[mm_type]
                    ref_base = stranded_mm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

                    stranded_edited += region_data_d[COUNTS_KEY][NUM_OF_MM_FORMAT % stranded_mm_type] * 0.5
                    stranded_edited_sites += region_data_d[COUNTS_KEY][
                                                 NUM_OF_MM_SITES_FORMAT % stranded_mm_type] * 0.5
                    stranded_overall_sites += region_data_d[COUNTS_KEY][
                                                  NUM_OF_N_SITES_COVERED_FORMAT % ref_base] * 0.5
                    stranded_canonical += region_data_d[COUNTS_KEY][NUM_OF_CANONICAL_FORMAT % ref_base] * 0.5

                sample_index_dict[INDEXED_EDITED_FORMAT % mm_type] += stranded_edited
                sample_index_dict[STRANDED_INDEXED_EDITED_FORMAT % (mm_type, strand)] += stranded_edited
                sample_index_dict[REGIONED_INDEXED_EDITED_FORMAT % (mm_type, genomic_position)] += stranded_edited

                sample_index_dict[INDEXED_CANONICAL_FORMAT % mm_type] += stranded_canonical
                sample_index_dict[STRANDED_INDEXED_CANONICAL_FORMAT % (mm_type, strand)] += stranded_canonical
                sample_index_dict[
                    REGIONED_INDEXED_CANONICAL_FORMAT % (mm_type, genomic_position)] += stranded_canonical

                sample_index_dict[NUM_OF_INDEXED_MM_SITES_FORMAT % mm_type] += stranded_edited_sites
                sample_index_dict[
                    STRANDED_NUM_OF_INDEXED_MM_SITES_FORMAT % (mm_type, strand)] += stranded_edited_sites
                sample_index_dict[
                    REGIONED_NUM_OF_INDEXED_MM_SITES_FORMAT % (mm_type, genomic_position)] += stranded_edited_sites

                sample_index_dict[NUM_OF_INDEXED_OVERALL_SITES_FORMAT % mm_type] += stranded_overall_sites
                sample_index_dict[
                    STRANDED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT % (mm_type, strand)] += stranded_overall_sites
                sample_index_dict[
                    REGIONED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT % (
                        mm_type, genomic_position)] += stranded_overall_sites

                sample_index_dict[NUM_OF_REGIONS % mm_type] += 1
                sample_index_dict[STRANDED_NUM_OF_REGIONS % (mm_type, strand)] += 1
                sample_index_dict[REGIONED_NUM_OF_REGIONS % (mm_type, genomic_position)] += 1

    logging.info("Done Summing Data for %s" % sample.sample_name)

    for strand_deciding_method in STRAND_DECIDING_OPTIONS:
        if strand_deciding_method != MAIN_STRAND_DECIDING_METHOD and not verbose:
            continue
        sample_rec = OrderedDict()
        sample_rec[STRAND_DECIDING_METHOD] = strand_deciding_method

        sample_index_dict = all_methods_sample_index_dict[strand_deciding_method]

        sample_rec[GROUP] = group
        sample_rec[SAMPLE] = sample.sample_name
        sample_rec[SAMPLE_PATH] = sample.sample_path

        for mm_type in possible_mismatches:  # .UNSTRANDED_MISMATCHES:
            try:
                sample_rec[EDITING_INDEX_FORMAT % mm_type] = str(100 * float(sample_index_dict[
                                                                                 INDEXED_EDITED_FORMAT % mm_type]) / (
                                                                         sample_index_dict[
                                                                             INDEXED_EDITED_FORMAT % mm_type] +
                                                                         sample_index_dict[
                                                                             INDEXED_CANONICAL_FORMAT % mm_type]))
            except ZeroDivisionError:
                sample_rec[EDITING_INDEX_FORMAT % mm_type] = "0.0"

        if verbose:
            # cast to str all of the values
            for key in sorted(sample_index_dict):
                sample_rec[key] = str(sample_index_dict[key])

        overall_index_recs.append(sample_rec)

    logging.info("Writing Into Index Output For All Samples, For %s" % sample_name)
    overall_index_headers = overall_index_recs[0].keys()
    csv_outputer.output([os.path.join(summery_dir, EDITING_INDEX_OUTPUT_FILE_NAME), ], overall_index_headers,
                        overall_index_recs)


def get_complement_region_counts():
    # type: () -> dict
    """
    This function inits an transform dict for the antisense counts
    :return dict: the transform mapping
    """
    region_counts_trans_key = dict()
    for ref_pos in MismatchesAndRefsEnum.ALL_REFS:
        comp_ref = MismatchesAndRefsEnum.COMPLEMENTARY_REF[ref_pos]
        region_counts_trans_key[COVERAGE_AT_N_POSITIONS_FORMAT % ref_pos] = COVERAGE_AT_N_POSITIONS_FORMAT % comp_ref
        region_counts_trans_key[NUM_OF_N_SITES_COVERED_FORMAT % ref_pos] = NUM_OF_N_SITES_COVERED_FORMAT % comp_ref
        region_counts_trans_key[NUM_OF_CANONICAL_FORMAT % ref_pos] = NUM_OF_CANONICAL_FORMAT % comp_ref
    for mm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
        comp_mm = MismatchesAndRefsEnum.COMPLEMENTARY_MM[mm_type]
        region_counts_trans_key[NUM_OF_MM_FORMAT % mm_type] = NUM_OF_MM_FORMAT % comp_mm
        region_counts_trans_key[NUM_OF_MM_SITES_FORMAT % mm_type] = NUM_OF_MM_SITES_FORMAT % comp_mm
        region_counts_trans_key[SNPS_NUM_OF_MM_FORMAT % mm_type] = SNPS_NUM_OF_MM_FORMAT % comp_mm
        region_counts_trans_key[SNPS_NUM_OF_MM_SITES_FORMAT % mm_type] = SNPS_NUM_OF_MM_SITES_FORMAT % comp_mm

    return region_counts_trans_key


def output_sample_stranded(summery_dir, sample, refseq_annotations, regions_counts,
                           by_groups_dict, verbose, possible_mismatches):
    # type: (str, Sample, dict, dict, dict, bool, list) -> None
    """
    This methods calc the output for stranded data
    :param summery_dir: The dir to create summery outputs in.
    :param sample: The sample processed.
    :param refseq_annotations: per region refseq annotation
    :param regions_counts: The processed data of the sample (per region).
    :param by_groups_dict: The samples group dict.
    :param bool verbose: A flag. If set, will output a verbose output for the index, including all the counts
    :param list possible_mismatches: what mismatch types to output.
    :return: None
    """

    per_sample_regions_dict = OrderedDict()

    csv_outputer = CSVOutputer(append=True, override=False)

    """:type sample Sample"""
    assert isinstance(sample, Sample)
    group = NA
    for group_p in by_groups_dict.get(sample.sample_name).keys():
        if group_p[1] <= 1:
            group = group_p[0]
            break

    sample_name = sample.sample_name
    per_sample_regions_dict[GROUP] = group
    per_sample_regions_dict[SAMPLE] = sample_name
    per_sample_regions_dict[SAMPLE_PATH] = sample.sample_path

    sample_index_dict = init_empty_index_dict()

    logging.info("Summing Data for %s" % sample.sample_name)
    for region, region_data_d in regions_counts[SENSE_STRAND].iteritems():
        for count_type, count_val in region_data_d.iteritems():
            sample_index_dict[count_type] += count_val

        annotations = refseq_annotations[region]
        antisense_position, _ = annotations[ANTISENSE_STRAND]
        sense_position, _ = annotations[SENSE_STRAND]

        strand = SENSE_STRAND
        genomic_position = sense_position
        index_counts_stranded_data(genomic_position, region_data_d, sample_index_dict, strand)

    for region, region_data_d in regions_counts[ANTISENSE_STRAND].iteritems():
        for ref_pos in MismatchesAndRefsEnum.ALL_REFS:
            pos_ref = MismatchesAndRefsEnum.COMPLEMENTARY_REF[ref_pos]
            sample_index_dict[COVERAGE_AT_N_POSITIONS_FORMAT % pos_ref] += region_data_d[
                COVERAGE_AT_N_POSITIONS_FORMAT % pos_ref]
            sample_index_dict[NUM_OF_N_SITES_COVERED_FORMAT % pos_ref] += region_data_d[
                NUM_OF_N_SITES_COVERED_FORMAT % pos_ref]
            sample_index_dict[NUM_OF_CANONICAL_FORMAT % pos_ref] += region_data_d[NUM_OF_CANONICAL_FORMAT % pos_ref]
        for mm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
            pos_mm_type = MismatchesAndRefsEnum.COMPLEMENTARY_MM[mm_type]
            sample_index_dict[NUM_OF_MM_FORMAT % pos_mm_type] += region_data_d[NUM_OF_MM_FORMAT % pos_mm_type]
            sample_index_dict[NUM_OF_MM_SITES_FORMAT % pos_mm_type] += region_data_d[
                NUM_OF_MM_SITES_FORMAT % pos_mm_type]
            sample_index_dict[SNPS_NUM_OF_MM_FORMAT % pos_mm_type] += region_data_d[SNPS_NUM_OF_MM_FORMAT % pos_mm_type]
            sample_index_dict[SNPS_NUM_OF_MM_SITES_FORMAT % pos_mm_type] += region_data_d[
                SNPS_NUM_OF_MM_SITES_FORMAT % pos_mm_type]

        annotations = refseq_annotations[region]
        antisense_position, _ = annotations[ANTISENSE_STRAND]

        strand = ANTISENSE_STRAND
        genomic_position = antisense_position
        index_counts_stranded_data(genomic_position, region_data_d, sample_index_dict, strand)

    logging.info("Done Summing Data for %s" % sample.sample_name)

    sample_rec = OrderedDict()
    sample_rec[STRAND_DECIDING_METHOD] = STRAND_DECIDING_STRANDED_DATA
    sample_rec[GROUP] = group
    sample_rec[SAMPLE] = sample.sample_name
    sample_rec[SAMPLE_PATH] = sample.sample_path

    for mm_type in possible_mismatches:
        try:
            sample_rec[EDITING_INDEX_FORMAT % mm_type] = str(100 * float(sample_index_dict[
                                                                             INDEXED_EDITED_FORMAT % mm_type]) / (
                                                                     sample_index_dict[
                                                                         INDEXED_EDITED_FORMAT % mm_type] +
                                                                     sample_index_dict[
                                                                         INDEXED_CANONICAL_FORMAT % mm_type]))
        except ZeroDivisionError:
            sample_rec[EDITING_INDEX_FORMAT % mm_type] = "0.0"

    if verbose:
        # cast to str all of the values
        for key in sorted(sample_index_dict):
            sample_rec[key] = str(sample_index_dict[key])

    logging.info("Writing Into Index Output For All Samples, For %s" % sample_name)
    overall_index_headers = sample_rec.keys()
    csv_outputer.output([os.path.join(summery_dir, EDITING_INDEX_OUTPUT_FILE_NAME), ], overall_index_headers,
                        [sample_rec, ])


def index_counts_stranded_data(genomic_position, region_data_d, sample_index_dict, strand):
    for imm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
        mm_type = imm_type if strand == SENSE_STRAND else MismatchesAndRefsEnum.COMPLEMENTARY_MM[imm_type]
        ref_base = imm_type.split(MismatchesAndRefsEnum.MM_SEP)[0]

        stranded_edited = region_data_d[NUM_OF_MM_FORMAT % imm_type]
        stranded_edited_sites = region_data_d[NUM_OF_MM_SITES_FORMAT % imm_type]
        stranded_overall_sites = region_data_d[NUM_OF_N_SITES_COVERED_FORMAT % ref_base]
        stranded_canonical = region_data_d[NUM_OF_CANONICAL_FORMAT % ref_base]

        sample_index_dict[INDEXED_EDITED_FORMAT % mm_type] += stranded_edited
        sample_index_dict[STRANDED_INDEXED_EDITED_FORMAT % (mm_type, strand)] += stranded_edited
        sample_index_dict[REGIONED_INDEXED_EDITED_FORMAT % (mm_type, genomic_position)] += stranded_edited

        sample_index_dict[INDEXED_CANONICAL_FORMAT % mm_type] += stranded_canonical
        sample_index_dict[STRANDED_INDEXED_CANONICAL_FORMAT % (mm_type, strand)] += stranded_canonical
        sample_index_dict[
            REGIONED_INDEXED_CANONICAL_FORMAT % (mm_type, genomic_position)] += stranded_canonical

        sample_index_dict[NUM_OF_INDEXED_MM_SITES_FORMAT % mm_type] += stranded_edited_sites
        sample_index_dict[
            STRANDED_NUM_OF_INDEXED_MM_SITES_FORMAT % (mm_type, strand)] += stranded_edited_sites
        sample_index_dict[
            REGIONED_NUM_OF_INDEXED_MM_SITES_FORMAT % (mm_type, genomic_position)] += stranded_edited_sites

        sample_index_dict[NUM_OF_INDEXED_OVERALL_SITES_FORMAT % mm_type] += stranded_overall_sites
        sample_index_dict[
            STRANDED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT % (mm_type, strand)] += stranded_overall_sites
        sample_index_dict[
            REGIONED_NUM_OF_INDEXED_OVERALL_SITES_FORMAT % (
                mm_type, genomic_position)] += stranded_overall_sites

        sample_index_dict[NUM_OF_REGIONS % mm_type] += 1
        sample_index_dict[STRANDED_NUM_OF_REGIONS % (mm_type, strand)] += 1
        sample_index_dict[REGIONED_NUM_OF_REGIONS % (mm_type, genomic_position)] += 1


def load_and_process_sample(all_methods_sample_index_dict, config_defaults, get_regions_data_dict, refseq_annotations,
                            refseq_strands, samples_groups_dict, snps, per_sample_output, max_processes, num_of_samples,
                            pipeline_output_dir, summery_dir, sample, verbose, is_stranded, possible_mismatches):
    # type: (dict, dict, dict, dict, dict, dict, dict, bool, int, int, str, str, Sample, bool, bool, list) -> None
    """
    This function parses each sample's cmpileup data and decide the strands for each mm type per region and creates the
        output.
    :param dict[str, dict[str, float]] all_methods_sample_index_dict: a counter dict, sent as param to save run time
    :param dict[str, str] config_defaults: The dict containing the defaults part of the sample's configuration.
    :param dict[str, object] get_regions_data_dict: The dict containing the overall regions data.
    :param dict[Site,dict[str, RefSeqPosEnum]] refseq_annotations: A dictionary of all possible refeseq annotations per region.
    :param dict[Site, str] refseq_strands: The most likely strand of the regions according to the annotations (by region).
    :param dict[str, dict] samples_groups_dict: A dictionary of all "known" samples and groups.
    :param dict[str, dict[Site, bool]] snps:  A dictionary with the SNPs per mismatch type per position.
    :param bool per_sample_output:
    :param int max_processes: The maximal number of threads for strands decision
    :param int num_of_samples: The overall number of samples, used in per region data.
    :param str pipeline_output_dir: The base output dir of the pipeline part.
    :param str summery_dir: The base dir for summery files.
    :param Sample sample: The sample instance the data is relevant to.
    :param bool verbose: A flag. If set, will output a verbose output for the index, including all the counts
    :param bool is_stranded: A flag. If set, will treat the data as stranded.
    :param list possible_mismatches: what mismatch types to output.
    :return: None
    """

    if per_sample_output:
        per_sample_output_dir = config_defaults["output_dir"].replace(pipeline_output_dir, summery_dir)
    else:
        per_sample_output_dir = None

    try:
        if is_stranded:
            # get cmpileup path
            count_pileup_file_s_1 = config_defaults[REGIONS_PILEUP_COUNT_STRAND_1_OPTION]
            count_pileup_file_s_2 = config_defaults[REGIONS_PILEUP_COUNT_STRAND_2_OPTION]
            logging.info("Processing The Data Of %s as Stranded" % sample.sample_name)
            regions_counts_by_strand = get_stranded_strands_and_counts(refseq_strands=refseq_strands,
                                                                       sample=sample,
                                                                       snps=snps,
                                                                       count_pileup_file_s_1=count_pileup_file_s_1,
                                                                       count_pileup_file_s_2=count_pileup_file_s_2,
                                                                       )

            output_sample_stranded(summery_dir=summery_dir,
                                   sample=sample,
                                   refseq_annotations=refseq_annotations,
                                   regions_counts=regions_counts_by_strand,
                                   by_groups_dict=samples_groups_dict,
                                   verbose=verbose,
                                   possible_mismatches=possible_mismatches)
            if get_regions_data_dict is not None:
                generate_per_region_output_stranded(get_regions_data_dict=get_regions_data_dict,
                                                    num_of_samples=num_of_samples,
                                                    per_sample_output=per_sample_output,
                                                    per_sample_output_dir=per_sample_output_dir,
                                                    sample=sample,
                                                    refseq_annotations=refseq_annotations,
                                                    regions_counts=regions_counts_by_strand,
                                                    by_groups_dict=samples_groups_dict,
                                                    split_by_strand=False)
        else:
            # get cmpileup path
            count_pileup_file = config_defaults[REGIONS_PILEUP_COUNT_OPTION]
            logging.info("Processing The Data Of %s as Unstranded" % sample.sample_name)

            regions_strands = get_strands_and_counts(refseq_strands=refseq_strands,
                                                     sample=sample,
                                                     snps=snps,
                                                     max_processes=max_processes,
                                                     count_pileup_file=count_pileup_file,
                                                     verbose=verbose)

            output_sample_indexes(summery_dir=summery_dir,
                                  sample=sample,
                                  refseq_annotations=refseq_annotations,
                                  regions_d=regions_strands,
                                  by_groups_dict=samples_groups_dict,
                                  all_methods_sample_index_dict=all_methods_sample_index_dict,
                                  verbose=verbose,
                                  possible_mismatches=possible_mismatches)
            if get_regions_data_dict is not None:
                generate_per_region_output(get_regions_data_dict=get_regions_data_dict,
                                           num_of_samples=num_of_samples,
                                           per_sample_output=per_sample_output,
                                           per_sample_output_dir=per_sample_output_dir,
                                           sample=sample,
                                           regions_d=regions_strands,
                                           by_groups_dict=samples_groups_dict)
    except Exception, e:
        logging.exception("Failed Processing of %s!" % sample.sample_name)
    finally:
        if is_stranded:
            del regions_counts_by_strand
        else:
            del regions_strands


def get_strands_and_counts(refseq_strands, sample, snps, max_processes, count_pileup_file, verbose):
    # type: (dict, Sample, dict, int, str, bool) -> dict
    """
    This function parallels (using threading) the processing of the sample's data (counting coverages and mismatches, deciding strands per region)
    :param dict[Site, str] refseq_strands: The most likely strand of the regions according to the annotations (by region).
    :param Sample sample: The sample we're processing now.
    :param dict snps: A dictionary with the SNPs per mismatch type per position.
    :param int max_processes: The maximal number of threads for strands decision
    :param str count_pileup_file: The path to the current sample cmpileup to process
    :param bool verbose: A flag. If set, will output a verbose output for the index, including all the counts
    :return dict: The regions counts and strands
    """
    sample_name = sample.sample_name
    sema_decide = threading.Semaphore(max_processes)

    # process cmpileup
    logging.info("Started Processing The Data Of %s" % sample_name)
    logging.info("Loading Coverage File %s" % count_pileup_file)
    try:
        sample_coverages = parse_count_pileup(count_pileup_file=count_pileup_file, snps=snps)
        logging.info("Done Loading Coverage Data of %s" % sample_name)
    except Exception, e:
        logging.exception("Failed Loading Coverage Data of %s! (Won't Delete cmpileup)" % sample_name)
        return dict()
    # create the threads (run later) for deciding strands for all regions
    threads = list()
    covered_chroms = dict()
    regions_strands = dict()
    for region, region_coverage_counts in sample_coverages.iteritems():
        regions_strands[region] = dict()
        covered_chroms[region.region] = True
        tr = threading.Thread(target=decide_strands,
                              args=(
                                  region_coverage_counts, refseq_strands[region], sema_decide,
                                  regions_strands[region], verbose),
                              name="EditingIndexStrandDecideSubprocess%s" % sample_name)
        tr.daemon = True

        threads.append(tr)
    # decide all strands (paralleling)
    logging.info("Started Deciding Strands For %s" % sample_name)
    for thr in threads:
        sema_decide.acquire()
        thr.start()
    for thr in multiprocessing.active_children():
        thr.join(3600)  # if waited for more than an hour something is wrong....
    logging.info(
        "Done Deciding Regions Strands of %s, found regions in: %s" % (sample_name, pformat(covered_chroms.keys())))

    return regions_strands


def get_stranded_strands_and_counts(refseq_strands, sample, snps, count_pileup_file_s_1,
                                    count_pileup_file_s_2):
    # type: (dict, Sample, dict, str, str) -> dict
    """
    This function parallels (using threading) the processing of the sample's data (counting coverages and mismatches, deciding strands per region)
    :param dict[Site, str] refseq_strands: The most likely strand of the regions according to the annotations (by region).
    :param Sample sample: The sample we're processing now.
    :param dict snps: A dictionary with the SNPs per mismatch type per position.
    :param str count_pileup_file_s_1: The path to the current sample cmpileup to process of one strand
    :param str count_pileup_file_s_2: The path to the current sample cmpileup to process of the other strand
    :return dict: The regions counts and strands
    """
    sample_name = sample.sample_name

    # process cmpileup
    logging.info("Started Processing The Data Of %s" % sample_name)
    try:
        logging.info("Loading First Strand Coverage File %s" % count_pileup_file_s_1)
        sample_coverages_s_1 = parse_count_pileup(count_pileup_file=count_pileup_file_s_1, snps=snps)
        logging.info("Done Loading First Strand Coverage Data of %s" % sample_name)
    except Exception, e:
        logging.exception("Failed Loading Coverage Data of %s! (Won't Delete cmpileup)" % sample_name)
        return dict()
    try:
        logging.info("Loading Second Strand Coverage File %s" % count_pileup_file_s_1)
        sample_coverages_s_2 = parse_count_pileup(count_pileup_file=count_pileup_file_s_2, snps=snps)
        logging.info("Done Loading Second Strand Coverage Data of %s" % sample_name)
    except Exception, e:
        logging.exception("Failed Loading Coverage Data of %s! (Won't Delete cmpileup)" % sample_name)
        return dict()

    covered_chroms = dict()
    refseq_vote_s_1 = get_refseq_majority(refseq_strands, sample_coverages_s_1, covered_chroms)
    refseq_vote_s_2 = get_refseq_majority(refseq_strands, sample_coverages_s_2, covered_chroms)

    both_not_unknown = not (refseq_vote_s_1 == UNKNOWN_STRAND and refseq_vote_s_2 == UNKNOWN_STRAND)
    both_different = refseq_vote_s_1 != refseq_vote_s_2

    assert both_not_unknown and both_different, "Strand Calling In Stranded Library is Ambiguous! Skipping This Sample"
    regions_strands = dict()
    if refseq_vote_s_1 == SENSE_STRAND:
        regions_strands[SENSE_STRAND] = sample_coverages_s_1
        regions_strands[ANTISENSE_STRAND] = sample_coverages_s_2
    elif refseq_vote_s_1 == ANTISENSE_NAME:
        regions_strands[ANTISENSE_STRAND] = sample_coverages_s_1
        regions_strands[SENSE_STRAND] = sample_coverages_s_2
    elif refseq_vote_s_2 == SENSE_STRAND:
        regions_strands[SENSE_STRAND] = sample_coverages_s_2
        regions_strands[ANTISENSE_STRAND] = sample_coverages_s_1
    elif refseq_vote_s_2 == ANTISENSE_NAME:
        regions_strands[ANTISENSE_STRAND] = sample_coverages_s_2
        regions_strands[SENSE_STRAND] = sample_coverages_s_1

    logging.info(
        "Done Deciding Regions Strands of %s, found regions in: %s" % (sample_name, pformat(covered_chroms.keys())))

    return regions_strands


def get_refseq_majority(refseq_strands, sample_coverages, covered_chroms):
    # type: (dict, dict, dict) -> str
    """
    This function calcs the majority  strand of the refseq annotations to decide stranded strand.
    :param dict[Site, str] refseq_strands: The refseq strands
    :param dict[str, int] sample_coverages: The counts of coverages and mismatches, as parsed from the cmpileup.
    :param dict[str, bool] covered_chroms: A dict for logging covered regions.
    :return str: The strand.
    """
    refseq_annotation_counter = {SENSE_STRAND: 0, ANTISENSE_STRAND: 0, UNKNOWN_STRAND: 0}
    for region in sample_coverages:
        covered_chroms[region.region] = True
        refseq_annotation_counter[refseq_strands[region]] += 1
    #  ration of ~2/3 of one strand in library
    if refseq_annotation_counter[SENSE_STRAND] / refseq_annotation_counter[ANTISENSE_STRAND] >= 0.67:
        refseq_majority_vote = SENSE_STRAND
    elif refseq_annotation_counter[ANTISENSE_STRAND] / refseq_annotation_counter[SENSE_STRAND] >= 0.67:
        refseq_majority_vote = ANTISENSE_STRAND
    else:
        refseq_majority_vote = UNKNOWN_STRAND

    return refseq_majority_vote


def decide_strands(region_coverage_counts, refseq_strand, semaphore_decide, regions_strands_of_sample, verbose):
    # type: (dict, str, multiprocessing.Semaphore, dict, bool) -> None
    """
    Decides, paralleling (using threading), the strand (for each type of mismatch) for the region.
    :param dict[str, int] region_coverage_counts: The counts of coverages and mismatches, as parsed from the cmpileup.
    :param str refseq_strand: The most likely strand of the region according to the annotations.
    :param multiprocessing.Semaphore semaphore_decide: The semaphore used for th4e threading.
    :param dict[str, dict] regions_strands_of_sample: The output variable, containing the decisions and count for all of the samples.
    :param bool verbose: A flag. If set, will output a verbose output for the index, including all the counts
    :return: None
    """
    try:
        regions_strands_of_sample[COUNTS_KEY] = region_coverage_counts
        mm_to_strands_d = dict()
        for strand_deciding_method, operation in STRAND_DECIDING_OPTIONS.iteritems():
            if strand_deciding_method != MAIN_STRAND_DECIDING_METHOD and not verbose:
                continue
            mm_to_strands_d[strand_deciding_method] = dict()
            for mm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
                mm_to_strands_d[strand_deciding_method][mm_type] = operation(region_coverage_counts, refseq_strand,
                                                                             mm_type)

        regions_strands_of_sample[STRANDS_KEY] = mm_to_strands_d
    except Exception, e:
        logging.exception("Failed On Strand Decisions")
    finally:
        semaphore_decide.release()


def load_resources_paths(r_ini_path):
    logging.info("Loading resources paths from %s" % r_ini_path)
    r_config = ConfigParser()
    r_config.read(r_ini_path)

    resources = {
        RESOURCES_INI_BEDTOOLS_OPT: r_config.get("DEFAULT", RESOURCES_INI_BEDTOOLS_OPT),
        RESOURCES_INI_SAMTOOLS_OPT: r_config.get("DEFAULT", RESOURCES_INI_SAMTOOLS_OPT),
        RESOURCES_INI_JAVA_OPT: r_config.get("DEFAULT", RESOURCES_INI_JAVA_OPT),
        RESOURCES_INI_JAVA_UTILS_OPT: r_config.get("DEFAULT", RESOURCES_INI_JAVA_UTILS_OPT),
        RESOURCES_INI_BAM_UTILS_OPT: r_config.get("DEFAULT", RESOURCES_INI_BAM_UTILS_OPT),
    }

    for genome in AVAILABLE_BUILTIN_GENOMES:
        resources[genome] = {
            RESOURCES_INI_GENOME_OPT: r_config.get(genome, RESOURCES_INI_GENOME_OPT),
            RESOURCES_INI_REGIONS_OPT: r_config.get(genome, RESOURCES_INI_REGIONS_OPT),
            RESOURCES_INI_CDS_REGIONS_OPT: r_config.get(genome, RESOURCES_INI_CDS_REGIONS_OPT),
            RESOURCES_INI_REFSEQ_OPT: r_config.get(genome, RESOURCES_INI_REFSEQ_OPT),
            RESOURCES_INI_SNPS_OPT: r_config.get(genome, RESOURCES_INI_SNPS_OPT),
            RESOURCES_INI_GENES_EXPRESSION_OPT: r_config.get(genome, RESOURCES_INI_GENES_EXPRESSION_OPT),
        }

    return resources


def save_compressed_pickle(resource, filename):
    with gzip.GzipFile(filename=filename, mode='w') as rzip:
        cPickle.dump(resource, rzip)


def load_compressed_pickle(filename):
    return cPickle.load(gzip.GzipFile(filename=filename, mode='rb'))


def load_resources(gene_expression_file_good, gene_expression_full_path, max_processes_strand_decision,
                   refseq_full_path, regions_bed_full_path, snps_file_good, snps_full_path, bedtools_path,
                   genome_name, regions_name, load_from_pickle=True, save_to_pickle=False):
    refseqs_pickle_file = PICKLE_FILE_FORMAT % dict(resource_file_name=refseq_full_path,
                                                    genome_name=genome_name,
                                                    regions_name=regions_name)
    refseqs_strands_pickle_file = PICKLE_FILE_FORMAT % dict(resource_file_name=refseq_full_path + ".FixedStrands",
                                                            genome_name=genome_name,
                                                            regions_name=regions_name)
    snps_pickle_file = PICKLE_FILE_FORMAT % dict(resource_file_name=snps_full_path,
                                                 genome_name=genome_name,
                                                 regions_name=regions_name)

    outof_temp_env()

    if load_from_pickle and os.path.exists(refseqs_pickle_file):
        logging.info("Loading RefSeq File From Pickled Data %s" % refseqs_pickle_file)
        refseqs = load_compressed_pickle(refseqs_pickle_file)
    else:
        logging.info("Loading RefSeq File %s" % refseq_full_path)
        refseqs = load_refseq_bed(refseq_file=refseq_full_path,
                                  regions_file=regions_bed_full_path,
                                  bedtools_path=bedtools_path)

    if load_from_pickle and os.path.exists(snps_pickle_file):
        logging.info("Loading SNPs File From Pickled Data %s" % snps_pickle_file)
        snps = load_compressed_pickle(snps_pickle_file)
    else:
        if not snps_file_good:
            logging.warn("No SNPs File Was Loaded!")
            snps = dict()
        else:
            logging.info("Loading SNPs File %s" % snps_full_path)
            snps = load_snps_file(snps_file=snps_full_path,
                                  intersect_first_with=regions_bed_full_path,
                                  bedtools_path=bedtools_path)

    if load_from_pickle and os.path.exists(refseqs_strands_pickle_file):
        logging.info("Loading RefSeq Strands File From Pickled Data %s" % refseqs_strands_pickle_file)
        refseqs_strands = load_compressed_pickle(refseqs_strands_pickle_file)
    else:
        logging.info("Loading Gene Expression File")
        if not gene_expression_file_good:
            logging.warn("No Gene Expression File Was Loaded!")
            _ = GeneRegularExpressionData(source_name=EDITING_INDEX_SOURCE)
        else:
            load_gene_expression_bed(gene_expression_file=gene_expression_full_path,
                                     source_name=EDITING_INDEX_SOURCE,
                                     intersect_first_with=refseq_full_path,
                                     bedtools_path=bedtools_path)
        logging.info("Calculating Strands For Regions According to RefSeq")
        refseqs_strands = dict()
        sema_decide = threading.Semaphore(max_processes_strand_decision)
        threads = list()
        for region in refseqs:
            tr = threading.Thread(target=derive_strand_from_refseq_p,
                                  args=(GeneRegularExpressionData(source_name=EDITING_INDEX_SOURCE), refseqs, region,
                                        sema_decide, refseqs_strands),
                                  name="EditingIndexStrandDecideSubprocessRefSeq")
            tr.daemon = True

            threads.append(tr)
        # decide all strands (paralleling)
        for thr in threads:
            sema_decide.acquire()
            thr.start()
        for thr in multiprocessing.active_children():
            thr.join(3600)  # if waited for more than an hour something is wrong....
        logging.info("Done Deciding Regions Strands For Regions According to RefSeq")

    if save_to_pickle:
        # save refseqs
        save_compressed_pickle(resource=refseqs, filename=refseqs_pickle_file)

        # save refseqs_strands
        save_compressed_pickle(resource=refseqs_strands, filename=refseqs_strands_pickle_file)

        # save snps
        save_compressed_pickle(resource=snps, filename=snps_pickle_file)

    return refseqs, refseqs_strands, snps


def main(root_dir, output_dir, output_dir_summery, genome_resources, log_path, bam_files_suffix,
         groups_file, include_paths, include_paths_operator, exclude_paths, exclude_paths_operator, recursion_depth,
         follow_links, defaults_override_conf, defaults_override_args, genome_path, regions_bed, snps_file, refseq_file,
         gene_expression_file, reload_pickled_resoures, use_refseq_all, get_regions_metadata, just_get_configs,
         per_sample_output, delete_cmpileups, max_processes_sample, max_processes_strand_decision, verbose, is_stranded,
         is_paired_end, possible_mismatches, q_threshold):
    # type: (str, str, str, str, str, str, str, list,  operator, list, operator, int, bool, str, dict, str, str, str, str, str, bool, bool, bool,bool,bool, bool, bool, int, int, bool, bool, bool, str, str) -> None
    """
    :param str root_dir: The directory in which the BAM files from the aligner are. Aligner logs should be there to if availabe.
    :param str output_dir: The output directory for the run outputs.
    :param str output_dir_summery: The output directory for the summary output.
    :param str genome_resources: one of the builtin genomes resources groups
    :param str groups_file: The path of the groups file if given, otherwise all sample are considered if the same group.
    :param str log_path: The path to the log directory.
    :param str bam_files_suffix: The suffixes of the BAM files to run on.
    :param list[str] include_paths: List of path fragments that must be in path for it to be included.
    :param operator include_paths_operator: The operator used to check for fragments (AND or OR)
    :param list[str] exclude_paths: List of path fragments that mustn't be in path for it to be included.
    :param operator exclude_paths_operator: The operator used to check for fragments (AND or OR)
    :param int recursion_depth: The recursion depth to look for files from root directory
    :param bool follow_links: If set will follow links (a dir command param)
    :param str defaults_override_conf: A path to a cofig file for overriding default setting conf.
    :param dict[str, str] defaults_override_args: A dictionary of options to override in the defaults,
     will also override defaults_override_conf.
    :param str regions_bed: The path of the edited regions bed file to use instead of the builtin file
    :param str genome_path: The path of the genome fasta to use instead of the builtin file
    :param str snps_file: A path to a user given SNPs file to use instead of the builtin file.
    :param str refseq_file: A path to a user given RefSeq file to use instead of the builtin file.
    :param str gene_expression_file: A path to a user given gene expression file to use instead of the builtin file.
    :param bool reload_pickled_resoures: If set, will try to load pre-processed resources saved by PreloadAEIResources.
    :param bool use_refseq_all: If set, will use the RefSeq All Table instead of the RefSeq Curated table
    :param bool get_regions_metadata: A flag. If set, will also print per region data (big output)
    :param bool just_get_configs: A flag. If set, will not run pipeline part (i.e. will not generate the cmpileups and assume they exist)
    :param bool per_sample_output: A flag. If set,  will also print per region *per sample* data (*very* big output)
    :param bool delete_cmpileups: A flag. If set, will delete the cmpileups after conversion, maintaining relatively low memory profile.
    :param int max_processes_sample: The maximal number of samples to process in parallel.
    :param int max_processes_strand_decision:  The maximal strand decisions per sample to process in parallel.
    :param bool verbose: A flag. If set, will output a verbose output for the index, including all the counts
    :param bool is_stranded: A flag. If set, will treat the data as stranded.
    :param bool is_paired_end: A flag. If set, will treat the data as paired end.
    :param str possible_mismatches: the key for POSSIBLE_MM, deciing what mismatch types to output.
    :return:None
    """
    script_dir = path_locator.module_path()
    timestamp = datetime.today().isoformat()
    init_logging_dict(os.path.join(log_path, ".".join(["EditingIndex", timestamp, "log"])))

    resources_paths = load_resources_paths(RESOURCES_INI_FORMAT % dict(script_dir=script_dir))

    if genome_resources != NOT_BUILTIN_GENOME:
        using_builtins_resources = True
        builtins_resources = resources_paths[genome_resources]
    else:
        using_builtins_resources = False
        builtins_resources = {}

    assert (genome_path == "" and using_builtins_resources) or os.path.exists(genome_path), PATH_UNAVAILABLE_ERR % (
        "Genome File", genome_path)
    assert (regions_bed == "" and using_builtins_resources) or os.path.exists(regions_bed), PATH_UNAVAILABLE_ERR % (
        "Regions BED File", regions_bed)
    assert (refseq_file == "" and using_builtins_resources) or os.path.exists(refseq_file), PATH_UNAVAILABLE_ERR % (
        "Refseq Annotations File", refseq_file)

    snps_file_good = (snps_file == "" and using_builtins_resources) or os.path.exists(snps_file)
    if not snps_file_good and snps_file:
        logging.warn(PATH_UNAVAILABLE_ERR % ("SNPs File", snps_file))

    gene_expression_file_good = (gene_expression_file == "" and using_builtins_resources) or os.path.exists(
        gene_expression_file)
    if not gene_expression_file and gene_expression_file:
        logging.warn(PATH_UNAVAILABLE_ERR % ("Genes Expression File", gene_expression_file))

    genome_full_path = builtins_resources[RESOURCES_INI_GENOME_OPT] if (
            genome_path == "" and using_builtins_resources) else genome_path

    refseq_r_opt = RESOURCES_INI_REFSEQ_ALL_OPT if use_refseq_all else RESOURCES_INI_REFSEQ_OPT
    cds_reg_opt = RESOURCES_INI_CDS_REGIONS_REFSEQ_ALL_OPT if use_refseq_all else RESOURCES_INI_CDS_REGIONS_OPT

    regions_bed_full_path = builtins_resources[cds_reg_opt] if (
            regions_bed == "" and using_builtins_resources) else regions_bed
    snps_full_path = builtins_resources[RESOURCES_INI_SNPS_OPT] if (
            snps_file == "" and using_builtins_resources) else snps_file
    refseq_full_path = builtins_resources[refseq_r_opt] if (
            refseq_file == "" and using_builtins_resources) else refseq_file
    gene_expression_full_path = builtins_resources[RESOURCES_INI_GENES_EXPRESSION_OPT] if (
            gene_expression_file == "" and using_builtins_resources) else gene_expression_file

    #  Old - genome_index_path = os.path.join(output_dir, os.path.basename(genome_full_path) + ".GenomeIndex.jsd")
    # Cloud Optimized - if created once it remains for all runs of the instance.
    regions_dir = os.path.dirname(regions_bed_full_path)
    regions_name = os.path.basename(regions_bed_full_path)
    genome_name = os.path.basename(genome_full_path)
    genome_index_name = "%(genome_name)s.%(regions_name)s.GenomeIndex.jsd" % dict(regions_name=regions_name,
                                                                                  genome_name=genome_name)
    genome_index_path = os.path.join(regions_dir, genome_index_name)

    resources_d = {GENOME_FASTA_OPTION: genome_full_path,
                   REGIONS_BED_OPTION: regions_bed_full_path,
                   BEDTOOLS_PATH_OPTION: resources_paths[RESOURCES_INI_BEDTOOLS_OPT],
                   SAMTOOLS_PATH_OPTION: resources_paths[RESOURCES_INI_SAMTOOLS_OPT],
                   JAVA_HOME_OPTION: resources_paths[RESOURCES_INI_JAVA_OPT],
                   JAVA_UTILS_HOME_OPTION: resources_paths[RESOURCES_INI_JAVA_UTILS_OPT],
                   BAM_UTILS_PATH_OPTION: resources_paths[RESOURCES_INI_BAM_UTILS_OPT],
                   BAM_FILE_SUFFIX_OPTION: bam_files_suffix,
                   GENOME_INDEX_PATH_OPTION: genome_index_path,
                   REGIONS_NAME_OPTION: regions_name,
                   OVERALL_MAX_PROCESSES_NUM_OPTION: str(max_processes_sample),
                   Q_THRESHOLD_OPTION: q_threshold}
    resources_d.update(defaults_override_args)

    logging.info("Running with Resource Files and Run Paths: %s", pformat(resources_d))
    predicates_dict = {ALL: lambda path: path.endswith(bam_files_suffix)}
    bam_files = get_paths(root_path=root_dir,
                          must_include_paths=include_paths,
                          must_include_operator=include_paths_operator,
                          exclude_paths=exclude_paths,
                          exclude_operator=exclude_paths_operator,
                          follow_links=follow_links,
                          recursion_depth=recursion_depth,
                          predicates_dict=predicates_dict)[ALL]

    if groups_file:
        load_groups_and_samples_from_file(groups_file=groups_file)

    if is_stranded:
        full_conf_path = STRANDED_FULL_CONFIG_PATH_FORMAT % dict(script_dir=script_dir)
    else:
        full_conf_path = FULL_CONFIG_PATH_FORMAT % dict(script_dir=script_dir)

    if is_paired_end:
        resources_d[IS_PE_OPTION] = "True"

    samples_confs = run_pipeline(root_path=root_dir,
                                 output_dir=output_dir,
                                 bam_files_suffix=bam_files_suffix,
                                 full_conf_path=full_conf_path,
                                 defaults_conf_path=DEFAULTS_CONFIG_PATH_FORMAT % dict(script_dir=script_dir),
                                 defaults_override_conf_path=defaults_override_conf,
                                 defaults_override_args=resources_d,
                                 files=bam_files,
                                 log_dir=log_path,
                                 create_na_group=None is groups_file,
                                 just_get_configs=just_get_configs)

    if len(samples_confs) == 0:
        logging.warn("No Samples Were Found To Run On! Exiting...")
        return

    bedtools_path = samples_confs.values()[0].defaults(formatted=True)[BEDTOOLS_PATH_OPTION]

    refseqs, refseqs_strands, snps = load_resources(gene_expression_file_good=gene_expression_file_good,
                                                    gene_expression_full_path=gene_expression_full_path,
                                                    max_processes_strand_decision=max_processes_strand_decision,
                                                    refseq_full_path=refseq_full_path,
                                                    regions_bed_full_path=regions_bed_full_path,
                                                    snps_file_good=snps_file_good,
                                                    snps_full_path=snps_full_path,
                                                    bedtools_path=bedtools_path,
                                                    genome_name=genome_name,
                                                    regions_name=regions_name,
                                                    load_from_pickle=reload_pickled_resoures)

    get_index(refseq_strands=refseqs_strands,
              refseq_annotations=refseqs,
              samples_configs=samples_confs,
              summery_dir=output_dir_summery,
              pipeline_output_dir=output_dir,
              snps=snps,
              get_regions_data=get_regions_metadata,
              per_sample_output=per_sample_output,
              delete_cmpileups=delete_cmpileups,
              verbose=verbose,
              max_processes_sample=max_processes_sample,
              max_processes_strand_decision=max_processes_strand_decision,
              is_stranded=is_stranded,
              possible_mismatches=POSSIBLE_MM[possible_mismatches])
    logging.info("Done Writing Index Data")


if __name__ == "__main__":
    desc = """Run on each file in a given directory a set of steps."""
    script_dir = path_locator.module_path()

    parser = argparse.ArgumentParser(prog='Editing Index Runner', description=desc,
                                     formatter_class=SortingHelpFormatter)

    inputs_g = parser.add_argument_group(title="Input Files Options:")
    add_get_paths_function_to_argparser(parser=inputs_g)
    inputs_g.add_argument('-f', '--bam_files_suffix', metavar="bam files suffix", dest='files_suffix', nargs='?',
                          required=False,
                          default="Aligned.sortedByCoord.out.bam",
                          help="A suffix of the BAM files to run on (e.g. Sorted.By.Coord.bam). Should be the *full* "
                               "suffix")
    inputs_g.add_argument('--stranded', dest='is_stranded', action='store_true', required=False,
                          help="A flag. If set, will treat the data as stranded")
    inputs_g.add_argument('--paired_end', dest='is_paired_end', action='store_true', required=False,
                          help="A flag. If set, will treat the data as paired end. *currently affecting only stranded"
                               " runs*")
    inputs_g.add_argument('--q_threshold', dest='q_threshold', nargs="?", required=False, default="30",
                          help="The minimal Phred quality score to count.")

    outputs_g = parser.add_argument_group(title="Output Files Options:")
    outputs_g.add_argument('-o', '--output_dir', metavar="output_dir", dest='output_dir', nargs='?', required=False,
                           default=".", help="The root directory for the cmpileup creation outputs"
                                             " (will create sub-dirs per sample). If keep_cmpileup is not set, all "
                                             "content here will be deleted.")
    outputs_g.add_argument('--per_region_output', dest='per_region_output', action='store_true', required=False,
                           help="If set, will output regions data (raw counts, big output)")
    outputs_g.add_argument('--keep_cmpileup', dest='keep_cmpileup', action='store_true', required=False,
                           help="If set, will not delete the cmpileups")
    outputs_g.add_argument('-os', '--output_dir_summery', metavar="output_dir_summery", dest='output_dir_summery',
                           nargs='?', required=False, default=".", help="The directory for the summary output.")
    outputs_g.add_argument('--per_sample_output', dest='per_sample_output', action='store_true', required=False,
                           help="If set,  will output regions metadata per sample too (very big output)")
    outputs_g.add_argument('--verbose', dest='verbose', action='store_true', required=False,
                           help="If set, will give a verbose output including all counts, otherwise will output only "
                                "indexes")
    outputs_g.add_argument("-mm", "--mismatches", metavar="wanted mismatches", dest="possible_mismatches",
                           required=False, choices=POSSIBLE_MM, default=UNSTRANDED_MM,
                           help="One of %s - Which mismatch types should be present in the output -"
                                " only one for each complementary pair (A2C, A2G, A2T, C2A, C2G, C2T), or all types"
                                "." % POSSIBLE_MM.keys())

    resources_g = parser.add_argument_group(title="Resources Options:",
                                            description="These specify the files to use for the analysis."
                                                        " Grouped builtin files can be chosen with the --genome option"
                                                        " and specific files can be overridden with the other options.")
    resources_g.add_argument('--genome', metavar="genome resources", dest='genome_resources', nargs='?', required=True,
                             help="one of the builtin genomes resources groups -(%s). If '%s' is provided, all"
                                  " resources files must be supplied manually" % (
                                      ", ".join(AVAILABLE_BUILTIN_GENOMES), NOT_BUILTIN_GENOME)
                             , choices=AVAILABLE_GENOMES)
    resources_g.add_argument('-rb', '--regions', metavar="regions bed", dest='regions_bed', nargs='?', default="",
                             required=False,
                             help="A path to an edited regions bed file to use instead of the builtin file")
    resources_g.add_argument('--snps', metavar="SNPs file", dest='snps_file', nargs='?', required=False,
                             default="", help="A path to a SNPs file to use instead of the builtin file")
    resources_g.add_argument('--refseq', metavar="refseq file", dest='refseq_file', nargs='?', required=False,
                             default="", help="A path to a refseq file to use instead of the builtin file")
    resources_g.add_argument('--genes_expression', metavar="genes expression file", dest='genes_expression_file',
                             default="", nargs='?', required=False,
                             help="A path to a genes expression file to use instead of the builtin file")
    resources_g.add_argument('-gf', '--genome_fasta', metavar="genome fasta", dest='genome_path', nargs='?',
                             default="", required=False,
                             help="The path of the genome fasta to use instead of the builtin file")
    resources_g.add_argument('--reload_pickled_resoures', dest='reload_pickled_resoures', action='store_true',
                             required=False, help="If set, will try to load pre-processed resources saved by "
                                                  "PreloadAEIResources.")
    resources_g.add_argument('--refseq_all', dest='use_refseq_all', action='store_true',
                             required=False,
                             help="If set, will use the RefSeq All Table instead of the RefSeq Curated table")
    add_groups_file_to_argparser(parser=resources_g)

    run_configs = parser.add_argument_group(title="Run Configuration Options:")
    run_configs.add_argument('-l', '--log_path', metavar="log_path", dest='log_path', nargs='?', required=False,
                             default=".",
                             help="The path where the logs (and flags) will be written.")
    run_configs.add_argument('-c', '--config_override_file', metavar="config_override_file",
                             dest='config_override_file',
                             nargs='?', required=False,
                             help="A path to a configuration file (in INI format) to override "
                                  "the values in the editing index config file found in %s." %
                                  os.path.join(script_dir, "Resources"))
    run_configs.add_argument('-a', '--args', metavar="override config extra args", dest='args', required=False,
                             default={}, nargs='*',
                             help='named args (in the format of <var>=\"<val>\",<var>=\"<val>\") to '
                                  'overrid the defaults config.')
    run_configs.add_argument('--recalc', dest='just_get_configs', action='store_true', required=False,
                             help="If set, will not the run pipeline part (i.e. will not generate the cmpileups,"
                                  " instead the script assumes they exist). Notice - this option ignores the flags")

    sys_opts = parser.add_argument_group(title="System Options:")
    sys_opts.add_argument('--ts', metavar="sample threads", dest='max_processes_sample', nargs='?', required=False,
                          default=10, type=int, help="The number of samples to process in parallel")
    sys_opts.add_argument('--tsd', metavar="sample strands threads", dest='max_processes_strand_decision', nargs='?',
                          required=False, default=50, type=int,
                          help="The maximal strand decisions per sample to process in parallel")

    options = parser.parse_args()
    args = convert_args_to_dict(options.args)

    if options.per_sample_output and not options.per_region_output:
        parser.error("ERROR: per_sample_output is an extension of per_region_output,"
                     " per_region_output must be set for per_sample_output")

    main(root_dir=options.root_dir,
         output_dir=options.output_dir,
         output_dir_summery=options.output_dir_summery,
         genome_resources=options.genome_resources,
         log_path=options.log_path,
         bam_files_suffix=options.files_suffix,
         groups_file=options.groups_file,
         include_paths=options.include_prefixes,
         include_paths_operator=options.include_operator,
         exclude_paths=options.exclude_prefixes,
         exclude_paths_operator=options.exclude_operator,
         recursion_depth=options.recursion_depth,
         follow_links=options.follow_links,
         defaults_override_conf=options.config_override_file,
         defaults_override_args=args,
         genome_path=options.genome_path,
         regions_bed=options.regions_bed,
         snps_file=options.snps_file,
         refseq_file=options.refseq_file,
         gene_expression_file=options.genes_expression_file,
         reload_pickled_resoures=options.reload_pickled_resoures,
         use_refseq_all=options.use_refseq_all,
         get_regions_metadata=options.per_region_output,
         just_get_configs=options.just_get_configs,
         per_sample_output=options.per_sample_output,
         delete_cmpileups=not options.keep_cmpileup,
         max_processes_sample=options.max_processes_sample,
         max_processes_strand_decision=options.max_processes_strand_decision,
         verbose=options.verbose,
         is_stranded=options.is_stranded,
         is_paired_end=options.is_paired_end,
         possible_mismatches=options.possible_mismatches,
         q_threshold=options.q_threshold
         )
