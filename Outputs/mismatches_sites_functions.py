from Commons.general_functions import get_sample_paths_to_unique_name

__author__ = 'Hillel'
"""
This file contains output records generating function for output of mismatches sites.
"""
# =====================imports====================#

from DataObjects.KnownProperties.MismatchesSites import MismatchesSites, GroupMMRate
from DataObjects.KnownProperties.Group import Group
from DataObjects.BaseClasses.EntityRecord import EntityRecord
from DataObjects.BaseClasses.Sample import Sample

from Commons.consts import GENOMIC_REGION, SITE_START, SITE_END, SITE_STRAND, MISMATCH_TYPE, RECORD_HEADER_FORMAT, \
    MISMATCH_PREVALENCE, MISMATCH_RATE, TOTAL_READS, MISMATCHED_READS, ALL, NA, DATA_SOURCE, CANONICAL_READS, \
    MULTI_VALUE_OUTPUT_SEP
from Commons.data_structs import Site

# =====================consts=====================#
GENERAL_HEADERS_KEY = "general_headers_for_output"
GROUPS_PREVALENCE_HEADERS_KEY = "groups_prevalence_headers_for_output"
GROUPS_RATIO_HEADERS_KEY = "groups_ratio_headers_for_output"

SEP = MULTI_VALUE_OUTPUT_SEP  # for readabilty


# =====================functions==================#

def get_mismatches_sites_quantification_per_site(source_name, coverage_threshold=1):
    # type: (str, int) -> (list,dict)
    """
    This fucntion generates output records containing all data from a MismatchesSites
    :param str source_name: The name of the source
    :param int coverage_threshold: The threshold beneath which to ignore data
    :return: A list of all the output recs dicts.
    """
    res_recs = []
    headers_per_sample = dict()
    headers_per_sample[GENERAL_HEADERS_KEY] = [DATA_SOURCE, GENOMIC_REGION, SITE_START, SITE_END, SITE_STRAND,
                                               MISMATCH_TYPE]
    samples_unique = Sample.sample_names_unique()
    if not samples_unique:
        samples_paths_to_unique_names = get_sample_paths_to_unique_name()
    mismatches = MismatchesSites(source=source_name)
    overall_site_prevalence_h = RECORD_HEADER_FORMAT % dict(record_name=ALL, data_type=MISMATCH_PREVALENCE)
    headers_per_sample.setdefault(ALL, []).append(overall_site_prevalence_h)
    prevalence_per_group = mismatches.get_mismatches_group_prevalence(coverage_threshold=coverage_threshold,
                                                                      by_name=True)
    prevalence_per_site = mismatches.get_mismatches_sample_prevalence(coverage_threshold=coverage_threshold)
    sites_rate_per_group = mismatches.get_mismatches_group_rates(by_name=True)
    for region, start_d in mismatches.get_sites().iteritems():
        for start, end_d in start_d.iteritems():
            for end, strand_d in end_d.iteritems():
                for strand, mm_dict in strand_d.iteritems():
                    site = Site(region, start, end, strand)
                    for mm, samples_dict in mm_dict.iteritems():
                        found_samples = dict()
                        rec = {DATA_SOURCE: source_name,
                               GENOMIC_REGION: str(region),
                               SITE_START: str(start),
                               SITE_END: str(end),
                               SITE_STRAND: str(strand),
                               MISMATCH_TYPE: mm,
                               overall_site_prevalence_h: str(prevalence_per_site[site][mm])}
                        for group in EntityRecord.get_all_records(of_types=(Group,)).values():
                            group_prevalence = str(prevalence_per_group[site].get(mm, {}).get(group.group_name, 0))

                            group_prevalence_h = RECORD_HEADER_FORMAT % dict(record_name=group.group_name,
                                                                             data_type=MISMATCH_PREVALENCE)
                            headers_per_sample.setdefault(GROUPS_PREVALENCE_HEADERS_KEY, set()).add(group_prevalence_h)
                            rec[group_prevalence_h] = group_prevalence

                            group_rates = sites_rate_per_group[site].get(mm, {}).get(group.group_name,
                                                                                     GroupMMRate(group, NA, NA, NA, NA,
                                                                                                 NA))
                            group_rate_h = RECORD_HEADER_FORMAT % dict(record_name=group.group_name,
                                                                       data_type=MISMATCH_RATE)
                            group_total_h = RECORD_HEADER_FORMAT % dict(record_name=group.group_name,
                                                                        data_type=TOTAL_READS)
                            group_mismatched_h = RECORD_HEADER_FORMAT % dict(record_name=group.group_name,
                                                                             data_type=MISMATCHED_READS)
                            group_canonical_h = RECORD_HEADER_FORMAT % dict(record_name=group.group_name,
                                                                            data_type=CANONICAL_READS)
                            headers_per_sample.setdefault(GROUPS_RATIO_HEADERS_KEY, []).extend([group_total_h,
                                                                                                group_canonical_h,
                                                                                                group_mismatched_h,
                                                                                                group_rate_h])
                            if group_rates.num_of_total_reads < coverage_threshold:
                                rec[group_rate_h] = NA
                                rec[group_canonical_h] = NA
                                rec[group_mismatched_h] = NA
                            else:
                                rec[group_rate_h] = str(group_rates.ratio)
                                rec[group_canonical_h] = str(group_rates.num_of_canonical_reads)
                                rec[group_total_h] = str(group_rates.num_of_total_reads)
                                rec[group_mismatched_h] = str(group_rates.num_of_mm_reads)

                        for sample_id, sample_mm_rate in samples_dict.iteritems():
                            sample_name = sample_mm_rate.sample.sample_name if samples_unique else \
                                samples_paths_to_unique_names[sample_mm_rate.sample.sample_path]
                            found_samples[sample_id] = True
                            sample_rate_h = RECORD_HEADER_FORMAT % dict(record_name=sample_name,
                                                                        data_type=MISMATCH_RATE)
                            sample_total_h = RECORD_HEADER_FORMAT % dict(record_name=sample_name,
                                                                         data_type=TOTAL_READS)
                            sample_mismatched_h = RECORD_HEADER_FORMAT % dict(record_name=sample_name,
                                                                              data_type=MISMATCHED_READS)
                            sample_canonical_h = RECORD_HEADER_FORMAT % dict(record_name=sample_name,
                                                                             data_type=CANONICAL_READS)
                            headers_per_sample[sample_name] = [sample_total_h, sample_canonical_h, sample_mismatched_h,
                                                               sample_rate_h]
                            if sample_mm_rate.num_of_total_reads < coverage_threshold:
                                rec[sample_rate_h] = NA
                                rec[sample_canonical_h] = NA
                                rec[sample_mismatched_h] = NA
                            else:
                                rec[sample_rate_h] = str(sample_mm_rate.ratio)
                                rec[sample_canonical_h] = str(sample_mm_rate.num_of_canonical_reads)
                                rec[sample_total_h] = str(sample_mm_rate.num_of_total_reads)
                                rec[sample_mismatched_h] = str(sample_mm_rate.num_of_mm_reads)

                        for sample in EntityRecord.get_all_records(of_types=(Sample,)).values():
                            if not found_samples.get(sample.record_id, False):
                                sample_name = sample.sample_name if samples_unique else \
                                    samples_paths_to_unique_names[sample_mm_rate.sample.sample_path]

                                sample_rate_h = RECORD_HEADER_FORMAT % dict(record_name=sample_name,
                                                                            data_type=MISMATCH_RATE)
                                sample_total_h = RECORD_HEADER_FORMAT % dict(record_name=sample_name,
                                                                             data_type=TOTAL_READS)
                                sample_mismatched_h = RECORD_HEADER_FORMAT % dict(record_name=sample_name,
                                                                                  data_type=MISMATCHED_READS)
                                sample_canonical_h = RECORD_HEADER_FORMAT % dict(record_name=sample_name,
                                                                                 data_type=CANONICAL_READS)
                                headers_per_sample[sample_name] = [sample_total_h, sample_canonical_h,
                                                                   sample_mismatched_h, sample_rate_h]
                                rec[sample_rate_h] = NA
                                rec[sample_canonical_h] = NA
                                rec[sample_total_h] = NA
                                rec[sample_mismatched_h] = NA

                        res_recs.append(rec)
    return res_recs, headers_per_sample
