
__author__ = 'Hillel'
"""
This file contains output records generating function for output of mhc alleles.
"""
# =====================imports====================#

from DataObjects.KnownProperties.MHCAlleles import MHCAlleles
from Commons.data_structs import MHCAllele
from DataObjects.KnownProperties.Group import Group
from DataObjects.BaseClasses.EntityRecord import EntityRecord
from DataObjects.BaseClasses.Sample import Sample

from Commons.consts import MHC_CLASS, MHC_LOCUS, MHC_ALLELE, GROUP, SAMPLE, SAMPLE_PATH, RECORD_HEADER_FORMAT, \
    MHC_PREVALENCE, ALL, NA, DATA_SOURCE, STATISTICS_DATA, MULTI_VALUE_OUTPUT_SEP, KEY_VAL_OUTPUT_SEP
from Commons.general_functions import get_groups_and_sample_names_dict
# =====================consts=====================#
GENERAL_HEADERS_KEY = "general_headers_for_output"
GROUPS_PREVALENCE_HEADERS_KEY = "groups_prevalence_headers_for_output"
GROUPS_RATIO_HEADERS_KEY = "groups_ratio_headers_for_output"


# =====================functions==================#

def mhc_alleles_prevalence_summery(source_name):
    # type: (str) -> (list,dict)
    """
    This fucntion generates output records containing all data from a MismatchesSites
    :param str source_name: The name of the source
    :return: A list of all the output recs dicts.
    """
    res_recs = []
    headers_per_sample = dict()
    headers_per_sample[ALL] = [DATA_SOURCE, MHC_CLASS, MHC_LOCUS, MHC_ALLELE]

    mhc_alleles = MHCAlleles.get_mhc_alleles_by_source(source=source_name)
    overall_allele_prevalence_h = RECORD_HEADER_FORMAT % dict(record_name=ALL, data_type=MHC_PREVALENCE)
    headers_per_sample[ALL].append(overall_allele_prevalence_h)
    prevalence_per_group = mhc_alleles.get_mhc_alleles_group_prevalence(by_name=True)
    prevalence_per_allele = mhc_alleles.get_mhc_alleles_sample_prevalence()

    for mhc_allele, samples_dict in mhc_alleles.get_mhc_alleles().iteritems():
        assert isinstance(mhc_allele, MHCAllele)
        rec = {DATA_SOURCE: source_name,
               MHC_CLASS: str(mhc_allele.mhc_class),
               MHC_LOCUS: str(mhc_allele.locus),
               MHC_ALLELE: str(mhc_allele.allele),
               overall_allele_prevalence_h: str(prevalence_per_allele[mhc_allele])}
        for group in EntityRecord.get_all_records(of_types=(Group,)).values():
            group_prevalence = str(prevalence_per_group[mhc_allele].get(group.group_name, 0))
            group_prevalence_h = RECORD_HEADER_FORMAT % dict(record_name=group.group_name, data_type=MHC_PREVALENCE)
            headers_per_sample[ALL].append(group_prevalence_h)
            rec[group_prevalence_h] = group_prevalence

        res_recs.append(rec)

    return res_recs, headers_per_sample


def mhc_alleles_per_sample_summery(source_name):
    # type: (str) -> (list,dict)
    """
    This fucntion generates output records containing all data from a MismatchesSites
    :param str source_name: The name of the source
    :return: A list of all the output recs dicts.
    """
    res_recs = []
    headers_per_sample = dict()
    headers_per_sample[ALL] = [DATA_SOURCE, GROUP, SAMPLE, SAMPLE_PATH, MHC_CLASS, MHC_LOCUS, MHC_ALLELE,
                               STATISTICS_DATA]
    mhc_alleles = MHCAlleles.get_mhc_alleles_by_source(source=source_name)
    groups_and_sample_names_dict = get_groups_and_sample_names_dict(by_sample=True)
    found_samples = dict()
    for mhc_allele, samples_dict in mhc_alleles.get_mhc_alleles().iteritems():
        assert isinstance(mhc_allele, MHCAllele)

        for sample_id, p_vals in samples_dict.iteritems():
            rec = {DATA_SOURCE: source_name,
                   MHC_CLASS: str(mhc_allele.mhc_class),
                   MHC_LOCUS: str(mhc_allele.locus),
                   MHC_ALLELE: str(mhc_allele.allele)}
            sample = Sample.entity_record_by_id(sample_id)
            """:type sample Sample"""
            sample_name = sample.sample_name
            rec[SAMPLE_PATH] = sample.sample_path
            rec[SAMPLE] = sample_name
            sample_groups = groups_and_sample_names_dict[sample.sample_name].keys()
            sample_groups = [group[0] for group in sorted(sample_groups, key=lambda x: x[1])]
            rec[GROUP] = MULTI_VALUE_OUTPUT_SEP.join(sample_groups)
            rec[STATISTICS_DATA] = MULTI_VALUE_OUTPUT_SEP.join(
                [KEY_VAL_OUTPUT_SEP.join([test, str(pval)]) for test, pval in p_vals.iteritems()]
                                                            )
            found_samples[sample_id] = True
            res_recs.append(rec)
    for sample in EntityRecord.get_all_records(of_types=(Sample,)).values():
        if not found_samples.get(sample.record_id, False):
            """:type sample Sample"""
            sample_groups = groups_and_sample_names_dict[sample.sample_name].keys()
            sample_groups = [group[0] for group in sorted(sample_groups, key=lambda x: x[1])]

            rec = {DATA_SOURCE: source_name,
                   GROUP: MULTI_VALUE_OUTPUT_SEP.join(sample_groups),
                   SAMPLE: sample.sample_name,
                   SAMPLE_PATH: sample.sample_path,
                   MHC_CLASS: str(NA),
                   MHC_LOCUS: str(NA),
                   MHC_ALLELE: str(NA),
                   STATISTICS_DATA: str(NA)
                   }

            res_recs.append(rec)

    return res_recs, headers_per_sample

