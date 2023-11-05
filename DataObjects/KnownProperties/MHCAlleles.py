__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Imports
# endregion

# region Internal Imports
from Commons.consts import NA
from Commons.data_structs import StatisticalTestPValue, MHCAllele
from Commons.general_functions import deep_copy_dict, convert_params_to_bool_dict
from DataObjects.BaseClasses.KnownProperty import KnownProperty
from DataObjects.BaseClasses.Sample import Sample
from DataObjects.KnownProperties.Group import Group


# endregion
# =====================constants===================#

# =====================classes=====================#

class MHCAlleles(KnownProperty):
    """
    This class holds the Predicted MHC alleles per sample.
    :ivar __mhc_alleles_by_alleles: (get - mismatch_sites, to add sites or remove use - add/remove_sites) the sites dict.
        <mismatch type> ><Chr/Transcript name>:<start>:<end>:[sample, mismatched reads, conincal ones, totoal]
    :type __mhc_alleles_by_alleles: C{dict} of C{str} to
    C{dict} of C{str} to C{dict} of C{str} to C{dict} of C{str} to  C{list} of L{SampleMMRate}
    :ivar __mhc_alleles_by_sample: same as __regions_data_by_sites except mismatchs type is the root for the
     mapping
    :type __mhc_alleles_by_sample: C{dict} of C{str} to
    C{dict} of C{str} to C{dict} of C{str} to C{dict} of C{str} to  C{list} of L{SampleMMRate}
    """

    def dump_me(self):
        pass

    def __init__(self, source):
        """
        :param source: The (data-wise) origin of the data.
        :type source: C{str}
        """
        if hasattr(self, "source"):
            return
        self.source = source
        self.__mhc_alleles_by_alleles = {}
        self.__mhc_alleles_by_sample = {}

    def add_allele(self, sample, mhc_class, locus, allele, p_val=NA):
        # type: (Sample, str, str, str, StatisticalTestPValue) -> object
        """
        This function updates the alleles dictionary with new alleles
        :param Sample sample: The instance of the sample
        :param str mhc_class: The class of the MHC.
        :param str locus: the locus (name) of the allele.
        :param str allele: the full allele name (e.g. A*01:21).
        :param StatisticalTestPValue p_val: p-value or ranking for this prediction.
        :return: None
        """

        assert isinstance(sample, Sample)
        mhc_allele = MHCAllele(mhc_class, locus, allele)
        s_id = sample.record_id
        self.__mhc_alleles_by_alleles.setdefault(mhc_allele, {})[s_id] = p_val
        self.__mhc_alleles_by_sample.setdefault(s_id, {})[mhc_allele] = p_val
        # "relate" the instances
        sample.add_related_entity(self)
        self.add_related_entity(sample)

    def remove_allele(self,  mhc_class, allele, locus, sample):
        """
        This function updates the sites dictionary removing the mismatch sites
        :param Sample sample: The instance of the sample
        :param str mhc_class: The class of the MHC.
        :param str locus: the locus (name) of the allele.
        :param str allele: the full allele name (e.g. A*01:21).
        :param allele: The allele (or transcript name)
        """
        assert isinstance(sample, Sample)
        mhc_allele = MHCAllele(mhc_class, locus, allele)

        if mhc_allele not in self.__mhc_alleles_by_alleles:
            return  # nothing to delete...
        _ = self.__mhc_alleles_by_alleles[mhc_allele].pop(sample.record_id)
        if not self.__mhc_alleles_by_alleles[mhc_allele].keys():
            _ = self.__mhc_alleles_by_alleles.pop(mhc_allele)
        _ = self.__mhc_alleles_by_sample[sample.record_id].pop(mhc_allele)
        if not self.__mhc_alleles_by_sample[sample.record_id].keys():
            self.__mhc_alleles_by_sample.pop(sample.record_id)

    def get_mhc_alleles(self, dont_filter=True, by_sample=False, mhc_classes=None, alleles=None, sample_ids=None):
        """
        Retrieves the sites.
        :param dont_filter: A flag. if set will retrieving a simple copy without filtration. Use to save run-time.
        :param by_sample: A flag. if set will return sites by sa,ple. Otherwise by alleles (see class doc)
        :param mhc_classes: if given will return only alleles that match this mismatch type.
        :param alleles: if given will return only alleles that belong to this sample_id.
        :param sample_ids: if given will return only alleles that belong to this sample_id.
        :return: the alleles relevant in the same structure as the calls dict.
        """
        mhc_classes_d = convert_params_to_bool_dict(mhc_classes)
        alleles_d = convert_params_to_bool_dict(alleles)
        sample_ids_d = convert_params_to_bool_dict(sample_ids)
        if by_sample:
            filtered_copy = deep_copy_dict(self.__mhc_alleles_by_sample)
            if dont_filter:
                return filtered_copy
            for s_id, sample_alleles_dict in filtered_copy.iteritems():
                if sample_ids and not sample_ids_d.get(s_id, False):
                    _ = filtered_copy.pop(s_id)
                    continue
                for mhc_class, alleles_dict in sample_alleles_dict.iteritems():
                    if mhc_classes and not mhc_classes_d.get(mhc_class, False):
                        _ = sample_alleles_dict.pop(mhc_class)
                        continue
                    for allele, pval in alleles_dict.iteritems():
                        if alleles and not alleles_d.get(allele, False):
                            _ = alleles_dict.pop(allele)
                            continue

        else:
            filtered_copy = deep_copy_dict(self.__mhc_alleles_by_alleles)
            if dont_filter:
                return filtered_copy

            for mhc_class, alleles_dict in filtered_copy.iteritems():
                if mhc_classes and not mhc_classes_d.get(mhc_class, False):
                    _ = filtered_copy.pop(mhc_class)
                    continue
                for allele, samples_dict in alleles_dict.iteritems():
                    if alleles and not alleles_dict.get(allele, False):
                        _ = alleles_dict.pop(allele)
                        continue
                    for sample_id, pval in samples_dict.iteritems():
                        if sample_ids and not sample_ids_d.get(sample_id, False):
                            _ = samples_dict.pop(sample_id)
                            continue
        return filtered_copy

    def get_mhc_alleles_sample_prevalence(self):
        """
        This function calculates the frequency of samples per site.
        :return: A dict of [MHCAllele, frequencies].
        """
        s_freq = {}
        num_of_samples = float(len(self.get_related_entities(of_types=(Sample)).get(Sample, [])))
        for mhc_allele, samples in self.__mhc_alleles_by_alleles.iteritems():
            s_freq[mhc_allele] = len(samples) / num_of_samples

        return s_freq

    def get_mhc_alleles_group_prevalence(self, by_name=False):
        # type: (int, bool) -> dict[dict[float]]
        """
        This function calculates the prevalence of groups per allele.
        :param bool by_name: If set will retrieve the frequencies mapped to names rather than record_id.
        :return: A dict of [MHCAllele, frequencies].
        """
        s_freq = {}
        for mhc_allele, sample_ids in self.__mhc_alleles_by_alleles.iteritems():
            s_freq[mhc_allele] = dict()
            samples = KnownProperty.entity_records_by_ids(sample_ids)
            for sample in samples:
                for group in sample.get_related_entities(of_types=(Group,))[Group]:
                    group_key = group.group_name if by_name else group.record_id
                    s_freq.setdefault(mhc_allele, {}).setdefault(group_key, 0)
                    s_freq[mhc_allele][group_key] += (1.0 / group.group_size)
        return s_freq

    @staticmethod
    def get_mhc_alleles_by_source(source):
        # type: (str) -> MHCAlleles
        """
        This function retrieve a mismatch of a specified source.
        :param source: The source to look for.
        :return: The instance or None if not found.
        """
        mhc_alleles = MHCAlleles.get_all_records(of_types=(MHCAlleles,)).values()
        for mhc_allele in mhc_alleles:
            if mhc_allele.source == source:
                return mhc_allele
        return None




