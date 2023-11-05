__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Imports
# endregion

# region Internal Imports
from Commons.consts import NA
from Commons.data_structs import Site, SampleCoverageData
from Commons.general_functions import deep_copy_dict, convert_params_to_bool_dict

from DataObjects.BaseClasses.KnownProperty import KnownProperty
from DataObjects.BaseClasses.Sample import Sample
# endregion

# =====================constants===================#


# =====================classes=====================#



class CoverageData(KnownProperty):
    """
    This class holds coverage data from a data source_name.
    :ivar __coverage_by_sites: (get - mismatch_sites, to add sites or remove use - add/remove_sites) the sites dict.
        <mismatch type> ><Chr/Transcript name>:<start>:<end>:<strand>:[SampleMMRate(s)]
    :type __coverage_by_sites: C{dict} of C{str} to
    C{dict} of C{int} to C{dict} of C{int} to C{dict} of C{str} to C{dict} of C{str} to C{list} of L{SampleMMRate}
    :ivar __coverage_by_sample: same as __regions_data_by_sites except mismatchs type is the root for the
     mapping
    :type __coverage_by_sites: C{dict} of C{str} to
    C{dict} of C{str} to C{dict} of C{int} to C{dict} of C{int} to  C{list} of L{SampleMMRate}
    :ivar source_name: The source of the data.
    :type source_name: C{str}
    :ivar __sites_to_samples_coverage_data: A dictionary of all sites (L{SampleCoverageData})to samples instances.
    """

    def dump_me(self):
        pass

    def __init__(self, source_name, genome):
        """
        :param str source_name: The (data-wise) origin of the data.
        :param str genome: The genome.
        """
        if hasattr(self, "source_name"):
            return
        self.source_name = source_name
        self.genome = genome
        self.__coverage_by_sites = {}
        self.__coverage_by_sample = {}
        self.__sites_to_samples_coverage_data = {}
        self.__regions_to_samples_coverage_data = {}

    @property
    def regions_to_samples_coverage_data(self):
        # type: () -> dict
        return deep_copy_dict(self.__regions_to_samples_coverage_data)

    def add_site(self, sample, region, start, end, strand, reference, adenosines, cytosines, guanosines, thymines,
                 containing_region_region_name=NA, containing_region_start=0, containing_region_end=0,
                 containing_region_strand=NA):
        # type: (Sample, str, int, int, str, str, int, int, int, int, str, int, int, str) -> None
        """
        This function updates the sites dictionary with new coverage values *overriding prev values*
        :param Sample sample: The instance of the sample
        :param str region: The chromosome (or transcript name)
        :param int start: start position
        :param int end: end position
        :param str strand: the strand of the site
        :param str reference: the reference base in the genome.
        :param int adenosines: the number of adenosines.
        :param int cytosines: the number of cytosines.
        :param int guanosines: the number of guanosines.
        :param int thymines: the number of thymines.
        :param str containing_region_region_name: The chromosome (or transcript name) of the genomic region containing this position.
        :param int containing_region_start: The start position of the genomic region containing this position.
        :param int containing_region_end: The end position of the genomic region containing this position.
        :param str containing_region_strand: The end position of the genomic region containing this position.
        :return: None
        """
        assert isinstance(sample, Sample)
        sample.add_related_entity(self)

        start = int(start)
        end = int(end)
        site = Site(region, start, end, strand)

        adenosines = int(adenosines)
        cytosines = int(cytosines)
        guanosines = int(guanosines)
        thymines = int(thymines)

        reg_start = int(containing_region_start)
        reg_end = int(containing_region_end)
        reg_site = Site(containing_region_region_name, reg_start, reg_end,
                    containing_region_strand)

        sample_coverage_data = SampleCoverageData(sample=sample, site=site, reference=reference,
                                                  adenosines=adenosines, cytosines=cytosines, guanosines=guanosines,
                                                  thymines=thymines, containing_region=reg_site)
        # self.__coverage_by_sites.setdefault(region, {}).setdefault(start, {}).setdefault(end, {}). \
        #     setdefault(strand, {})[sample.record_id] = sample_coverage_data
        # self.__coverage_by_sample.setdefault(sample.record_id, {}).setdefault(region, {}). \
        #     setdefault(start, {}).setdefault(end, {})[strand] = sample_coverage_data
        self.__regions_to_samples_coverage_data.setdefault(reg_site, {}).setdefault(sample.record_id, []).append(sample_coverage_data)

    # def get_sites(self, dont_filter=True, by_sample=False, sample_ids=None, chromosomes=None, starts=None, ends=None,
    #               strands=None):
    #     """
    #     Retrieves the sites.
    #     :param dont_filter: A flag. if set will retrieving a simple copy without filtration. Use to save run-time.
    #     :param by_sample: A flag. if set will return sites by mismathc. Otherwise by location (see class doc)
    #     :param sample_ids: if given will return only mismatches that belong to these samples.
    #     :param chromosomes: if given will return only mismatches that are inside this chromosomes.
    #     :param starts: if given will return only mismatches that match starts.
    #     :param ends: if given will return only mismatches that match ends.
    #     :param strands: if given will return only mismatches that match strands.
    #     :return: the sites relevant in the same structure as the calls dict.
    #     """
    #     sample_ids_d = convert_params_to_bool_dict(sample_ids)
    #     chromosomes_d = convert_params_to_bool_dict(chromosomes)
    #     starts_d = convert_params_to_bool_dict(starts)
    #     ends_d = convert_params_to_bool_dict(ends)
    #     strands_d = convert_params_to_bool_dict(strands)
    #     if by_sample:
    #         filtered_copy = deep_copy_dict(self.__coverage_by_sample)
    #         if dont_filter:
    #             return filtered_copy
    #         for sample_id, sample_dict in self.__coverage_by_sample.iteritems():
    #             if sample_ids and not sample_ids_d.get(sample_id, False):
    #                 _ = filtered_copy.pop(sample_id)
    #                 continue
    #             for chrom, chrom_d in sample_dict.iteritems():
    #                 if chromosomes and not chromosomes_d.get(chrom, False):
    #                     _ = sample_dict.pop(chrom)
    #                     continue
    #                 for st, st_d in chrom_d.iteritems():
    #                     if starts and not starts_d.get(st, False):
    #                         _ = chrom_d.pop(st)
    #                         continue
    #                     for e, e_d in st_d.iteritems():
    #                         if ends and not ends_d.get(e, False):
    #                             _ = st_d.pop(e)
    #                             continue
    #                         for strand, samples_dict in e_d.iteritems():
    #                             if strands and not strands_d.get(strand, False):
    #                                 _ = strand.pop(strand)
    #                                 continue
    #
    #     else:
    #         filtered_copy = deep_copy_dict(self.__coverage_by_sites)
    #         if dont_filter:
    #             return filtered_copy
    #
    #         for chrom, chrom_d in self.__coverage_by_sites.iteritems():
    #             if chromosomes and not chromosomes_d.get(chrom, False):
    #                 _ = filtered_copy.pop(chrom)
    #                 continue
    #             for st, st_d in chrom_d.iteritems():
    #                 if starts and not starts_d.get(st, False):
    #                     _ = chrom_d.pop(st)
    #                     continue
    #                 for e, e_d in st_d.iteritems():
    #                     if ends and not ends_d.get(e, False):
    #                         _ = st_d.pop(e)
    #                         continue
    #                     for strand, sample_dict in st_d.iteritems():
    #                         if strands and not strands_d.get(strand, False):
    #                             _ = st_d.pop(e)
    #                             continue
    #                         for s_id, sample_coverage in sample_dict.iteritems():
    #                             if sample_ids and not sample_ids_d.get(s_id, False):
    #                                 _ = sample_coverage.pop(s_id)
    #
    #     return filtered_copy

    def filter_coverage_by_sites(self, sites):
        # type: (list) -> dict
        """
        This function return a collection of the coverage data filtered by the given sites
        :param list[Site] sites: A list to filter according to.
        :return dict[Site, dict: The filtered site.
        """
        res_dict = dict()
        sites_d = convert_params_to_bool_dict(sites)
        for region, start_d in self.__coverage_by_sites.iteritems():
            for start, end_d in start_d.iteritems():
                for end, strand_d in end_d.iteritems():
                    for strand, cov_dict in strand_d.iteritems():
                        site = Site(region, start, end, strand)
                        if sites_d.get(site, False):
                            res_dict[site] = deep_copy_dict(cov_dict)
        return res_dict
