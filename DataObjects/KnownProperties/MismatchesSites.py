__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Imports
from collections import namedtuple
# endregion

# region Internal Imports
from Commons.consts import MismatchesAndRefsEnum, NA
from Commons.data_structs import Site, KeepOnlyEnum
from Commons.general_functions import deep_copy_dict, convert_params_to_bool_dict, add_na

from DataObjects.BaseClasses.KnownProperty import KnownProperty
from DataObjects.BaseClasses.Sample import Sample
from DataObjects.KnownProperties.Group import Group

# endregion

# =====================constants===================#
# for saving runtime
POSSIBLE_MISMATCHES = set(MismatchesAndRefsEnum.MISMATCH_TYPE_PARSE.values() + [MismatchesAndRefsEnum.CANONICAL, ])

ANNOVAR_FORMAT = "annovar"
BED_FORMAT = "bed"
AVAILABLE_FORMATS = [ANNOVAR_FORMAT, BED_FORMAT]


# =====================classes=====================#

# noinspection PyClassHasNoInit
class SampleMMRate(namedtuple("SampleMMRate", "sample num_of_mm_reads num_of_canonical_reads num_of_total_reads "
                                              "ratio mismatch_type site")):
    """
    This holds the mismatch rate for a sample for a sites.
    :ivar Sample sample: the instance of the L{Sample} this is relating to.
    :ivar int num_of_mm_reads: the number of reads supporting the mismatch.
    :ivar int num_of_canonical_reads: the number of reads supporting the canonical sequence.
    :ivar int num_of_total_reads: the total number of reads.
    :ivar float ratio: the ratio of reads supporting the mismatch.
    :ivar dict extra_data: extra_data about this position.

    """
    __slots__ = ()


# noinspection PyClassHasNoInit
class GroupMMRate(namedtuple("GroupMMRate", "group num_of_mm_reads num_of_canonical_reads num_of_total_reads "
                                            "ratio site")):
    """
    This holds the mismatch rate for a sample for a sites.
    :ivar str group: the instance of the L{Group} this is relating to.
    :ivar int num_of_mm_reads: the number of reads supporting the mismatch.
    :ivar int num_of_canonical_reads: the number of reads supporting the canonical sequence.
    :ivar int num_of_total_reads: the total number of reads.
    :ivar float ratio: the ratio of reads supporting the mismatch.
    """
    __slots__ = ()


class MismatchesSites(KnownProperty):
    """
    This class holds the mismatches from a data source.
    :ivar __mismatch_sites_by_sites: (get - mismatch_sites, to add sites or remove use - add/remove_sites) the sites dict.
        <mismatch type> ><Chr/Transcript name>:<start>:<end>:<strand>:[SampleMMRate(s)]
    :type __mismatch_sites_by_sites: C{dict} of C{str} to
    C{dict} of C{int} to C{dict} of C{int} to C{dict} of C{str} to C{dict} of C{str} to C{list} of L{SampleMMRate}
    :ivar __mismatch_sites_by_mismatch: same as __regions_data_by_sites except mismatchs type is the root for the
     mapping
    :type __mismatch_sites_by_sites: C{dict} of C{str} to
    C{dict} of C{str} to C{dict} of C{int} to C{dict} of C{int} to  C{list} of L{SampleMMRate}
    :ivar source: The sourec of the data.
    :type source: C{str}
    :ivar __sites_to_samples: A dictionary of all sites (L{MismatchSite})to samples instances.
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
        self.__mismatch_sites_by_sites = {}
        self.__mismatch_sites_by_mismatch = {}
        self.__sites_to_samples = {}

    @property
    def sites_to_samples(self):
        # type: () -> dict
        return deep_copy_dict(self.__sites_to_samples)

    def add_site(self, mismatch_type, sample, region, start, end, strand, num_of_mismatched_reads,
                 num_of_canonical_reads=NA, ratio=None, total_coverage=NA):
        """
        This function updates the sites dictionary with new mismatch sites *overriding prev values*
        :param str mismatch_type: The type of the mismatch. *only from from MismatchesTypesEnum!*
        :param Sample sample: The instance of the sample
        :param str region: The chromosome (or transcript name)
        :param int start: start position
        :param int end: end position
        :param str strand: the strand of the site
        :param int num_of_mismatched_reads: how many reads support the mismatch
        :param int num_of_canonical_reads: how many reads support the canonical sequence.
        :param float ratio: The editing ratio. If not provided will automatically be calculated.
        :param int total_coverage: how many reads in total cover the site.
        :return: None
        """

        assert isinstance(sample, Sample)
        assert mismatch_type in POSSIBLE_MISMATCHES
        sample.add_related_entity(self)

        start = int(start)
        end = int(end)
        mm_reads = int(num_of_mismatched_reads)
        canonical_reads = int(num_of_canonical_reads) if num_of_canonical_reads != NA else NA
        if total_coverage != NA:
            total = total_coverage
        else:
            total = mm_reads + canonical_reads if num_of_canonical_reads != NA else NA
        if None is ratio:
            rate = float(mm_reads) / total if total not in (NA, 0) else NA
        else:
            rate = ratio
        site = Site(region, start, end, strand)
        sample_rate = SampleMMRate(sample, mm_reads, canonical_reads, total, rate, mismatch_type, site)
        self.__mismatch_sites_by_sites.setdefault(region, {}).setdefault(start, {}).setdefault(end, {}). \
            setdefault(strand, {}).setdefault(mismatch_type, {})[sample.record_id] = sample_rate
        self.__mismatch_sites_by_mismatch.setdefault(mismatch_type, {}).setdefault(region, {}). \
            setdefault(start, {}).setdefault(end, {}).setdefault(strand, {})[sample.record_id] = sample_rate
        self.__sites_to_samples.setdefault(site, {}).setdefault(mismatch_type, {})[sample.record_id] = sample

    def get_sites(self, dont_filter=True, by_mismatch=False, mismatch_types=None, sample_ids=None,
                  chromosomes=None, starts=None, ends=None, strands=None):
        """
        Retrieves the sites.
        :param dont_filter: A flag. if set will retrieving a simple copy without filtration. Use to save run-time.
        :param by_mismatch: A flag. if set will return sites by mismathc. Otherwise by location (see class doc)
        :param mismatch_types: if given will return only mismatches that match this mismatch type.
        :param sample_ids: if given will return only mismatches that belong to these samples.
        :param chromosomes: if given will return only mismatches that are inside this chromosomes.
        :param starts: if given will return only mismatches that match starts.
        :param ends: if given will return only mismatches that match ends.
        :param strands: if given will return only mismatches that match strands.
        :return: the sites relevant in the same structure as the calls dict.
        """
        mismatch_types_d = convert_params_to_bool_dict(mismatch_types)
        sample_ids_d = convert_params_to_bool_dict(sample_ids)
        chromosomes_d = convert_params_to_bool_dict(chromosomes)
        starts_d = convert_params_to_bool_dict(starts)
        ends_d = convert_params_to_bool_dict(ends)
        strands_d = convert_params_to_bool_dict(strands)
        if by_mismatch:
            filtered_copy = deep_copy_dict(self.__mismatch_sites_by_mismatch)
            if dont_filter:
                return filtered_copy
            for mm, mm_dict in self.__mismatch_sites_by_mismatch.iteritems():
                if mismatch_types and not mismatch_types_d.get(mm, False):
                    _ = filtered_copy.pop(mm)
                    continue
                for chrom, chrom_d in mm_dict.iteritems():
                    if chromosomes and not chromosomes_d.get(chrom, False):
                        _ = mm_dict.pop(chrom)
                        continue
                    for st, st_d in chrom_d.iteritems():
                        if starts and not starts_d.get(st, False):
                            _ = chrom_d.pop(st)
                            continue
                        for e, e_d in st_d.iteritems():
                            if ends and not ends_d.get(e, False):
                                _ = st_d.pop(e)
                                continue
                            for strand, samples_dict in e_d.iteritems():
                                if strands and not strands_d.get(strand, False):
                                    _ = strand.pop(strand)
                                    continue
                                for s_id, sample_rate in samples_dict.iteritems():
                                    if sample_ids and not sample_ids_d.get(s_id, False):
                                        _ = samples_dict.pop(s_id)
                                        continue

        else:
            filtered_copy = deep_copy_dict(self.__mismatch_sites_by_sites)
            if dont_filter:
                return filtered_copy

            for chrom, chrom_d in self.__mismatch_sites_by_sites.iteritems():
                if chromosomes and not chromosomes_d.get(chrom, False):
                    _ = filtered_copy.pop(chrom)
                    continue
                for st, st_d in chrom_d.iteritems():
                    if starts and not starts_d.get(st, False):
                        _ = chrom_d.pop(st)
                        continue
                    for e, e_d in st_d.iteritems():
                        if ends and not ends_d.get(e, False):
                            _ = st_d.pop(e)
                            continue
                        for strand, mm_dict in st_d.iteritems():
                            if strands and not strands_d.get(strand, False):
                                _ = st_d.pop(e)
                                continue
                            for mm, sample_rate_d in mm_dict.iteritems():
                                if mismatch_types and not mismatch_types_d.get(mm, False):
                                    _ = e_d.pop(mm)
                                    continue
                                for s_id, sample_rate in sample_rate_d.iteritems():
                                    if sample_ids and not sample_ids_d.get(mm, False):
                                        _ = sample_rate_d.pop(s_id)

        return filtered_copy

    def merge_with(self, mismatches_sites_by_sites, combine_rates=KeepOnlyEnum.BOTH, na_value=None):
        # type: (dict, KeepOnlyEnum, int) -> None
        """
        This function merges data from MismatchesSites.get_sites with the current instance.
        make sure the data is by sites!
        :param dict[Site, dict[str,SampleMMRate]] mismatches_sites_by_sites: A dict of site to mismatches.
        :param KeepOnlyEnum combine_rates: will use the choice of keep this (LEFT), other (RIGHT) or combine (BOTH) operands
        :param int na_value: If provided will use this value for atts with value=NA.
        :return: None
        """
        # TODO: add specific exceptions and logging
        for site, mismatch_dict in mismatches_sites_by_sites.iteritems():
            """:type site Site"""
            region = site.region
            start = site.start
            end = site.end
            strand = site.strand
            for mm, samples_dict in mismatch_dict.iteritems():
                for sample_id, other_sample_mm_rate in samples_dict.iteritems():
                    if self.__sites_to_samples.get(site, {}).get(mm, {}).get(sample_id, None):  # mismatch exist here
                        if combine_rates == KeepOnlyEnum.BOTH:
                            num_of_mismatched_reads = other_sample_mm_rate.num_of_mm_reads
                            num_of_canonical_reads = other_sample_mm_rate.num_of_canonical_reads
                            sample_mm_rates = self.__mismatch_sites_by_sites[site.region][site.start][site.end][
                                site.strand][mm][sample_id]
                            """:type sample_mm_rates SampleMMRate"""
                            num_of_mismatched_reads = add_na(num_of_mismatched_reads, sample_mm_rates.num_of_mm_reads,
                                                             na_value)
                            num_of_canonical_reads = add_na(num_of_canonical_reads,
                                                            other_sample_mm_rate.num_of_canonical_reads, na_value)
                            ratio = None  # will be re-calculated at add_site
                        elif combine_rates == KeepOnlyEnum.LEFT:
                            continue  # no need to update site
                        else:
                            num_of_mismatched_reads = other_sample_mm_rate.num_of_mm_reads
                            num_of_canonical_reads = other_sample_mm_rate.num_of_canonical_reads
                            ratio = other_sample_mm_rate.ratio
                    else:  # mismatch doesn't exist.
                        num_of_mismatched_reads = other_sample_mm_rate.num_of_mm_reads
                        num_of_canonical_reads = other_sample_mm_rate.num_of_canonical_reads
                        ratio = other_sample_mm_rate.ratio

                    self.add_site(mismatch_type=mm,
                                  sample=other_sample_mm_rate.sample,
                                  region=region,
                                  start=start,
                                  end=end,
                                  strand=strand,
                                  num_of_mismatched_reads=num_of_mismatched_reads,
                                  num_of_canonical_reads=num_of_canonical_reads,
                                  ratio=ratio
                                  )

    def filter_mismatches_by_sites(self, sites):
        # type: (list) -> dict
        """
        This function return a collection of the mismatches filtered by the given sites
        :param list[MismatchSite] sites: A list to filter according to.
        :return dict[MismatchSite, dict: The filtered site.
        """
        res_dict = dict()
        sites_d = convert_params_to_bool_dict(sites)
        for region, start_d in self.__mismatch_sites_by_sites.iteritems():
            for start, end_d in start_d.iteritems():
                for end, strand_d in end_d.iteritems():
                    for strand, mm_dict in strand_d.iteritems():
                        site = Site(region, start, end, strand)
                        if sites_d.get(site, False):
                            res_dict[site] = deep_copy_dict(mm_dict)
        return res_dict

    def contained(self, mismatch_types=None, sample_ids=None, chromosome=None,
                  start=None, end=None, strands=None):
        """
        Returns True if the site coordinate sent is contained in this instance's collection of sites.
        *Note* params left as default will be ignored in the predicate!*
        :param mismatch_types: the mismatch types
        :param sample_ids: the samples ids
        :param chromosome: the chromosome.
        :param start: the start position
        :param end: the end position
        :param strands: the possible strands
        :return: True if found False otherwise.
        :rtype C{bool}
        """
        if mismatch_types or sample_ids:
            run_on = self.get_sites(dont_filter=False, mismatch_types=mismatch_types, sample_ids=sample_ids,
                                    strands=strands)

        if chromosome not in run_on:
            return False
        chrom_d = run_on[chromosome]
        if start in chrom_d:
            start_d = chrom_d[start]
            if end in start_d:
                return True
            else:
                for e in start_d:
                    if e > end:
                        return True
        else:
            for s in chrom_d:
                if s < start:
                    start_d = chrom_d[s]
                    if end in start_d:
                        return True
                    else:
                        for e in start_d:
                            if e > end:
                                return True
        return False

    def get_mismatches_sample_prevalence(self, coverage_threshold=1):
        """
        This function calculates the frequency of samples per site.
        :param int coverage_threshold: A threshold for counting an event as existing.
        :return: A dict of frequencies.
        """
        s_freq = {}
        num_of_sample = float(len(self.get_all_records(of_types=(Sample,))))
        for site, mismatch_type_d in self.__sites_to_samples.iteritems():
            for mismatch_type, samples in mismatch_type_d.iteritems():
                for sample_id in samples:
                    s_freq.setdefault(site, {}).setdefault(mismatch_type, 0)
                    coverage = self.__mismatch_sites_by_sites[site.region][site.start][site.end][site.strand][
                        mismatch_type][sample_id].num_of_total_reads
                    rate = self.__mismatch_sites_by_sites[site.region][site.start][site.end][site.strand][
                        mismatch_type][sample_id].ratio
                    if coverage < coverage_threshold:
                        continue
                    if rate == 0:
                        continue
                    s_freq[site][mismatch_type] += 1 / num_of_sample
        return s_freq

    def get_mismatches_group_prevalence(self, coverage_threshold=1, by_name=False):
        # type: (int, bool) -> dict
        """
        This function calculates the frequency of groups per site.
        :param int coverage_threshold: A threshold for counting an event as existing.
        :param bool by_name: If set will retrieve the frequencies mapped to names rather than record_id.
        :return: A dict of frequencies.
        """
        s_freq = {}

        for site, mismatch_type_d in self.__sites_to_samples.iteritems():
            s_freq.setdefault(site, {})
            for mismatch_type, samples in mismatch_type_d.iteritems():
                for sample_id, sample in samples.iteritems():
                    coverage = self.__mismatch_sites_by_sites[site.region][site.start][site.end][site.strand][
                        mismatch_type][sample_id].num_of_total_reads
                    rate = self.__mismatch_sites_by_sites[site.region][site.start][site.end][site.strand][
                        mismatch_type][sample_id].ratio
                    if coverage < coverage_threshold:
                        continue
                    if rate == 0:
                        continue
                    for group in sample.get_related_entities(of_types=(Group,))[Group]:
                        group_key = group.group_name if by_name else group.record_id
                        s_freq.setdefault(site, {}).setdefault(mismatch_type, {}).setdefault(group_key, 0)
                        s_freq[site][mismatch_type][group_key] += (1.0 / group.group_size)

        return s_freq

    def get_mismatches_group_rates(self, by_name=False):
        """
        This function calculates the frequency of groups per site.
        :param bool by_name: If set will retrieve the frequencies mapped to names rather than record_id.
        :return: A dict of GroupMMRate
        """
        s_freq = {}
        for site, mismatch_type_d in self.__sites_to_samples.iteritems():
            s_freq.setdefault(site, {})
            for mismatch_type, samples in mismatch_type_d.iteritems():
                for sample_id, sample in samples.iteritems():
                    sample_mm_rate = self.__mismatch_sites_by_sites[site.region][site.start][site.end][site.strand][
                        mismatch_type][sample_id]
                    if sample_mm_rate.num_of_total_reads == 0:
                        continue
                    total_reads = sample_mm_rate.num_of_total_reads if sample_mm_rate.num_of_total_reads != NA else 0
                    canonical_reads = sample_mm_rate.num_of_canonical_reads if sample_mm_rate.num_of_canonical_reads != \
                                                                               NA else 0
                    mismatched_reads = sample_mm_rate.num_of_mm_reads

                    for group in sample.get_related_entities(of_types=(Group,))[Group]:
                        group_key = group.group_name if by_name else group.record_id
                        s_freq.setdefault(site, {}).setdefault(mismatch_type, {}).setdefault(group_key,
                                                                                             GroupMMRate(group, 0, 0, 0,
                                                                                                         0.0, site))
                        tmp = s_freq[site][mismatch_type][group_key]
                        g_total_reads = tmp.num_of_total_reads + total_reads
                        g_canonical_reads = tmp.num_of_canonical_reads + canonical_reads
                        g_mismatched_reads = tmp.num_of_mm_reads + mismatched_reads
                        ratio = float(g_mismatched_reads) / g_total_reads if g_total_reads != 0 else NA
                        s_freq[site][mismatch_type][group_key] = GroupMMRate(group, g_mismatched_reads,
                                                                             g_canonical_reads, g_total_reads,
                                                                             ratio, site)
        return s_freq

    def get_site_data(self, region, start, end, strand, mismatch_type, sample_id, predicate=False):
        # type: (str, int, int, str, bool) -> object
        """
        This function retrieves data of a specific mismatch site for a sample.
        :param str region: The chromosome (or transcript name)
        :param int start: start position
        :param int end: end position
        :param str strand: the strand of the site
        :param str mismatch_type: the mismatch type
        :param str sample_id: the sample record id.
        :param bool predicate: If set will just return True or False
        :return object: The region data if found, else - None
        """
        try:
            data = self.__mismatch_sites_by_sites[region][start][end][strand][mismatch_type][sample_id]
            if predicate:
                data = True
        except KeyError:
            data = SampleMMRate(MismatchesSites.entity_record_by_id(sample_id), 0, 0, 0, 0.0, mismatch_type,
                                Site(region, start, end, strand))
            if predicate:
                data = False

        return data
