from itertools import starmap

from Commons.general_functions import convert_params_to_bool_dict, deep_copy_dict

__author__ = 'Hillel'

# =====================imports=====================#
# region Builtin Imports
from collections import namedtuple
from csv import reader
# endregion

# region Internal Imports
from Commons.data_structs import MHCAllele, Site, AffinitiesEnum
from Commons.consts import NA
from DataObjects.BaseClasses.DataRecord import DataRecord
from Outputs.Outputers.CSVOutputer import CSVOutputer


# endregion

# =====================consts======================#

# =====================classes=====================#

# noinspection PyClassHasNoInit
class EpitopesNumber(namedtuple("EpitopesNumber", "peptide_site mhc_allele mismatches_sites num_of_epitopes "
                                                  "peptide_type")):
    """
    This hold the coordinates for a mismatch site.
    :ivar Site peptide_site: The site of the peptide.
    :ivar MHCAllele mhc_allele: The allele the affinity is to.
    :ivar list[Site] mismatches_sites: The sites of mismatches from the "canonical" *DNA* seq.
    :ivar int num_of_epitopes: The predicted number of epitopes.
    :ivar str peptide_type: The type of the peptide (e.g. Random, Real, Edited).
    """
    __slots__ = ()


class EpitopesNumbers(DataRecord):
    """
    This entity holds the data of  peptides affinity to MHC alleles
    :ivar source_name: The name of the source.
    :type source_name: C{str}
    """

    def __init__(self, source_name):
        # type: (str) -> None
        if hasattr(self, "source_name"):
            return
        self.source_name = source_name

        self.__num_of_epitopes_by_site = dict()
        self.__num_of_epitopes_by_mhc = dict()

    def dump_me(self, dump_path):
        pass
        # recs = []
        # headers = ["source_name", "region", "start", "end", "strand", "mhc_class", "locus_name", "allele_name",
        #            "offset",
        #            "ligand_length", "rank", "affinity_rank", "affinity", "binding_core", "mismatches_sites"]
        # for site, alleles_dict in self.__num_of_epitopes_by_site.iteritems():
        #     """:type site Site"""
        #     for allele, aff_rank_d in alleles_dict.iteritems():
        #         """"type allele MHCAllele"""
        #         for aff_rank, peptides in aff_rank_d.iteritems():
        #             for peptide in peptides:
        #                 rec = {"source_name": self.source_name,
        #                        "region": site.region,
        #                        "start": str(site.start),
        #                        "end": str(site.end),
        #                        "strand": str(site.strand),
        #                        "mhc_class": allele.mhc_class,
        #                        "locus_name": allele.locus,
        #                        "allele_name": allele.allele,
        #                        "offset": str(peptide.offset),
        #                        "ligand_length": str(peptide.ligand_length),
        #                        "rank": str(peptide.rank),
        #                        "affinity_rank": peptide.affinity_rank,
        #                        "affinity": str(peptide.affinity),
        #                        "binding_core": peptide.binding_core,
        #                        # "mismatches_sites": NA  # TODO: add...
        #                        }
        #                 recs.append(rec)
        # CSVOutputer().output([dump_path, ], headers, recs)

    @staticmethod
    def load_me(dump_path):
        pass
        # with open(dump_path, "rb") as dp:
        #     instance_data = [line for line in reader(dp)]
        #
        # headers_line = instance_data[0]
        #
        # source_name_i = headers_line.index("source_name")
        # region_i = headers_line.index("region")
        # start_i = headers_line.index("start")
        # end_i = headers_line.index("end")
        # strand_i = headers_line.index("strand")
        # mhc_class_i = headers_line.index("mhc_class")
        # locus_name_i = headers_line.index("locus_name")
        # allele_name_i = headers_line.index("allele_name")
        # offset_i = headers_line.index("offset")
        # ligand_length_i = headers_line.index("ligand_length")
        # rank_i = headers_line.index("rank")
        # affinity_rank_i = headers_line.index("affinity_rank")
        # affinity_i = headers_line.index("affinity")
        # binding_core_i = headers_line.index("binding_core")
        #
        # copy = PeptidesAffinities(instance_data[1][source_name_i])
        # """:type instance PeptidesAffinities"""
        #
        # for line in instance_data[1:]:
        #     region = line[region_i]
        #     start = int(line[start_i])
        #     end = int(line[end_i])
        #     strand = line[strand_i]
        #     mhc_class = line[mhc_class_i]
        #     locus_name = line[locus_name_i]
        #     allele_name = line[allele_name_i]
        #     offset = int(line[offset_i])
        #     ligand_length = int(line[ligand_length_i])
        #     rank = float(line[rank_i])
        #     affinity_rank = line[affinity_rank_i]
        #     affinity = float(line[affinity_i])
        #     binding_core = line[binding_core_i]
        #     copy.add_epitopes_number(region, start, end, strand, mhc_class, locus_name, allele_name, offset, ligand_length,
        #                          rank, affinity_rank, affinity, binding_core)
        #
        # del instance_data
        #
        # return copy

    def add_epitopes_number(self, region, start, end, strand, mhc_class, locus_name, allele_name, num_of_epitopes,
                            peptide_type, mismatches_sites=NA):
        # type: (str, int, int, str, str, str, str, int, str, list) -> None
        """
        This function adds to the affinity dictionary (overriding existing values)
        :param str region: The chromosome (or transcript name)
        :param int start: start position
        :param int end: end position
        :param str strand: the strand of the site
        :param str mhc_class: The class of the MHC.
        :param str locus_name: the locus (name) of the allele.
        :param str allele_name: the full allele name (e.g. A*01:21).
        :param int num_of_epitopes: The predicted number of epitopes.
        :param str peptide_type: The type of the peptide (e.g. Random, Real, Edited).
        :param list mismatches_sites: A list of mismatches site contributing to the current peptide sequence
        :return: None
        """

        mhc_allele = MHCAllele(mhc_class, locus_name, allele_name)
        peptide_site = Site(region, start, end, strand)
        epis_n = EpitopesNumber(peptide_site, mhc_allele, mismatches_sites, num_of_epitopes, peptide_type)
        self.__num_of_epitopes_by_site.setdefault(peptide_site, {}).setdefault(mhc_allele, {})[peptide_type] = epis_n
        self.__num_of_epitopes_by_mhc.setdefault(mhc_allele, {}).setdefault(peptide_site, {})[peptide_type] = epis_n

    def get_epitopes_numbers(self, dont_filter=True, by_allele=False, mhc_alleles=None, sites=None, peptide_types=None):
        """
        Retrieves the number of epitopes.
        :param bool dont_filter: A flag. if set will retrieving a simple copy without filtration. Use to save run-time.
        :param bool by_allele: A flag. if set will return numbers by the MHC alleles otherwise by site.
        :param list[MHCAllele] mhc_alleles: if given will return only epitopes numbers to alleles that match this.
        :param list[Site] sites: if given will return only epitopes numbers to these sites.
        :param list[str] peptide_types: if given will return only epitopes numbers to these types.
        :return: the alleles relevant in the same structure as the calls dict.
        """
        mhc_alleles_d = convert_params_to_bool_dict(mhc_alleles)
        sites_d = convert_params_to_bool_dict(sites)
        peptide_types_d = convert_params_to_bool_dict(peptide_types)
        if by_allele:
            filtered_copy = deep_copy_dict(self.__num_of_epitopes_by_mhc)
            if dont_filter:
                return filtered_copy
            for allele, sites_dict in self.__num_of_epitopes_by_mhc.iteritems():
                for site, epitopes_numbers_d in sites_dict.iteritems():
                    if sites and not mhc_alleles_d.get(site, False):
                        _ = filtered_copy[allele].pop(site)
                        continue
                    for peptide_type in epitopes_numbers_d:
                        if peptide_types and not peptide_types_d.get(peptide_type, False):
                            _ = filtered_copy[allele][site].pop(peptide_type)
                            continue

        else:
            filtered_copy = deep_copy_dict(self.__num_of_epitopes_by_site)
            if dont_filter:
                return filtered_copy

            for site, alleles_dict in self.__num_of_epitopes_by_site.iteritems():
                if sites and not sites_d.get(site, False):
                    _ = filtered_copy.pop(site)
                    continue
                for allele, epitopes_numbers_d in alleles_dict.iteritems():
                    if mhc_alleles and not mhc_alleles_d.get(allele, False):
                        _ = filtered_copy[site].pop(allele)
                    for peptide_type in epitopes_numbers_d:
                        if peptide_types and not peptide_types_d.get(peptide_type, False):
                            _ = filtered_copy[allele][site].pop(peptide_type)
                            continue
        return filtered_copy

