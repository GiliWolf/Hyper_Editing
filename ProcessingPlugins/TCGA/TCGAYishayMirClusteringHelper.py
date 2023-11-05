__author__ = 'Hillel'
# =====================imports=====================#
from datetime import datetime
import os

from ProcessingPlugins.RNAEditingQuantificationClusterer import RNAQuantificationClustererHelper, S_ID
from Commons.help_functions import PARTICIPANT_NUM, SAMPLE_NAME, SAMPLE_TYPE, TISSUE_SOURCE, parse_TCGA_sample_name,\
    TCGA_HEADERS_LIST
from DataConverters.RNAEditingQuantificationConverter import YishayKnownAnalysisParser, PRE_RNA
from DataObjects.BaseClasses.Sample import Sample
from DataObjects.Loci.PositionInRNA import PositionInRNA

# =====================constants===================#

FILTERED = ".editing"
UNFILTERED = "test.txt"

# =====================classes=====================#

class TCGAYishayMirClusteringHelper(RNAQuantificationClustererHelper):
    def __init__(self, filtered_suffix=FILTERED, unfiltered_suffix=UNFILTERED, groups_file=None):
        self.filtered_found = {}
        self.unfiltered_found = {}
        if groups_file:
            self.init_group_data_from_groups_file(groups_file)
        self.filtered_suffix = filtered_suffix
        self.unfiltered_suffix = unfiltered_suffix
        self.filtered_converter = YishayKnownAnalysisParser
        self.unfiltered_converter = None

    def sort_file_type(self, path, filtered, unfiltered):
        sample_params = parse_TCGA_sample_name(path)
        if sample_params is None and self.groups_file_data == {}:
            return
        group,s_id = self.get_group(path, sample_params)

        if path.endswith(self.filtered_suffix):
            if self.filtered_found.get(sample_params[PARTICIPANT_NUM] + sample_params[SAMPLE_TYPE] + group, True):
                filtered.append(path)
                self.filtered_found[sample_params[PARTICIPANT_NUM] + sample_params[SAMPLE_TYPE] + group] = False

    def create_sample(self, path):
        sample_params = parse_TCGA_sample_name(os.path.basename(path))
        group,s_id = self.get_group(path, sample_params)
        if s_id:
            sample_params[S_ID] = s_id
        return Sample(group=group, sample_name=sample_params[SAMPLE_NAME], extra_data=sample_params)

    def get_group(self, path, sample_params):
        return sample_params[TISSUE_SOURCE] + "_" + sample_params[SAMPLE_TYPE],None if self.groups_file_data == {} else \
            self.get_group_from_groups_file(path)

    def filter_loci(self, loci):
        for locus in loci:
            if isinstance(locus, PositionInRNA) and locus.annotation == PRE_RNA:
                yield locus
            else:
                continue

    def get_rna_att(self, locus):
        if isinstance(locus, PositionInRNA):
            return locus.rna_name

    @staticmethod
    def print_atts_help():
        """
        This function returns a help string representing the created L{Sample} params for user's help.
        """
        atts_help = "(extra_data is extracted from TCGA codon in sample name, to see meaning check on TCGA site)\n"
        atts_help += "<att> [=> key(s)]\n\n1. group\n2. sample_name\n"
        atts_help += "3. extra_data =>\n" + "\t[(if groups file is given) sample_source_id (third column in the file)]\n"
        atts_help += "".join(map(lambda s: "\t%s\n" % s, TCGA_HEADERS_LIST))
        return atts_help

