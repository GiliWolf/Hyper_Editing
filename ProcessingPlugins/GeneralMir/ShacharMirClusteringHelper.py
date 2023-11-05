__author__ = 'Hillel'
# =====================imports=====================#
from datetime import datetime
import os

from ProcessingPlugins.RNAEditingQuantificationClusterer import RNAQuantificationClustererHelper, S_ID
from Commons.help_functions import PARTICIPANT_NUM, SAMPLE_NAME, SAMPLE_TYPE, TISSUE_SOURCE, parse_TCGA_sample_name
from DataConverters.RNAEditingQuantificationConverter import ShacharFilteredResultParser, ShacharUnFilteredResultParser,\
    PRE_RNA
from DataObjects.BaseClasses.Sample import Sample
from DataObjects.Loci.PositionInRNA import PositionInRNA

# =====================constants===================#

FILTERED = "editing.txt"
UNFILTERED = "test.txt"

# =====================classes=====================#


class ShacharMirClusteringHelper(RNAQuantificationClustererHelper):
    def __init__(self, filtered_suffix=FILTERED, unfiltered_suffix=UNFILTERED, groups_file=None):
        self.filtered_found = {}
        self.unfiltered_found = {}
        if groups_file:
            self.init_group_data_from_groups_file(groups_file)
        else:
            raise AssertionError("General Mir Clusterer Requires A Group File!")
        self.filtered_suffix = filtered_suffix
        self.unfiltered_suffix = unfiltered_suffix
        self.filtered_converter = ShacharFilteredResultParser
        self.unfiltered_converter = ShacharUnFilteredResultParser


    def sort_file_type(self, path, filtered, unfiltered):
        group,s_id = self.get_group(path)

        if path.endswith(self.filtered_suffix):
            if self.filtered_found.get(s_id + group, True):
                filtered.append(path)
                self.filtered_found[s_id + group] = False

        elif path.endswith(self.unfiltered_suffix):
            if self.unfiltered_found.get(s_id + group, True):
                unfiltered.append(path)
                self.unfiltered_found[s_id + group] = False

    def create_sample(self, path):
        group,s_id = self.get_group(path)
        sample_params = {S_ID: s_id, "group": group}
        return Sample(group=group, sample_name=s_id, extra_data=sample_params)

    def get_group(self, path):
        return self.get_group_from_groups_file(path)

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
        atts_help = "(groups file must be given) \n<att> [=> key(s)]\n" \
                    "\n1. group (second column in the file)\n2. sample_name (third column in the file)\n"
        atts_help += "3. extra_data =>\n" + \
                     "\tsample_source_id (third column in the file)\n" \
                     "\tgroup (second column in the file)"
        return atts_help
