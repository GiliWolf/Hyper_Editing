__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Imports
import os
# endregion


# region Internal Imports
from Commons.consts import UNKNOWN

from DataConverters.MismatchesQuantification.RNAEditingQuantificationConverter import REDIToolsKnownParser

from DataObjects.BaseClasses.Sample import Sample
from DataObjects.KnownProperties.Group import Group

from ProcessingPlugins.MismatchesQuantification.RNAEditingQuantificationClustererHelper import RNAEditingQuantificationClustererHelper, \
    GROUP_ORDER_SAMPLE_EXTRA_DATA_KEY
# endregion
# =====================constants===================#

REDITOOLKNOWN_PIPELINE_SUFFIX = ".REDItoolKnown.out.tab"

REDIToolsKnownClusteringHelper_FILES_GROUP = "All"
REDIToolsKnownClusteringHelper_FILTERED_SOURCE_NAME = "REDIToolKnownFiltered"
REDIToolsKnownClusteringHelper_COMBINED_SOURCE_NAME = "REDIToolKnown"

# =====================classes=====================#


class REDIToolsKnownClusteringHelper(RNAEditingQuantificationClustererHelper):
    def __init__(self, files_suffix=REDITOOLKNOWN_PIPELINE_SUFFIX, groups_order=[0,],
                 filtered_source_name=REDIToolsKnownClusteringHelper_FILTERED_SOURCE_NAME,
                 combined_source_name=REDIToolsKnownClusteringHelper_COMBINED_SOURCE_NAME, filtered_converter_kwargs={}):
        self.files_suffix = files_suffix
        self.filtered_converter = REDIToolsKnownParser
        self.filtered_converter_kwargs = eval(filtered_converter_kwargs) if filtered_converter_kwargs else dict()
        self.unfiltered_converter = REDIToolsKnownParser
        self.unfiltered_converter_kwargs = {}
        self.filtered_files_group = REDIToolsKnownClusteringHelper_FILES_GROUP
        self.groups_order = groups_order
        self.filtered_source_name = filtered_source_name
        self.combined_source_name = combined_source_name

    @staticmethod
    def print_atts_help():
        """
        This function returns a help string representing the created L{Sample} params for user's help.
        """
        atts_help = "(Info folders are folders relevant to grouping of samples," \
                    " given to the clustering helper at init as indexes of folders in the samples' paths.)\n"
        atts_help += "<att> [=> key(s)]\n\n1. group\n2. sample_name\n"
        atts_help += "3. extra_data =>\n" + "\t".join(["main_group (usually the parent folder)\n",
                                                       "sub_groups (all following 'info' folders)\n",
                                                       "[(if groups file is given) sample_source_id (third column in the file)]\n"])
        return ""#atts_help

    def get_files_predicates_dict(self):
        return {self.filtered_files_group: lambda x: x.endswith(self.files_suffix)}

    def create_sample(self, path):
        sample_name = os.path.basename(path).replace(self.files_suffix, "")
        sample_o_dir = os.path.dirname(path)
        sample = Sample(sample_name=sample_name, sample_path=sample_o_dir + '*')
        sample.extra_data[GROUP_ORDER_SAMPLE_EXTRA_DATA_KEY] = self.groups_order
        if not sample.get_related_entities(of_types=(Group,)):
            group = Group(group_name=UNKNOWN)
            group.add_related_entity(sample)
            sample.add_related_entity(group)

        return sample