__author__ = 'Hillel'
# =====================imports=====================#
import abc
import importlib
import logging
from _csv import reader


from Commons.data_structs import KeepOnlyEnum

# =====================consts======================#
GROUP_ORDER_SAMPLE_EXTRA_DATA_KEY = "GroupsOrder"

# =====================classes=====================#


class RNAEditingQuantificationClustererHelper(object):
    __metaclass__ = abc.ABCMeta
    """
    This class is a base class for helpers, which are dependant on the sample exact type (e.g. tcga, REDITools).
    The following function should be implemented so that the general code can run.
    Note: Each subclass must have a filtered_converter and unfiltered_converter fields, holding the relevant
    RNAEditingQuantificationConverter converter type, and the args as well.
    """
    unfiltered_converter = None
    unfiltered_converter_kwargs = {}
    groups_file_data = None
    filtered_converter = None
    filtered_converter_kwargs = {}
    filtered_files_group = ""
    unfiltered_files_group = ""
    sum_all_mm_to_coverage = False
    filtered_source_name = "RNAEditingQuantificationFiltered"
    unfiltered_source_name = "RNAEditingQuantificationUnfiltered"
    combined_source_name = "RNAEditingQuantificationCombined"
    keep_filtered_or_unfiltered = KeepOnlyEnum.BOTH

    @abc.abstractmethod
    def get_files_predicates_dict(self, *args, **kwargs):
        """
        This function append the path to either the filtered or unfiltered (depending on implementation)
        :return: None
        """
        pass

    @abc.abstractmethod
    def create_sample(self, path, samples_created=False, *args, **kwargs):
        """
        This function creates a L{Sample} instance depending on implementation.
        :param path: The sample's file path.
        :param bool samples_created: If set will not create samples and will look for them by name.
        :return: The relevant Sample instance.
        :rtype: L{Sample}
        """
        pass

    @staticmethod
    def print_atts_help():
        """
        This function returns a help string representing the created L{Sample} params for user's help.
        """
        raise NotImplementedError()

    def init_group_data_from_groups_file(self, groups_file):
        """
        This function retrun the group of a path if groups file was given.
        :param groups_file: The C{str} path to the groups file.
        :return: The correlating group from the grops file.
        """

        with open(groups_file, 'rb') as gf:
            g_data = [line for line in reader(gf, delimiter="\t")]
        for line in g_data:
            if line:
                self.groups_file_data[line[0]] = (line[1], line[2]) if len(line) == 3 else (line[1], None)


class RNAQuantificationClustererHelperFactory(object):
    """
    This class is a factory for RNAQuantificationClustererHelper instances.
    """
    CLUSTERING_HELPERS = {
     "REDIToolsKnownClusteringHelper": "ProcessingPlugins.MismatchesQuantification.REDITools.REDIToolsKnownClusteringHelper",
     "TCGAShacharMirClusteringHelper": "ProcessingPlugins.MismatchesQuantification.TCGA.TCGAShacharMirClusteringHelper",
     "TCGAYishayMirClusteringHelper": "ProcessingPlugins.MismatchesQuantification.TCGA.TCGAYishayMirClusteringHelper",
     "ShacharMirClusteringHelper": "ProcessingPlugins.MismatchesQuantification.GeneralMir.ShacharMirClusteringHelper"
    }

    @staticmethod
    def get_clustering_helper(clustering_helper, clustering_helper_kwargs):
        # type: (str, dict) -> RNAEditingQuantificationClustererHelper
        """
        The factory function.
        :param clustering_helper: The name of the clustering helper.
        :param clustering_helper_kwargs: any kwargs for its init.
        :return: An instance of RNAQuantificationClustererHelper.
        :rtype L{RNAQuantificationClustererHelper}
        """
        try:
            cl_hm = importlib.import_module(RNAQuantificationClustererHelperFactory.CLUSTERING_HELPERS[clustering_helper])
        except KeyError:
            logging.error("Tried To Initialize A Clustering Helper That Doesn't Exist! (%s)" % clustering_helper)
            return None
        except ImportError:
            logging.exception("Couldn't Import Clustering Helper! (%s)" % clustering_helper)

        clustering_helper = getattr(cl_hm, clustering_helper)(**clustering_helper_kwargs)
        return clustering_helper
