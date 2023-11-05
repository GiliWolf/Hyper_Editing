__author__ = 'Hillel'
# =====================imports=====================#
# region builtin Imports
import argparse
import logging
import os
import operator
from datetime import datetime

# endregion

# region Internal Imports
if __name__ == "__main__":
    from sys import path

    path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from Commons.consts import ALL, NO_OUTPUT_DETECTED
from Commons.general_functions import init_logging_dict, get_paths, add_get_paths_function_to_argparser, \
    load_groups_and_samples_from_file, add_groups_file_to_argparser

from DataConverters.EpitopeAnalysis.OptitypeConverter import OptitypeConverter, OPTITYPE_FILENAME_RE

from ProcessingPlugins.WhiteBoxPlugin import WhiteBoxPlugin
from ProcessingPlugins.EpitopeAnalysis.EpitopeAnalysisConsts import OUTPUT_FILENAME

from Outputs.Outputers.CSVOutputer import CSVOutputer
from Outputs import mhc_alleles_functions

# endregion


# =====================constants===================#
OPTITYPE_TOOL_NAME = "Optitype"

OPTITYPE_PREVALENCE_OUTPUT_FILENAME = OUTPUT_FILENAME % dict(prediction_tool=OPTITYPE_TOOL_NAME,
                                                             summary_type="AllelesPrevalence")
OPTITYPE_SAMPLES_OUTPUT_FILENAME = OUTPUT_FILENAME % dict(prediction_tool=OPTITYPE_TOOL_NAME,
                                                          summary_type="SamplesAlleles")

# ----logging consts-----
OPTITYPE_LOG_FILE = r'%s_optitype_analyser.log'


# =====================functions=====================#

OPTITYPE_PREDICATE = lambda x: OPTITYPE_FILENAME_RE.search(x)


# =====================classes=====================#

class OptitypeAnalyser(WhiteBoxPlugin):
    def process_input(self, root_path, output_path=None, include_paths=[], include_paths_operator=operator.or_,
                      exclude_paths=[], exclude_paths_operator=operator.and_, recursion_depth=100, follow_links=False,
                      sample_name_path_i=-1, top_group_name_path_i=-2, groups_file=None):
        # type: (dict, int, int, str) -> list
        """
        This function converts outputs from Optitype and loads the connections
        :param str root_path: the root path to start the dir.
        :param list include_paths: paths that must appear if the path to be included
        :param operator include_paths_operator: The operator used when checking for include_paths
        :param list exclude_paths: paths that mustn't appear if the path to be included
        :param operator exclude_paths_operator: The operator used when checking for exclude_paths
        :param int recursion_depth: The recursion depth.
        :param bool follow_links: If set will follow links.
        :param str output_path: The dir for the summaries output. *NOTE* if None will not generate output.
        :param sample_name_path_i: The index of the sample name in the path (negative index is recommended)
        :param top_group_name_path_i:  The index of the top group name in the path (negative index is recommended)
        :param str groups_file: The path to the groups and samples file.
        :return: list of outputs paths
        """
        if groups_file:
            load_groups_and_samples_from_file(groups_file=groups_file)

        predicates_dict = {ALL: OPTITYPE_PREDICATE}  # update the predicate
        inputs = get_paths(root_path=root_path,
                           must_include_paths=include_paths,
                           must_include_operator=include_paths_operator,
                           exclude_paths=exclude_paths,
                           exclude_operator=exclude_paths_operator,
                           follow_links=follow_links,
                           recursion_depth=recursion_depth,
                           predicates_dict=predicates_dict
                           )
        if [] is inputs:
            logging.warn(NO_OUTPUT_DETECTED % "OptitypeAnalyser")
        samples = OptitypeConverter.convert_optitype_output(file_names=inputs[ALL], sample_name_path_i=sample_name_path_i,
                                                            top_group_name_path_i=top_group_name_path_i)
        if output_path:
            prevalence_recs, prevalence_headers_d = \
                mhc_alleles_functions.mhc_alleles_prevalence_summery(OptitypeConverter.OPTITYPE_SOURCE_NAME)
            sample_alleles_recs, sample_alleles_headers_d = \
                mhc_alleles_functions.mhc_alleles_per_sample_summery(OptitypeConverter.OPTITYPE_SOURCE_NAME)

            outputer = CSVOutputer()
            outputer.output(paths=[os.path.join(output_path, OPTITYPE_PREVALENCE_OUTPUT_FILENAME)],
                            headers=prevalence_headers_d[ALL], records=prevalence_recs)
            outputer.output(paths=[os.path.join(output_path, OPTITYPE_SAMPLES_OUTPUT_FILENAME)],
                            headers=sample_alleles_headers_d[ALL], records=sample_alleles_recs)


if __name__ == '__main__':
    desc = 'A Script for extracting the predicted hla alleles from optitype output.' \
           ' Summaries are written into a given dir.'


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='OptitypeAnalyser', description=desc, formatter_class=MyFormatter)

    parser.add_argument('-o', '--output_dir', metavar="output dir", dest='output_dir', nargs='?', required=True,
                        help='Outputer dir for the output')
    parser.add_argument('-sni', '--sample_name_path_i', metavar="sample name path index", dest='sample_name_path_i',
                        nargs='?', required=True, type=int, default=None,
                        help="The index ofthe sampe name in the path."
                             " e.g. /drive/folder/group/sample_name/file, sample name i will be either -2 or 3")
    parser.add_argument('-gni', '--group_name_path_i', metavar="group name path index", dest='group_name_path_i',
                        nargs='?', required=True, type=int, default=None,
                        help="The index ofthe group name in the path."
                             " e.g. /drive/folder/group/sample_name/file, group name i will be either -3 or 2")
    add_get_paths_function_to_argparser(parser)
    add_groups_file_to_argparser(parser)
    args = parser.parse_args()

    init_logging_dict(log_file=os.path.join(args.output_dir, OPTITYPE_LOG_FILE % datetime.today().isoformat()))

    analyser = OptitypeAnalyser()

    analyser.process_input(root_path=args.root_dir,
                           output_path=args.output_dir,
                           include_paths=args.include_prefixes,
                           include_paths_operator=args.include_operator,
                           exclude_paths=args.exclude_prefixes,
                           exclude_paths_operator=args.exclude_operator,
                           recursion_depth=args.recursion_depth,
                           follow_links=args.follow_links,
                           sample_name_path_i=args.sample_name_path_i,
                           top_group_name_path_i=args.group_name_path_i,
                           groups_file=args.groups_file)
