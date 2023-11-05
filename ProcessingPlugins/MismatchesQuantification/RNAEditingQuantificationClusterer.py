__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Imports
import argparse
import logging
import operator
import os
from csv import reader
from datetime import datetime

# endregion

# region Internal Imports
if __name__ == "__main__":
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from ProcessingPlugins.MismatchesQuantification.RNAEditingQuantificationClustererHelper import \
    RNAQuantificationClustererHelperFactory, RNAEditingQuantificationClustererHelper
from ProcessingPlugins.WhiteBoxPlugin import WhiteBoxPlugin

from Outputs.Outputers.CSVOutputer import CSVOutputer
from Outputs.mismatches_sites_functions import *

from Commons.general_functions import add_get_paths_function_to_argparser, get_paths, load_groups_and_samples_from_file, \
    SORT_SAMPLES_FUNCTIONS_DICT, SORT_BY_GROUPS_DEC, SORT_BY_GROUPS_INC, get_sorted_unique_sample_name, \
    convert_params_to_bool_dict, init_logging_dict, add_groups_file_to_argparser, convert_args_to_dict
from Commons.data_structs import Site

# endregion

# =====================constants===================#
LOG_FILE = "%s/%sRNAEditingQuantification.log"
PREKNOWN_INTERESTING_SITES_FILE_SEP = "\t"

KNOWN_SITE_DOESNT_EXIST_WARN = "The Given PreKnown Sites %s Doesn't Exist!"

FILTERED_WIDTH_ALL_FORMAT = r"MismatchesRatesAbove%sLooseCutoff.csv"
FILTERED_WIDTH_STRICT_FORMAT = r"MismatchesRatesAbove%sStrictCutoff.csv"
FILTERED_WIDTH_ALL_FORMAT_GROUPS = r"GroupsMismatchesRatesAbove%sLooseCutoff.csv"
FILTERED_WIDTH_STRICT_FORMAT_GROUPS = r"GroupsMismatchesRatesAbove%sStrictCutoff.csv"

BED_CHR_I = 0
BED_START_I = 1
BED_END_I = 2
BED_STRND_I = -1
BED_NAME_I = -3
BED_SCORE_I = -2

FILTERED_TYPE = "Filtered"
UNFILTERED_TYPE = "Unfiltered"
# =====================functions=====================#
def get_coupled_atts(atts):
    if not isinstance(atts, (list, tuple)):
        atts = [atts]
    val = []
    for att in atts:
        if not isinstance(att, (list, tuple)):
            att = att
            key = None
        else:
            key = att[1]
            att = att[0]
        val.append((att, key))
    return val

def bind_samples(samples, binder_atts):
    binding_dict = {}

    binder_atts = get_coupled_atts(binder_atts)
    for sample in samples:
        binding_val = []
        for binder_att in binder_atts:
            tmp = sample.__getattribute__(binder_att[0])
            if binder_att[1]:
                tmp = tmp.get(binder_att[1])
            if None is tmp:
                continue
            binding_val.append(tmp)
        binding_val = tuple(binding_val)
        if binding_val not in binding_dict:
            binding_dict[binding_val] = []

        binding_dict[binding_val].append(sample)
    for samples in binding_dict.itervalues():
        Sample.bind_samples([sample.record_id for sample in samples])

    return binding_dict


# =====================classes=====================#


class RNAQuantificationClusterer(WhiteBoxPlugin):
    """
    This plugin process outputs that contain quantification of a mismatch\editing per locus for several samples.
    """

    csv_output = CSVOutputer()

    @staticmethod
    def convert_dir(clustering_helper, run_on_paths, preknown_sites_file=None):
        # type: (RNAEditingQuantificationClustererHelper, dict, str) -> list
        """

        :param RNAEditingQuantificationClustererHelper clustering_helper: The clustering helper to use
        :param dict[str, list[str]] run_on_paths: The paths to convert.
        :param str preknown_sites_file: A pre-known interesting sites containing file.
        :return list[Sample]: The samples created.
        """
        filtered_paths = run_on_paths.get(clustering_helper.filtered_files_group, [])
        unfiltered_paths = run_on_paths.get(clustering_helper.unfiltered_files_group, [])

        filtered_samples = RNAQuantificationClusterer.covert_and_load_input_files(clustering_helper=clustering_helper,
                                                                                  input_paths=filtered_paths,
                                                                                  filter_type=FILTERED_TYPE)
        interesting_loci = RNAQuantificationClusterer.find_interesting_sites(clustering_helper, preknown_sites_file)
        interesting_loci_dict = convert_params_to_bool_dict(interesting_loci)

        combined_mismatches_sites = MismatchesSites(clustering_helper.combined_source_name)
        filtered_mismatches_sites = MismatchesSites(clustering_helper.filtered_source_name)
        combined_mismatches_sites.merge_with(filtered_mismatches_sites.filter_mismatches_by_sites(interesting_loci))

        unfiltered_samples = []
        if unfiltered_paths:
            unfiltered_samples = RNAQuantificationClusterer.covert_and_load_input_files(
                clustering_helper=clustering_helper,
                input_paths=unfiltered_paths,
                filter_type=UNFILTERED_TYPE,
                keep_only=interesting_loci_dict)
            unfiltered_mismatches_sites = MismatchesSites(clustering_helper.unfiltered_source_name)
            combine = clustering_helper.keep_filtered_or_unfiltered
            combined_mismatches_sites.merge_with(unfiltered_mismatches_sites.filter_mismatches_by_sites(interesting_loci),
                                             combine_rates=combine, na_value=0)

        return set(filtered_samples + unfiltered_samples)

    @staticmethod
    def covert_and_load_input_files(clustering_helper, input_paths, filter_type, keep_only=None):
        # type: (RNAEditingQuantificationClustererHelper, list, str, dict) -> list
        """
        This function converts paths using RNAEditingQuantificationClustererHelper.converter
        :param RNAEditingQuantificationClustererHelper clustering_helper: The clustering helper to use.
        :param list[str] input_paths:
        :param str filter_type: To use args for filtered or unfiltered
        :param dict[Site, bool] keep_only: If given will retain only these sites.
        :return list[Sample]: A list of the created samples.
        """
        created_samples = []
        debug_counter = len(input_paths)
        if filter_type == FILTERED_TYPE:
            source_name = clustering_helper.filtered_source_name
            kwargs = clustering_helper.filtered_converter_kwargs
            converter = clustering_helper.filtered_converter
        else:
            source_name = clustering_helper.unfiltered_source_name
            kwargs = clustering_helper.unfiltered_converter_kwargs
            converter = clustering_helper.unfiltered_converter
        for sfile in input_paths:
            logging.debug("Converting (using %s) %s; %s left [%s]\n" % (source_name, sfile, str(debug_counter),
                                                                        datetime.now().strftime("%c")))
            sample = clustering_helper.create_sample(sfile)
            created_samples.append(sample)
            with open(sfile, 'rb') as ifile:
                converter.parse_known_sites_editing_file(
                    data=ifile.read(),
                    sample=sample,
                    source_name=source_name,
                    keep_only=keep_only,
                    **kwargs)
            debug_counter -= 1

        return created_samples

    @staticmethod
    def find_interesting_sites(clustering_helper, min_common_samples, preknown_interesting_sites_file=None):
        # type: (RNAEditingQuantificationClustererHelper, int, str) -> list
        """
        This function find interesting locations (as defined by the prevalence limit) for further usage
        :param RNAEditingQuantificationClustererHelper clustering_helper: The clustering helper to use.
        :param int min_common_samples: The minimal prevalence of a mismatch to retain.
        :param str preknown_interesting_sites_file: A path to an interesting sites file
        :return: a list of the interesting sites
        """
        interesting_sites = list()
        mismatches = MismatchesSites(clustering_helper.filtered_source_name)
        assert isinstance(mismatches, MismatchesSites)
        freqs_per_group = mismatches.get_mismatches_group_prevalence()
        for site, mismatch_type_d in freqs_per_group.iteritems():
            assert isinstance(site, Site)
            for mismatch_type, groups_dict in mismatch_type_d.iteritems():
                for group, frequency in groups_dict.iteritems():
                    if min_common_samples <= 100 * frequency:
                        interesting_sites.append(Site(site.region, site.start, site.end, site.strand))

        if preknown_interesting_sites_file:
            with open(preknown_interesting_sites_file) as ksf:
                for line in reader(ksf, delimiter=PREKNOWN_INTERESTING_SITES_FILE_SEP):
                    interesting_sites.append(
                        Site(line[BED_CHR_I], line[BED_START_I], line[BED_END_I], line[BED_STRND_I]))

        return interesting_sites

    @staticmethod
    def process_input(min_common_samples, clustering_helper_class_name, loose_min_expression, strict_min_expression,
                      root_path, output_path, include_paths=[], exclude_paths=[], follow_links=False,
                      include_paths_operator=operator.or_, exclude_paths_operator=operator.or_, recursion_depth=100,
                      groups_file=None,
                      clustering_helper_args={},
                      matched_samples_binder_atts=None, matched_sort_atts=SORT_BY_GROUPS_DEC,
                      all_sort_atts=SORT_BY_GROUPS_DEC, essential_groups_suffixes=[],
                      min_groups_bound=2,
                      preknown_sites_file="", known_sites_sanity_file=""):
        # type: (int, str, int, int, str, str, list, list, bool, operator, operator, int, str, dict, list, list, list, list, int, str, str) -> object

        clustering_helper = RNAQuantificationClustererHelperFactory.get_clustering_helper(
            clustering_helper_class_name,
            clustering_helper_args)
        run_on_paths = get_paths(root_path=root_path,
                                 must_include_paths=include_paths,
                                 must_include_operator=include_paths_operator,
                                 exclude_paths=exclude_paths,
                                 exclude_operator=exclude_paths_operator,
                                 follow_links=follow_links,
                                 recursion_depth=recursion_depth,
                                 predicates_dict=clustering_helper.get_files_predicates_dict()
                                 )
        if groups_file:
            load_groups_and_samples_from_file(groups_file=groups_file)

        samples = RNAQuantificationClusterer.convert_dir(clustering_helper, run_on_paths,
                                                         preknown_sites_file)

        strict_output_recs, strict_headers = get_mismatches_sites_quantification_per_site(
            source_name=clustering_helper.combined_source_name,
            coverage_threshold=strict_min_expression)
        loose_output_recs, _ = get_mismatches_sites_quantification_per_site(
            source_name=clustering_helper.combined_source_name,
            coverage_threshold=loose_min_expression)

        RNAQuantificationClusterer.create_output(output_path, samples, all_sort_atts, loose_output_recs, strict_headers,
                                                 strict_output_recs, strict_min_expression, loose_min_expression)

        if matched_samples_binder_atts:
            matched_result_dict, matched_samples = RNAQuantificationClusterer.calc_matched(samples, min_groups_bound,
                                                                                           essential_groups_suffixes,
                                                                                           matched_samples_binder_atts)
            RNAQuantificationClusterer.create_output(os.path.join(output_path, "MatchedSamples"), matched_samples,
                                                     matched_sort_atts,
                                                     loose_output_recs, strict_headers, strict_output_recs)

        logging.info("RNAEditingQuantification Done!")

        return True

    @staticmethod
    def calc_matched(samples, min_groups_bound, essential_groups_suffixes, binder_atts):
        bound_samples = bind_samples(samples, binder_atts)
        matched_groups = {}
        matched_samples = []
        for binder, samples in bound_samples.iteritems():
            if RNAQuantificationClusterer.is_matched(essential_groups_suffixes, min_groups_bound, samples):
                matched_groups[binder] = samples
                matched_samples.extend(samples)

        return matched_samples

    @staticmethod
    def is_matched(essential_groups_suffixes, min_groups_bound, samples):
        if len(samples) < min_groups_bound:
            return False
        for essential_ending in essential_groups_suffixes:
            if not any([sample.group.endswith(essential_ending) for sample in samples]):
                return False
        return True

    @staticmethod
    def create_output(output_path, samples, sort_way, loose_output_recs, output_headers, strict_output_recs,
                      strict_cutoff, loose_cutoff):
        # :type: (str, list[DataObjects.BaseClasses.Sample], str, list, dict, list, int, int ) -> None
        """
        This function generates the outputs of this script.
        :param strict_output_recs: The output recs for the strict cutoff
        :param output_headers: The output headers
        :param loose_output_recs: The output recs for the loose cutoff
        :param strict_cutoff: The strict cutoff
        :param loose_cutoff: The loose cutoff
        :param output_path: The output path.
        :param samples: The samples converted here.
        :param sort_way: How to sort the output (one of general_functions.SORT_SAMPLES_FUNCTIONS_DICT keys)
        :return: None.
        """
        logging.debug("Building Outputs [in %s]\n" % output_path)

        sample_names = convert_params_to_bool_dict(get_sorted_unique_sample_name(samples, sort_way))

        samples_headers = []
        samples_headers.extend(output_headers.pop(GENERAL_HEADERS_KEY))

        groups_headers = samples_headers[:]
        groups_headers.extend(sorted(output_headers.pop(GROUPS_PREVALENCE_HEADERS_KEY)))
        groups_headers.extend(sorted(output_headers.pop(GROUPS_RATIO_HEADERS_KEY)))

        for sample in sample_names:
            samples_headers.extend(sorted(output_headers.pop(sample)))
        for sample in output_headers:  # to add any samples missing from the analysis (only NAs will be added)
            samples_headers.extend(sorted(output_headers[sample]))

        RNAQuantificationClusterer.csv_output.output([os.path.join(output_path, FILTERED_WIDTH_ALL_FORMAT %
                                                                   str(loose_cutoff))], samples_headers,
                                                     loose_output_recs)
        RNAQuantificationClusterer.csv_output.output([os.path.join(output_path, FILTERED_WIDTH_STRICT_FORMAT %
                                                                   str(strict_cutoff))], samples_headers,
                                                     strict_output_recs)

        RNAQuantificationClusterer.csv_output.output([os.path.join(output_path, FILTERED_WIDTH_ALL_FORMAT_GROUPS %
                                                                   str(loose_cutoff))], groups_headers,
                                                     loose_output_recs)
        RNAQuantificationClusterer.csv_output.output([os.path.join(output_path, FILTERED_WIDTH_STRICT_FORMAT_GROUPS %
                                                                   str(strict_cutoff))], groups_headers,
                                                     strict_output_recs)


def main(min_common_samples, clustering_helper, minimal_expression, strict_min_expression, root_path, output_path,
         include_paths=[], exclude_paths=[], follow_links=False, include_paths_operator=operator.or_,
         exclude_paths_operator=operator.or_, recursion_depth=100, groups_file=None, clustering_helper_args={},
         matched_samples_binder_atts=None, sort_matched=SORT_BY_GROUPS_INC,
         all_sort_atts=SORT_BY_GROUPS_DEC, essential_groups_suffixes=[], min_groups_bound=2,
         known_sites_file="", known_sites_sanity_file=""):
    init_logging_dict(LOG_FILE % (output_path, datetime.today().isoformat()))

    RNAQuantificationClusterer.process_input(root_path=root_path,
                                             output_path=output_path,
                                             include_paths=include_paths,
                                             include_paths_operator=include_paths_operator,
                                             exclude_paths=exclude_paths,
                                             exclude_paths_operator=exclude_paths_operator,
                                             recursion_depth=recursion_depth,
                                             follow_links=follow_links,

                                             clustering_helper_class_name=clustering_helper,
                                             clustering_helper_args=clustering_helper_args,

                                             loose_min_expression=minimal_expression,
                                             strict_min_expression=strict_min_expression,

                                             groups_file=groups_file,

                                             min_common_samples=min_common_samples,
                                             min_groups_bound=min_groups_bound,
                                             matched_samples_binder_atts=matched_samples_binder_atts,
                                             matched_sort_atts=sort_matched,
                                             all_sort_atts=all_sort_atts,
                                             essential_groups_suffixes=essential_groups_suffixes,

                                             preknown_sites_file=known_sites_file,
                                             known_sites_sanity_file=known_sites_sanity_file)


if __name__ == "__main__":
    # don't ask why...
    # RNAQuantificationClustererHelper = \
    #     getattr(importlib.import_module('ProcessingPlugins.RNAEditingQuantificationClusterer'),
    #             'RNAQuantificationClustererHelper')
    desc = "Quantify editing level across a list of given positions, divided into groups"


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='RNAEditingQuantification', description=desc, formatter_class=MyFormatter)
    add_get_paths_function_to_argparser(parser=parser)

    parser.add_argument('-o', '--output_dir', metavar="output_path", dest='out_path', nargs='?', required=True,
                        help='Outputer dir for the output')

    parser.add_argument('-m', '--min_common_samples', dest='min_common_sample', nargs='?', required=True, type=int,
                        help='The minimal frequency of a mismatch in a group to retain (in percentages)')
    parser.add_argument('-lc', '--loose_expression_cutoff', metavar="loose_cutoff", dest='min_expression', nargs='?',
                        required=True, type=int, help='The loose cutoff to count a sample as expressing a locus')
    parser.add_argument('-sc', '--strict_expression_cutoff', metavar="strict_cutoff", dest='strict_min_expression',
                        nargs='?', type=int, required=True,
                        help='The strict cutoff to count a sample as expressing a locus')
    parser.add_argument('-ch', '--clustring_helper', dest='cl_h', nargs='?', required=True,
                        choices=RNAQuantificationClustererHelperFactory.CLUSTERING_HELPERS.keys(),
                        help='The clustring helper class (A class for each input type).\n')# + print_atts_help())
    parser.add_argument('-cha', '--clustering_helper_args', dest='clustering_helper_args', nargs='*', required=False,
                        default={},
                        help="""Args to init the clustering helper with. The should be in te format <arg>,<val>
                        '""")

    # parser.add_argument('-b', '--binding_atts', dest='bindings', nargs='*', required=False, default=None,
    #                     help="""The attributes to bind (i.e. "couple") the samples (for matched processing).
    #                Should be provided in the format of either a list of bare atts (e.g. extra_data group),
    #                 or a list of atts and their sub keys (e.g. (extra_data, participant_num)))""",
    #                     type=convert_args_to_dict)
    parser.add_argument('-bm', '--bound_min', dest='min_bound_lim', nargs='?', required=False, type=int,
                        help='The minimal amount of groups to classify samples as matched', default=2)
    parser.add_argument('-e', '--essential_groups_suff', dest='essential_groups', nargs='*', required=False,
                        help='A list of essential groups suffixes to classify samples as matched', default=[])

    parser.add_argument('-sa', '--sort_all', dest='sort_all', nargs='?', required=False, default=SORT_BY_GROUPS_DEC,
                        choices=SORT_SAMPLES_FUNCTIONS_DICT.keys(), help="""How to sort the output.""")
    parser.add_argument('-sm', '--sort_atts_matched', dest='sort_matched', nargs='*', required=False,
                        choices=SORT_SAMPLES_FUNCTIONS_DICT.keys(), default=SORT_BY_GROUPS_DEC,
                        help="""How to sort the output for matched samples""")

    parser.add_argument('-kn', '--known_sites', dest='known_sites', nargs=2, required=False, default=[None, None],
                        help=("Two params:\n"
                              "     - A path to BED file contiaing pre-known interesting sites to check.\n"
                              "     - A path to a BED file containing all possible sites in your run (e.g."
                              " using Shachar's script a file contatains all mirs coordinates and names)"))

    add_groups_file_to_argparser(parser)
    args = parser.parse_args()

    root_path = args.root_dir

    main(min_common_samples=args.min_common_sample,
         clustering_helper=args.cl_h,
         minimal_expression=args.min_expression,
         strict_min_expression=args.strict_min_expression,
         root_path=args.root_dir,
         output_path=args.out_path,
         include_paths=args.include_prefixes,
         include_paths_operator=args.include_operator,
         exclude_paths=args.exclude_prefixes,
         exclude_paths_operator=args.exclude_operator,
         recursion_depth=args.recursion_depth,
         follow_links=args.follow_links,
         groups_file=args.groups_file,
         clustering_helper_args=convert_args_to_dict(args.clustering_helper_args),
         #matched_samples_binder_atts=args.bindings,
         sort_matched=args.sort_matched,
         all_sort_atts=args.sort_all,
         min_groups_bound=args.min_bound_lim,
         essential_groups_suffixes=args.essential_groups,
         known_sites_file=args.known_sites[0],
         known_sites_sanity_file=args.known_sites[1])
