__author__ = 'Hillel'
"""
This module analyzes output from the run of Hagit's HE tool.
"""
# region Imports
# region Builtin Imports
from collections import namedtuple, OrderedDict
from csv import reader
from datetime import date
from operator import and_, or_
from pprint import pformat

import argparse
import glob
import inspect
import logging
import multiprocessing
import os
import subprocess

# endregion

# region Imports From Internal Modules
if __name__ == "__main__":
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

from Outputs.Outputers.RawOutputer import RawOutputer
from Outputs.Outputers.CSVOutputer import CSVOutputer

from DataConverters.HyperEditingDetectConverter import HyperEditingBEDConverter

from Session.BlackBoxesWrapperSctipts.STAR242_functios import get_STAR_total_reads, get_STAR_unmapped, \
    get_STAR_uniquely_aligned_reads

from Commons.data_structs import SortingHelpFormatter
from Commons.general_functions import init_logging_dict, convert_args_to_dict, remove_files, \
    add_get_paths_function_to_argparser, get_paths, get_motif, make_recursive_output_dir, format_to_re, get_path_depth

from Commons.consts import *

from Tools.HyperEditing.HyperEditingConsts import *
from Tools.HyperEditing.HyperEditingAnalysis.Consts import SAMPLE_NAME_REGEX, POSITION, MOTIF_OUTPUT_FORMAT

# endregion
# endregion

# region Consts
INTERSECT_SENSE_FORMAT = "%(bedtools)s intersect -a '%(he_bed_file)s' -b '%(intersect_file)s' -wo -s"
HE_BED_NUM_OF_RECS = 6
INTERSECTED_OUTPUT_FORMAT = "IntersectedFiles/%(category)s/%(strand_match)s/%(sample)s/%(orig_filename)s"


# endregion


def intersect_with_regions(he_bed_file, sample_name, output_dir, intersect_file, category_column_index,
                           bedtools_path="bedtools"):
    summary_dict = dict()
    outputer = RawOutputer(pprint_flag=False)
    he_intersect_same_strand_fh = os.popen(
        INTERSECT_SENSE_FORMAT % dict(bedtools=bedtools_path, he_bed_file=he_bed_file,
                                      intersect_file=intersect_file))
    categories_d = dict()
    for line in he_intersect_same_strand_fh:
        recs = line.strip("\n").split("\t")
        categories_d.setdefault(recs[category_column_index], list()).append("\t".join(recs[:HE_BED_NUM_OF_RECS]))

    for category, lines in categories_d.iteritems():
        summary_dict[category] = str(len(lines))
        opath = os.path.join(output_dir, INTERSECTED_OUTPUT_FORMAT % dict(category=category, strand_match="Sense",
                                                                          sample=sample_name,
                                                                          orig_filename=os.path.basename(he_bed_file)))

        outputer.output(paths=[opath, ], raw_data="\n".join(lines))

    return summary_dict


# region Main

def analyze_by_region_type(root_dir,
                           include_paths,
                           include_paths_operator,
                           exclude_paths,
                           exclude_paths_operator,
                           recursion_depth,
                           follow_links,
                           unique_he_sites_files_suffix,
                           # Intersect Options
                           regions_intersect_file,
                           bedtools,
                           split_category_index,

                           # Output Options
                           output_dir,
                           per_sample_output):
    regex = format_to_re(unique_he_sites_files_suffix)
    predicates_dict = {ALL: lambda path: len(regex.findall(path)) >= 1}
    he_sites_bed_files = get_paths(root_path=root_dir if per_sample_output else os.path.join(root_dir,
                                                                                             HE_DETECT_DIR_PER_SAMPLE_BEDS_FORMAT % dict(
                                                                                                 sample_name="all",
                                                                                                 data_type=HE_SITES_VAL)),
                                   must_include_paths=include_paths,
                                   must_include_operator=include_paths_operator,
                                   exclude_paths=exclude_paths,
                                   exclude_operator=exclude_paths_operator,
                                   follow_links=follow_links,
                                   recursion_depth=recursion_depth,
                                   predicates_dict=predicates_dict)[ALL]
    summary_recs = list()
    for bed_file in he_sites_bed_files:
        sample_dir = os.path.basename(os.path.dirname(bed_file))
        # new_o_dir = make_recursive_output_dir(outdir=output_dir,
        #                                       root=root_dir,
        #                                       extension=sample_name,
        #                                       depth=get_path_depth(os.path.dirname(bed_file)) - get_path_depth(
        #                                           root_dir))
        summary_dict = intersect_with_regions(he_bed_file=bed_file,
                                              sample_name=sample_dir,
                                              output_dir=output_dir,  # new_o_dir,
                                              intersect_file=regions_intersect_file,
                                              category_column_index=split_category_index,
                                              bedtools_path=bedtools)
        match = regex.findall(bed_file)[0]
        summary_dict["File"] = bed_file
        summary_dict["PathMatch"] = match
        summary_recs.append(summary_dict)

    outputer = CSVOutputer()
    headers = OrderedDict({"File": 1, "PathMatch": 1})
    for rec in summary_recs:
        headers.update(OrderedDict({key: rec[key] for key in sorted(rec)}))
    outputer.output([os.path.join(output_dir, "SitesSummary.csv")], headers, summary_recs)


if __name__ == '__main__':
    desc = "This script analyses the basic ADAR motif from Hyper Editing analysis"

    parser = argparse.ArgumentParser(prog='HEMotifAnalyzer', description=desc, formatter_class=SortingHelpFormatter)
    # region CMD Line - Parser Options
    # region Input Options
    inputs_g = parser.add_argument_group(title="Input Files Options")
    add_get_paths_function_to_argparser(parser=inputs_g)
    inputs_g.add_argument('-f', '--files_name', metavar="he sites file name", dest='he_sites_file_name', nargs='?',
                          required=False, help="The hyper editing sites BED file name to run on.",
                          default="ESuniqS.A2G.bed")
    # endregion

    # region Intersect Options
    intersect_g = parser.add_argument_group(title="Intersect Options")
    intersect_g.add_argument('-if', '--intersect_file', metavar="regions intersect file", dest='regions_intersect_file',
                             nargs='?', required=True, help="The file to intersect sites BED with.")
    intersect_g.add_argument('-si', '--category_index', metavar="split category index", dest='split_category_index',
                             nargs='?', required=True, help="The index of the category column after intersection")
    intersect_g.add_argument('-b', '--bedtools', dest='bedtools',
                             nargs='?', required=False, help="The bedtools path.", default="bedtools")
    # endregion

    # region Output Options
    outputs_g = parser.add_argument_group(title="Outputs Options")
    outputs_g.add_argument('-o', '--output_dir', metavar="output dir", dest='output_dir', nargs='?',
                           required=True, help="The path to the output dir.")
    outputs_g.add_argument('-ps', '--per_sample_output', dest='per_sample_output', required=False,
                           action='store_true',
                           help="If set, the motif will be calculated per sample rather than for all samples combined.")
    # endregion

    # endregion
    options = parser.parse_args()
    analyze_by_region_type(
        # Inputs options
        root_dir=options.root_dir,
        include_paths=options.include_prefixes,
        include_paths_operator=options.include_operator,
        exclude_paths=options.exclude_prefixes,
        exclude_paths_operator=options.exclude_operator,
        recursion_depth=options.recursion_depth,
        follow_links=options.follow_links,
        unique_he_sites_files_suffix=options.he_sites_file_name,

        # Intersect Options
        regions_intersect_file=options.regions_intersect_file,
        bedtools=options.bedtools,
        split_category_index=int(options.split_category_index),

        # Output Options
        output_dir=options.output_dir,
        per_sample_output=options.per_sample_output)

    # endregion
