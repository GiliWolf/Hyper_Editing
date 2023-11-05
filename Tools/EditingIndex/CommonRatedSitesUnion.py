__author__ = 'Hillel'
"""
This module analyzes output from naive de novo pipeline
"""
# region Imports
# region Builtin Imports
import operator
import argparse
import os
from collections import Counter

# endregion

# region Imports From Internal Modules
if __name__ == "__main__":
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from Commons.general_functions import path_included

from Outputs.Outputers.RawOutputer import RawOutputer
from Commons.data_structs import SortingHelpFormatter

# endregion
# endregion

INTERSECT_SENSE_FORMAT = "%(bedtools)s intersect -a '%(he_bed_file)s' -b '%(intersect_file)s' -wo -s"


def get_common_sites(bed_files, min_fraction_edited=1, min_coverage=10, unrated_beds=[]):
    sites_c = Counter()
    mm_c = Counter()
    for bfile in unrated_beds:
        with open(bfile) as bed_recs:
            for bed in bed_recs:
                chr, start, end, coverage, mismatched, _ = bed.split("\t")
                if int(coverage) < min_coverage:
                    continue
                sites_c[(chr, start, end)] += 1

    for bfile in bed_files:
        with open(bfile) as bed_recs:
            for bed in bed_recs:
                chr, start, end, coverage, mismatched, _ = bed.split("\t")
                if int(coverage) < min_coverage or int(mismatched) == 0:
                    continue
                mm_c[(chr, start, end)] += 1
    c_sites = dict()
    nsamples = len(bed_files)
    for site in sites_c:
        if sites_c[site] == nsamples and mm_c[site] >= nsamples * min_fraction_edited:
            c_sites[site] = True
    return c_sites


def filter_file(bed_path, opath, common_sites_dict):
    res = ""
    outputer = RawOutputer(pprint_flag=False)
    with open(bed_path) as bfile:
        for brec in bfile:
            chr, start, end, _, _, _ = brec.split("\t")
            if common_sites_dict.get((chr, start, end), False):
                res += brec

    outputer.output([opath, ], res)


def divide_paths(path_groups_csv, path_i, category_is, must_include_paths, u_paths_i):
    categories_d = dict()
    categories_du = dict()
    with open(path_groups_csv) as path_groups:
        _ = path_groups.readline()
        for line in path_groups:
            if line == "":
                continue
            recs = line.strip().split(",")
            if must_include_paths and not path_included(must_include_paths, recs[path_i], operator.or_):
                continue
            categories_d.setdefault(tuple([recs[category_i] for category_i in category_is]), list()).append(
                recs[path_i])
            categories_du.setdefault(tuple([recs[category_i] for category_i in category_is]), list()).append(
                recs[u_paths_i])

    return categories_d, categories_du


def create_common_rated_sites(path_groups_csv, path_i, categories_is, min_common_fraction, output_dir,
                              must_include_paths, min_coverage, u_paths_i):
    categories_d, categories_du = divide_paths(path_groups_csv, path_i, categories_is, must_include_paths, u_paths_i)
    for category, paths in categories_d.iteritems():
        unrated_beds = categories_du[category]
        common_sites = get_common_sites(bed_files=paths, min_fraction_edited=min_common_fraction,
                                        min_coverage=min_coverage, unrated_beds=unrated_beds)
        scategories = [str(c) for c in category]
        for bedf in unrated_beds:
            opath = os.path.join(
                *([output_dir] + scategories + [os.path.basename(bedf) + "CommonIn%s.bed" % min_common_fraction, ]))
            filter_file(bed_path=bedf, opath=opath, common_sites_dict=common_sites)


if __name__ == '__main__':
    desc = "This script intersect rates sites "

    parser = argparse.ArgumentParser(prog='CommonRatedSites', description=desc, formatter_class=SortingHelpFormatter)
    # region CMD Line - Parser Options
    # region Input Options
    inputs_g = parser.add_argument_group(title="Input Files Options")
    inputs_g.add_argument('-p', '--files_csv', metavar="naive rated sites summery file", dest='path_groups_csv',
                          nargs='?',
                          required=True, help="")
    inputs_g.add_argument('-s', '--subdirs_prefixes', metavar="subdirs_prefixes", dest='must_include_paths', nargs='*',
                          required=False, help="Fragments of paths that has to be present for the file to be included"
                                               " (OR operator is the default).",
                          default=[])
    # endregion

    # region Intersect Options
    intersect_g = parser.add_argument_group(title="Intersect Options")
    intersect_g.add_argument('-si', '--category_indexes', metavar="split category index", dest='categories_is',
                             nargs='*', required=True, help="The index(es) of the category(ies) column", type=int)
    intersect_g.add_argument('-pi', '--path_index', metavar="path index", dest='paths_i',
                             nargs='?', required=True, help="The index of the path column", type=int)
    intersect_g.add_argument('-upi', '--unrated_path_index', metavar="unrated path index", dest='u_paths_i',
                             nargs='?', required=True, help="The index of the path column", type=int)
    intersect_g.add_argument('-f', '--fraction', metavar="minimal common fraction", dest='min_common_fraction',
                             nargs='?', required=False, default=1,
                             help="The minimal fraction of samples to include site", type=float)
    intersect_g.add_argument('-c', '--coverage', metavar="minimal coverage", dest='min_coverage',
                             nargs='?', required=False, default=10,
                             help="The minimal fraction of samples to include site", type=int)
    # endregion

    # region Output Options
    outputs_g = parser.add_argument_group(title="Outputs Options")
    outputs_g.add_argument('-o', '--output_dir', metavar="output dir", dest='output_dir', nargs='?',
                           required=True, help="The path to the output dir.")
    # endregion

    # endregion
    options = parser.parse_args()
    create_common_rated_sites(path_groups_csv=options.path_groups_csv,
                              path_i=options.paths_i,
                              categories_is=options.categories_is,
                              min_common_fraction=options.min_common_fraction,
                              output_dir=options.output_dir,
                              must_include_paths=options.must_include_paths,
                              min_coverage=options.min_coverage,
                              u_paths_i=options.u_paths_i)

    # endregion
