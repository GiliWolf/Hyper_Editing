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
from Outputs.Outputers.CSVOutputer import CSVOutputer
from Commons.data_structs import SortingHelpFormatter

# endregion
# endregion

MERGE_SITES_FORMAT = "cat %(beds)s|%(bedtools)s sort -i stdin| %(bedtools)s merge -i stdin > %(output_file)s"
INTERSECT_FORMAT = "%(bedtools)s intersect -a '%(bed_file)s' -b '%(union_sites)s' -wa  -u > %(output_file)s"


def get_sites_u(bed_files, output_file, bedtools_path="bedtools"):
    os.system(MERGE_SITES_FORMAT % dict(bedtools=bedtools_path, beds=" ".join(bed_files), output_file=output_file))


def filter_file(bed_path, opath, union_sites_file, bedtools_path="bedtools"):
    os.system(INTERSECT_FORMAT % dict(bedtools=bedtools_path, bed_file=bed_path, output_file=opath,
                                      union_sites=union_sites_file))


def divide_paths(path_groups_csv, path_i, unrated_path_i, category_is, must_include_paths):
    categories_d = dict()
    with open(path_groups_csv) as path_groups:
        _ = path_groups.readline()
        for line in path_groups:
            if line == "":
                continue
            recs = line.strip().split(",")
            if must_include_paths and not path_included(must_include_paths, recs[path_i], operator.or_):
                continue
            categories_d.setdefault(tuple([recs[category_i] for category_i in category_is]), list()).append(
                (recs[path_i], recs[unrated_path_i]))

    return categories_d


def create_common_rated_sites(path_groups_csv, path_i, categories_is, output_dir,
                              must_include_paths, u_paths_i, bedtools):
    categories_d = divide_paths(path_groups_csv, path_i, u_paths_i, categories_is, must_include_paths)
    orecs = list()

    for category, paths in categories_d.iteritems():
        scategories = os.path.join(*[str(c) for c in category])
        u_sites_p = os.path.join(output_dir, scategories, "UnionSites.bed")
        if not os.path.isdir(os.path.dirname(u_sites_p)):
            os.makedirs(os.path.dirname(u_sites_p))
        get_sites_u(bed_files=[p[0] for p in paths],
                    output_file=u_sites_p,
                    bedtools_path=bedtools)
        for bedf in paths:
            opath = os.path.join(output_dir, scategories, os.path.basename(bedf[1]) + "UnionSites.bed")
            filter_file(bed_path=bedf[1],
                        opath=opath,
                        union_sites_file=u_sites_p,
                        bedtools_path=bedtools)


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
    intersect_g.add_argument('-b', '--bedtools', dest='bedtools',
                             nargs='?', required=False, help="The bedtools path.", default="bedtools")

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
                              u_paths_i=options.u_paths_i,
                              output_dir=options.output_dir,
                              must_include_paths=options.must_include_paths,
                              bedtools=options.bedtools)

    # endregion
