__author__ = 'Hillel'
"""
This module analyzes output from naive de novo pipeline
"""
# region Imports
# region Builtin Imports
import operator
import argparse
import os
from numpy import histogram

# endregion

# region Imports From Internal Modules
if __name__ == "__main__":
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from Commons.general_functions import path_included, add_get_paths_function_to_argparser, format_to_re, get_paths

from Outputs.Outputers.CSVOutputer import CSVOutputer
from Commons.data_structs import SortingHelpFormatter

# endregion
# endregion

INTERSECT_FORMAT = "%(bedtools)s intersect -a '%(bed_file)s' -b '%(intersect_file)s' -wa -s -u"


def get_sites_histogram(bed_file, fraction_bins, min_coverage, refseq, bedtools):
    freqs_c = list()
    intersect_rated_fh = os.popen(INTERSECT_FORMAT % dict(bedtools=bedtools, bed_file=bed_file,
                                                          intersect_file=refseq))
    for line in intersect_rated_fh:
        chr, start, end, coverage, mismatched, _ = line.strip("\n").split("\t")
        if int(coverage) < min_coverage:
            continue
        freqs_c.append(float(mismatched) / float(coverage))

    hist = histogram(freqs_c, bins=sorted(fraction_bins) + [1.0,])
    return {str(key): str(val) for key, val in zip(hist[1], hist[0])}


def rated_sites_histogram(root_dir, include_paths, include_paths_operator, exclude_paths, exclude_paths_operator,
                          recursion_depth, follow_links, rated_sites_files_suffix, min_coverage, fractions, refseq,
                          bedtools, output_dir):
    regex = format_to_re(rated_sites_files_suffix)
    predicates_dict = {"all": lambda path: len(regex.findall(path)) >= 1}
    rated_sites_bed_files = get_paths(root_path=root_dir,
                                      must_include_paths=include_paths,
                                      must_include_operator=include_paths_operator,
                                      exclude_paths=exclude_paths,
                                      exclude_operator=exclude_paths_operator,
                                      follow_links=follow_links,
                                      recursion_depth=recursion_depth,
                                      predicates_dict=predicates_dict)["all"]
    recs = list()
    for rated_sites_bed_file in rated_sites_bed_files:
        sites_freqs = get_sites_histogram(bed_file=rated_sites_bed_file, fraction_bins=fractions,
                                           min_coverage=min_coverage, refseq=refseq, bedtools=bedtools)
        sites_freqs["SamplePath"] = rated_sites_bed_file
        recs.append(sites_freqs)

    outputer = CSVOutputer()
    headers = ["SamplePath",] + ["%s" % freq for freq in sorted(fractions)]
    outputer.output([os.path.join(output_dir, "RatedStesFreqs.csv"),], headers, recs)


if __name__ == '__main__':
    desc = "This script intersect rates sites "

    parser = argparse.ArgumentParser(prog='CommonRatedSites', description=desc, formatter_class=SortingHelpFormatter)
    # region CMD Line - Parser Options
    # region Input Options
    inputs_g = parser.add_argument_group(title="Input Files Options")
    add_get_paths_function_to_argparser(parser=inputs_g)
    inputs_g.add_argument('-f', '--files_name', metavar="he sites file name", dest='rated_sites_file_name', nargs='?',
                          required=False, help="The hyper editing sites BED file name to run on.",
                          default=".OnlyA2GNoStrand.SigOnly.bed")
    # endregion

    # region Intersect Options
    intersect_g = parser.add_argument_group(title="Intersect Options")
    intersect_g.add_argument('-c', '--coverage', metavar="minimal coverage", dest='min_coverage',
                             nargs='?', required=False, default=10,
                             help="The minimal fraction of samples to include site", type=int)
    intersect_g.add_argument('--refeseq', dest='refseq', nargs='?', required=False,
                             default="/home/alu/hillelr/scripts/Releases/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz",
                             help="Refseq intersect file")
    intersect_g.add_argument('-pc', '--fractions', metavar="count fraction", dest='fractions',
                             nargs='*', required=False, default=[0.1, 0.2, 0.3, 0.4, 0.5],
                             help="The fractions to count in between", type=float)
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
    rated_sites_histogram(root_dir=options.root_dir,
                          include_paths=options.include_prefixes,
                          include_paths_operator=options.include_operator,
                          exclude_paths=options.exclude_prefixes,
                          exclude_paths_operator=options.exclude_operator,
                          recursion_depth=options.recursion_depth,
                          follow_links=options.follow_links,
                          rated_sites_files_suffix=options.rated_sites_file_name,

                          min_coverage=options.min_coverage,
                          fractions=options.fractions,
                          refseq=options.refseq,
                          bedtools=options.bedtools,
                          output_dir=options.output_dir)

    # endregion
