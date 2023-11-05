__author__ = 'Hillel'
"""
This module analyzes output from the run of Hagit's HE tool.
"""
# region Imports
# region Builtin Imports

import argparse
import os

# endregion

# region Imports From Internal Modules
if __name__ == "__main__":
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

from Outputs.Outputers.CSVOutputer import CSVOutputer

from DataConverters.HyperEditingDetectConverter import HyperEditingBEDConverter

from Commons.data_structs import SortingHelpFormatter
from Commons.general_functions import add_get_paths_function_to_argparser, get_paths, get_motif, \
    make_recursive_output_dir, \
    get_path_depth

from Commons.consts import *
from Tools.HyperEditing.HyperEditingAnalysis.Consts import SAMPLE_NAME_REGEX, POSITION, MOTIF_OUTPUT_FORMAT
from Tools.HyperEditing.HyperEditingConsts import *


# endregion
# endregion



def generate_sample_motif(sequences, sample_output_dir, sample_name):
    # type: (list, str, str) -> None
    """
    This function gets the immediate neighbors motif form a hyper editing detect file.
    :param list sequences: The path to the hyper editing unique sites file
    :param str sample_output_dir: The path of the output dir.
    :param str sample_name: The name of the sample
    :return: None
    """

    motif_recs = list()
    motif_headers = [POSITION,
                     MismatchesAndRefsEnum.A,
                     MismatchesAndRefsEnum.C,
                     MismatchesAndRefsEnum.G,
                     MismatchesAndRefsEnum.T,
                     ]

    motif = get_motif(mismatch_sense=MismatchesAndRefsEnum.A, sites_seqs=sequences, window_size=1)
    for position in sorted(motif):
        freqs_d = {base: str(freq) for base, freq in motif[position].iteritems()}
        freqs_d[POSITION] = str(position)
        motif_recs.append(freqs_d)

    outputer = CSVOutputer(delim="\t")
    sample_full_o_p = os.path.join(sample_output_dir, MOTIF_OUTPUT_FORMAT % dict(sample_name=sample_name))
    outputer.output(paths=[sample_full_o_p, ], headers=motif_headers, records=motif_recs)


# region Main
def analyze_motif(
        # Inputs options
        root_dir,
        include_paths,
        include_paths_operator,
        exclude_paths,
        exclude_paths_operator,
        recursion_depth,
        follow_links,
        unique_he_sites_files_suffix,

        # Output Options
        output_dir,
        per_sample_output):
    predicates_dict = {ALL: lambda path: path.endswith(unique_he_sites_files_suffix)}
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

    for bed_file in he_sites_bed_files:
        sample_dir = os.path.basename(os.path.dirname(bed_file))
        sample_name = SAMPLE_NAME_REGEX.findall(sample_dir)[0][0]
        new_o_dir = make_recursive_output_dir(outdir=output_dir,
                                              root=os.path.dirname(bed_file),
                                              extension=sample_name,
                                              depth=get_path_depth(bed_file) - get_path_depth(root_dir))
        sequences = [bed_rec.neighbors for bed_rec in HyperEditingBEDConverter.convert(bed_file)]
        generate_sample_motif(sequences=sequences,
                              sample_output_dir=new_o_dir,
                              sample_name=sample_name)


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
    analyze_motif(
        # Inputs options
        root_dir=options.root_dir,
        include_paths=options.include_prefixes,
        include_paths_operator=options.include_operator,
        exclude_paths=options.exclude_prefixes,
        exclude_paths_operator=options.exclude_operator,
        recursion_depth=options.recursion_depth,
        follow_links=options.follow_links,
        unique_he_sites_files_suffix=options.he_sites_file_name,
        # Output Options
        output_dir=options.output_dir,
        per_sample_output=options.per_sample_output
    )

# endregion
