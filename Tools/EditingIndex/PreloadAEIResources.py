# <editor-fold desc="Imports">
import argparse
import os

if __name__ == '__main__':
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from Commons.data_structs import SortingHelpFormatter
from Tools.EditingIndex.A2IEditingIndex import load_resources

# </editor-fold>


if __name__ == "__main__":
    desc = """Run on each file in a given directory a set of steps."""
    parser = argparse.ArgumentParser(prog='Editing Index Runner', description=desc,
                                     formatter_class=SortingHelpFormatter)
    parser.add_argument('--genome', metavar="genome name", dest='genome_name', nargs='?', required=True,
                        help="The file name of the genome assembly (the one used with the index)")
    parser.add_argument('-rb', '--regions', metavar="regions bed", dest='regions_bed', nargs='?', default="",
                        required=False, help="A path to an edited regions bed file to use (the one used with the index)")
    parser.add_argument('--snps', metavar="SNPs file", dest='snps_file', nargs='?', required=False,
                        default="", help="A path to a SNPs file to use (the one used with the index).")
    parser.add_argument('--refseq', metavar="refseq file", dest='refseq_file', nargs='?', required=False,
                        default="", help="A path to a refseq file to use (the one used with the index).")
    parser.add_argument('--genes_expression', metavar="genes expression file", dest='genes_expression_file',
                        default="",
                        nargs='?',
                        required=False, help="A path to a genes expression file to use.")
    parser.add_argument('--tsd', metavar="strands threads", dest='max_processes_strand_decision', nargs='?',
                        required=False, default=50, type=int,
                        help="The maximal strand decisions per sample to process in parallel")
    parser.add_argument('--bedtools_path', metavar="bedtools path", dest='bedtools_path', nargs='?',
                        required=False, default="bedtools", help="The command to invoke bedtools")

    options = parser.parse_args()

    _, _, _ = load_resources(gene_expression_file_good=os.path.exists(options.genes_expression_file),
                             gene_expression_full_path=options.genes_expression_file,
                             max_processes_strand_decision=options.max_processes_strand_decision,
                             refseq_full_path=options.refseq_file,
                             regions_bed_full_path=options.regions_bed,
                             snps_file_good=os.path.exists(options.snps_file),
                             snps_full_path=options.snps_file,
                             bedtools_path=options.bedtools_path,
                             genome_name=options.genome_name,
                             regions_name=os.path.basename(options.regions_bed),
                             load_from_pickle=False,
                             save_to_pickle=True)
