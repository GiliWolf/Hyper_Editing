import argparse
from datetime import datetime
from collections import OrderedDict
import multiprocessing
import logging
import os

if __name__ == '__main__':
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from scipy.stats import binom_test

from Commons.consts import UNKNOWN_STRAND, MismatchesAndRefsEnum, SAMPLE
from Tools.EditingIndex.EditingIndexConsts import NUM_OF_MM_FORMAT, NUM_OF_MM_SITES_FORMAT, SNPS_NUM_OF_MM_FORMAT, \
    SNPS_NUM_OF_MM_SITES_FORMAT
from Commons.data_structs import Site
from Commons.general_functions import add_get_paths_function_to_argparser, get_paths, get_sample_name, init_logging_dict
from Outputs.Outputers.CSVOutputer import CSVOutputer

from Tools.EditingIndex.DataConverters.SNPsDBCoverter import load_snps_file


def get_num_of_sites(cmpileup_path, snps, output_path, sample_name, binom_p, pval_threshold, semaphore):
    # type: (str, dict, str, str, float, float, multiprocessing.Semaphore) -> None

    outputer = CSVOutputer(override=False, append=True)
    mm_counts = OrderedDict({SAMPLE: sample_name})
    for mm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
        mm_counts[NUM_OF_MM_FORMAT % mm_type] = 0
        mm_counts[NUM_OF_MM_SITES_FORMAT % mm_type] = 0
        mm_counts[SNPS_NUM_OF_MM_FORMAT % mm_type] = 0
        mm_counts[SNPS_NUM_OF_MM_SITES_FORMAT % mm_type] = 0

    try:
        logging.info("Loading Coverage File %s" % cmpileup_path)
        with open(cmpileup_path) as count_pileup_fh:
            for line in count_pileup_fh:
                chrom, start, end, position, ref_pos, total_coverage, adenosines, cytosines, guanosines, thymines, unreconized, low_q = line.strip(
                    "\n").split("\t")

                position = int(position)

                site = Site(chrom, position - 1, position, UNKNOWN_STRAND)

                coverage_dict = {
                    MismatchesAndRefsEnum.A: int(adenosines),
                    MismatchesAndRefsEnum.C: int(cytosines),
                    MismatchesAndRefsEnum.G: int(guanosines),
                    MismatchesAndRefsEnum.T: int(thymines),
                }

                ref_pos = MismatchesAndRefsEnum.REFS_TRANSLATE.get(ref_pos)
                if not ref_pos:
                    continue  # this is an unknown position
                for base, coverage in coverage_dict.iteritems():
                    if coverage == 0:
                        continue
                    if ref_pos == base:
                        continue

                    mm_type = MismatchesAndRefsEnum.MISMATCH_TYPE_PARSE[ref_pos + base]
                    has_snp = snps.get(mm_type, {}).get(site, False)

                    if not binom_test(coverage, total_coverage, binom_p) <= pval_threshold:
                        continue

                    if has_snp:
                        mm_counts[SNPS_NUM_OF_MM_FORMAT % mm_type] += coverage
                        mm_counts[SNPS_NUM_OF_MM_SITES_FORMAT % mm_type] += 1
                    else:
                        mm_counts[NUM_OF_MM_FORMAT % mm_type] += coverage
                        mm_counts[NUM_OF_MM_SITES_FORMAT % mm_type] += 1

        outputer.output([output_path, ], mm_counts.keys(), [{key: str(val) for key, val in mm_counts.iteritems()}, ])
        logging.info("Done Loading Coverage Data of %s" % sample_name)
    except Exception, e:
        logging.exception("Failed To Process Sample %s!" % sample_name)

    finally:
        semaphore.release()


def translate_phred_score_to_seq_err_chance(score):
    return 10 ** (-score / 10)


def main(root_dir, output_dir, regions_bed, bedtools_path, include_paths, include_paths_operator, exclude_paths,
         exclude_paths_operator, recursion_depth, follow_links, bam_files_suffix, snps_file, max_processes_sample,
         min_phred_score, pval_threshold, truncs):
    timestamp = datetime.today().isoformat()
    init_logging_dict(os.path.join(output_dir, ".".join(["BinomialRatedIndex", timestamp, "log"])))

    semaphore = multiprocessing.Semaphore(max_processes_sample)
    predicates_dict = {"All": lambda path: path.endswith(bam_files_suffix)}
    seq_err_p = translate_phred_score_to_seq_err_chance(min_phred_score)
    cmpileup_files = get_paths(root_path=root_dir, must_include_paths=include_paths,
                               must_include_operator=include_paths_operator, exclude_paths=exclude_paths,
                               exclude_operator=exclude_paths_operator, follow_links=follow_links,
                               recursion_depth=recursion_depth, predicates_dict=predicates_dict)["All"]
    print cmpileup_files
    logging.info("Loading SNPs File")
    snps = load_snps_file(snps_file=snps_file, intersect_first_with=regions_bed, bedtools_path=bedtools_path)
    threads = list()
    out_path = os.path.join(output_dir, "BinomialRatedIndex.csv")
    for cfile in cmpileup_files:
        sample_name = get_sample_name(os.path.basename(cfile), truncs[0], truncs[1])
        tr = multiprocessing.Process(target=get_num_of_sites,
                                     args=(cfile,
                                           snps,
                                           out_path,
                                           sample_name,
                                           seq_err_p,
                                           pval_threshold,
                                           semaphore),
                                     name="IndexByRatedSites%s" % sample_name)
        tr.daemon = True
        threads.append(tr)

    for thr in threads:
        semaphore.acquire()
        thr.run()
    for thr in multiprocessing.active_children():
        thr.join(3600)

    logging.info("Done!")


if __name__ == "__main__":
    script_dir = os.path.dirname(__file__)
    desc = """Rate and sums mismatches from cmpileup according to the binomial distribution"""


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='Editing Index Runner', description=desc, formatter_class=MyFormatter)
    add_get_paths_function_to_argparser(parser=parser)

    parser.add_argument('-o', '--outpath', metavar="outpath", dest='outpath', nargs='?', required=True,
                        help="The path of the out to write to (appending)")
    parser.add_argument('-f', '--files_suffix', metavar="bam files suffix", dest='files_suffix', nargs='?',
                        required=True, help="A suffix of the cmpileup files to run on")
    parser.add_argument('--snps', metavar="SNPs file", dest='snps_file', nargs='?', required=True,
                        help="The path of the SNPs db file")
    parser.add_argument('-q', '--qScore', metavar="minimal Phred score", dest='minQ', nargs='?', required=False,
                        default=0.001, help="The minimal Phred score used for the mismatch calling", type=float)
    parser.add_argument('-pt', '--pvalThreshold', metavar="p-value threshold", dest='pval', nargs='?', required=False,
                        default=0.05, help="The maximal p-value for the binomial test", type=float)
    parser.add_argument('--ts', metavar="sample threads", dest='max_processes_sample', nargs='?', required=False,
                        default=10, type=int, help="The number of samples to process in parallel")
    parser.add_argument('-rb', '--regions', metavar="regions bed", dest='regions_bed', nargs='?', required=False,
                        help="The path of a regions-of-interest bed file (improves run time and memory print)",
                        default='')
    parser.add_argument('-b', '--bedtools', metavar="bedtools path", dest='bedtools', nargs='?', required=False,
                        help="The path of bedtools to use (in the case of usage with regions-of-interest file)",
                        default='bedtools')
    parser.add_argument('-t', '--truncs', metavar="truncs", dest='truncs', nargs='?', required=False,
                        default="1,0", help="The truncation of the file name "
                                            "<suffixes to retain>,<chars to drop from filename>")

    options = parser.parse_args()
    truncs = [int(options.truncs.split(",")[0]), int(options.truncs.split(",")[1])]

    main(root_dir=options.root_dir,
         output_dir=options.outpath,
         regions_bed=options.regions_bed,
         include_paths=options.include_prefixes,
         include_paths_operator=options.include_operator,
         exclude_paths=options.exclude_prefixes,
         exclude_paths_operator=options.exclude_operator,
         recursion_depth=options.recursion_depth,
         follow_links=options.follow_links,
         bam_files_suffix=options.files_suffix,
         snps_file=options.snps_file,
         bedtools_path=options.bedtools,
         max_processes_sample=options.max_processes_sample,
         min_phred_score=options.minQ,
         pval_threshold=options.pval,
         truncs=truncs
         )
