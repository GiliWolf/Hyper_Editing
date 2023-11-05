__author__ = 'Hillel'
# =====================imports=====================#
import argparse
from datetime import datetime
import os
import logging
import re
from collections import OrderedDict

if __name__ == '__main__':
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from Commons.data_structs import SortingHelpFormatter
from Commons.general_functions import add_get_paths_function_to_argparser, get_paths, init_logging_dict

from Outputs.Outputers.CSVOutputer import CSVOutputer

# =====================constants===================#

UNCMP_ADD = "_unompressed_tmp_"

FASTQC_DATA_SUFFIX = r'fastqc_data.txt'
# e.g. /home/alu/hillelr/melanoma/fastqc/SRR3083584_2_fastqc/fastqc_data.txt
FASTQC_SUMMARY_SUFFIX = r'summary.txt'
# e.g. /home/alu/hillelr/melanoma/fastqc/SRR3083584_2_fastqc/summary.txt
SAMPLE_RE = re.compile(r"Filename\t(.*)\t")

PHRED_TEST_RE = re.compile(r">>Per base sequence quality\t(.+?)\r?\n(.+?)>>END_MODULE", re.DOTALL)
PHRED_LINE_RE = re.compile(r"(?:\r?\n|#)(.*?)\t(.*?)\t(.*?)\t")

BASE_CONTENT_RE = re.compile(r">>Per base sequence content\t(.+?)\r?\n(.+?)>>END_MODULE", re.DOTALL)
BASE_CONTENT_LINE_RE = re.compile(r"(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\r?\n")

SEQ_LENGTH_DISTRIBUTION_RE = re.compile(r">>Sequence Length Distribution\t(.+?)\r?\n(.+?)>>END_MODULE", re.DOTALL)
SEQ_LENGTH_DISTRIBUTION_LINE_RE = re.compile(r"(.*?)\t(.*?)\r?\n")

SEQ_DUP_RE = re.compile(r">>Sequence Duplication Levels\t(.+?)\r?\n(.+?)>>END_MODULE", re.DOTALL)
SEQ_DUP_TOT_LINE_RE = re.compile(r"#Total Duplicate Percentage\t(.*?)\r?\n")
SEQ_DUP_LINE_RE = re.compile(r"(\d.*?)\t(.*?)\r?\n")

SUMMARY_LINE_RE = re.compile(r"(.+)\t(.+)\t(.+)\r?\n")

PHRED_SUMMARY_FILENAME = "FastqCPhredScores.csv"
BASE_CONTENT_FILENAME = "FastqPerBaseContent.csv"
SEQ_LENGTH_DISTRIBUTION_FILENAME = "SeqLengthDistribution.csv"
SEQ_DUP_FILENAME = "SeqDuplications.csv"

SUMMARY_FILENAME = "TestResSummary.csv"

BASIC_HEADER_FORMAT = r"%s_For_%s"

DATA_SUMMARY_FILE_HEADERS = ["Path", "Sample", "Result"]
SUMMARY_FILE_HEADERS = ["Path", "Sample"]

TEST_FAIL_ERR_MSG = "Failed To Summarize %s Data for %s(%s)!"


# =====================functions===================#


def summarize_spcified_tests(tests_files, output_path):
    outputer = CSVOutputer()
    q_score_recs = list()
    q_score_headers = OrderedDict({h: True for h in DATA_SUMMARY_FILE_HEADERS})

    base_content_recs = list()
    base_content_headers = OrderedDict({h: True for h in DATA_SUMMARY_FILE_HEADERS})

    seq_length_recs = list()
    seq_length_headers = OrderedDict({h: True for h in DATA_SUMMARY_FILE_HEADERS})

    seq_dup_recs = list()
    seq_dup_headers = OrderedDict({h: True for h in DATA_SUMMARY_FILE_HEADERS})

    for f in tests_files:
        with open(os.path.join(f)) as log:
            data = log.read()
            sample = SAMPLE_RE.findall(data)[0].replace(UNCMP_ADD, "")
            try:
                rec = analyse_q_score_test(data, f, sample)
                q_score_headers.update(OrderedDict([(h, True) for h in rec]))
                q_score_recs.append(rec)
            except:
                logging.exception(TEST_FAIL_ERR_MSG % ("PHRED Scores", sample, f))

            try:
                rec = analyse_base_content(data, f, sample)
                base_content_headers.update(OrderedDict([(h, True) for h in rec]))
                base_content_recs.append(rec)
            except:
                logging.exception(TEST_FAIL_ERR_MSG % ("Per Base Content", sample, f))

            try:
                rec = analyse_seq_len_dist(data, f, sample)
                seq_length_headers.update(OrderedDict([(h, True) for h in rec]))
                seq_length_recs.append(rec)
            except:
                logging.exception(TEST_FAIL_ERR_MSG % ("Sequenced Length Distribution", sample, f))

            try:
                rec = analyse_seq_dup(data, f, sample)
                seq_dup_headers.update(OrderedDict([(h, True) for h in rec]))
                seq_dup_recs.append(rec)
            except:
                logging.exception(TEST_FAIL_ERR_MSG % ("Duplicates", sample, f))

    outputer.output([os.path.join(output_path, SEQ_DUP_FILENAME)], seq_dup_headers, seq_dup_recs)
    outputer.output([os.path.join(output_path, SEQ_LENGTH_DISTRIBUTION_FILENAME)], seq_length_headers, seq_length_recs)
    outputer.output([os.path.join(output_path, BASE_CONTENT_FILENAME)], base_content_headers, base_content_recs)
    outputer.output([os.path.join(output_path, PHRED_SUMMARY_FILENAME)], q_score_headers, q_score_recs)


def analyse_seq_len_dist(data, f, sample):
    rec = OrderedDict()
    rec["Path"] = f
    rec["Sample"] = sample
    rec["Result"], length_data = SEQ_LENGTH_DISTRIBUTION_RE.findall(data)[0]
    seq_len_dist_data = SEQ_LENGTH_DISTRIBUTION_LINE_RE.findall(length_data)
    for line in seq_len_dist_data[1:]:
        header = BASIC_HEADER_FORMAT % ("#Reads", line[0])
        rec[header] = line[1]
    return rec


def analyse_seq_dup(data, f, sample):
    rec = OrderedDict()
    rec["Path"] = f
    rec["Sample"] = sample
    rec["Result"], dup_data = SEQ_DUP_RE.findall(data)[0]
    tot_dup = SEQ_DUP_TOT_LINE_RE.findall(dup_data)[0]
    seq_dup_data = SEQ_DUP_LINE_RE.findall(dup_data)
    header = BASIC_HEADER_FORMAT % ("DuplicatesPC", "Total")
    rec[header] = tot_dup
    for line in seq_dup_data:
        header = BASIC_HEADER_FORMAT % ("DuplicatesPC", line[0])
        rec[header] = line[1]
    return rec


def analyse_base_content(data, f, sample):
    rec = OrderedDict()
    rec["Path"] = f
    rec["Sample"] = sample
    rec["Result"], base_data = BASE_CONTENT_RE.findall(data)[0]
    base_content_data = BASE_CONTENT_LINE_RE.findall(base_data)
    bases = "|".join(base_content_data[0][1:])
    for line in base_content_data[1:]:
        header = BASIC_HEADER_FORMAT % (bases, line[0])
        rec[header] = "|".join(line[1:])
    return rec


def analyse_q_score_test(data, f, sample):
    rec = OrderedDict()
    rec["Path"] = f
    rec["Sample"] = sample
    rec["Result"], phred_data = PHRED_TEST_RE.findall(data)[0]
    for line in PHRED_LINE_RE.findall(phred_data)[1:]:
        rec[BASIC_HEADER_FORMAT % ("Mean", line[0])] = line[1]
        rec[BASIC_HEADER_FORMAT % ("Median", line[0])] = line[2]
    return rec


def summarize_tests_pass_fail(sum_files, output_path):
    for f in sum_files:
        if f.endswith(FASTQC_SUMMARY_SUFFIX):
            with open(f) as log:
                data = log.read()
                outputer = CSVOutputer(override=False, append=True)
                headers = SUMMARY_FILE_HEADERS[:]
                rec = {}
                rec["Path"] = f
                for line in SUMMARY_LINE_RE.findall(data):
                    rec["Sample"] = line[-1]
                    rec[line[1]] = line[0]
                    headers.append(line[1])

                outputer.output([os.path.join(output_path, SUMMARY_FILENAME)], headers, [rec])


def get_fastq_summary(root_dir, output_dir, include_paths, include_paths_operator, exclude_paths,
                      exclude_paths_operator, recursion_depth, follow_links):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    timestamp = datetime.now().isoformat()
    init_logging_dict(os.path.join(output_dir, "FastqCSummarizer" + "." + timestamp + ".log"))

    predicates_dict = {
        FASTQC_DATA_SUFFIX: lambda path: path.endswith(FASTQC_DATA_SUFFIX),
        FASTQC_SUMMARY_SUFFIX: lambda path: path.endswith(FASTQC_SUMMARY_SUFFIX)
    }

    all_files = get_paths(root_path=root_dir, must_include_paths=include_paths,
                          must_include_operator=include_paths_operator, exclude_paths=exclude_paths,
                          exclude_operator=exclude_paths_operator, follow_links=follow_links,
                          recursion_depth=recursion_depth, predicates_dict=predicates_dict)

    summarize_spcified_tests(all_files[FASTQC_DATA_SUFFIX], output_dir)
    summarize_tests_pass_fail(all_files[FASTQC_SUMMARY_SUFFIX], output_dir)


if __name__ == "__main__":
    desc = """Get summary for fastqc analysis"""

    parser = argparse.ArgumentParser(prog='FastqC Summarizer', description=desc,
                                     formatter_class=SortingHelpFormatter)
    add_get_paths_function_to_argparser(parser=parser)
    parser.add_argument('-o', '--output_dir', metavar="output dir", dest='output_dir', nargs='?', required=True,
                        help='The path of the output dir to create.')

    options = parser.parse_args()
    get_fastq_summary(root_dir=options.root_dir,
                      output_dir=options.output_dir,
                      include_paths=options.include_prefixes,
                      include_paths_operator=options.include_operator,
                      exclude_paths=options.exclude_prefixes,
                      exclude_paths_operator=options.exclude_operator,
                      recursion_depth=options.recursion_depth,
                      follow_links=options.follow_links)
