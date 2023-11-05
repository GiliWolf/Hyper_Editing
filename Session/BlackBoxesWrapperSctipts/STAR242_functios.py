__author__ = 'Hillel'
# =====================imports=====================#
import argparse
import os
import re

# =====================constants===================#
# TODO: safe delete this copy
STAR_LOG_SUFFIX = r'Log.final.out'
# e.g. /private5/Projects/AutoImmune/MS/RRMS/RRMS/SRR1839798/SRR1839798Log.final.out
TOTAL_READS_RE = re.compile(r"Number of input reads .*?(\d+)")
UNIQUE_READS_RE = re.compile(r"Uniquely mapped reads number .*?(\d+)")

SUMMERY_FILE_HEADERS = ["SamplePath", "OverallReads", "UniquelyMappedReads", "UnmappedReads"]

# =====================functions===================#


def get_STAR_total_reads(output_dir):
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            if f.endswith(STAR_LOG_SUFFIX):
                with open(os.path.join(root, f)) as log:
                    data = log.read()
                    return TOTAL_READS_RE.findall(data)[0], os.path.join(root, f)


def get_STAR_unmapped(output_dir):
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            if f.endswith(STAR_LOG_SUFFIX):
                with open(os.path.join(root, f)) as log:
                    data = log.read()
                    tot = int(TOTAL_READS_RE.findall(data)[0])
                    uniq = int(UNIQUE_READS_RE.findall(data)[0])
                    return str(tot - uniq), os.path.join(root, f)


def get_STAR_uniquely_aligned_reads(output_dir):
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            if f.endswith(STAR_LOG_SUFFIX):
                with open(os.path.join(root, f)) as log:
                    data = log.read()
                    return UNIQUE_READS_RE.findall(data)[0], os.path.join(root, f)


def write_summery_file(output_dir, output_path):
    from Outputs.Outputers.CSVOutputer import CSVOutputer
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            if f.endswith(STAR_LOG_SUFFIX):
                with open(os.path.join(root, f)) as log:
                    data = log.read()
                    path = root
                    total_reads = TOTAL_READS_RE.findall(data)[0]
                    uniq_mapped = UNIQUE_READS_RE.findall(data)[0]
                    unmapped = str(int(total_reads) - int(uniq_mapped))
                    outputer = CSVOutputer(override=False, append=True)
                    outputer.output([output_path],["SamplePath", "OverallReads", "UniquelyMappedReads", "UnmappedReads"],
                                    [{"SamplePath": path, "OverallReads": total_reads, "UniquelyMappedReads": uniq_mapped,
                                      "UnmappedReads": unmapped}])
                    return


if __name__ == "__main__":
    ACTIONS_FUNC_DICT = {
        "total_reads_count": get_STAR_total_reads,
        "uniquely_aligned_reads_count": get_STAR_uniquely_aligned_reads,
        "unmmaped_reads_count":get_STAR_unmapped
    }
    desc = """Get metadata from STAR 242 runs."""

    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    parser = argparse.ArgumentParser(prog='STAR242 Help funcs', description=desc, formatter_class=MyFormatter)

    parser.add_argument('-s', '--star_output_dir', metavar="star output dir",dest='s_output_dir', nargs='?', required=True,
                        help="The folder containing STAR output from the run.")

    parser.add_argument('-o', '--output_path', metavar="output path",dest='output_path', nargs='?', required=False,
                        help="A path to output the results.", default="./star_meta_out.csv")

    parser.add_argument('-a', '--action', metavar="action",dest='action', nargs='?', required=False,
                        default="total_reads_count", help="The wanted data from STAR's run",
                        choices=ACTIONS_FUNC_DICT.keys())

    parser.add_argument('-sf', '--summery_file', dest='summery_file', required=False,
                        help="If set, will append to output file summery all statistics.", action='store_true')


    options = parser.parse_args()
    if options.summery_file:
        write_summery_file(options.s_output_dir, options.output_path)
    else:
        res, path = ACTIONS_FUNC_DICT[options.action](options.s_output_dir)
        with open(options.output_path, 'wb') as o:
            o.write(" ".join([str(res), path]))
