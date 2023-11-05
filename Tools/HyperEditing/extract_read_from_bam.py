#! /private/common/Software/anaconda/anaconda3/envs/python3/bin/python
# https://timoast.github.io/blog/2015-10-12-extractreads/
"""
extract_reads.py
Created by Tim Stuart
"""

import pysam
import argparse
import os

def get_read_names(bam):
    name_set = set()
    bam_file = pysam.AlignmentFile(bam, 'rb')
    [name_set.add(read.to_dict()["name"]) for read in bam_file]
    bam_file.close()
    print(len(name_set))

    return name_set


def extract_reads(bam, output, read_names):
    temp_unsorted_bam_path = output.replace(".bam", ".unsortedTemp.bam")
    # open BAM
    bam_file = pysam.AlignmentFile(bam, "rb")
    # index by read name
    index = pysam.IndexedReads(bam_file)
    index.build()
    # open an output BAM file
    output_bam = pysam.AlignmentFile(temp_unsorted_bam_path, "wb", template=bam_file)
    # write all names found to file
    for name in read_names:
        for read in index.find(name):
            output_bam.write(read)
    # close files
    output_bam.close()
    bam_file.close()
    # sort output
    pysam.sort("-o", output, temp_unsorted_bam_path)
    # Remove unsorted BAM file
    os.remove(temp_unsorted_bam_path)


def main(input_bam, other_bam_with_names, output_bam_path):
    unwanted_read_names = get_read_names(other_bam_with_names)
    input_read_names = get_read_names(input_bam)
    wanted_reads = input_read_names-unwanted_read_names
    print(len(wanted_reads))

    extract_reads(bam=input_bam, read_names=wanted_reads, output=output_bam_path)


if __name__ == "__main__":
    # Create parser
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description='Create sorted BAM file without reads that appear in other given BAM')
    parser.add_argument('-i', '--input_bam', metavar="input_bam", help='Input BAM file', required=True)
    parser.add_argument('-b', '--names_bam', metavar="names_bam", help='BAM with names to be excluded', required=True)
    parser.add_argument('-o', '--out_bam', metavar="out_bam", help='Extracted alignments output', required=True)
    options = parser.parse_args()

    main(options.input_bam, options.names_bam, options.out_bam)