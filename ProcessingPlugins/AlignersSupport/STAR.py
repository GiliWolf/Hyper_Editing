__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Imports
import os
import re
# endregion

# region Internal Imports
from Commons.data_structs import AlignmentData
# endregion

# =====================constants===================#
STAR_LOG_FORMAT = r'%(aligner_output_format)sLog.final.out'

TOTAL_READS_RE = re.compile(r"Number of input reads .*?(\d+)")
UNIQUE_READS_RE = re.compile(r"Uniquely mapped reads number .*?(\d+)")
MULTI_MAPPED_READS_RE = re.compile(r"Number of reads mapped to multiple loci .*?(\d+)")

# =====================functions===================#


def get_star_alignment_data(aligner_output_format, sequencing_method=AlignmentData.SEQ_METHOD_SE):
    log_path = STAR_LOG_FORMAT % dict(aligner_output_format=aligner_output_format)
    pe = sequencing_method == AlignmentData.SEQ_METHOD_PE
    with open(log_path) as log:
        data = log.read()
        total_reads = int(TOTAL_READS_RE.findall(data)[0])
        uniquely_mapped = int(UNIQUE_READS_RE.findall(data)[0])
        multi_mapped = int(MULTI_MAPPED_READS_RE.findall(data)[0])

        if pe:
            total_reads = total_reads * 2
            uniquely_mapped = uniquely_mapped * 2

        align_data = AlignmentData(total_fragments=total_reads,
                                   uniquely_mapped=uniquely_mapped,
                                   multi_mapped=multi_mapped,
                                   unmapped=total_reads - uniquely_mapped - multi_mapped,
                                   sequencing_method=sequencing_method)
    return align_data
