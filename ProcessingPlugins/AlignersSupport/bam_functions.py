__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Imports
import os
# endregion

# region Internal Imports
# endregion

# =====================constants===================#
SAMTOOLS_COUNT_UNQ_MAPPED_CMD_FORMAT = "%(samtools)s view -c -F 260 %(bam)s"
SAMTOOLS_COUNT_MULTI_MAPPED_CMD_FORMAT = "%(samtools)s view -f 256 %(bam)s| awk -F \"\t\" '!a[$1]++ {count++} END {" \
                                         "print count}' "
SAMTOOLS_COUNT_PE_CMD_FORMAT = "%(samtools)s view -c -f 1 %(bam)s"
# =====================functions===================#


def count_uniquely_mapped_reads_in_bam(samtools_path, bam_path):
    """
    This function retrieves the number of uniquely mapped reads from a bam file.
    :param str bam_path: The path of the bam to count reads in.
    :param str samtools_path: The run path of samtools.
    :return int: The  number of mapped reads in the bam file.
    """
    return int(os.popen(SAMTOOLS_COUNT_UNQ_MAPPED_CMD_FORMAT % dict(samtools=samtools_path, bam=bam_path)).read())


def count_multi_mapped_reads_in_bam(samtools_path, bam_path):
    """
    This function retrieves the number of multi mapped reads from a bam file.
    :param str bam_path: The path of the bam to count reads in.
    :param str samtools_path: The run path of samtools.
    :return int: The  number of mapped reads in the bam file.
    """
    return int(os.popen(SAMTOOLS_COUNT_MULTI_MAPPED_CMD_FORMAT % dict(samtools=samtools_path, bam=bam_path)).read())


def check_for_pe_reads_in_bam(samtools_path, bam_path):
    """
    This function checks if a bam file is derived from paired end reads.
    :param str bam_path: The path of the bam to count reads in.
    :param str samtools_path: The run path of samtools.
    :return bool: True if the data was derived from paired end reads, False if from single end.
    """
    return 0 < int(os.popen(SAMTOOLS_COUNT_PE_CMD_FORMAT % dict(samtools=samtools_path, bam=bam_path)).read())
