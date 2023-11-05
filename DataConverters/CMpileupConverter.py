from collections import OrderedDict

__author__ = 'Hillel'
"""
This module loads a gzipped DB RefSeqData file of the tab format (the format is as found in the UCSC Genome Browser,
 a header is assumed to be present):
1.	Chromosome
2.	Start (0-based)
3.	End
4.	RefSeq Name (e.g. NR_075077.1)
5.	Gene Common name (e.g. C1orf141)
6.	Strand
7.	Exons Start Positions (0-based) separated by ","
8.	Exons End Positions separated by ","

"""
# =====================imports=====================#
# region Builtin Import
import logging
from pprint import pformat
# endregion

# region InternalImports
from Commons.consts import UNKNOWN_STRAND, MismatchesAndRefsEnum, GENOMIC_REGION, SITE_START, SITE_END
from Tools.EditingIndex.EditingIndexConsts import COVERAGE_AT_N_POSITIONS_FORMAT, NUM_OF_MM_FORMAT, \
    NUM_OF_CANONICAL_FORMAT, NUM_OF_MM_SITES_FORMAT, NUM_OF_N_SITES_COVERED_FORMAT, SNPS_NUM_OF_MM_FORMAT, \
    SNPS_NUM_OF_MM_SITES_FORMAT
from Commons.data_structs import Site

# endregion

# =====================constants===================#
TABLE_SEP = "\t"
REGION_SEP = "_"


# =====================functions===================#
def get_empty_region_counts():
    """
    This function inits an empty counts dict.
    :return dict: an initialized count dict.
    """
    empty_region_counts = dict()
    for ref_pos in MismatchesAndRefsEnum.ALL_REFS:
        empty_region_counts[COVERAGE_AT_N_POSITIONS_FORMAT % ref_pos] = 0
        empty_region_counts[NUM_OF_N_SITES_COVERED_FORMAT % ref_pos] = 0
    for mm_type in MismatchesAndRefsEnum.ALL_MISMATCHES:
        empty_region_counts[NUM_OF_MM_FORMAT % mm_type] = 0
        empty_region_counts[NUM_OF_MM_SITES_FORMAT % mm_type] = 0
        empty_region_counts[SNPS_NUM_OF_MM_FORMAT % mm_type] = 0
        empty_region_counts[SNPS_NUM_OF_MM_SITES_FORMAT % mm_type] = 0
    for base in MismatchesAndRefsEnum.ALL_REFS:
        empty_region_counts[NUM_OF_CANONICAL_FORMAT % base] = 0
    return empty_region_counts


def parse_count_pileup(count_pileup_file, snps={}):
    # type: (str, dict) -> list
    """
    This function loads a RefSeqData DB file
    :param str count_pileup_file: The file in the described format.
    :param dict snps: A dictionary index of the known SNPs (to bool for ease of use).
    :return dict of the records
    :rtype dict[Site, list[SampleCoverageData]]
    """
    records = list()
    num_of_coverages = 0
    chroms = dict()
    empty_region_counts = get_empty_region_counts()

    with open(count_pileup_file) as count_pileup_fh:
        for line in count_pileup_fh:
            chrom, start, end, ref_pos, total_coverage, adenosines, cytosines, guanosines, thymines, unreconized, low_q = line.strip(
                "\n").split(TABLE_SEP)
            site_counts = OrderedDict()
            site_counts[GENOMIC_REGION] = chrom
            site_counts[SITE_START] = start
            site_counts[SITE_END] = end
            site_counts["RefBase"] = ref_pos
            site_counts.update(empty_region_counts.copy())

            chroms[chrom] = True
            start = int(start)
            end = int(end)

            site = Site(chrom, start, end, UNKNOWN_STRAND)

            coverage_dict = {
                MismatchesAndRefsEnum.A: int(adenosines),
                MismatchesAndRefsEnum.C: int(cytosines),
                MismatchesAndRefsEnum.G: int(guanosines),
                MismatchesAndRefsEnum.T: int(thymines),
            }

            num_of_coverages += 1

            ref_pos = MismatchesAndRefsEnum.REFS_TRANSLATE.get(ref_pos)
            if not ref_pos:
                continue  # this is an unknown position
            total_wo_snps = sum(coverage_dict.itervalues())
            for base, coverage in coverage_dict.iteritems():
                if coverage == 0:
                    continue
                if ref_pos == base:
                    site_counts[NUM_OF_CANONICAL_FORMAT % base] += coverage
                    continue

                mm_type = MismatchesAndRefsEnum.MISMATCH_TYPE_PARSE[ref_pos + base]
                has_snp = snps.get(mm_type, {}).get(site, False)
                if has_snp:
                    site_counts[SNPS_NUM_OF_MM_FORMAT % mm_type] += coverage
                    site_counts[SNPS_NUM_OF_MM_SITES_FORMAT % mm_type] += 1
                    total_wo_snps -= coverage
                else:
                    site_counts[NUM_OF_MM_FORMAT % mm_type] += coverage
                    site_counts[NUM_OF_MM_SITES_FORMAT % mm_type] += 1

            if total_wo_snps > 0:
                site_counts[COVERAGE_AT_N_POSITIONS_FORMAT % ref_pos] += total_wo_snps
                site_counts[NUM_OF_N_SITES_COVERED_FORMAT % ref_pos] += 1

            site_counts["Total"] = str(total_wo_snps)
            records.append(site_counts)

    logging_msg = "Found %s regions, %s coverages in chromosomes: %s" % (str(len(records)), str(num_of_coverages),
                                                                         pformat(chroms.keys()))
    logging.info(logging_msg)
    return records
