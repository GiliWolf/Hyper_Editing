__author__ = 'Hillel'
# =====================imports=====================#
# region Builtin Imports
import logging
# endregion

from Commons.consts import GENERAL_TEST, BONFERRONI_P_TEST, BH_P_TEST, MATURE_MIR_POS, PRE_MIR_POS, STATISTICS_DATA, \
    MIR_NAME, NA
from Commons.data_structs import StatisticalTestPValue, Site
from Commons.consts import SENSE_STRAND, ANTISENSE_STRAND
from DataObjects.KnownProperties.MismatchesSites import *

# =====================constants===================#
# ----error messages----
PARSER_PROVIDED_ERROR = "RNAEditingQuantificationConverter.convert Must Receive a subclass of" \
                        " EditingQuantificationParser Type As The Second Argument!"

FAILED_TO_CREATE_PARSER = "Failed To Create The Parser!\n%s"

MATURE_RNA = "MatureRNA"
PRE_RNA = "PreRNA"

# **********Shachar's scripts************#

# ---- "readability" consts ----
ShacharUnFilteredResultParser_LINE_SPLITTER = "\n"
ShacharUnFilteredResultParser_COL_SPLITTER = "\t"
ShacharUnFilteredResultParser_P_VALUE = "P-val="

# ---- splicers ----
ShacharFilteredResultParser_LINE_SPLITTER = "\n"
ShacharFilteredResultParser_COL_SPLITTER = "\t"

# region REDItoolKnown

# ---- "readability" consts ----
REDIToolsKnownParser_CHR_I = 0
REDIToolsKnownParser_POSITION_I = 1
REDIToolsKnownParser_REFERENCE_I = 2
REDIToolsKnownParser_STRAND_I = 3
REDIToolsKnownParser_COVERAGE_I = 4
REDIToolsKnownParser_RATIO_I = 6
REDIToolsKnownParser_ALLSUBS_I = 7

# ---- splicers ----
REDIToolsKnownParser_CSV_DELIM = "\t"
REDIToolsKnownParser_RATIO_STRIP_CHARS = "[]"
REDIToolsKnownParser_RATIO_DELIM = ","
REDIToolsKnownParser_SUBS_DELIM = " "
REDIToolsKnownParser_NO_SUBS = '-'

# ---- logging ----
SKIPPING_SITE_DBG_MSG = "MismatchSite was not found in keep_only for sample %s, skipping... (site chr: %s, start: %s, end: %s)"
NO_COVERAGE_SITE_DBG_MSG = "MismatchSite has not coverage for sample %s, skipping... (site chr: %s, start: %s, end: %s)"
# endregion

# **********YishayKnownAnalysis***************#
# ---- "readability" consts ----
YishayKnownAnalysisParser_CHR_I = 0
YishayKnownAnalysisParser_START_I = 1
YishayKnownAnalysisParser_END_I = 2
YishayKnownAnalysisParser_STRAND_I = 3
YishayKnownAnalysisParser_MIR_NAME_I = 4
YishayKnownAnalysisParser_PRE_MIR_POS_I = 5
YishayKnownAnalysisParser_MATURE_MIR_POS_I = 6
YishayKnownAnalysisParser_MISMATCH_I = 7
YishayKnownAnalysisParser_REFERENCE_COUNT_I = 8
YishayKnownAnalysisParser_MISMATCHED_COUNT_I = 9
YishayKnownAnalysisParser_RATIO_I = 10

# ---- splicers ----
YishayKnownAnalysisParser_CSV_DELIM = "\t"


# =====================classes=====================#


class RatesAtSitesParser(object):
    @staticmethod
    def parse_known_sites_editing_file(data, sample, source_name, *args, **kwargs):
        """
        This is the parsing function (that should be overridden differently for each input format.
        :param data: The file handler to read from.
        :return: the resulting objects.
        """
        raise NotImplementedError

    @staticmethod
    def contained(keep_only_sites, region, start, end, strand=None):
        if strand:
            return Site(region, start, end, strand) in keep_only_sites
        else:
            return (Site(region, start, end, SENSE_STRAND) in keep_only_sites or
                    Site(region, start, end, ANTISENSE_STRAND) in keep_only_sites)


class REDIToolsKnownParser(RatesAtSitesParser):
    @staticmethod
    def parse_known_sites_editing_file(data, sample, source_name, keep_only=None,
                                       default_mismatch=MismatchesAndRefsEnum.CANONICAL, *args, **kwargs):
        mismatches_sites = MismatchesSites(source=source_name)

        lines = [l.split(REDIToolsKnownParser_CSV_DELIM) for l in data.split("\n")]

        for line in lines[1:]:
            if line == [""]:
                continue
            pos = int(line[REDIToolsKnownParser_POSITION_I])
            pos_start = pos - 1
            pos_end = pos
            strand = SENSE_STRAND if int(line[REDIToolsKnownParser_STRAND_I]) else ANTISENSE_STRAND
            chromosome = line[REDIToolsKnownParser_CHR_I]

            if keep_only and not RatesAtSitesParser.contained(region=chromosome, start=pos_start, end=pos_end,
                                                              strand=strand):
                logging.debug(SKIPPING_SITE_DBG_MSG % (sample, chromosome, pos_start, pos_end))
                continue
            overall_count = int(line[REDIToolsKnownParser_COVERAGE_I])

            if overall_count <= 0:
                logging.debug(NO_COVERAGE_SITE_DBG_MSG % (sample, chromosome, pos_start, pos_end))
                continue

            a_count, c_count, g_count, t_count = \
                line[REDIToolsKnownParser_RATIO_I].strip(REDIToolsKnownParser_RATIO_STRIP_CHARS).split(
                    REDIToolsKnownParser_RATIO_DELIM)

            type_to_count = {"A": a_count, "C": c_count, "T": t_count, "G": g_count}

            reference = line[REDIToolsKnownParser_REFERENCE_I]
            subs = line[REDIToolsKnownParser_ALLSUBS_I].split(REDIToolsKnownParser_SUBS_DELIM)

            has_mismatch = any([REDIToolsKnownParser_NO_SUBS != sub for sub in subs])
            for sub in subs:
                if REDIToolsKnownParser_NO_SUBS == sub and has_mismatch:
                    continue
                elif REDIToolsKnownParser_NO_SUBS == sub:
                    mismatch_type = default_mismatch
                    edited = 0
                else:
                    mismatch_type = MismatchesAndRefsEnum.MISMATCH_TYPE_PARSE[sub]
                    edited = int(type_to_count[sub[1]])
                mismatches_sites.add_site(mismatch_type=mismatch_type, sample=sample, region=chromosome,
                                          start=pos_start, end=pos_end, strand=strand, num_of_mismatched_reads=edited,
                                          num_of_canonical_reads=overall_count - edited)

        return mismatches_sites  # to save search time.


class YishayKnownAnalysisParser(RatesAtSitesParser):
    @staticmethod
    def parse_known_sites_editing_file(data, sample, source_name, *args, **kwargs):
        mismatches_sites = MismatchesSites(source=source_name)

        lines = [l.split(REDIToolsKnownParser_CSV_DELIM) for l in data.split("\n")]

        for line in lines:
            if line == [""]:
                continue
            chrom = line[YishayKnownAnalysisParser_CHR_I]
            start = line[YishayKnownAnalysisParser_START_I]
            end = line[YishayKnownAnalysisParser_END_I]
            strand = line[YishayKnownAnalysisParser_STRAND_I]
            mir_name = line[YishayKnownAnalysisParser_MIR_NAME_I]
            pre_mir_pos = line[YishayKnownAnalysisParser_PRE_MIR_POS_I]
            mature_mir_pos = line[YishayKnownAnalysisParser_MATURE_MIR_POS_I]
            mismatch_type = line[YishayKnownAnalysisParser_MISMATCH_I]
            reference_count = int(line[YishayKnownAnalysisParser_REFERENCE_COUNT_I])
            mismatched_count = int(line[YishayKnownAnalysisParser_MISMATCHED_COUNT_I])
            ratio = float(line[YishayKnownAnalysisParser_RATIO_I]) if line[
                                                                          YishayKnownAnalysisParser_RATIO_I] != "NA" else -1
            overall_count = mismatched_count + reference_count
            if overall_count <= 0:
                logging.debug(NO_COVERAGE_SITE_DBG_MSG % (sample, mir_name, start, end))
            mismatches_sites.add_site(mismatch_type=mismatch_type, sample=sample, region=chrom, start=start, end=end,
                                      strand=strand, num_of_mismatched_reads=mismatched_count,
                                      num_of_canonical_reads=overall_count, ratio=ratio,
                                      extra_data={MATURE_MIR_POS: mature_mir_pos, PRE_MIR_POS: pre_mir_pos,
                                                  MIR_NAME: mir_name})

        return mismatches_sites  # to save search time.


class ShacharFilteredResultParser(RatesAtSitesParser):
    """
    This class is the parser for the filtered results file from Shachar's script (editing.txt)
    """
    @staticmethod
    def parse_known_sites_editing_file(data, sample, source_name, *args, **kwargs):
        lines = data.split(ShacharUnFilteredResultParser_LINE_SPLITTER)
        mismatches_sites = MismatchesSites(source=source_name)

        for line in lines[1:]:
            if '' == line:
                continue

            chromosome, start, end, strand, mir_name, location_in_pre_mir, mismatch_type, mismatched_count, \
            overall_count, location_in_mir, raw_p_value, bonferroni_p_value, bh_p_value = line.split(
                ShacharFilteredResultParser_COL_SPLITTER)
            mismatch_type = MismatchesAndRefsEnum.MISMATCH_TYPE_PARSE[mismatch_type]
            start = int(start)
            end = int(end)
            mismatched_count = int(mismatched_count)
            overall_count = int(overall_count)
            raw_p_value = float(raw_p_value)
            bonferroni_p_value = float(bonferroni_p_value)
            bh_p_value = float(bh_p_value)

            stats_d = StatisticalTestPValue()
            stats_d[GENERAL_TEST] = raw_p_value
            stats_d[BONFERRONI_P_TEST] = bonferroni_p_value
            stats_d[BH_P_TEST] = bh_p_value

            mismatches_sites.add_site(mismatch_type=mismatch_type, sample=sample, region=mir_name, start=start, end=end,
                                      strand=strand, num_of_mismatched_reads=mismatched_count,
                                      num_of_canonical_reads=overall_count, ratio=ratio,
                                      extra_data={MATURE_MIR_POS: location_in_mir, PRE_MIR_POS: location_in_pre_mir,
                                                  MIR_NAME: mir_name, STATISTICS_DATA: stats_d})

        return mismatches_sites  # to save search time.


class ShacharUnFilteredResultParser(RatesAtSitesParser):
    """
    This class is the parser for the unfiltered results file from Shachar's script (test.txt)
    """
    @staticmethod
    def parse_known_sites_editing_file(data, sample, source_name, keep_only=None, *args, **kwargs):
        assert isinstance(keep_only, (None, MismatchesSites))
        mismatches_sites = MismatchesSites(source_name)
        lines = data.split(ShacharUnFilteredResultParser_LINE_SPLITTER)

        # due to the need to iterate with different advancing rate (sometime one line sometime two),
        #  an external index is used to determine position.
        i = 0
        while i < len(lines):
            if lines[i] == "":
                i += 1
                continue
            mir_name, editing_pos, mismatch_type, edited_count, overall_count = lines[i].split(
                ShacharUnFilteredResultParser_COL_SPLITTER)

            editing_pos = int(editing_pos)
            start = editing_pos - 1
            end = editing_pos
            edited_count = int(edited_count)
            overall_count = int(overall_count)
            mismatch_type = MismatchesAndRefsEnum.MISMATCH_TYPE_PARSE[mismatch_type]
            if overall_count < 1:
                logging.debug(NO_COVERAGE_SITE_DBG_MSG % (sample, mir_name, start, end))
                i += 1
                continue

            if keep_only and not RatesAtSitesParser.contained(region=mir_name, start=start, end=end):
                    logging.debug(SKIPPING_SITE_DBG_MSG % (sample, mir_name, start, end))
                    i += 1
                    if edited_count > 0:
                        i += 1
                    continue

            if edited_count > 0:
                i += 1
                p_val = lines[i].split(ShacharUnFilteredResultParser_COL_SPLITTER)[-1].strip(
                    ShacharUnFilteredResultParser_P_VALUE)
                p_val = StatisticalRecordDescriptor(GENERAL_TEST, float(p_val))
                stats_d = StatisticalTestPValue()
                stats_d[GENERAL_TEST]= p_val

            mismatches_sites.add_site(mismatch_type=mismatch_type, sample=sample, region=mir_name, start=start, end=end,
                                      strand=NA, num_of_mismatched_reads=edited_count,
                                      num_of_canonical_reads=overall_count, ratio=ratio,
                                      extra_data={MATURE_MIR_POS: NA, PRE_MIR_POS: NA,
                                                  MIR_NAME: mir_name, STATISTICS_DATA: stats_d})

            i += 1

        return mismatches_sites

