__author__ = 'Hillel'
"""This file contain the parser for optitype output"""
# =====================imports=====================#
from csv import reader
import os
import re
from datetime import datetime
import logging

from Commons.consts import IFILE_CONVERT_ERR, IFILE_READ_ERR, IFILE_UNKNOWN_GROUP_WRN_MSG, IFILE_UNKNOWN_SAMPLE_WRN_MSG,\
    MHC_CLASS_I, OBJECTIVE_VALUE
from Commons.data_structs import StatisticalTestPValue
from Commons.general_functions import get_groups_and_sample_names_dict
from DataObjects.KnownProperties.MHCAlleles import MHCAlleles

# =====================consts=====================#
# path format
OPTITYPE_FILENAME_RE = re.compile(r"(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_result.tsv")

# output format
TOPMOST_PREDICTION_LINE_I = 2
OPTITYPE_SEP = "\t"
# alleles indexes in the line
A1_I = 1
A2_I = 2
B1_I = 3
B2_I = 4
C1_I = 5
C2_I = 6
OBJECTIVE_VAL_I = 7

HLA_LOCUS_A = "HLA_A"
HLA_LOCUS_B = "HLA_B"
HLA_LOCUS_C = "HLA_C"

# =====================classes=====================#

class OptitypeConverter:

    OPTITYPE_SOURCE_NAME = "OptitypePredictions"

    @staticmethod
    def convert_optitype_output(file_names, sample_name_path_i, top_group_name_path_i):
        # type: (list, int, int) -> bool
        """
        This function converts optitype outputs.
        :param file_names: The paths of the file to convert
        :param sample_name_path_i: The index of the sample name in the path
        :param top_group_name_path_i: The index of the sample name in the path
        :return: samples found
        """
        samples = []
        mhc_allels = MHCAlleles(OptitypeConverter.OPTITYPE_SOURCE_NAME)
        group_and_sample_names_dict = get_groups_and_sample_names_dict()

        samples_to_timestamp = dict()
        for ifile in file_names:
            timestamp = OPTITYPE_FILENAME_RE.findall(os.path.basename(ifile))
            if not timestamp:  # shouldn't ever reach here
                logging.exception(IFILE_CONVERT_ERR % ifile)
                continue
            # current version of optitype output foramt. change if needed
            timestamp = [int(val) for val in timestamp[0]]
            timestamp = datetime(timestamp[0], timestamp[1], timestamp[2], timestamp[3], timestamp[4], timestamp[5])

            try:
                sample_name = ifile.split(os.path.sep)[sample_name_path_i]
                group_name = ifile.split(os.path.sep)[top_group_name_path_i]
                group_d = group_and_sample_names_dict.get(group_name, None)
                if None is group_d:
                    logging.warn(IFILE_UNKNOWN_GROUP_WRN_MSG % (ifile, group_name))
                    continue

                sample = group_d.get(sample_name, None)
                if None is sample:
                    logging.warn(IFILE_UNKNOWN_SAMPLE_WRN_MSG % (ifile, sample))
                    continue

                if sample.record_id in samples_to_timestamp and samples_to_timestamp[sample.record_id] > timestamp:
                    continue  # this is an older run
                samples_to_timestamp[sample.record_id] = timestamp
                samples.append(sample)
                ifile_handler = open(ifile, 'rb')
                data = [line for line in reader(ifile_handler, delimiter=OPTITYPE_SEP)][TOPMOST_PREDICTION_LINE_I]
                obj = StatisticalTestPValue()
                obj[OBJECTIVE_VALUE] = float(data[OBJECTIVE_VAL_I])
                mhc_allels.add_allele(sample, MHC_CLASS_I, HLA_LOCUS_A, data[A1_I], obj)
                mhc_allels.add_allele(sample, MHC_CLASS_I, HLA_LOCUS_A, data[A2_I], obj)
                mhc_allels.add_allele(sample, MHC_CLASS_I, HLA_LOCUS_B, data[B1_I], obj)
                mhc_allels.add_allele(sample, MHC_CLASS_I, HLA_LOCUS_B, data[B2_I], obj)
                mhc_allels.add_allele(sample, MHC_CLASS_I, HLA_LOCUS_C, data[C1_I], obj)
                mhc_allels.add_allele(sample, MHC_CLASS_I, HLA_LOCUS_C, data[C2_I], obj)
            except IOError:
                logging.exception(IFILE_READ_ERR % ifile)
            except:
                logging.exception(IFILE_CONVERT_ERR % ifile)
        return samples

