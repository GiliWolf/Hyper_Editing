__author__ = 'Hillel'
# =====================imports=====================#
from csv import reader
from collections import namedtuple
import logging

# =====================constants===================#
HyperEditingBEDRec = namedtuple("HyperEditingBEDRec", "chr start end neighbors reads_count strand")

# error messages
FILE_OPEN_ERR_MSG = "Failed to Open Hyper Editing Detect File %s"


# =====================classes=====================#


class HyperEditingBEDConverter(object):
    @staticmethod
    def convert(he_bed_path):
        recs = list()
        try:
            he_bed_file = open(he_bed_path, "rb")
        except IOError:
            logging.exception(FILE_OPEN_ERR_MSG % he_bed_path)
            return None

        for line in reader(he_bed_file, delimiter="\t"):
            if '' == line:
                continue
            recs.append(HyperEditingBEDRec(*line))

        return recs
