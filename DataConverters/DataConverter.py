__author__ = 'Hillel'
# =====================imports=====================#
import abc
from collections import namedtuple

# =====================constants===================#
# error messages

# =====================classes=====================#


ConvertResults = namedtuple("ConvertResults", "loci known_properties genomic_records sample")


class DataConverter(object):
    """
    This is the base class for the data converters.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def convert(self, file_handler, *args, **kwargs):
        """
        This function is the converting function.
        :param file_handler: The file handler of the input to convert.
        :param sample_id: The ID of the sample instance the file is related to.
        :param args: optional args.
        :param kwargs: optional key word args.
        :return: The corresponding data object to the input.
        """