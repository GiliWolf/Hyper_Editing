__author__ = 'Hillel'
# =====================imports=====================#
from Outputer import Outputer
from json import dump
import logging
import os

# =====================constants===================#
EMPTY_VAL = "-"
NOT_AVAILABLE = "NA"

OVERRIDE_APPEND_CONFLICT = "Cannot Init CSVOutputer with both override and append on True!"


# =====================classes=====================#


class JSONOutputer(Outputer):
    """
    This class implements creating XML output
    """

    def output(self, paths, output_dict):
        """
        This function creates the output.
        :param list[str] paths: The paths to output to.
        :param dict output_dict: The records to write, header name => value
        :return: None
        """

        mode = 'wb'
        for path in paths:
            if os.path.isfile(path):
                pass

            self.make_dir(path)

            with open(path, mode) as o:
                dump(obj=output_dict, fp=o)

            self.chmod_to_all(path)
