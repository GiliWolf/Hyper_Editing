__author__ = 'Hillel'
# =====================imports=====================#
from Outputer import Outputer

import os
from pprint import pprint
# =====================classes=====================#


class RawOutputer(Outputer):
    """
    This class implement creating raw output (i.e. writing raw data into a file)
    """
    def __init__(self, override=True, pprint_flag=True):
        self.override = override
        self.pprint = pprint_flag

    def output(self, paths, raw_data):
        """
        This function creates the output.
        :return: None
        """
        for path in paths:
            if os.path.isfile(path) and not self.override:
                continue

            self.make_dir(path)

            with open(path, 'wb') as o:
                if self.pprint:
                    pprint(raw_data, o)
                else:
                    o.write(raw_data)

            self.chmod_to_all(path)
