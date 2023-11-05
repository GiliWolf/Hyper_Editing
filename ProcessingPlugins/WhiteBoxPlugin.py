__author__ = 'Hillel'
# =====================imports=====================#
import abc

# =====================constants===================#
# error messages
OUTPUTS_TYPE_ERR = "Outputs Must Be Instances Of Outputs.Outputer!"

# =====================classes=====================#


class WhiteBoxPlugin(object):
    """
    This is the base class for all "white box" plugins (and is abstract) - an interface class.
    It contains only one virtual method - process_input, to allow polymorphic activating of the plugins.
    """
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def process_input(self, *args, **kwargs):
        raise NotImplementedError()

