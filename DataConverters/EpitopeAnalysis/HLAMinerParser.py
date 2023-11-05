__author__ = 'Hillel'
# =====================imports=====================#
import re
from numpy import mean
from collections import namedtuple

from DataConverters.DataConverter import DataConverter


# =====================constants===================#
HLA_RE = re.compile(r"((HLA-\w*).*?)(?=$|(?:HLA)|--)", re.DOTALL)
PREDICTION_RE = re.compile(r"Prediction #\d+ - (.*?)\n(.*?)(?=(?:\tPrediction)|$|(?:HLA)|--)", re.DOTALL)
PREDICTIONS_RE = re.compile(r"\t\t.*?,.*?,(.*?),.*?(?:\n|$)")

SAME_SCORE = " (same score as above)"

CLASS_I = "ClassI"
CLASS_II = "ClassII"

CLASS_I_HLAS = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", ]
# error messages


# =====================classes=====================#
HLAPrediction = namedtuple("HLAPrediction", "allele scores mean_score")
HLAClassPrediction = namedtuple("HLAPrediction", "hla_class predictions")


class HLAMinerConverter(DataConverter):
    def convert(self, file_handler, *args, **kwargs):
        """
        This conversion converts an HLAMiner output returning all predictions by thier scores.
        :param file_handler: The file hanlder of the output.
        :param args: Currently not used.
        :param kwargs: Currently not used.
        :return: A dictionary of the prediction by their class and scores: <class I/II> => <hla class> =>
        <predictions by score>
        :rtype: C{dict} of C{str} to C{dict} of C{str} to L{HLAClassPrediction}
        """
        res = {}
        for hla_predictions, hla_class in HLA_RE.findall(file_handler.read()):
            class_predictions = {}
            i = 1
            for allele, predictions in PREDICTION_RE.findall(hla_predictions):
                scores = []
                prediction = class_predictions.setdefault(str(i), [])
                for prediction_score in PREDICTIONS_RE.findall(predictions):
                    scores.append(float(prediction_score))

                mean_score = str(mean(scores))
                prediction.append(HLAPrediction(allele=allele.replace(SAME_SCORE, ""), scores=[str(s) for s in scores],
                                                mean_score=mean_score))
                i += 1 if SAME_SCORE not in allele else 0
            if class_predictions == {}:
                continue
            class_type = CLASS_I if hla_class in CLASS_I_HLAS else CLASS_II
            res.setdefault(class_type, {})[hla_class] = HLAClassPrediction(hla_class=hla_class,
                                                                           predictions=class_predictions)
        return res



