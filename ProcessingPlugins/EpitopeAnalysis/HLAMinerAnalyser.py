__author__ = 'Hillel'
# =====================imports=====================#
import argparse
import logging
import os

if __name__ == "__main__":
    from sys import path
    path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from Outputs.Outputers.CSVOutputer import CSVOutputer
from Commons.general_functions import init_logging_dict
from Commons.help_functions import get_input_files
from Commons.consts import GROUP, SAMPLE, MHC_CLASS, MHC_LOCUS, MHC_ALLELE
from ProcessingPlugins.WhiteBoxPlugin import WhiteBoxPlugin
from DataConverters.EpitopeAnalysis.HLAMinerParser import HLAPrediction, HLAMinerConverter


# =====================constants===================#
DEFAULT_HLAMINER_OUTPUT_NAME = "HLAminer_HPTASR.csv"

# ----output consts-----

CLASS_TYPE = MHC_CLASS
ALLELE = MHC_ALLELE
CLASS = MHC_LOCUS
PREDICTION_RANK = "PredictionRank"
SCORES = "Scores"
MEAN_SCORE = "MeanScore"
OUTPUT_HEADERS = [GROUP, SAMPLE, CLASS_TYPE, CLASS, PREDICTION_RANK, ALLELE, MEAN_SCORE, SCORES]

SCORES_SEP = ";"
OUTPUT_FILENAME = "HLAMinerAnslysisSummery.csv"
# ----logging consts-----
LOG_FILE = r'hla_miner_analyser.log'

# error messages
IFILE_READ_ERR = "Could Not Open Input File %s!"
IFILE_CONVERT_ERR = "Could Not Convert Input File %s!"

# =====================classes=====================#


class HLAMinerAnalyser(WhiteBoxPlugin):
    def __init__(self, instance_nickname="HLAMinerAnalyser", log_file=None):
        super(HLAMinerAnalyser, self).__init__(instance_nickname)
        self.converter = HLAMinerConverter()
        init_logging_dict(log_file=os.path.join(args.output_dir, log_file if log_file else LOG_FILE),
                          instance_nickname=instance_nickname)

    def process_input(self, input_dir, hla_files_names=[DEFAULT_HLAMINER_OUTPUT_NAME], output_dir=".",
                      sample_name_path_i=None, group_name_path_i=None):
        inputs = get_input_files(root_dir=input_dir, fingerprints=hla_files_names)

        c_predictions = {}
        for ifile in inputs:
            try:
                ifile_handler = open(ifile, 'rb')
                c_predictions[ifile] = self.converter.convert(ifile_handler)
            except IOError:
                logging.exception(IFILE_READ_ERR % ifile)
            except:
                logging.exception(IFILE_CONVERT_ERR % ifile)

        out_recs = []
        rec = {}
        for s_path, c_predict in c_predictions.iteritems():
            path_parts = s_path.split(os.path.sep)
            if group_name_path_i:
                rec[GROUP] = path_parts[group_name_path_i]
            rec[SAMPLE] = path_parts[sample_name_path_i] if sample_name_path_i else s_path

            for class_type in c_predict:
                rec[CLASS_TYPE] = class_type
                for hla_class in c_predict[class_type]:
                    rec[CLASS] = hla_class

                    for rank, predictions in c_predict[class_type][hla_class].predictions.iteritems():
                        rec[PREDICTION_RANK] = rank
                        for prediction in predictions:
                            if not isinstance(prediction, HLAPrediction):
                                continue
                            rec[ALLELE] = prediction.allele
                            rec[MEAN_SCORE] = prediction.mean_score
                            rec[SCORES] = SCORES_SEP.join(prediction.scores)

                            out_recs.append(rec.copy())

        outputer = CSVOutputer()
        outputer.output(paths=[os.path.join(output_dir, OUTPUT_FILENAME)], headers=OUTPUT_HEADERS, records=out_recs)



if __name__ == '__main__':
    desc = 'A Script for extracting the predicted hla alleles from HLMiner output'


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='HLAMinerAnalyser', description=desc, formatter_class=MyFormatter)
    parser.add_argument('-i', '--input_dir', metavar="input_path",dest='input_path', nargs='?', required=True,
                        help='The dir containing the HLAMiner files to run on.')
    parser.add_argument('-f', '--fingerprints', metavar="fingerprints", dest='fingerprints', nargs='*', required=True,
                    help='name parts to recognize files for input.')
    parser.add_argument('-o', '--output_dir', metavar="outputdir", dest='output_dir', nargs='?', required=True,
                    help='The path of the output dir to create.')
    parser.add_argument('-sni', '--sample_name_path_i', metavar="sample name path index", dest='sample_name_path_i',
                        nargs='?', required=False, type=int, default=None, help="The index ofthe sampe name in the path (counting from 0)."
                                                    " e.g. /drive/folder/group/sample_name/file, sample name i will be either -2 or 3")
    parser.add_argument('-gni', '--group_name_path_i', metavar="group name path index", dest='group_name_path_i',
                        nargs='?', required=False, type=int, default=None, help="The index ofthe group name in the path (counting from 0)."
                                                    " e.g. /drive/folder/group/sample_name/file, group name i will be either -3 or 2")

    args = parser.parse_args()

    analyser = HLAMinerAnalyser()
    analyser.process_input(input_dir=args.input_path, hla_files_names=args.fingerprints, output_dir=args.output_dir,
                      sample_name_path_i=args.sample_name_path_i, group_name_path_i=args.group_name_path_i)
