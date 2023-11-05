__author__ = 'Hillel'
# =====================imports=====================#
import re
import os
from csv import reader
from collections import OrderedDict

import TCGA_consts
from consts import *

# =====================constants=====================#

TCGA_HEADER_RE = re.compile(r"((.*?)-(.*?)-(.*?)-(\d{2})(\w)-(\d{2})(\w)-(.*?)-(.*?))(?:_|$)")

SAMPLE_NAME = 'sample_name'
PROJECT = 'project'
TISSUE_SOURCE = 'tissue_source'
PARTICIPANT_NUM = 'participant_num'
SAMPLE_TYPE = 'sample_type'
VIAL = 'vial'
PORTION = 'portion'
ANALYTE = 'analyte'
PLATE_NAME = 'plate_name'
SEQ_CENTER = 'seq_center'

TCGA_HEADERS_LIST = [SAMPLE_NAME, PROJECT, TISSUE_SOURCE, PARTICIPANT_NUM, SAMPLE_TYPE, VIAL,
                     PORTION, ANALYTE, PLATE_NAME, SEQ_CENTER]


REVERSE_BASES_DICT = dict(A='T', C='G', G='C', T='A', a='t', c='g', g='c', t='a')

POSSIBLE_SEPERATOR = ["\t", ',']
# =====================functions=====================#


def convert_args_to_dict(arg):
    assert all([KW_ARGS_RE.match(a) for a in arg.split(",")])
    return eval("dict(" + arg + ")")


def parse_TCGA_sample_name(sample_name):
    """
    This function parse the name of TCGA samples and returns the represented fields.
    :param sample_name: The name of the sample, e.g. TCGA-A3-3358-01A-01R-1540-13.
    :return: A dictionary containing the fields,
     according to what written in the site (wiki.nci.nih.gov/display/TCGA/TCGA+Barcode)
    """
    search_obj = TCGA_HEADER_RE.search(sample_name)

    if search_obj is None:
        return None

    sample_name, project, tissue_source, participant_num, sample_type, vial, \
    portion, analyte, plate_name, seq_center = search_obj.groups()

    return {
        SAMPLE_NAME: sample_name,
        PROJECT: project,
        TISSUE_SOURCE: TCGA_consts.SOURCE_SITE_TRANSLATE[tissue_source],
        PARTICIPANT_NUM: participant_num,
        SAMPLE_TYPE: TCGA_consts.SAMPLE_D_TYPE[sample_type],
        VIAL: vial,
        PORTION: portion,
        ANALYTE: TCGA_consts.ANALYTE[analyte],
        PLATE_NAME: plate_name,
        SEQ_CENTER: seq_center
    }


def get_sra_sample_data(sra_path, path_fragment_indexes):
    fragments = sra_path.split(os.sep)[:-1]

    res = []

    for i in path_fragment_indexes:
        res.append(fragments[i])

    res.reverse()
    return res


def reverse_strand(seq, complementary=True, ignore_non_bases=False):
    rev = ""
    for base in seq:
        def_res = base if ignore_non_bases else "?"
        rev += REVERSE_BASES_DICT.get(base, def_res)

    if complementary:
        rev = "".join(reversed(rev))
    return rev


def get_input_files(root_dir, fingerprints, excludes=[], require_all=False):
    """
    This function returns a list of file names containing one (or all) of the fingerprints provided
    :param root_dir: The root dir to run on.
    :param fingerprints: a C{list} of C{str}, strings wanted in the path.
    :param excludes: a C{list} of C{str}, strings *not* wanted in the path. (checked with OR operand)
    :param require_all: If True will require all <fingerprints> to be in the path, otherwise any will suffice.
    :return: the paths.
    :rtype C{list} of C{str}
    """
    inputs = []
    predicate = all if require_all else any
    for root, dirs, files in os.walk(root_dir):
        for ifile in files:
            if any([ex in root for ex in excludes]):
                continue
            if any([ex in ifile for ex in excludes]):
                continue
            if predicate([fingerprint in ifile for fingerprint in fingerprints ]):
                inputs.append(os.path.join(root, ifile))
    return inputs


def determine_sep(sample_str, known_length=None):
    """
    This function determines if the seperator in a file is tab or comma.
    :param sample_str: a string sample from the file (e.g. a line)
    :param known_length: if given will be used as sanity check to make sure number of recs got is identical.
    :return:  the seperator (either \t or ,)
    """
    for sep in POSSIBLE_SEPERATOR:
        if sep in sample_str:
            if known_length:
                if len(sample_str.split(sep)) == known_length:
                    return sep
            else:
                return sep


