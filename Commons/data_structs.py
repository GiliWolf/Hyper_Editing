__author__ = 'Hillel'
"""
This module contains simple DS commonly used across various classes.
recommended usages are found in each class' docs.
"""
# =====================imports=====================#
from collections import namedtuple
import argparse
from operator import attrgetter


# =====================classes=====================#


# noinspection PyClassHasNoInit
class StatisticalTestPValue(dict):
    """
    This class holds the data about a statistical test related to an information.
    Keys should be statistical test name as found in GGPSResources.consts, and the values - the p values.
    """


# noinspection PyClassHasNoInit
class Site(namedtuple("Site", "region start end strand")):
    """
    This hold the coordinates for a mismatch site.
    :ivar str mismatch_type: The type of the mismatch. *only from from MismatchesTypesEnum!*
    :ivar str chromosome: The chromosome (or transcript name)
    :ivar int start: start position
    :ivar int end: end position
    :ivar str strand: the strand of the site
    """
    __slots__ = ()

    def __str__(self):
        return "Coordinate: %s:%s-%s; Strand: %s" % (self.region, self.start, self.end, self.strand)


# noinspection PyClassHasNoInit
class AlignmentData(namedtuple("AlignmentData", "total_reads uniquely_mapped multi_mapped unmapped "
                                                "sequencing_method")):
    """
    This hold the coordinates for a mismatch site.
    :ivar int total_reads: The number of reads  in the sample.
    :ivar int uniquely_mapped: The number of uniquely mapped reads in the sample.
    :ivar int multi_mapped: The number of multi mapped reads in the sample.
    :ivar int unmapped: The number of unmapped reads in the sample.
    :ivar str sequencing_method: (paired end)PE/(single end)SE
    """
    __slots__ = ()
    SEQ_METHOD_PE = "PairedEnd"
    SEQ_METHOD_SE = "SingleEnd"


# noinspection PyClassHasNoInit
class MHCAllele(namedtuple("MHCAllele", "mhc_class locus allele")):
    """
    This holds an allele of MHC
    :ivar str mhc_class: the class of the mhc.
    :ivar str locus: the locus (name) of the allele (e.g. HLA-A).
    :ivar str allele: the full allele name (e.g. A*01:21).
    """
    __slots__ = ()


# noinspection PyClassHasNoInit
class AffinitiesEnum:
    """
    Contains the values for the affinities ranking - weak, intermediate, and strong.
    """
    WEAK = "Weak"
    INTERMEDIATE = "Intermediate"
    STRONG = "Strong"


# noinspection PyClassHasNoInit
class KeepOnlyEnum(object):
    """
    This class contain 3 values for combining values - keep left operand (LEFT), right operand (RIGHT), or combine them (BOTH)
    """
    LEFT = 1
    RIGHT = 2
    BOTH = 3


# noinspection PyClassHasNoInit
class PredicateAnswersEnum(object):
    """
    This class contain 3 values for combining values - keep left operand (LEFT), right operand (RIGHT), or combine them (BOTH)
    """
    YES = 1
    NO = 2
    UNKNOWN = 3


# noinspection PyClassHasNoInit
class SampleCoverageData(namedtuple("SampleCoverageData", "sample site reference total adenosines cytosines guanosines "
                                                          "thymines containing_region")):
    """
    This holds the mismatch rate for a sample for a sites.
    :ivar Sample sample: the instance of the L{Sample} this is relating to.
    :ivar Site site: the site in the genome.
    :ivar str reference: the reference in the genome.
    :ivar int total: the total coverage not including Ns and low quality reads.
    :ivar int adenosines: the number of adenosines.
    :ivar int cytosines: the number of cytosines.
    :ivar int guanosines: the number of guanosines.
    :ivar int thymines: the number of thymines.
    :ivar Site containing_region: the genomic region containing this position.
    """
    __slots__ = ()


# noinspection PyClassHasNoInit
class RefSeqPosEnum(str):
    """
    Enum for position in genome.
    """
    __slots__ = ()

    EXONIC = "Exonic"
    INTERGENIC = "Intergenic"
    INTRONIC = "Intronic"

    ALL_ANNOTATIONS = [EXONIC, INTRONIC, INTERGENIC]

    __HIERARCHY_D = {EXONIC: 3, INTRONIC: 2, INTERGENIC: 1}
    __HIERARCHY_D_REVERSE = {3: EXONIC, 2: INTRONIC, 1: INTERGENIC}

    @staticmethod
    def stronger_annotation(this, other):
        return RefSeqPosEnum.__HIERARCHY_D_REVERSE[
            max(RefSeqPosEnum.__HIERARCHY_D[this], RefSeqPosEnum.__HIERARCHY_D[other])]


# noinspection PyClassHasNoInit
class RefSeq(namedtuple("RefSeq", "site exons refseq_name common_name")):
    """
    This holds the mismatch rate for a sample for a sites.
    :ivar Site site: the site of the SNP.
    :ivar dict[int, int] exons: The exons indexes (start=>end, relative to '+' strand).
    :ivar str refseq_name: The gene RefSeq name (e.g. NR_075077.1).
    :ivar str common_name: The gene common name (e.g. C1orf141).
    """
    __slots__ = ()


# noinspection PyClassHasNoInit
class SNP(namedtuple("SNP", "site prevalence function reference_seq snp_seq")):
    """
    This holds the mismatch rate for a sample for a sites.
    :ivar Site site: the site of the SNP.
    :ivar float prevalence: The prevalence of the reference.
    :ivar str function: the function of this region (e.g. intergenic)
    :ivar str reference_seq: the reference sequence
    :ivar str snp_seq: the SNP sequence
    """
    __slots__ = ()


class SortingHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)
