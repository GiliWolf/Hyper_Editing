# =====================imports=====================#
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from numpy import average, mean, std, var, median
from csv import reader
import re
import operator
from collections import namedtuple
from Outputs.Outputers.CSVOutputer import CSVOutputer
from ProcessingPlugins.EpitopeAnalysis.AffinityAssement import get_hla_allels_per_sample, CHR, POSITION, \
    SUM_ORIG_BIND_HEADER_S, SUM_ORIG_BIND_HEADER_I, SUM_ORIG_BIND_HEADER_W, SUM_NEO_BIND_HEADER_S, BIND_DIFF_HEADER, \
    SUM_NEO_BIND_HEADER_I, SUM_NEO_BIND_HEADER_W, transform_allele
# from ProcessingPlugins.EpitopeAnalysis.AffinityAssement import get_hla_allels_per_sample, CHR, \
#     POSITION, BIND_DIFF_HEADER, SUM_REC_SEP, HLA_HEADER, BIND_DIFF_PC_HEADER, SUM_NEO_BIND_HEADER_S, \
#     SUM_NEO_BIND_HEADER_I, SUM_NEO_BIND_HEADER_W, SUM_ORIG_BIND_HEADER_S, SUM_ORIG_BIND_HEADER_I, SUM_ORIG_BIND_HEADER_W
from Commons.data_structs import Site, MHCAllele
from Commons.consts import NA, MHC_CLASS_I
from Commons.general_functions import convert_params_to_bool_dict

from Commons.general_functions import load_groups_and_samples_from_file, get_groups_and_sample_names_dict
from DataObjects.KnownProperties.MismatchesSites import MismatchesSites, SampleMMRate
from DataObjects.BaseClasses.Sample import Sample
from DataObjects.DataRecords.EpitopesNumbers import EpitopesNumbers, EpitopesNumber

from Outputs.mismatches_sites_functions import get_mismatches_sites_quantification_per_site, GENERAL_HEADERS_KEY

# =====================constants=====================#


EPITOPE_NUMBER_SOURCE_NAME = "EpitopesNumberSourceName"
ORIG_REAL = "OrigReal"
NEO_REAL = "NeoReal"
ORIG_RAND = "OrigRand"
NEO_RAND = "NeoRand"

GROUP_SIZE = "GroupSize"
NUMOF = "NumOf"
ALLELES_COUNT = "NumOfAlleles"
NORM_NUMOF = "NormNumOf"
NON_ZERO_EFPMS = "NonZero"
ALL_EFPMS = "All"
AVG = "Average"
MEAN = "Mean"
MED = "Median"
DEV = "Deviation"
VAR = "Variance"
GROUP = "Group"
SAMPLE = "Sample"
EFMP_FILTER = "EFMPSType"
EFPM = "EFPMs"
PEP_TYPE = "PeptideType"
AFF_THRESHOLD = "AffinityThreshold"
EPI_NUM_HEADERS = [PEP_TYPE, GROUP, SAMPLE, NUMOF, ALLELES_COUNT]
EPI_NUM_RAW_HEADERS = [PEP_TYPE, GROUP, SAMPLE, NUMOF, CHR, POSITION]
SAMPLE_FILENAME = "PerSamplesNumOFEpitopesStats.csv"
GROUP_FILENAME = "PerGroupNumOFEpitopesStats.csv"
RAW_FILENAME = "NumOFEpitopesRaw.csv"
EFMPsStats = namedtuple("EFMPsStats", "average mean median deviation variance numof norm_numof")


# =====================functions=====================#
def get_sites(sites_bed):
    sites = {}
    with open(sites_bed) as bed:
        for line in reader(bed):
            chrom = line[0]
            pos = line[1]
            sites[chrom + pos] = True
    return sites


def sum_epitopes_nums(hla_per_samples):
    # type: (dict, dict, dict) -> dict
    epi_nums = {}
    samples_allele_count = {}
    hlas_per_pos = EpitopesNumbers(source_name=EPITOPE_NUMBER_SOURCE_NAME).get_epitopes_numbers(by_allele=True)

    for mhc_allele, sites_d in hlas_per_pos.iteritems():
        relevant_samples_d = hla_per_samples.get(mhc_allele, {})
        for site, epitopes_numbers_d in sites_d.iteritems():
            for pep_type in epitopes_numbers_d:
                for group, samples in relevant_samples_d.iteritems():
                    for sample_name in samples:
                        samples_allele_count.setdefault(sample_name, set()).add(mhc_allele)
                        epi_nums.setdefault(pep_type, {}).setdefault(group, {}).setdefault(sample_name, 0)
                        epi_nums[pep_type][group][sample_name] = epi_nums[pep_type][group][sample_name] +\
                                                                 epitopes_numbers_d[pep_type].num_of_epitopes
    return epi_nums, samples_allele_count

def parse_num_of_epitopes(num_of_epi_file):
    epitopes_numbers = EpitopesNumbers(source_name=EPITOPE_NUMBER_SOURCE_NAME)
    with open(num_of_epi_file) as epi_file:
        data = [line for line in reader(epi_file)]

        chr_i = data[0].index("Chr")
        pos_i = data[0].index("Position")
        allele_i = data[0].index("Allele")
        orig_real_i = data[0].index("Eps_Orig_Real")
        neo_real_i = data[0].index("Eps_Neo_Real")
        orig_rand_i = data[0].index("Eps_Orig_Rand")
        neo_rand_i = data[0].index("Eps_Neo_Rand")

        for line in data[1:]:
            chrom = line[chr_i]
            pos = int(line[pos_i])
            start = pos - 1
            end = pos
            allele_str = line[allele_i]
            allele_l, allele_num = allele_str.split("_")
            allele_locus = "HLA_" + allele_l
            allele = allele_l + "*" + allele_num[:2] + ":" + allele_num[2:]
            orig_real_n = float(line[orig_real_i])
            neo_real_n = float(line[neo_real_i])
            orig_rand_n = float(line[orig_rand_i])
            neo_rand_n = float(line[neo_rand_i])

            epitopes_numbers.add_epitopes_number(chrom, start, end, "+", MHC_CLASS_I, allele_locus, allele,
                                                 orig_real_n, ORIG_REAL)
            epitopes_numbers.add_epitopes_number(chrom, start, end, "+", MHC_CLASS_I, allele_locus, allele,
                                                 orig_rand_n, ORIG_RAND)
            epitopes_numbers.add_epitopes_number(chrom, start, end, "+", MHC_CLASS_I, allele_locus, allele,
                                                 neo_real_n, NEO_REAL)
            epitopes_numbers.add_epitopes_number(chrom, start, end, "+", MHC_CLASS_I, allele_locus, allele,
                                                 neo_rand_n, NEO_RAND)


def sum_epitopes(num_of_epi_file, hla_to_samples_file, out_dir, sites_bed=None):
    groups_to_samples = get_groups_and_sample_names_dict()

    if sites_bed:
        limit_to_sites = get_sites(sites_bed)

    parse_num_of_epitopes(num_of_epi_file, )

    hla_per_samples = get_hla_allels_per_sample(hla_to_samples_file)

    epi_nums, samples_allele_count = sum_epitopes_nums(hla_per_samples)

    create_files(epi_nums, groups_to_samples, samples_allele_count, out_dir)


def get_groups_for_out(groups_to_samples, sample):
    groups = []
    for g, samples in groups_to_samples.iteritems():
        if sample in samples:
            groups.append(g)
    return groups


def create_files(epi_nums, groups_to_samples, samples_allele_count, out_dir):
    outpter = CSVOutputer()

    epi_n_headers = EPI_NUM_HEADERS[:]
    epi_n_headers_g = epi_n_headers[:] + [GROUP_SIZE,]
    epi_n_headers_g.pop(epi_n_headers_g.index(SAMPLE))  # pop sample out

    s_recs = []
    g_recs = []

    for pep_type, group_dict in epi_nums.iteritems():
        for group, samples_d in group_dict.iteritems():
            try:
                group_samples = str(len(groups_to_samples[group]))
            except:
                import pdb;pdb.set_trace()
            group_n_epi = 0
            group_allele_c = 0
            for sample in samples_d:
                group_n_epi += samples_d[sample]
                group_allele_c += len(samples_allele_count[sample])
                rec = {PEP_TYPE: pep_type,
                       SAMPLE: sample,
                       GROUP: group,
                       NUMOF: str(samples_d[sample]),
                       ALLELES_COUNT: str(len(samples_allele_count[sample]))
                       }
                s_recs.append(rec)
            rec = {PEP_TYPE: pep_type,
               GROUP: group,
               NUMOF: str(group_n_epi),
               ALLELES_COUNT: str(group_allele_c),
               GROUP_SIZE: group_samples
               }
            g_recs.append(rec)

    outpter.output([os.path.join(out_dir, SAMPLE_FILENAME)], epi_n_headers, s_recs)

    outpter.output([os.path.join(out_dir, GROUP_FILENAME)], epi_n_headers_g, g_recs)



if __name__ == "__main__":
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

    num_of_epi_file, hla_to_samples_file, out_dir, groups_to_samples_file = sys.argv[1:]
    load_groups_and_samples_from_file(groups_file=groups_to_samples_file)

    sum_epitopes(num_of_epi_file, hla_to_samples_file, out_dir)

"""
CMD examples:
D:\Hillel\Dropbox\Genomics\Scripts\GGPS\ProcessingPlugins\EpitopeAnalysis\pep_aff_exp_script.py
 "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\reads_counts\HESummery.csv"
 "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\reads_counts\filtered_csv_reads_coun.AGOnly.NonSyn.csv"
 "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\reads_counts\HESpecNonSyn160119.csv"
 "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\NonSynSitesA2gWONoisySamplesNoAntisense170116.csv"
 "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\BindStrengthPerAllelePerGroupSum.BindingAtLeastStrong.csv.Sense.csv"
 "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\BindStrengthPerAllelePerGroupSum.BindingAtLeastIntermediate.csv.Sense.csv"
 "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\BindStrengthPerAllelePerGroupSum.BindingAtLeastWeak.csv.Sense.csv"
 "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\groups_to_samples.tab" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\HLAMinerAnslysisSummery.csv"
  "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119"



D:\Hillel\Dropbox\Genomics\Scripts\GGPS\ProcessingPlugins\EpitopeAnalysis\pep_aff_exp_script.py "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\reads_counts\HESummery.csv" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\reads_counts\filtered_csv_reads_coun.AGOnly.NonSyn.csv" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\reads_counts\HESpecNonAntisense160119.csv" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\AntisenseSites160119.csv" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\BindStrengthPerAllelePerGroupSum.BindingAtLeastStrong.csv.Antisense.csv" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\BindStrengthPerAllelePerGroupSum.BindingAtLeastIntermediate.csv.Antisense.csv" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\BindStrengthPerAllelePerGroupSum.BindingAtLeastWeak.csv.Antisense.csv" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\groups_to_samples.tab" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Inputs\HLAMinerAnslysisSummery.csv" "D:\Hillel\Dropbox\Genomics\AutoImmun\Lupus\PeptidesAnalysis\170119\Antisense"
"""

'''
Parse R stats output:

import re

path  = "/home/alu/hillelr/AutoImmun/Lupus/Proteins/PeptidesAnalysis/170301_MakeSure_170511_noC0702/R/170511_EFPMS_Raw_statistics_R.txt"
data = open(path).read()
props = re.findall(r'("([\w\-\.]*?)".*?)#########################################################', data, re.DOTALL)

recs = [["Threshold", "CoreAffDelta","Test", "Group1", "Group2", 'P-Value'],]
for p in props:
    for res in re.findall(r"\t(\w+?) .*?\n\ndata:  (\w+?) and (\w+?)\n.*?p-value = (.+?)\n", p[0]):
        rec = p[1].split("-")
        rec.extend(res)
        recs.append(rec)


out_path  = "/home/alu/hillelr/AutoImmun/Lupus/Proteins/PeptidesAnalysis/170301_MakeSure_170511_noC0702/R/EFPMS_Raw_statistics_R_Summery.csv"

with open(out_path, 'wb') as o:
  for r in recs:
    o.write(",".join(r) + "\n")

path  = "/home/alu/hillelr/AutoImmun/Lupus/Proteins/PeptidesAnalysis/170301_MakeSure_170511_noC0702/R/170511_EFPMS_stats_summery_R.txt"
data = open(path).read()
props = re.findall(r'("([\w\-\.]*?)".*?)#########################################################', data, re.DOTALL)

recs = [["Att", "Threshold", "CoreAffDelta", "Filtert", "Test", "Group1", "Group2", 'P-Value'],]
for p in props:
    for res in re.findall(r"\t(\w+?) .*?\n\ndata:  (\w+?) and (\w+?)\n.*?p-value = (.+?)\n", p[0]):
        rec = p[1].split("-")
        rec.extend(res)
        recs.append(rec)


out_path  = "/home/alu/hillelr/AutoImmun/Lupus/Proteins/PeptidesAnalysis/170301_MakeSure_170511_noC0702/R/EFPMS_stats_summery_R_Summery.csv"

with open(out_path, 'wb') as o:
  for r in recs:
    o.write(",".join(r) + "\n")
'''
