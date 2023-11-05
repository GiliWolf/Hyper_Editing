__author__ = 'Hillel'
# =====================imports=====================#

import optparse
import os.path

import numpy
from scipy import stats

from Outputs.Outputers.CSVOutputer import CSVOutputer, EMPTY_VAL

# =====================constants===================#
USAGE = "<script.py> -i <input file> -o <output path> [-g]"

MANN_WHITNEY_U_INDEPENDANT = "Mann-WhitneyUTest-IndependatSamples"
TTEST_IND = "TTest-IndependanSamples"

SIGNIFICANCE = 0.055

PVAL_HEADER = "P-Value For"
REC_COUNT_HEADER_FORMAT = "Num of Records of %s"
REC_MEAN_R_HEADER_FORMAT = "Means Ratio for %s"
COMPARE_TO_HEADER = "Comparing To"
SIG_IN_GROUP_FORMAT = "Sig. In %s"


# =====================functions===================#

def get_arrays(path, group_col, exclude_cols=[], selected_cols=[]):
    with open(path) as input_f:
        res = {}
        array = numpy.genfromtxt(input_f, names=True, delimiter=",", dtype=None)
        run_on_cols = selected_cols if selected_cols else array.dtype.names
        for vals_col in array.dtype.names:
            if vals_col in exclude_cols or vals_col == group_col:
                continue
            groups = {}
            for group in set(array[group_col]):
                groups[group] = [c[0] for c in zip(array[vals_col].tolist(),array[group_col]) if c[1] == group]
            res[vals_col] = groups
    return res


def divied_to_groups(indexes, row):
    res = {}
    for i,val in enumerate(row):
        had_val = False
        if indexes[i] != EMPTY_VAL:
            group = res.setdefault(indexes[i],[])
            if val != "NA":
                had_val = True
                group.append(float(val))

    return res if had_val else None


def get_arrays_from_rna_qunt(path, exclude_cols=[], selected_cols=[]):
    with open(path) as input_f:
        res = {}
        indexes = {}
        array = numpy.genfromtxt(input_f, names=True, delimiter=",", dtype=None)
        for i, val in enumerate(array[-1]):
            indexes[i] = val
        for i,row in enumerate(array[:-1]):
            l_row = row.tolist()
            group = divied_to_groups(indexes, row)
            if group:
                res["-".join(l_row[:3])] = group

        return res



def ttest(a, b):
    try:
        var_equality = stats.levene(a, b).pvalue > SIGNIFICANCE

        return stats.ttest_ind(a,b,equal_var=var_equality).pvalue
    except:
        return 0.0

def mann_whitney_wilcoxon_u_test(a,b):
    return stats.mannwhitneyu(a,b).pvalue
        
TESTS_DICT = {
            MANN_WHITNEY_U_INDEPENDANT:mann_whitney_wilcoxon_u_test,
            TTEST_IND: ttest
            }


def calc_test(groups_arrs, compare_to_group, tests, exclude_groups=[], selected_groups=[]):
    all_tests = {}
    run_on_groups = selected_groups if selected_groups else groups_arrs
    for test in tests:
        ttests = {}
        for val_groups in run_on_groups:
            ttests_res = {}
            for group in groups_arrs[val_groups]:
                if group == compare_to_group or group in exclude_groups:
                    continue
                try:
                    ttests_res[group] = TESTS_DICT[test](groups_arrs[val_groups].get(compare_to_group, []),
                                    groups_arrs[val_groups].get(group, []))
                except ValueError:
                    ttests_res[group] = "NA"
            ttests[val_groups] = ttests_res
            
        all_tests[test] = ttests

    return all_tests


def get_statitics(ipath, opath, groups_col, compare_to_group, tests, exclude_cols=[], exclude_groups=[],
                  selected_groups=[], selected_cols=[], parse_rna_qunt=False):
    if parse_rna_qunt:
        groups = get_arrays_from_rna_qunt(ipath, exclude_cols, selected_cols)
    else:
        groups = get_arrays(ipath, groups_col, exclude_cols, selected_cols)

    test_results = calc_test(groups, compare_to_group, tests, exclude_groups,selected_groups)
    output = CSVOutputer()

    for test in tests:
        headers = [PVAL_HEADER, COMPARE_TO_HEADER]
        recs = []
        sig_recs = []
        group_count = {}
        for val_col in test_results[test]:
            rec = {COMPARE_TO_HEADER: compare_to_group, PVAL_HEADER: val_col}
            any_sig = False
            for group in sorted(test_results[test][val_col]):
                curr_count = REC_COUNT_HEADER_FORMAT % group
                curr_mean_r = REC_MEAN_R_HEADER_FORMAT % group
                curr_sig_in_group = SIG_IN_GROUP_FORMAT % group
                if group not in headers:
                    headers.insert(2, curr_sig_in_group)
                    headers.append(group)
                    headers.append(curr_count)
                    headers.append(curr_mean_r)
                test_res = test_results[test][val_col][group]
                rec[group] = str(test_res)
                if test_res != "NA":
                    if test_res < SIGNIFICANCE and test_res != 0:
                        group_count[group] = group_count.get(group,0) + 1
                        rec[curr_sig_in_group] = "Yes"
                        any_sig = True
                    try:
                        mean_ratio = numpy.nanmean(groups[val_col][group])/numpy.nanmean(groups[val_col][compare_to_group])
                    except Exception, e:
                        mean_ratio = "NA"
                else:
                    mean_ratio = "NA"
                rec[curr_count] = str(len(groups[val_col][group]))
                rec[curr_mean_r] = str(mean_ratio)
            recs.append(rec)
            if any_sig:
                sig_recs.append(rec)
        recs.append(group_count)
        sig_recs.append(group_count)
        output.output([os.path.join(opath, "Statistics", test+ ".csv")],headers,recs)
        output.output([os.path.join(opath, "Statistics","SignificantsOnly", test+ ".csv")],headers,sig_recs)

def main():
    pass
if __name__ == "__main__":
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-c", dest="config_file", help="The configuration file's path")
    parser.add_option("-d", dest="root_dir", default=None, help="The input directory")
    parser.add_option("-s", dest = "subdirs_prefix", default="", help="A prefix of sub dir to run on")
    parser.add_option("-f", dest = "files_suffix", default="", help="A prefix of sub dir to run on")
    parser.add_option("-l", dest = "log_path", default=None, help="The path where the logs will be written")
    parser.add_option("-k", dest = "kill_flag", default=None, help="The path of the kill flag")
    parser.add_option("-r", dest = "recursion_depth", default=100, help="The depth of of recursion into the folders")
    parser.add_option("-t", dest = "truncs", default="0,2", help="The truncation of the file name "
                                                               "<suffixes to retain>,<chars to drop from filename>")

    options = parser.parse_args()[0]
    if not main(options):
        parser.print_usage()
        parser.print_help()
'''
s.get_statitics(r"C:\Users\User\Documents\dropbox\Dropbox\Genomics\AutoImmun\Lupus\RedIToolsKnown\filtered_csv_width_above_strict_min_repNonAlu.csv",
              r"C:\Users\User\desktop\width\repNonAlu", "RNA",
              "PAXgeneBlood_healthy_control", [s.MANN_WHITNEY_U_INDEPENDANT,s.TTEST_IND], ["sample",]
              ,["GM12878Ro60KO_wildtype-IFNa", "GM12878Ro60KO_wildtype-untreated", "GM12878Ro60KO_Ro60_null-IFNa",
                "GM12878Ro60KO_Ro60_null-untreated", "HealthyPBMC_untreated","HealthyPBMC_IFNa" ],parse_rna_qunt=True)
'''

