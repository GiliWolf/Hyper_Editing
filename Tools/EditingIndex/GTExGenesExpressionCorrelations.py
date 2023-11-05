# =====================imports=====================#
# region Builtin Import
import argparse
import logging
import os
import sys
from collections import OrderedDict
import multiprocessing

import datetime
from scipy.stats import spearmanr, pearsonr, rankdata, kendalltau
from scipy import nan
from numpy import std, average, nanmean, nanstd

# endregion

# region InternalImports
if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from Commons.general_functions import get_file_handle, init_logging_dict
from Commons.consts import GROUP, SAMPLE, MismatchesAndRefsEnum, MISMATCH_TYPE
from Outputs.Outputers.CSVOutputer import CSVOutputer
from Tools.EditingIndex.EditingIndexConsts import EDITING_INDEX_FORMAT, STRAND_DECIDING_METHOD
from Tools.EditingIndex.A2IEditingIndex import STRAND_DECIDING_BY_REFSEQ_AND_UNFILTERED_MM_SITES as BY_REFSEQ

# endregion

# =====================constants===================#
MISSING_SAMPLES_DICT = dict()

TABLE_SEP = "\t"

ENSEMBLE_ID_HEADER = "Name"
COMMON_NAME_HEADER = "Description"
COMBINED_GROUP = "AllCombined"
COMBINED_FROM_SAMPLES = "FromSamples"

NAN_STR = str(nan)


# =====================functions===================#
def create_files_for_r_analyses(index_summery_path, gene_expression_file, output_path, groups_header, max_processes,
                                name_trans_file):
    timestamp = datetime.datetime.today().isoformat()
    init_logging_dict(os.path.join(os.path.dirname(output_path), ".".join(["GTExCorrPreProcessing", timestamp, "log"])))
    trans_dict = dict()
    with open(name_trans_file) as ntf:
        indexes = {header: i for i, header in enumerate(ntf.readline().strip().split("\t"))}
        for line in ntf:
            recs = line.strip().split("\t")
            sample = recs[indexes[SAMPLE]]
            combined_from = [s for s in recs[indexes[COMBINED_FROM_SAMPLES]].split(";") if s]
            trans_dict[sample] = combined_from

    csv_outputer = CSVOutputer(append=True, override=False)
    samples_data = dict()

    with open(index_summery_path) as isp:
        indexes = {header: i for i, header in enumerate(isp.readline().strip().split(","))}
        for line in isp:
            recs = line.strip().split(",")
            if not recs[indexes[STRAND_DECIDING_METHOD]] == BY_REFSEQ:
                continue
            index_sample_data = dict()
            sample = recs[indexes[SAMPLE]]

            index_sample_data[SAMPLE] = sample
            try:
                index_sample_data[GROUP] = recs[indexes[groups_header]]
            except KeyError:
                print "Groups Header Provided is Wrong!"
            index_sample_data[EDITING_INDEX_FORMAT % MismatchesAndRefsEnum.A2G] = \
                recs[indexes[EDITING_INDEX_FORMAT % MismatchesAndRefsEnum.A2G]]
            if sample in trans_dict:
                samples_data[sample] = index_sample_data

    print "Loaded joined samples indexes: %s" % str(len(samples_data))

    samples_to_use = trans_dict.keys()
    gene_exp_fh = get_file_handle(gene_expression_file)
    ensemble_id_i = 0
    common_name_i = 1
    indexes = dict()
    for line in gene_exp_fh:
        records = line.strip("\n").split(TABLE_SEP)
        if not (ENSEMBLE_ID_HEADER in line and COMMON_NAME_HEADER):
            continue
        else:
            ensemble_id_i = records.index(ENSEMBLE_ID_HEADER)
            common_name_i = records.index(COMMON_NAME_HEADER)
            indexes = {rec: index for index, rec in enumerate(records)}
            _ = indexes.pop(ENSEMBLE_ID_HEADER)
            _ = indexes.pop(COMMON_NAME_HEADER)
            break
    sema = multiprocessing.Semaphore(max_processes)
    headers = [SAMPLE, GROUP, EDITING_INDEX_FORMAT % MismatchesAndRefsEnum.A2G, "GeneName", "EnsemblID",
               "GeneExpression"]
    ps = []
    output_recs = list()
    c_out_c = 0
    for line in gene_exp_fh:
        if c_out_c == 100:
            csv_outputer.output([output_path, ], headers, output_recs)
            output_recs = list()
            c_out_c = 0
        process_gene_line(common_name_i, csv_outputer, ensemble_id_i, headers, indexes, line, output_path, samples_data,
                          samples_to_use, sema, trans_dict, output_recs)
        c_out_c += 1

    #     thr = multiprocessing.Process(target=process_gene_line,
    #                                    args=(common_name_i, csv_outputer, ensemble_id_i, headers, indexes, line, output_path, samples_data,
    #                       samples_to_use, sema, trans_dict))
    #     thr.daemon = True
    #     ps.append(thr)
    #
    # for thr in ps:
    #     sema.acquire()
    #     thr.run()
    #
    # for thr in multiprocessing.active_children():
    #     thr.join()



def process_gene_line(common_name_i, csv_outputer, ensemble_id_i, headers, indexes, line, output_path, samples_data,
                      samples_to_use, sema, trans_dict, output_recs):
    index_counter = 0
    samples_counter = 0

    try:
        records = line.strip("\n").split(TABLE_SEP)
        logging.debug("Num of records: %s" % str(len(records)))
        ensembl_id = records[ensemble_id_i]
        common_name = records[common_name_i]
        # output_recs = list()
        for sample, merged_from in trans_dict.iteritems():
            if sample not in samples_to_use:
                continue
            output_rec = dict()
            output_rec.update(samples_data[sample])
            output_rec["GeneName"] = common_name
            output_rec["EnsemblID"] = ensembl_id
            samples_counter += len(merged_from)
            index_counter += 1
            exps = list()
            for s in trans_dict[sample]:
                try:
                    exps.append(float(records[indexes[s]]))
                except KeyError:
                    if s not in MISSING_SAMPLES_DICT:
                        MISSING_SAMPLES_DICT[s] = True
                        logging.debug("Warning: Sample %s Wasn't Found In Genes Expression Data!" % s)
                    continue  # if there's no gene expression data don't use this sample index values
            output_rec["GeneExpression"] = str(average(exps))
            output_recs.append(output_rec)
        # csv_outputer.output([output_path, ], headers, output_recs)
    except Exception,e:
        logging.exception("Failed in the middle")
    finally:
        logging.debug("%s, %s samples, %s merged" % (common_name, samples_counter, index_counter))
        sema.release()


def get_correlations(index_summery_path, gene_expression_file, output_path, groups_header, max_processes):
    csv_outputer = CSVOutputer(append=True, override=False)
    index_data = dict()
    index_sample_groups = {COMBINED_GROUP: list()}

    with open(index_summery_path) as isp:
        indexes = {header: i for i, header in enumerate(isp.readline().strip().split(","))}
        for line in isp:
            recs = line.strip().split(",")
            if not recs[indexes[STRAND_DECIDING_METHOD]] == BY_REFSEQ:
                continue
            sample = recs[indexes[SAMPLE]]
            try:
                index_sample_groups.setdefault(recs[indexes[groups_header]], list()).append(sample)
                index_sample_groups[COMBINED_GROUP].append(sample)
            except KeyError:
                print "Groups Header Provided is Wrong!"
            index_data[sample] = dict()

            for mm_type in MismatchesAndRefsEnum.UNSTRANDED_MISMATCHES:
                index_data.setdefault(mm_type, dict())
                index_data[mm_type][sample] = recs[indexes[EDITING_INDEX_FORMAT % mm_type]]

    samples_to_use = index_sample_groups[COMBINED_GROUP]
    gene_exp_fh = get_file_handle(gene_expression_file)
    ensemble_id_i = 0
    common_name_i = 1
    indexes = dict()
    for line in gene_exp_fh:
        records = line.strip("\n").split(TABLE_SEP)
        if not (ENSEMBLE_ID_HEADER in line and COMMON_NAME_HEADER):
            continue
        else:
            ensemble_id_i = records.index(ENSEMBLE_ID_HEADER)
            common_name_i = records.index(COMMON_NAME_HEADER)
            indexes = {rec: index for index, rec in enumerate(records)}
            _ = indexes.pop(ENSEMBLE_ID_HEADER)
            _ = indexes.pop(COMMON_NAME_HEADER)
            break
    ps = list()
    sema = multiprocessing.Semaphore(max_processes)
    for line in gene_exp_fh:
        records = line.strip("\n").split(TABLE_SEP)
        ensembl_id = records[ensemble_id_i]
        common_name = records[common_name_i]
        genes_expression = {sample: float(records[sample_i]) for sample, sample_i in indexes.iteritems() if
                            sample in samples_to_use}
        tr = multiprocessing.Process(target=calc_correlations,
                                     args=(csv_outputer, genes_expression, index_data, index_sample_groups,
                                           output_path, sema, ensembl_id, common_name),
                                     name="EditingIndexCorrSubprocess")
        tr.daemon = True
        ps.append(tr)

    for thr in ps:
        sema.acquire()
        thr.run()
    for thr in multiprocessing.active_children():
        thr.join()


def calc_correlations(csv_outputer, genes_expression, index_data, index_sample_groups, output_path, sema,
                      ensembl_id, common_name):
    for mm_type in MismatchesAndRefsEnum.UNSTRANDED_MISMATCHES:
        output_path_formatted = output_path % mm_type
        tr = multiprocessing.Process(target=corr_per_mm,
                                     args=(csv_outputer, genes_expression, index_data, index_sample_groups,
                                           mm_type, output_path_formatted, ensembl_id, common_name))
        tr.daemon = True
        tr.start()
        for thr in multiprocessing.active_children():
            thr.join()

    sema.release()


def corr_per_mm(csv_outputer, genes_expression, index_data, index_sample_groups, mm_type, output_path_formatted,
                ensembl_id, common_name):
    results = list()
    ps = list()
    sema = multiprocessing.Semaphore(10)
    for group in index_sample_groups:
        tr = multiprocessing.Process(target=corr_per_group,
                                     args=(genes_expression, group, index_data, index_sample_groups, mm_type,
                                           results, ensembl_id, common_name, sema))
        tr.daemon = True
        ps.append(tr)

    for thr in ps:
        sema.acquire()
        thr.run()
    for thr in multiprocessing.active_children():
        thr.join()

    csv_outputer.output([output_path_formatted, ], results[0].keys(), results)


def corr_per_group(genes_expression, group, index_data, index_sample_groups, mm_type, results, ensembl_id,
                   common_name, sema):
    global MISSING_SAMPLES_DICT
    record = OrderedDict()
    record[GROUP] = group
    record[MISMATCH_TYPE] = mm_type
    record["GeneName"] = common_name
    record["EnsemblID"] = ensembl_id
    samples = sorted(index_sample_groups[group])
    index_vals = list()
    gene_exp = list()
    for sample in samples:
        try:
            sample_ge = genes_expression[sample]
            if str(sample_ge) == NAN_STR:
                continue
            else:
                gene_exp.append(float(sample_ge))
        except KeyError:
            if sample not in MISSING_SAMPLES_DICT:
                MISSING_SAMPLES_DICT[sample] = True
                print "Warning: Sample %s Wasn't Found In Genes Expression Data!" % sample
            continue  # if there's no gene expression data don't use this sample index values

        index_vals.append(float(index_data[mm_type][sample]))
    ge_avg = mean(gene_exp)
    index_avg = average(index_vals)
    record["AvgExpression"] = str(ge_avg)
    record["ExpressionSTD"] = str(std(gene_exp))
    record["IndexAvg"] = str(index_avg)
    record["IndexSTD"] = str(std(index_vals))
    if ge_avg < 1:
        spearman_rho = spearman_pval = pearson_rho = pearson_pval = kendall_tau = kendall_pval = nan
    else:
        index_ranks = rankdata(index_vals, method='ordinal')
        gene_exp_ranks = rankdata(gene_exp, method='ordinal')
        genes_mean = nanmean(gene_exp)
        genes_std = nanstd(gene_exp)
        norm_genes_exp = [(exp - genes_mean) / genes_std for exp in gene_exp]
        try:
            spearman_rho, spearman_pval = spearmanr(index_vals, norm_genes_exp)
            pearson_rho, pearson_pval = pearsonr(index_vals, norm_genes_exp)
            kendall_tau, kendall_pval = kendalltau(index_ranks, gene_exp_ranks)
        except ValueError:
            spearman_rho = nan
            spearman_pval = nan
            pearson_rho = nan
            pearson_pval = nan
            kendall_tau = nan
            kendall_pval = nan

    record["SpearmanRho"] = str(spearman_rho)
    record["SpearmanPVal"] = str(spearman_pval)
    record["PearsonRho"] = str(pearson_rho)
    record["PearsonPVal"] = str(pearson_pval)
    record["KendallTau"] = str(kendall_tau)
    record["KendallPVal"] = str(kendall_pval)
    results.append(record)

    sema.release()


def trans_gtex_tpms(gtex_tpm_path, out_path, name_trans_file):
    trans_dict = dict()
    with open(name_trans_file) as ntf:
        indexes = {header: i for i, header in enumerate(ntf.readline().strip().split("\t"))}
        for line in ntf:
            recs = line.strip().split("\t")
            sample = recs[indexes[SAMPLE]]
            combined_from = [s for s in recs[indexes[COMBINED_FROM_SAMPLES]].split(";") if s]
            trans_dict[sample] = combined_from

    print "Num of joined Samples: %s" % len(trans_dict)
    indexes = dict()
    gene_exp_fh = get_file_handle(gtex_tpm_path)
    for line in gene_exp_fh:
        records = line.strip("\n").split(TABLE_SEP)
        if not (ENSEMBLE_ID_HEADER in line and COMMON_NAME_HEADER):
            continue
        else:
            ensemble_id_i = records.index(ENSEMBLE_ID_HEADER)
            common_name_i = records.index(COMMON_NAME_HEADER)
            indexes = {rec: index for index, rec in enumerate(records)}
            _ = indexes.pop(ENSEMBLE_ID_HEADER)
            _ = indexes.pop(COMMON_NAME_HEADER)
            break

    res = "\t".join([ENSEMBLE_ID_HEADER, COMMON_NAME_HEADER, ] + sorted(trans_dict)) + "\n"
    for line in gene_exp_fh:
        records = line.strip("\n").split(TABLE_SEP)
        res += records[ensemble_id_i] + "\t"
        res += records[common_name_i] + "\t"
        for sample in sorted(trans_dict):
            exps = list()
            for s in trans_dict[sample]:
                try:
                    exps.append(float(records[indexes[s]]))
                except KeyError:
                    if s not in MISSING_SAMPLES_DICT:
                        MISSING_SAMPLES_DICT[s] = True
                        print "Warning: Sample %s Wasn't Found In Genes Expression Data!" % s
                    continue  # if there's no gene expression data don't use this sample index values
            res += str(average(exps)) + "\t"
        res += "\n"

        with open(out_path, 'ab') as out:
            out.write(res)
            res = ""


def translate_samples_rin(rins_file, out_path, name_trans_file):
    trans_dict = dict()
    values_dict = dict()
    with open(name_trans_file) as ntf:
        indexes = {header: i for i, header in enumerate(ntf.readline().strip().split("\t"))}
        for line in ntf:
            recs = line.strip().split("\t")
            sample = recs[indexes[SAMPLE]]
            values_dict[sample] = list()
            for s in recs[indexes[COMBINED_FROM_SAMPLES]].split(";"):
                if s:
                    trans_dict[s] = sample
    samples_c = 0
    sid_c = 0
    with open(rins_file) as rinf:
        indexes = {header: i for i, header in enumerate(rinf.readline().strip().split(","))}
        for line in rinf:
            recs = line.strip().split(",")
            sid = recs[indexes["SAMPID"]]
            rin = float(recs[indexes["SMRIN"]])

            try:
                sample = trans_dict[sid]
            except KeyError:
                print "Sample %s not in indexes!" % sid
                continue
            sid_c += 1
            values_dict[sample].append(rin)

    res_recs = list()

    for sample, rins in values_dict.iteritems():
        res_recs.append({SAMPLE: sample, "RIN": str(nanmean(rins))})
        samples_c += 1

    outputter = CSVOutputer()
    outputter.output([out_path, ], res_recs[0].keys(), res_recs)
    print "Done! %s SIDs, %s Samples" % (sid_c, samples_c)



if __name__ == "__main__":
    script_dir = os.path.dirname(__file__)
    desc = """Creates a correlation tests output between the indexes and genes expression data"""


    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(prog='Editing Index Correlation', description=desc, formatter_class=MyFormatter)

    parser.add_argument('-i', '--indexFile', metavar="A2IEditingIndex with counts output", dest='index_file', nargs='?',
                        required=True, help="The path of the A2IEditingIndex output.")

    parser.add_argument('-g', '--genesExpresion', metavar="genes expression table", dest='gene_expression_file',
                        nargs='?',
                        required=True, help="The path of the genes expressions file (tab file with the headers "
                                            "<'Name' (Ensemble ID)><'Description' (Common Name)><Sample name 1>"
                                            "<Sample name2> ....  .")
    parser.add_argument('-gh', '--groupsHeader', metavar="groups header", dest='groups_header', nargs='?',
                        required=True,
                        help="The header of the groups column, to split the data according to the correct"
                             " groups")

    parser.add_argument('-o', '--output_dir', metavar="output_dir", dest='output_dir', nargs='?', required=False,
                        default=".", help="The root directory for the output")
    parser.add_argument('-n', '--names', metavar="trans_file", dest='trans_file', nargs='?', required=False,
                        default=".", help="")

    parser.add_argument('-t', '--threads', metavar="threads", dest='max_processes', nargs='?', required=False,
                        default=15, type=int, help="The number of threads to use")

    options = parser.parse_args()
    create_files_for_r_analyses(index_summery_path=options.index_file,
                                gene_expression_file=options.gene_expression_file,
                                output_path=os.path.join(options.output_dir, "EditingIndexGenesExpressions.csv"),
                                # "EditingIndexGenesCorrelations%s.csv"),
                                name_trans_file=options.trans_file,
                                groups_header=options.groups_header,
                                max_processes=options.max_processes)
