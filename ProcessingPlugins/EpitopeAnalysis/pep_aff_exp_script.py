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
from Commons.data_structs import AffinitiesEnum, PredicateAnswersEnum
from Commons.consts import NA
from Commons.general_functions import convert_params_to_bool_dict

from Commons.general_functions import load_groups_and_samples_from_file, get_groups_and_sample_names_dict
from DataObjects.KnownProperties.MismatchesSites import MismatchesSites, SampleMMRate
from DataObjects.BaseClasses.Sample import Sample
from DataObjects.DataRecords.PeptidesAffinities import PeptidesAffinities, PeptideAffinity

from Outputs.mismatches_sites_functions import get_mismatches_sites_quantification_per_site, GENERAL_HEADERS_KEY
# =====================constants=====================#
ANY_CA = "AnyCoreAffDeltaNeg"
ANY_CA_EQ = "AnyCoreAffDeltaNegOrZero"
ALL_CA = "AllCoreAffDeltaNeg"
ALL_CA_EQ = "AllCoreAffDeltaNegOrZero"
AVG_CA = "AvgCoreAffDeltaNeg"
AVG_CA_EQ = "AvgCoreAffDeltaNegOrZero"
ANY_CA_POS = "AnyCoreAffDeltaPos"
ANY_CA_POS_EQ= "AnyCoreAffDeltaPosOrZero"
ALL_CA_POS = "AllCoreAffDeltaPos"
ALL_CA_POS_EQ = "AllCoreAffDeltaPosOrZero"
AVG_CA_POS = "AvgCoreAffDeltaPos"
AVG_CA_ALL = "AvgCoreAffDeltaAll"
ALL_CA_ALL = "AnyCoreAffDeltaAll"
ANY_CA_ALL = "AllCoreAffDeltaAll"
AVG_CA_POS_EQ = "AvgCoreAffDeltaPosOrZero"
# CA_PREDICATES = {
#     ANY_CA: lambda rec: any([float(aff) < 0 for aff in rec.split(SUM_REC_SEP)]),
#     ALL_CA: lambda rec: all([float(aff) < 0 for aff in rec.split(SUM_REC_SEP)]),
#     AVG_CA: lambda rec: average([float(aff) for aff in rec.split(SUM_REC_SEP)]) < 0,
#     ANY_CA_POS: lambda rec: any([float(aff) > 0 for aff in rec.split(SUM_REC_SEP)]),
#     ALL_CA_POS: lambda rec: all([float(aff) > 0 for aff in rec.split(SUM_REC_SEP)]),
#     AVG_CA_POS: lambda rec: average([float(aff) for aff in rec.split(SUM_REC_SEP)]) > 0,
# }

READS_COUNT_SEP = ":"
NUMOF = "NumOf"
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
CA_TYPE = "CoreAffinityDelta"
AFF_THRESHOLD = "AffinityThreshold"
EFPMS_HEADERS = [EFMP_FILTER, AFF_THRESHOLD, CA_TYPE, GROUP, SAMPLE, MEAN, MED, AVG, DEV, VAR, NUMOF, NORM_NUMOF]
EFPMS_RAW_HEADERS = [AFF_THRESHOLD, CA_TYPE, GROUP, SAMPLE, EFPM, CHR, POSITION]
SAMPLE_FILENAME = "PerSamplesEFPMsStats.csv"
GROUP_FILENAME = "PerGroupEFPMsStats.csv"
RAW_FILENAME = "EFPMsRaw.csv"
EFMPsStats = namedtuple("EFMPsStats", "average mean median deviation variance numof norm_numof")
# =====================functions=====================#

def get_deltas(affs_orig, affs_neo,):
    deltas = [aff_orig - aff_noe for aff_orig, aff_noe in zip(affs_orig, affs_neo)]
    deltas_pc = [100 * abs(delta) / float(aff_orig) for delta, aff_orig in zip(deltas, affs_orig)]
    return deltas, deltas_pc


def aff_diff(affs_orig, affs_neo, pc_diff, op, boolean_oprator):
    deltas, deltas_pcs = get_deltas(affs_orig, affs_neo)
    avg_aff_diff_pc = average(deltas_pcs)
    if avg_aff_diff_pc < pc_diff:
        return PredicateAnswersEnum.UNKNOWN
    if boolean_oprator([op(delta, 0) for delta in deltas]):
        return PredicateAnswersEnum.YES
    return PredicateAnswersEnum.NO

def aff_diff_avg(affs_orig, affs_neo, pc_diff, op):
    deltas, deltas_pcs = get_deltas(affs_orig, affs_neo)
    avg_aff_diff_pc = average(deltas_pcs)
    if avg_aff_diff_pc < pc_diff:
        return PredicateAnswersEnum.UNKNOWN
    if op(average(deltas), 0):
        return PredicateAnswersEnum.YES
    return PredicateAnswersEnum.NO


CA_PREDICATES = {
    ANY_CA: lambda affs_orig, affs_neo, pc_diff: aff_diff(affs_orig, affs_neo, pc_diff, operator.lt, any),
    ANY_CA_EQ: lambda affs_orig, affs_neo, pc_diff: aff_diff(affs_orig, affs_neo, pc_diff, operator.le, any),
    ALL_CA: lambda affs_orig, affs_neo, pc_diff: aff_diff(affs_orig, affs_neo, pc_diff, operator.lt, all),
    ALL_CA_EQ: lambda affs_orig, affs_neo, pc_diff: aff_diff(affs_orig, affs_neo, pc_diff, operator.le, all),
    AVG_CA: lambda affs_orig, affs_neo, pc_diff: aff_diff_avg(affs_orig, affs_neo, pc_diff, operator.lt),

    ANY_CA_POS: lambda affs_orig, affs_neo, pc_diff: aff_diff(affs_orig, affs_neo, pc_diff, operator.gt, any),
    ANY_CA_POS_EQ: lambda affs_orig, affs_neo, pc_diff: aff_diff(affs_orig, affs_neo, pc_diff, operator.ge, any),
    ALL_CA_POS: lambda affs_orig, affs_neo, pc_diff: aff_diff(affs_orig, affs_neo, pc_diff, operator.gt, all),
    ALL_CA_POS_EQ: lambda affs_orig, affs_neo, pc_diff: aff_diff(affs_orig, affs_neo, pc_diff, operator.ge, all),
    AVG_CA_POS: lambda affs_orig, affs_neo, pc_diff: aff_diff_avg(affs_orig, affs_neo, pc_diff, operator.gt),

    AVG_CA_ALL: lambda affs_orig, affs_neo, pc_diff: (aff_diff_avg(affs_orig, affs_neo, pc_diff, operator.gt) or
                                                     aff_diff_avg(affs_orig, affs_neo, pc_diff, operator.lt)),
    ANY_CA_ALL: lambda affs_orig, affs_neo, pc_diff: (aff_diff(affs_orig, affs_neo, pc_diff, operator.lt, any) or
                                                     aff_diff(affs_orig, affs_neo, pc_diff, operator.gt, any)),
    ALL_CA_ALL: lambda affs_orig, affs_neo, pc_diff: (aff_diff(affs_orig, affs_neo, pc_diff, operator.lt, all) or
                                                     aff_diff(affs_orig, affs_neo, pc_diff, operator.gt, all)),
    #AVG_CA_POS_EQ: lambda affs_orig, affs_neo, pc_diff: average([aff_diff(aff_orig, aff_neo, pc_diff) for aff_orig, aff_neo in zip(affs_orig, affs_neo)]) >= 0
}

def get_sites(sites_bed):
    sites = {}
    with open(sites_bed) as bed:
        for line in reader(bed):
            chrom = line[0]
            pos = line[1]
            sites[chrom + pos] = True
    return sites


def get_exp_per_sites_per_group(STAR_stats_file, he_known_reads_count_file, he_spec_file, groups_to_samples_file,
                                sites_bed, out_dir):
    load_groups_and_samples_from_file(groups_file=groups_to_samples_file)
    sample_and_groups_dict = get_groups_and_sample_names_dict(by_sample=True)
    # groups_to_samples, samples_to_groups = get_groups(groups_to_samples_file)

    mapped_fragments_count = {}
    he_summery_data = [line for line in reader(open(STAR_stats_file))]
    for line in he_summery_data[1:]:
        sample = os.path.basename(line[0])
        mapped_reads_count = line[2]
        mapped_fragments_count.setdefault(sample, 0)
        mapped_fragments_count[sample] += int(mapped_reads_count)

    sites = get_sites(sites_bed)
    rate_per_s = {}

    spec_samples_i = {}
    he_spec_data = [line for line in reader(open(he_spec_file))]
    for i, header in enumerate(he_spec_data[0]):
        if 'UniqueSites' in header:
            sample = re.findall(r"(SRR.+?)(?:\.|$)", header)[0]
            spec_samples_i[i] = sample_and_groups_dict[sample].values()[0]

    he_mismatches = MismatchesSites(source="HE_for_aff")
    for line in he_spec_data[1:]:
        chrom = line[0]
        pos = int(line[1])
        for i, sample in spec_samples_i.iteritems():
            # rate_per_s.setdefault(chrom, {}).setdefault(pos, {}).setdefault(s, 0)
            if line[i] != '-':
                he_mismatches.add_site("A2G", sample, chrom, pos -1, pos, "+", int(line[i]))
                    # rate_per_s[chrom][pos][sample_name] += int(line[i])



    known_mismatches = MismatchesSites(source="known_for_aff")
    knwon_samples_i = {}
    he_known_data = [line for line in reader(open(he_known_reads_count_file))]
    for i, header in enumerate(he_known_data[0]):
        if '.REDItoolKnown.out.tab' in header:

            try:
                knwon_samples_i[i] = sample_and_groups_dict[header.replace('.REDItoolKnown.out.tab', '')].values()[0]
            except KeyError:
                continue

    for line in he_known_data[1:]:
        chrom = line[0]
        pos = int(line[1])
        if not sites.get(line[0] + line[1], False):
            continue
        for i,s in knwon_samples_i.iteritems():
                if line[i] != '-':
                    try:
                        edited, total, ratio = line[i].split(READS_COUNT_SEP)
                    except:
                        print "ejhre"
                    try:
                        if total == '0':
                            # rate_per_s[chrom][pos][s] = NA
                            continue
                        else:
                            known_mismatches.add_site("A2G", s, chrom, pos-1, pos, "+", int(edited), int(total) - int(edited), float(ratio))
                            # if s in rate_per_s[chrom][pos]:
                            #     rate_per_s[chrom][pos][s] = (rate_per_s[chrom][pos][s]+float(edited))/total#+= float(edited)
                    except Exception, e:
                        import pdb;pdb.set_trace()
    # recs = []
    # headers = ["Chr", "Pos", ]
    # headers.extend(["Edited-" + s for s in samples_to_groups.keys()])
    # headers.extend(["EFPM-" + s for s in samples_to_groups.keys()])
    #
    #
    # for chrom, positions in rate_per_s.iteritems():
    #     for pos, count in positions.iteritems():
    #         rec = {"Chr": chrom, "Pos": pos}
    #         for s, c in count.iteritems():
    #             rec["Edited-" + s] = str(c)
    #             if c != NA:
    #                 rec["EFPM-" + s] = str(c)#str(1000000 * float(c) / mapped_fragments_count[s])
    #             else:
    #                 rec["EFPM-" + s] = NA
    #         recs.append(rec)

    combined_mismatches = MismatchesSites("combined_for_aff")
    combined_mismatches.merge_with(known_mismatches.filter_mismatches_by_sites(known_mismatches.sites_to_samples.keys()))
    MismatchesSites.remove_entity(known_mismatches.record_id)
    combined_mismatches.merge_with(he_mismatches.filter_mismatches_by_sites(he_mismatches.sites_to_samples.keys()))
    MismatchesSites.remove_entity(he_mismatches.record_id)

    '''records, headers_d = get_mismatches_sites_quantification_per_site("combined_for_aff", 1)

    headers = []
    headers.extend(headers_d.pop(GENERAL_HEADERS_KEY))
    for key, val in headers_d.iteritems():
        headers.extend(val)


    outper = CSVOutputer()
    out_path = os.path.join(out_dir, "SamplesEditingRates.csv")
    outper.output([out_path], headers, records)'''

    # recs = []
    # rate_per_g = {}
    # headers = ["Chr", "Pos", ]
    # headers.extend(["Edited-" + g for g in groups_to_samples.keys()])
    # headers.extend(["EFPM-" + g for g in groups_to_samples.keys()])
    # for chrom, positions in rate_per_s.iteritems():
    #     for pos, count in positions.iteritems():
    #         for g in groups_to_samples:
    #             rate_per_g.setdefault(chrom, {}).setdefault(pos, {}).setdefault(g, 0.0)
    #         for s, c in count.iteritems():
    #             if c != NA:
    #                 rate_per_g[chrom][pos][samples_to_groups[s]] += c
    # for chrom, positions in rate_per_g.iteritems():
    #     for pos, count in positions.iteritems():
    #         rec = {"Chr": chrom, "Pos": pos}
    #         for g, c in count.iteritems():
    #             rec["Edited-" + g] = str(c)
    #             rec["EFPM-" + g] = str(1000000 * float(c) / mapped_fragments_count[s])
    #         recs.append(rec)
    #
    # outper = CSVOutputer()
    # outper.output([os.path.join(out_dir, "GroupsEditingRates.csv")], headers, recs)
    return mapped_fragments_count


def get_groups(groups_to_samples_file):
    groups_data = [line for line in reader(open(groups_to_samples_file), delimiter="\t")]
    groups_to_samples = {}
    samples_to_groups = {}
    for line in groups_data[1:]:
        if [] == line:
            continue
        groups_to_samples.setdefault(line[0], []).append(line[1])
        samples_to_groups[line[1]] = line[0]
    return groups_to_samples, samples_to_groups


def get_fgh(groups_to_samples, hla_to_samples, merge_sub_strs_dict, merge_ignores ):
    fgh = {}
    for hla, groups in get_hla_allels_per_sample(hla_to_samples, merge_sub_strs_dict, merge_ignores).iteritems():
        for group, samples in groups.iteritems():
            group_size = float(len(groups_to_samples[group]))
            s_num = float(len(samples))
            fgh.setdefault(hla, {})[group] = s_num / group_size

    return fgh


def get_alleles_per_site(orig_affs, neo_affs, limit_to_sites_dict=None, pc_cutoff=10):
    hlas_pos = {}
    # all_binding_d_i = affinity_assessment_summed_data[0].index(BIND_DIFF_HEADER)
    # # all_binding_d_pc_i = affinity_assessment_summed_data[0].index(BIND_DIFF_PC_HEADER)
    # chrom_i = affinity_assessment_summed_data[0].index(CHR)
    # pos_i = affinity_assessment_summed_data[0].index(POSITION)
    # allele_i = affinity_assessment_summed_data[0].index(HLA_HEADER)
    #
    # # counter
    # num_of_n_strong_i = affinity_assessment_summed_data[0].index(SUM_NEO_BIND_HEADER_S)
    # num_of_n_inter_i = affinity_assessment_summed_data[0].index(SUM_NEO_BIND_HEADER_I)
    # num_of_n_weak_i = affinity_assessment_summed_data[0].index(SUM_NEO_BIND_HEADER_W)
    #
    # num_of_o_strong_i = affinity_assessment_summed_data[0].index(SUM_ORIG_BIND_HEADER_S)
    # num_of_o_inter_i = affinity_assessment_summed_data[0].index(SUM_ORIG_BIND_HEADER_I)
    # num_of_o_weak_i = affinity_assessment_summed_data[0].index(SUM_ORIG_BIND_HEADER_W)
    for site, mhc_dict in orig_affs.iteritems():
        if site not in neo_affs:
            continue

        if limit_to_sites_dict and not limit_to_sites_dict.get(site.region + str(site.end), False):
            continue

        for mhc_allele, ranks_dict in mhc_dict.iteritems():
            orig_peptides_d = dict()
            neo_peptides_d = dict()
            neo_aff_pos = neo_affs[site][mhc_allele]

            for rank, peps in ranks_dict.iteritems():
                for peptide in peps:
                    """:type peptide PeptideAffinity"""
                    orig_peptides_d[peptide.offset] = peptide.affinity
            for rank, peps in neo_aff_pos.iteritems():
                for peptide in peps:
                    neo_peptides_d[peptide.offset] = peptide.affinity

            sorted_orig_affs = [orig_peptides_d[offset] for offset in sorted(orig_peptides_d)]
            sorted_neo_affs = [neo_peptides_d[offset] for offset in sorted(neo_peptides_d)]
            num_of_strong = len(ranks_dict.get(AffinitiesEnum.STRONG, [])) + len(neo_aff_pos.get(AffinitiesEnum.STRONG, []))
            num_of_inter = len(ranks_dict.get(AffinitiesEnum.INTERMEDIATE, [])) + len(neo_aff_pos.get(AffinitiesEnum.INTERMEDIATE, []))
            num_of_weak = len(ranks_dict.get(AffinitiesEnum.WEAK, [])) + len(neo_aff_pos.get(AffinitiesEnum.WEAK, []))

            for core_affinity_type, prediacte in CA_PREDICATES.iteritems():
                if prediacte(sorted_orig_affs, sorted_neo_affs, pc_cutoff) == PredicateAnswersEnum.YES:
                    if core_affinity_type == AVG_CA:
                        import pdb;pdb.set_trace()
                    #print core_affinity_type, "succeeded"
                    res_dict = {SUM_NEO_BIND_HEADER_S: num_of_strong ,
                                SUM_NEO_BIND_HEADER_I: num_of_inter,
                                SUM_NEO_BIND_HEADER_W: num_of_weak}
                    allele = mhc_allele#transform_allele(mhc_allele.allele)
                    hlas_pos.setdefault(site, {}).setdefault(core_affinity_type, {})[allele] = res_dict
                else:
                    pass#print core_affinity_type, "Failed"


    '''prev version:
    avg_binding_d_i = affinity_assessment_summed_data[0].index(BIND_DIFF_AVG_HEADER)
    all_binding_d_i = affinity_assessment_summed_data[0].index(BIND_DIFF_HEADER)
    chrom_i = affinity_assessment_summed_data[0].index(CHR)
    pos_i = affinity_assessment_summed_data[0].index(POSITION)
    allele_i = affinity_assessment_summed_data[0].index(HLA_HEADER)

    for line in affinity_assessment_summed_data[1:]:
        keep = float(line[avg_binding_d_i]) <= 0 if use_average \
            else any([float(aff) <= 0 for aff in line[all_binding_d_i].split(SUM_REC_SEP)])
        if not keep:
            continue

        hlas_pos.setdefault(line[chrom_i], {}).setdefault(line[pos_i], {})[line[allele_i]] = line[
            avg_binding_d_i] if use_average \
            else line[all_binding_d_i]

    '''

    return hlas_pos


def calc_avg_efpm_h_per_s(hlas_per_pos, groups_to_samples, fgh, efpm_data, samples_i, chrom_i, pos_i):
    avg_efpm_h_per_s = {}

    for line in efpm_data[1:]:
        chrom = line[chrom_i]
        pos = line[pos_i]

        if not (chrom in hlas_per_pos and pos in hlas_per_pos[chrom]):
            continue

        avg_efpm_h_per_s.setdefault(chrom, {}).setdefault(pos, {})
        for allele in hlas_per_pos[chrom][pos]:
            num_of_alleles = len(hlas_per_pos[chrom][pos])
            for group, samples in groups_to_samples.iteritems():
                avg_efpm_h_per_s[chrom][pos].setdefault(group, {})
                if fgh is None:
                    fgh_g = 1
                else:
                    fgh_g = fgh[allele].get(group, 0)
                for sample in samples:
                    avg_efpm_h_per_s[chrom][pos][group].setdefault(sample, 0)
                    sample_i = samples_i[sample]
                    efpm = float(line[sample_i])
                    avg_efpm_h_per_s[chrom][pos][group][sample] += fgh_g * efpm * num_of_alleles

    return avg_efpm_h_per_s


def calc_avg_efpm_h_per_group(avg_efpm_h_per_s):
    avg_efpm_h_per_group = {}

    for chrom, positions in avg_efpm_h_per_s.iteritems():
        avg_efpm_h_per_group.setdefault(chrom, {})
        for pos, groups in positions.iteritems():
            avg_efpm_h_per_group[chrom].setdefault(pos, {})
            for group, samples_efpms in groups.iteritems():
                group_size = len(samples_efpms)
                avg_efpm_h_per_group[chrom][pos][group] = sum(samples_efpms.values()) / group_size

    return avg_efpm_h_per_group


def get_sampels_i(headers_line, samples):
    samples_i = {}
    for sample in samples:
        samples_i[sample] = headers_line.index("EFPM-%s" % sample)

    return samples_i


def create_merged_groups(merge_sub_strs_dict, merge_ignores, groups_to_samples):
    merged_groups = {}
    for group, samples in groups_to_samples.iteritems():
        if group in merge_ignores:
            continue
        for sub_str in merge_sub_strs_dict:
            if sub_str in group:
                merged_group = merge_sub_strs_dict[sub_str]
                merged_groups.setdefault(merged_group, {sample: samples[sample] for sample in samples})
    groups_to_samples.update(merged_groups)

def get_rank_specification(num_of_b_s, num_of_b_i, num_of_b_w):
    res = ["All", ]
    if num_of_b_s > num_of_b_i:
        if num_of_b_i >= num_of_b_w:
            res.append("DomStrong")
        else:
            if num_of_b_s >= num_of_b_w:
                res.append("DomStrong")
            else:
                res.append("DomWeak")
    else:
        if num_of_b_i >= num_of_b_w:
            res.append("DomIntermediate")
        else:
            res.append("DomWeak")
    if num_of_b_s:
        res.append("AnyStrong")
    if num_of_b_i:
        res.append("AnyIntermediate")
    return res

def create_efpms_groups(orig_affs, neo_affs, groups_to_samples, efpm_data, hla_per_samples, hlas_per_pos):
    # type: (dict, dict, dict, MismatchesSites, dict, dict) -> dict
    efpms = {}

    sites_to_samples = efpm_data.sites_to_samples
    mismatches_dict = efpm_data.get_sites()
    for site, efpm_samples in sites_to_samples.iteritems():
        if site not in orig_affs and site not in neo_affs:  # position didn't have any alleles relevant.
            continue
        if site not in hlas_per_pos:
            continue

        for ca_type in hlas_per_pos[site]:
            for allele in hlas_per_pos[site][ca_type]:
                for group, samples in groups_to_samples.iteritems():
                    for sample_name, sample in samples.iteritems():
                        if sample_name not in hla_per_samples.get(allele, {}).get(group,{}):
                            continue

                        sample_id = sample.record_id
                        try:
                            efpm = mismatches_dict[site.region][site.start][site.end][site.strand]["A2G"][sample_id]#float(line[sample_i]) if line[sample_i] != NA else NA
                            """:type efpm SampleMMRate"""
                            rate = efpm.ratio
                        except KeyError:
                            rate = NA
                        num_of_b_s = hlas_per_pos[site][ca_type][allele][SUM_NEO_BIND_HEADER_S]
                        num_of_b_i = hlas_per_pos[site][ca_type][allele][SUM_NEO_BIND_HEADER_I]
                        num_of_b_w = hlas_per_pos[site][ca_type][allele][SUM_NEO_BIND_HEADER_W]
                        for b_rank in get_rank_specification(num_of_b_s, num_of_b_i, num_of_b_w):
                            efpms.setdefault(b_rank, {}).setdefault(ca_type, {}).setdefault(group, {}).\
                                setdefault(sample_name, []).append([rate, site.region, site.end])


    return efpms


def compute_stats(all_efpms, mapped_fragments_per_sample):
    stats_groups = {}
    stats_samples = {}

    for threshold, efpms in all_efpms.iteritems():
        for ca_type, groups in efpms.iteritems():
            stats_groups.setdefault(threshold, {}).setdefault(ca_type, {})
            stats_samples.setdefault(threshold, {}).setdefault(ca_type, {})
            for group, samples in groups.iteritems():
                g_efpms = []
                g_nz_efpms = []
                g_mapped_frags = 0
                g_sites_d = {}
                if not samples.keys():
                    continue
                for sample, efpms in samples.iteritems():
                    clean_efpms = filter(lambda x: x[0] != NA, efpms)
                    clean_sites = map(lambda x: (x[1], x[2]), clean_efpms)
                    clean_sites_d = convert_params_to_bool_dict(clean_sites)
                    clean_rates = map(lambda x: x[0], clean_efpms)
                    mapped_fragments = float(mapped_fragments_per_sample[sample])
                    nz_efpms = [e for e in clean_rates if e != 0]
                    g_efpms.extend(clean_rates)
                    g_nz_efpms.extend(nz_efpms)
                    g_mapped_frags += mapped_fragments
                    g_sites_d.update(clean_sites_d)
                    stats_samples[threshold][ca_type][sample] = {ALL_EFPMS: EFMPsStats(str(average(clean_rates)), str(mean(clean_rates)),
                                                                            str(median(clean_rates)),str(std(clean_rates)),
                                                                                                  str(var(clean_rates)),
                                                                                       str(len(clean_sites_d)),
                                                                                       str(len(clean_sites_d)/ mapped_fragments)),
                                                                 NON_ZERO_EFPMS:
                                                                     EFMPsStats(str(average(nz_efpms)), str(mean(nz_efpms)),
                                                                                str(median(nz_efpms)), str(std(nz_efpms)),
                                                                                str(var(nz_efpms)), str(len(nz_efpms)),
                                                                                       str(len(nz_efpms)/ mapped_fragments)),
                                                                 }

                try:
                    stats_groups[threshold][ca_type][group] = {ALL_EFPMS: EFMPsStats(str(average(g_efpms)), str(mean(g_efpms)),
                                                                                     str(median(g_efpms)),str(std(g_efpms)),
                                                                                     str(var(g_efpms)),str(len(g_efpms)),
                                                                                     str(len(g_efpms)/ g_mapped_frags)),
                                                               NON_ZERO_EFPMS: EFMPsStats(str(average(g_nz_efpms)),
                                                                                          str(mean(g_nz_efpms)),
                                                                                          str(median(g_nz_efpms)),
                                                                                          str(std(g_nz_efpms)),
                                                                                          str(var(g_nz_efpms)),
                                                                                          str(len(g_sites_d)),
                                                                                     str(len(g_sites_d)/ g_mapped_frags))}
                except:
                    import pdb;pdb.set_trace()
    return stats_groups, stats_samples



IGNORES = ['high_anti_Ro-ISM_low', 'med_anti_Ro-ISM_low', 'Systemic_lupus_erythematosus', 'control']
S_DICT = {'ISM_low': 'ISM_low', 'ISM_high': 'ISM_high'}


def sum_peps(affinity_assessment_orig_dump, affinity_assessment_neo_dump, mapped_fragments_per_sample, hla_to_samples_file,
             out_dir, merge_sub_strs_dict={}, merge_ignores=[], sites_bed=None, pc_cutoff=10):


    groups_to_samples = get_groups_and_sample_names_dict()

    if sites_bed:
        limit_to_sites = get_sites(sites_bed)
    # groups_to_samples, samples_to_groups = get_groups(groups_to_samples_file)
    create_merged_groups(merge_sub_strs_dict, merge_ignores, groups_to_samples)
    # affinity_assessment_summed_data_s = [line for line in reader(open(affinity_assessment_summed_file_strong))]
    # hlas_per_pos_strong = get_alleles_per_site(affinity_assessment_summed_data_s, limit_to_sites)
    # del affinity_assessment_summed_data_s
    #
    # affinity_assessment_summed_data_i = [line for line in reader(open(affinity_assessment_summed_file_inter))]
    # hlas_per_pos_inter = get_alleles_per_site(affinity_assessment_summed_data_i, limit_to_sites)
    # del affinity_assessment_summed_data_i

    # affinity_assessment_summed_data_a = [line for line in reader(open(affinity_assessment_summed_file_all))]
    # hlas_per_pos_all = get_alleles_per_site(affinity_assessment_summed_data_a, limit_to_sites)
    # del affinity_assessment_summed_data_a
    orig_affs = PeptidesAffinities.load_me(affinity_assessment_orig_dump)
    orig_affs_d = orig_affs.get_peptides_affinities()
    neo_affs = PeptidesAffinities.load_me(affinity_assessment_neo_dump)
    neo_affs_d = neo_affs.get_peptides_affinities()
    efpm_data = MismatchesSites("combined_for_aff")
    # samples_i = get_sampels_i(efpm_data[0], samples_to_groups.keys())
    # chrom_i = efpm_data[0].index(CHR)
    # pos_i = efpm_data[0].index("Pos")

    hla_per_samples = get_hla_allels_per_sample(hla_to_samples_file, merge_sub_strs_dict, merge_ignores)
    hlas_per_pos = get_alleles_per_site(orig_affs_d, neo_affs_d, limit_to_sites)
    efpms = create_efpms_groups(orig_affs_d, neo_affs_d, groups_to_samples, efpm_data, hla_per_samples, hlas_per_pos)
    # efpms["Strong"] = sum_epitopes_nums(hlas_per_pos_strong, groups_to_samples, efpm_data, samples_i, chrom_i, pos_i, hla_per_samples, )
    # efpms["Intermediate"] = sum_epitopes_nums(hlas_per_pos_inter, groups_to_samples, efpm_data, samples_i, chrom_i, pos_i, hla_per_samples)
    # efpms["All"] = sum_epitopes_nums(hlas_per_pos_all, groups_to_samples, efpm_data, samples_i, chrom_i, pos_i, hla_per_samples)

    stats_groups, stats_samples = compute_stats(efpms, mapped_fragments_per_sample)

    create_files(stats_groups, stats_samples, groups_to_samples,efpms, out_dir)

    ##############################

    # hlas_per_pos_avg = get_alleles_per_site(affinity_assessment_summed_data, use_average=True)
    # hlas_per_pos_any = get_alleles_per_site(affinity_assessment_summed_data, use_average=False)
    # fgh = get_fgh(groups_to_samples, hla_to_samples_file, merge_sub_strs_dict, merge_ignores)
    #avg_efpm_h_per_s_avg = calc_avg_efpm_h_per_s(hlas_per_pos_avg, groups_to_samples, fgh,
    #                                              efpm_data, samples_i, chrom_i,
    #                                              pos_i, merge_sub_strs_dict, merge_ignores)
    #
    # avg_efpm_h_per_g_avg = calc_avg_efpm_h_per_group(avg_efpm_h_per_s_avg)
    #
    # create_files(avg_efpm_h_per_g_avg, avg_efpm_h_per_s_avg, hlas_per_pos_avg, out_dir, "AvgEFPMH", "deltaAvg")
    #
    # avg_efpm_per_s_avg = calc_avg_efpm_h_per_s(hlas_per_pos_avg, groups_to_samples, None,
    #                                              efpm_data, samples_i, chrom_i,
    #                                              pos_i, merge_sub_strs_dict, merge_ignores)
    #
    # avg_efpm_per_g_avg = calc_avg_efpm_h_per_group(avg_efpm_h_per_s_avg)
    #
    # create_files(avg_efpm_per_g_avg, avg_efpm_per_s_avg, hlas_per_pos_avg, out_dir, "AvgEFPM", "deltaAvg")
    # ##############################
    # avg_efpm_h_per_s_any = calc_avg_efpm_h_per_s(hlas_per_pos_any, groups_to_samples, fgh,
    #                                              efpm_data, samples_i, chrom_i,
    #                                              pos_i, merge_sub_strs_dict, merge_ignores)
    #
    # avg_efpm_h_per_g_any = calc_avg_efpm_h_per_group(avg_efpm_h_per_s_any)
    #
    # create_files(avg_efpm_h_per_g_any, avg_efpm_h_per_s_any, hlas_per_pos_any, out_dir, "AvgEFPMH", "deltaAny")
    #
    # avg_efpm_per_s_any = calc_avg_efpm_h_per_s(hlas_per_pos_any, groups_to_samples, None,
    #                                              efpm_data, samples_i, chrom_i,
    #                                              pos_i, merge_sub_strs_dict, merge_ignores)
    #
    # avg_efpm_per_g_any = calc_avg_efpm_h_per_group(avg_efpm_h_per_s_any)
    #
    # create_files(avg_efpm_per_g_any, avg_efpm_per_s_any, hlas_per_pos_any, out_dir, "AvgEFPM", "deltaAny")


def get_groups_for_out(groups_to_samples, sample):
    groups = []
    for g, samples in groups_to_samples.iteritems():
        if sample in samples:
            groups.append(g)
    return groups


def create_files(stats_groups, stats_samples, groups_to_samples, raw_efpms, out_dir):#avg_efpm_h_per_g_avg, avg_efpm_h_per_s_avg, hlas_per_pos_avg, out_dir, data_type, filter_type):
    outpter = CSVOutputer()

    efpms_headers = EFPMS_HEADERS[:]
    efpms_headers_g = efpms_headers[:]
    efpms_headers_g.pop(efpms_headers_g.index(SAMPLE))  # pop sample out
    # stats_samples[threshold][ca_type][sample] = EFMPsStats(average(efpms), mean(efpms), median(efpms), std(efpms),
    #                                                        var(efpms))
    #
    # stats_groups[threshold][ca_type][group] = EFMPsStats(average(g_efpms), mean(g_efpms), median(g_efpms), std(g_efpms),
    #          deviation variance                                               var(g_efpms))
    s_recs = []
    for threshold, ca_types_data in stats_samples.iteritems():
        for ca_type, samples_efpms in ca_types_data.iteritems():
            for sample, efpms in samples_efpms.iteritems():
                for group in get_groups_for_out(groups_to_samples, sample):
                    rec = {EFMP_FILTER: NON_ZERO_EFPMS,
                            AFF_THRESHOLD:threshold,
                           CA_TYPE: ca_type,
                           SAMPLE:  sample,
                           GROUP:   group,
                           MEAN:    efpms[NON_ZERO_EFPMS].mean,
                           MED:     efpms[NON_ZERO_EFPMS].median,
                           AVG:     efpms[NON_ZERO_EFPMS].average,
                           DEV:     efpms[NON_ZERO_EFPMS].deviation,
                           VAR:     efpms[NON_ZERO_EFPMS].variance,
                           NUMOF:   efpms[NON_ZERO_EFPMS].numof,
                           NORM_NUMOF:   efpms[NON_ZERO_EFPMS].norm_numof }
                    s_recs.append(rec)
                    rec = {EFMP_FILTER: ALL_EFPMS,
                            AFF_THRESHOLD:  threshold,
                           CA_TYPE: ca_type,
                           SAMPLE:  sample,
                           GROUP:   group,
                           MEAN:    efpms[ALL_EFPMS].mean,
                           MED:     efpms[ALL_EFPMS].median,
                           AVG:     efpms[ALL_EFPMS].average,
                           DEV:     efpms[ALL_EFPMS].deviation,
                           VAR:     efpms[ALL_EFPMS].variance,
                           NUMOF:   efpms[ALL_EFPMS].numof,
                           NORM_NUMOF: efpms[ALL_EFPMS].norm_numof}
                    s_recs.append(rec)
    outpter.output([os.path.join(out_dir, SAMPLE_FILENAME)], efpms_headers, s_recs)
                    # outpter.output([avg_g_out_path], avg_g_headers, avg_g_recs)
    g_recs = []
    for threshold, ca_types_data in stats_groups.iteritems():
        for ca_type, groups_efpms in ca_types_data.iteritems():
            for group, efpms in groups_efpms.iteritems():
                rec = {EFMP_FILTER: NON_ZERO_EFPMS,
                    AFF_THRESHOLD: threshold,
                       CA_TYPE: ca_type,
                       GROUP: group,
                       MEAN: efpms[NON_ZERO_EFPMS].mean,
                       MED: efpms[NON_ZERO_EFPMS].median,
                       AVG: efpms[NON_ZERO_EFPMS].average,
                       DEV: efpms[NON_ZERO_EFPMS].deviation,
                       VAR: efpms[NON_ZERO_EFPMS].variance,
                       NUMOF:   efpms[NON_ZERO_EFPMS].numof,
                       NORM_NUMOF: efpms[NON_ZERO_EFPMS].norm_numof}

                g_recs.append(rec)
                rec = {EFMP_FILTER: ALL_EFPMS,
                    AFF_THRESHOLD: threshold,
                       CA_TYPE: ca_type,
                       GROUP: group,
                       MEAN: efpms[ALL_EFPMS].mean,
                       MED: efpms[ALL_EFPMS].median,
                       AVG: efpms[ALL_EFPMS].average,
                       DEV: efpms[ALL_EFPMS].deviation,
                       VAR: efpms[ALL_EFPMS].variance,
                       NUMOF:   efpms[ALL_EFPMS].numof,
                       NORM_NUMOF: efpms[ALL_EFPMS].norm_numof}

                g_recs.append(rec)
    outpter.output([os.path.join(out_dir, GROUP_FILENAME)], efpms_headers_g, g_recs)
    #efpms[ca_type][group].setdefault(sample, []).append(efpm)
    raw_recs =[]
    for threshold in raw_efpms:
        for ca_type in raw_efpms[threshold]:
            for group in raw_efpms[threshold][ca_type]:
                for sample in raw_efpms[threshold][ca_type][group]:
                    for efpm in raw_efpms[threshold][ca_type][group][sample]:
                        raw_recs.append({AFF_THRESHOLD:threshold, CA_TYPE:ca_type, GROUP:group, SAMPLE:sample,
                                         EFPM:str(efpm[0]), CHR: str(efpm[1]), POSITION:str(efpm[2])})

    outpter.output([os.path.join(out_dir, RAW_FILENAME)], EFPMS_RAW_HEADERS, raw_recs)
                # avg_g_headers = [CHR, POSITION, "Group", data_type, BIND_DIFF_AVG_HEADER, "HLASPSize"]
    # avg_s_headers = [CHR, POSITION, "Group", "Sample", data_type, BIND_DIFF_AVG_HEADER, "HLASPSize"]
    # avg_s_out_path = os.path.join(out_dir, "%s_per_sample.%s.csv" % (data_type, filter_type))
    # avg_g_out_path = os.path.join(out_dir, "%s_per_group.%s.csv"% (data_type, filter_type))
    # avg_s_recs = []
    # avg_g_recs = []
    # for chrom, positions in avg_efpm_h_per_s_avg.iteritems():
    #     for pos, groups in positions.iteritems():
    #         for group, samples_efpms in groups.iteritems():
    #             avg_g_recs.append({CHR: chrom, POSITION: pos, "Group": group,
    #                                data_type: str(avg_efpm_h_per_g_avg[chrom][pos][group]),
    #                                BIND_DIFF_AVG_HEADER: SUM_REC_SEP.join(
    #                                    [allele + delta_aff for allele, delta_aff in
    #                                     hlas_per_pos_avg[chrom][pos].iteritems()]),
    #                                "HLASPSize": str(len(hlas_per_pos_avg[chrom][pos]))})
    #             for sample in samples_efpms:
    #                 avg_s_recs.append({CHR: chrom, POSITION: pos, "Group": group, "Sample": sample,
    #                                    data_type: str(avg_efpm_h_per_s_avg[chrom][pos][group][sample]),
    #                                    BIND_DIFF_AVG_HEADER: SUM_REC_SEP.join(
    #                                        [allele + delta_aff for allele, delta_aff in
    #                                         hlas_per_pos_avg[chrom][pos].iteritems()]),
    #                                    "HLASPSize": str(len(hlas_per_pos_avg[chrom][pos]))})
    # outpter = CSVOutputer()
    # outpter.output([avg_s_out_path], avg_s_headers, avg_s_recs)
    # outpter.output([avg_g_out_path], avg_g_headers, avg_g_recs)
if __name__ == "__main__":
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))



    # he_summery_file, he_known_reads_count_file, he_spec_file,sites_bed, affinity_assessment_summed_file_strong,\
    # affinity_assessment_summed_file_inter,\
    # affinity_assessment_summed_file_all, groups_to_samples_file, \
    # hla_to_samples_file, out_dir,  = sys.argv[1:]0
    he_summery_file, he_known_reads_count_file, he_spec_file, sites_bed, affinity_assessment_orig_dump, \
    affinity_assessment_neo_dump, groups_to_samples_file, hla_to_samples_file, out_dir, pc_cutoff =  sys.argv[1:]
    # affinity_assessment_summed_file_all, groups_to_samples_file, \
    # hla_to_samples_file, out_dir,  = sys.argv[1:]
    mapped_fragments_per_sample = get_exp_per_sites_per_group(he_summery_file, he_known_reads_count_file, he_spec_file, groups_to_samples_file,
                                sites_bed, out_dir)
    # efpm_data_file = os.path.join(out_dir, "SamplesEditingRates.csv")
    sum_peps(affinity_assessment_orig_dump, affinity_assessment_neo_dump, mapped_fragments_per_sample, hla_to_samples_file,
             out_dir, S_DICT, IGNORES, sites_bed, int(pc_cutoff))


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
