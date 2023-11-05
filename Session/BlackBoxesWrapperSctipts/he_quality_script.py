import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from Outputs.Outputers.CSVOutputer import CSVOutputer
from Commons.help_functions import reverse_strand
import re
from csv import reader 
SAMPLE_RE = re.compile(r"Analysis\(stranded\) of detected.*?/([^/]*?)-\d\n\n(Edit_type.*?)\n\n\n", re.DOTALL)
EU_MM_FORAMT = "ClustersOfMismatch%s"
ES_MM_FORAMT = r"HESitesOfMismatch%s"
HEADERS = ["Sample", ]
def quality_asessment(statistics_path = "statistic\analyse.SE_0.05_0.8_30_0.6_0.1_0.8_0.2", out_path = "."):
    data = open(statistics_path).read()
    recs = SAMPLE_RE.findall(data)
    outputer = CSVOutputer()
    samples = {}
    added_headers = {}
    orecs= []
    for rec in recs:
        sample = rec[0]
        lines = [l.split("\t") for l in rec[1].split("\n")]
        mm_i = lines[0].index("Edit_type")
        eu_s_i = lines[0].index(r"%UE_of_tot")
        es_s_i = lines[0].index(r"%ES_of_tot")
        for line in lines[1:]:
            if line[mm_i] not in added_headers:
                added_headers[line[mm_i]] = True
                HEADERS.append(EU_MM_FORAMT%line[mm_i])
                HEADERS.append(ES_MM_FORAMT%line[mm_i])
            samples.setdefault(sample, {}).setdefault(line[mm_i], {})[r"%UE_of_tot"] = line[eu_s_i]
            samples[sample][line[mm_i]][r"%ES_of_tot"] = line[es_s_i]
            
    HEADERS.sort()
    
    for sample in samples:
        orec = {"Sample" :sample}
        for mm in samples[sample]:
            orec[EU_MM_FORAMT%mm] = samples[sample][mm][r"%UE_of_tot"]
            orec[ES_MM_FORAMT%mm] = samples[sample][mm][r"%ES_of_tot"]
        orecs.append(orec)
    
    outputer.output([out_path], HEADERS,  orecs)
    

def get_sites_unique_to_sample(sample_sites_spec_file, problematic_samples, out_dir):
    data = [line for line in reader(open(sample_sites_spec_file), delimiter="\t")]
    problematic_lines = []
    non_problematic_lines = []
    
    for line in data[1:]:
        chrom = line[0]
        pos = line[2]
        mm = line[3]
        samples = line[4]
        other_samples = False
        for sample in samples.split(";"):
            if sample and sample not in problematic_samples:
                other_samples = True
                break
        if other_samples:
            non_problematic_lines.append(line)
        else:
            problematic_lines.append(line)
    with open(os.path.join(out_dir,'problematic_sites.csv'), 'wb') as out:
        out.write(",".join(data[0]) + "\n")
        out.write("\n".join([",".join(line) for line in problematic_lines]))
    with open(os.path.join(out_dir,'non_problematic_sites.csv'), 'wb') as out:
        out.write(",".join(data[0]) + "\n")
        out.write("\n".join([",".join(line) for line in non_problematic_lines]))
        

def filter_problematic_sites(problematic_sites_file, uniq_sites_spec_file, out_path):
    prob_sites = {}
    for line in reader(open(problematic_sites_file)):
        chrom = line[0]
        pos = line[2]
        mm = line[3]
        prob_sites.setdefault(chrom, {}).setdefault(pos, {})[mm] = False
    
    good_sites = []
    sites = [l for l in reader(open(uniq_sites_spec_file))]
    for site in sites[1:]:
        chrom = site[0]
        pos = site[2]
        mm = site[3]
        strand = site[4]
        if strand == "-":
            mm = reverse_strand(mm, False, True)
        
        if prob_sites.get(chrom, {}).get(pos, {}).get(mm, True):
            good_sites.append(site)
    
    with open(out_path, 'w') as out:
        out.write(",".join(sites[0]) + "\n")
        out.write("\n".join([",".join(line) for line in good_sites]))
        
    