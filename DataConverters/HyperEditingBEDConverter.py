__author__ = 'Hillel'
# =====================imports=====================#
from csv import reader
from collections import namedtuple
import re

from DataConverter import DataConverter

# =====================constants===================#
HyperEditingBEDRec = namedtuple("HyperEditingBEDRec", "chr start end sample reads_count neighbors strand")

SAMPLE_RE = re.compile(r"\)(.*?)\.")

CHR_I = 0
START_I = 1
END_I = 2
SAMPLE_I = 3
NEIGHBORS_I = 4
STRAND_I = 5

DICT_KEY_SEP = "#!@"
# error messages

# =====================classes=====================#


class HyperEditingBEDConverter(DataConverter):
    def __init__(self, use_folder_names=False):
        self.use_folder_names = use_folder_names

    def convert(self, he_bed, *args, **kwargs):

        recs = {}
        res = {}

        lines = [line for line in reader(he_bed, delimiter="\t")]

        for line in lines:
            if '' == line:
                continue
            chromosome = recs.setdefault(line[CHR_I], {})
            if "sample_name" not in kwargs and not self.use_folder_names:
                sample_name = SAMPLE_RE.findall(line[SAMPLE_I])
                sample_name = sample_name[0] if sample_name else "Unknown"
            else:
                if "sample_name" in kwargs:
                    sample_name = kwargs["sample_name"]
                else:
                    sample_name = "Unknown"
            sample = chromosome.setdefault(sample_name, {})

            s_key = DICT_KEY_SEP.join([line[START_I], line[END_I], line[STRAND_I]])

            if s_key in sample:
                count, neighbors = sample[s_key]
                count += 1
            else:
                neighbors = line[NEIGHBORS_I] if "no_neighbors" not in kwargs else "NA"
                count = 1

            sample[s_key] = [count, neighbors]
        if "unique" in kwargs:
            sample_double = {}
        for chromosome in recs:
            for sample in recs[chromosome]:
                for site in recs[chromosome][sample]:
                    start, end, strand = site.split(DICT_KEY_SEP)
                    count, neighbors = recs[chromosome][sample][site]

                    chr_res = res.setdefault(chromosome, {})
                    site_res = chr_res.setdefault(site, [])
                    if "unique" in kwargs:
                        if sample_double.get(site, {}).get(sample, False):
                            continue
                        else:
                            sample_double.setdefault(site, {})[sample] = True
                    site_res.append(HyperEditingBEDRec(chromosome, start, end, sample, count, neighbors, strand))#count would be wrong for merged clusters
        return res
