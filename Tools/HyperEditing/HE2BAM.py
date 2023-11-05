#! /private/common/Software/anaconda/anaconda3/envs/python3/bin/python
#/usr/bin/python3.6

import os
import re
import subprocess
import sys
import traceback
import gc
from datetime import datetime

from Bio.Seq import Seq


# RUNNING_COMMAND = "Running: %(cmd)s"
# ERROR_MSG = "Command %(cmd)s; Exited With Code: %(ret_code)s"

def execute(command, input_str=None):
    try:
        # info
        # logging.debug(RUNNING_COMMAND % {'cmd': command})
        sys.stdout.flush()
        if input_str:
            p = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               input=input_str, universal_newlines=True, text=True)
        else:
            p = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
        # return output
        return p.stdout
    except subprocess.CalledProcessError as e:
        # logging.error(ERROR_MSG % {'cmd': command, 'ret_code': e.returncode})
        # print("ERROR:",  p.stderr)
        print(traceback.format_exc())
        exit(1)

HE_dir = sys.argv[1]  # full path to HE results folder
parm = sys.argv[2]  # 0.05_0.6_30_0.6_0.1_0.8_0.2
PE = sys.argv[3]  # SE or PE
samples = sys.argv[4]  # list of samples names (w/o path) separated by "," or "all" for analysing all the samples at the folder, for PE only the name w/o the -1/-2
out_sam_name = sys.argv[5]  # not including .sam
sam_for_header = sys.argv[6]  # bam file for the header of the converted sam
samtools_prog = sys.argv[7] if len(sys.argv) >= 8 else "samtools-1.8"
samtools_prog_params = "-l 9 -m 1G -@ 40"

sam_str_output = ""

if samples == "all":  # take all samples names from HE_dir, if PE names w/o -1/-2
    samples_str = ""
    d = HE_dir + "/UEdetect." + PE + "_" + parm

    for f in os.listdir(d):
        if ".UE.list" in f:
            if PE == "PE":
                samples_str = samples_str + "," + f[:-10]  # remove ".UE.list" at the end and -1/2
            else:
                samples_str = samples_str + "," + f[:-8]  # remove ".UE.list" at the end
    samples = samples_str[1:]  # remove the last ","

samples_list = samples.split(",")
samples_list = list(set(samples_list))  # uniq elements in the list in the case of PE

# print (samples_list)
if PE == "SE":
    for sample in samples_list:
        ReadsDict = {}
        file_name = HE_dir + "/UEdetect." + PE + "_" + parm + "/" + sample + ".UE.list"
        # print (file_name)
        with open(file_name) as listF:
            # T2C	(-)SRR1467567.18054784	chrX	39647356	CCAAATCTGTTGACTCCGACCTTCACCTTCCCCATGGTGTCTCAGGGATGTGGTTTGGCTGACGATGCAAAAGAAG	CCAAATCCGTTGACTCCGACCTTCACCTTCCCCATGGTGTCTCAGCGATGTGGCTCGGCTGGCGACGCAAAAGAAG	76	59	4	2
            for line in listF:
                fields = line.rstrip().split("\t")
                ReadsDict[fields[1][3:]] = {}
                ReadsDict[fields[1][3:]]["line"] = line.rstrip()
                ReadsDict[fields[1][3:]]["orientation"] = fields[1][1:2]

        file_name = HE_dir + "/AnalyseMM/" + sample + ".analyseMM"
        # print(file_name)
        with open(file_name) as analyseF:
            for line in analyseF:
                if "hits" in line:
                    fields = line.rstrip().split("\t")
                    if fields[0] in ReadsDict:
                        list_fields = ReadsDict[fields[0]]["line"].split("\t")
                        if (int(fields[
                                    1]) >> 4 & 1):  # sam flag bit "read reverse strand" is on, it means that read and quality sequences are reversed comp
                            if ReadsDict[fields[0]]["orientation"] == "-":
                                flag = "16"
                                quality_seq = fields[10]
                            else:
                                flag = "0"
                                quality_seq = fields[10][::-1]

                        else:  # sam flag bit "read reverse strand" is off, it means that read and quality sequences are not reversed comp
                            if ReadsDict[fields[0]]["orientation"] == "-":
                                flag = "16"
                                quality_seq = fields[10][::-1]

                            else:
                                flag = "0"
                                quality_seq = fields[10]

                        read_len = len(list_fields[5])
                        sam_fields = [fields[0], flag, list_fields[2], list_fields[3], "255", str(read_len) + "M", "*",
                                      "0", "0", list_fields[5], quality_seq]

                        sam_str_output = sam_str_output + "\t".join(sam_fields) + "\n"
                        ##out_sam.write("\t".join(sam_fields)+"\n")

else:  # PE
    for s in samples_list:
        for mate in ["-1", "-2"]:
            if mate == "-1":
                mapSamp = s + "-2"  # HE read is found at the first file (_1), so the mate should be found at the mapping of the second file
            else:  # mate==-2
                mapSamp = s + "-1"

            mateMapDict = {}
            map_files = [HE_dir + "/unMap/" + mapSamp + ".aln.bam", HE_dir + "/unMap/" + mapSamp + ".mem.bam"]

            for mf in map_files:
                # os.system("samtools-1.8 view " + mf + " > " + mapSamp +".sam") # convert bam to sam in order to read it ### P-open!!!
                mateMapF = execute(samtools_prog + " view " + mf).strip().split("\n")
                # with open(mapSamp +".sam") as mateMapF:
                for line in mateMapF:
                    # A00604:202:HLYW3DSXY:4:2639:20175:29606 16      chr17   77498956        37      151M    *       0       0       GCAGCTTCGGTGTGCAGATCATCCGTCCGTGTGGGGTTCTCAGTGCCGGAGGTCTTGGGGTGGGGGCCAGGCCTCGCACTTGCAGAGGAGCCCAGTGGGCTGCACGCTCCCCTCCATCCCCATCGGCCCTGTCCCCTGGAGTGTGTCAGAG        FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF        XT:A:U  NM:i:2  X0:i:1  X1:i:0  XM:i:2  XO:i:0  XG:i:0  MD:Z:27T24C98
                    # A00604:202:HLYW3DSXY:4:2639:21486:29528 4       *       0       0       *       *       0       0       GCTCTGGTCATGCAAAATAAACAAATTAAAAGATACAAAAAGCCACAGGGTTTCCTTCAGACCATTATATAATATCAACAGTTTTGAAAATTTTTAATGATATTTACACAAATTATCCAGGCAAAACATAACCGTAAAATAATTCAACAAA        FFF,,,FF:,FFF::FF,FFF,F,F,F,FFF,F:F,FFFFF,F:,:,,FFFFF,:,F,F,F,:,,,F:F,FF,F,,FF,F,:FF,FF,:,:,F,F,,,:F,F,:F,,F,FF,:,,,,,F,,:,F,FF:,,F:,FFFF,:,,,,,,F,,,FF
                    fields = line.rstrip().split("\t")
                    if (int(fields[1]) >> 2 & 1):  # sam flag bit "read unmapped" is on
                        continue

                    if (int(fields[
                                1]) >> 4 & 1):  # sam flag bit "read reverse strand" is on, it means read was mapped to antisense strand
                        strand = "-"
                        if mate == "-1":  # means this is the mapping of the second mate
                            flag = "147"  # flag for second and antisense
                        else:  # mate is "-2", means this is the mapping of the first mate
                            flag = "83"  # flag for first and antisense

                    else:  # read was mapped to sense
                        strand = "+"
                        if mate == "-1":  # means this is the mapping of the second mate
                            flag = "163"  # flag for second and sense
                        else:  # mate is "-2", means this is the mapping of the first mate
                            flag = "99"  # flag for first and sense

                    if fields[0] in mateMapDict:
                        mateMapDict[fields[0]]["locations"].append(fields[2] + "," + fields[3] + "," + strand)
                        mateMapDict[fields[0]]["lines"].append("\t".join([fields[0], flag, fields[2], fields[3], "255",
                                                                          fields[5], fields[6], fields[7], fields[8],
                                                                          fields[9], fields[10]]))

                    else:  # this is the first mapping of this read
                        mateMapDict[fields[0]] = {}
                        mateMapDict[fields[0]]["lines"] = ["\t".join([fields[0], flag, fields[2], fields[3], "255",
                                                                      fields[5], fields[6], fields[7], fields[8],
                                                                      fields[9], fields[10]])]
                        mateMapDict[fields[0]]["locations"] = [fields[2] + "," + fields[3] + "," + strand]

                    # add secondary alignments
                    # A00604:198:HVCNLDMXX:2:1102:14271:22106 16      chrX    39093158        0       151M    *       0       0       CGGATCGCGAGGTCAGGAGATCGAGACTATCCCGGCTAAAACGGTGAAGCCCCGTCTCTACTAAAAATGCAAAAAAATTGGCCAGGCGTAGTGGCGGGCGCCTGTGGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCC :FFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFF:FFFFFFF:FFFFFFFFFFFFFF:F::FFFFFFFFFFFFFF:FFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF XT:A:R  NM:i:6  X0:i:5  X1:i:0  XM:i:6  XO:i:0  XG:i:0  MD:Z:6A20C20A19A10A25A45        XA:Z:chr15,+43551696,151M,6;chr17,-565978,151M,6;chr12,+72381489,151M,6;chr7,-87244234,151M,6;
                    pattern = re.compile('XA:Z:(.+?)(\t|$)')
                    for f in fields:
                        if pattern.match(f):
                            # print(f[5:-1].split(";"))
                            for map in f[5:-1].split(";"):
                                # chr17,-565978,151M,6
                                # chr12,+72381489,151M,6

                                if map.split(",")[1][
                                   0:1] != strand:  # the secondary map is to the opposite sense verses the first map, means the seq and quality should be revers
                                    seq_Seq = Seq(fields[9])
                                    seq = str(seq_Seq.reverse_complement())
                                    quality_seq = fields[10][::-1]
                                else:
                                    seq = fields[9]
                                    quality_seq = fields[10]

                                if map.split(",")[1][0:1] == "-":  # the secondary map is to antisense
                                    if mate == "-1":  # means this is the mapping of the second mate
                                        flag = "147"  # flag for second and antisense
                                    else:  # mate is "-2", means this is the mapping of the first mate
                                        flag = "83"  # flag for first and antisense
                                else:  # the secondary map is to sense
                                    if mate == "-1":  # means this is the mapping of the second mate
                                        flag = "163"  # flag for second and sense
                                    else:  # mate is "-2", means this is the mapping of the first mate
                                        flag = "99"  # flag for first and sense

                                mateMapDict[fields[0]]["lines"].append(
                                    "\t".join([fields[0], flag, map.split(",")[0], map.split(",")[1][1:], "255",
                                               map.split(",")[2], fields[6], fields[7], fields[8], seq, quality_seq]))
                                mateMapDict[fields[0]]["locations"].append(
                                    map.split(",")[0] + "," + map.split(",")[1][1:] + "," + map.split(",")[1][0:1])

                # os.system("rm " + mapSamp + ".sam")
                del(mf)
                gc.collect()

            # print (mateMapDict)
            sample = s + mate
            ReadsDict = {}
            file_name = HE_dir + "/UEdetect." + PE + "_" + parm + "/" + sample + ".UE.list"
            # print (file_name)
            with open(file_name) as listF:
                for line in listF:
                    fields = line.rstrip().split("\t")
                    ReadsDict[fields[1][3:]] = {}
                    ReadsDict[fields[1][3:]]["line"] = line.rstrip()
                    ReadsDict[fields[1][3:]]["orientation"] = fields[1][1:2]

                    if fields[1][3:] in mateMapDict:
                        min_dis = 500000
                        ReadsDict[fields[1][3:]]["mate_line"] = ""
                        i = 0
                        for loc in mateMapDict[fields[1][3:]]["locations"]:
                            if loc.split(",")[0] == fields[2]:  # same chr
                                if loc.split(",")[2] != fields[1][1:2]:  # opposite orientation
                                    if abs(int(loc.split(",")[1]) - int(fields[3])) < min_dis:  # new location is min
                                        min_dis = abs(int(loc.split(",")[1]) - int(fields[3]))
                                        ReadsDict[fields[1][3:]]["mate_line"] = mateMapDict[fields[1][3:]]["lines"][i]
                            i = i + 1

            file_name = HE_dir + "/AnalyseMM/" + sample + ".analyseMM"
            # print(file_name)
            with open(file_name) as analyseF:
                for line in analyseF:
                    if "hits" in line:
                        fields = line.rstrip().split("\t")
                        if fields[0] in ReadsDict:
                            list_fields = ReadsDict[fields[0]]["line"].split("\t")
                            if (int(fields[
                                        1]) >> 4 & 1):  # sam flag bit "read reverse strand" is on, it means that read and quality sequences are reversed comp
                                if ReadsDict[fields[0]]["orientation"] == "-":
                                    quality_seq = fields[10]
                                else:
                                    quality_seq = fields[10][::-1]
                            else:  # sam flag bit "read reverse strand" is off, it means that read and quality sequences are not reversed comp
                                if ReadsDict[fields[0]]["orientation"] == "-":
                                    quality_seq = fields[10][::-1]
                                else:
                                    quality_seq = fields[10]

                            # determine the flag according to mate and orientation of mapping:
                            if ReadsDict[fields[0]]["orientation"] == "-":  # HE was mapped to the antisense
                                if mate == "-1":  # means this is HE from the first file
                                    flag = "83"  # flag for first and antisense
                                else:  # mate is "-2", means this is HE from the second file
                                    flag = "147"  # flag for second and antisense
                            else:  # # HE was mapped to the sense
                                if mate == "-1":  # means this is HE from the first file
                                    flag = "99"  # flag for first and sense
                                else:  # mate is "-2", means this is HE from the second file
                                    flag = "163"  # flag for second and sense

                            read_len = len(list_fields[5])
                            sam_fields = [fields[0], flag, list_fields[2], list_fields[3], "255", str(read_len) + "M",
                                          "*", "0", "0", list_fields[5], quality_seq]

                            if ReadsDict[fields[0]]["mate_line"] == "":
                                ##print ("mate is empty")
                                ##raise Exception (fields[0]+":"+" mate was not found")
                                continue  # logging warning

                            sam_str_output = sam_str_output + "\t".join(sam_fields) + "\n" + ReadsDict[fields[0]][
                                "mate_line"] + "\n"
                            ##out_sam.write("\t".join(sam_fields) + "\n")
                            ##out_sam.write(ReadsDict[fields[0]]["mate_line"] + "\n")

# os.system("samtools-1.8 view -H " + sam_for_header + " > " + out_sam_name + ".sam")
header = execute(samtools_prog + " view -H " + sam_for_header)
with open(out_sam_name + ".sam", "w") as out_sam:
    out_sam.write(header + sam_str_output)
# convert the result sam to sorted bam
##print ("samtools-1.8 view -bS "+ out_sam_name+".sam > " + out_sam_name+".bam")
os.system(" ".join([samtools_prog, "sort", "-O bam", samtools_prog_params, "-o", out_sam_name + ".sortedByCoord.out.bam", out_sam_name + ".sam"]))
os.system ("rm "+ out_sam_name+".sam")


### conf 1. this 2. sort 3. merge with given bam to given suff
