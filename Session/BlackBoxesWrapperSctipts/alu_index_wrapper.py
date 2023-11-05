__author__ = 'Hillel'
# =====================imports=====================#
import os
import argparse
import glob
import subprocess

# =====================constants===================#
IO_FILE = "alu_index_io_file.txt"
OUTPUT_FILE = "%s_alu_index_output.csv"
OUTPUT_ERROR_FILE = "%s_alu_index_error_log.log"
# =====================functions===================#
def get_i_paths(root_path, fp):
    paths = []
    for root, dirs, files in os.walk(root_path):
        for path in files:
            f_path = os.path.join(root, path)
            if fp in f_path:
                paths.append(f_path)

    return paths


def get_in_out_str(i_paths, part_a_output_dir, input_dir, truncs, max_rec_depth=5):
    in_out_str = ""
    samples_encountered = {}
    for path in i_paths:
        sample = os.path.basename(path)
        dir = os.path.dirname(path)
        o_p = dir.replace(input_dir, part_a_output_dir)
        if truncs[1] != 0:
            short_file_name = ".".join(sample.split(".")[:truncs[0]+1])[:0-truncs[1]]
        else:
            short_file_name = ".".join(sample.split(".")[:truncs[0]+1])
        if samples_encountered.get(short_file_name, False):
            continue
        else:
            samples_encountered[short_file_name] = True
        for i in xrange(max_rec_depth):
            path_parts = [o_p] + i*["*"] + [short_file_name]
            o_path = os.path.join(*path_parts)
            g_path = glob.glob(o_path)
            if g_path:
                o_path = g_path[0].replace(short_file_name,"")
                break
            else:
                o_path = None
        if o_path:
            in_out_str += o_path + "\t" + short_file_name + "\n"

    return in_out_str


def run_alu_index_part2(input_dir, part_a_output_dir, input_file_fingerprint, output_dir, part_b_script, truncs):
    input_paths = get_i_paths(input_dir, input_file_fingerprint)
    io_str = get_in_out_str(input_paths, part_a_output_dir, input_dir, truncs)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        os.chmod(output_dir, 511)

    io_file = os.path.join(output_dir, IO_FILE)
    with open(io_file, 'wb') as of:
        of.write(io_str)
    output_file = os.path.join(output_dir, OUTPUT_FILE % os.path.basename(part_a_output_dir))
    error_file = os.path.join(output_dir, OUTPUT_ERROR_FILE % os.path.basename(part_a_output_dir))
    call_str = " ".join(["perl516", part_b_script, io_file, output_file, error_file])
    subprocess.call(call_str, shell=True)


def convert_truncs(arg):
    return [int(arg.split(",")[0]), int(arg.split(",")[1])]

if __name__ == "__main__":
    desc = "create input file fo and run alu index"

    parser = argparse.ArgumentParser(prog='Alu Index Wrapper',description=desc)

    parser.add_argument('-i', '--input_dir', metavar="input_path",dest='input_path', nargs='?', required=True,
                        help='The dir containing the fastq files to run on.')
    parser.add_argument('-aio', '--alu_index_output_dir', metavar="alu index output_path",dest='ai_output_path', nargs='?', required=True,
                        help='The output dir to run on (must be identical to the dir of runAlu_partA pipeline).')
    parser.add_argument('-f', '--fingerprint', metavar="fingerprint", dest='fingerprint', nargs='?', required=True,
                        help='name parts to recognize files for input.')
    parser.add_argument('-o', '--output_dir', metavar="outputdir", dest='output_dir', nargs='?', required=True,
                        help='The path of the output dir to create.')
    parser.add_argument('-t', '--truncs', metavar="truncs",dest='truncs', nargs='?', required=False, type=convert_truncs,
                        default=[0,2], help="The truncation of the file name when created earlier on the pipeline"
                        "<suffixes to retain>,<chars to drop from filename>")
    parser.add_argument('-pbs', '--part_b_script', metavar="part_b_script", dest='part_b_script', nargs='?', required=False,
                        help='The path of the part B script.',
                        default=r"/home/alu/hillelr/scripts/AluIndex/get_sample_stat_new_index_PE.pl")
    options = parser.parse_args()
    run_alu_index_part2(options.input_path, options.ai_output_path, options.fingerprint, options.output_dir,
                        options.part_b_script, options.truncs)