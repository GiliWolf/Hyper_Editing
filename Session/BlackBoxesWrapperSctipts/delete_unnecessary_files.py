__author__ = 'Hillel'
# =====================imports=====================#
import os
import shutil
import argparse


# =====================constants===================#
TMP_FOLDER_NAME = "set_aside"

# =====================functions===================#


def move_if_exists(dir, fingerprint, just_print, recursive=False):
    entries = []
    if recursive:
        for root, dirs, files in os.walk(dir):
            for f in files:
                entries.append(os.path.join(root,f))
    else:
        entries = [os.path.join(dir, e)for e in os.listdir(dir)]
    set_aside_dir = os.path.join(dir,TMP_FOLDER_NAME)
    if not os.path.isdir(set_aside_dir) and not just_print:
        os.mkdir(set_aside_dir)

    for entry in entries:
        if fingerprint in entry:
            if just_print:
                print os.path.join(set_aside_dir,entry)
                continue
            shutil.move(entry, os.path.join(set_aside_dir,os.path.basename(entry)))


def rm_everything_else(dir):
    entries = []
    for root, dirs, files in os.walk(dir):
            for f in files:
                entries.append(os.path.join(root,f))
        
    for entry in entries:
        if TMP_FOLDER_NAME not in entry:
            try:
                os.remove(os.path.join(dir,entry))
            except:
                shutil.rmtree(os.path.join(dir,entry), ignore_errors=True)


def delete_unnecessary_files(dir, fingerprints, just_print, dont_delete, recursive):
    
    for fingerprint in fingerprints:
        move_if_exists(dir, fingerprint, just_print, recursive)
    
    if just_print or dont_delete:
        return
    rm_everything_else(dir)

    set_aside_dir = os.path.join(dir,TMP_FOLDER_NAME)

    entries = os.listdir(set_aside_dir)
    for entry in entries:
        shutil.move(os.path.join(set_aside_dir,entry), os.path.join(dir,entry))

    os.rmdir(set_aside_dir)


if __name__ == "__main__":
    desc = "Remove any entry from from folder that doesn't match any of given names parts"


    parser = argparse.ArgumentParser(prog='RNAEditingQuantification',description=desc)

    parser.add_argument('-i', '--input_dir', metavar="input_path",dest='input_path', nargs='?', required=True,
                        help='The input dir to run on.')
    parser.add_argument('-f', '--fingerprints', metavar="fingerprints", dest='fingerprints', nargs='+', required=True,
                        help='name parts to recognize and kepp files according to.')
    parser.add_argument('-p', '--print_run', dest='just_print', required=False,
                        help='If set will only print the files not deleted', action='store_true', default=False)
    parser.add_argument('-d', '--dont_delete',  dest='dont_delete', required=False,
                        help='If set will only move the white list files to "<input_path>/%s"' % TMP_FOLDER_NAME, action='store_true', default=False)
    parser.add_argument('-r', '--recursive',  dest='recursive', required=False,
                        help='If set will run recursively into <input_path>' , action='store_true', default=False)
                        
    options = parser.parse_args()
    delete_unnecessary_files(options.input_path, options.fingerprints, options.just_print, options.dont_delete, options.recursive)