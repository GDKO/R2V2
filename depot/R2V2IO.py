"""Functions that have to do with IO operations."""

import sys
import os
import shutil
import gzip

from Bio import SeqIO

def check_programs(*programs):
    """Check if program is in path."""
    error_list = []
    for program in programs:
        if not shutil.which(program):
            error_list.append(f"\t[!] {program} not found! Please install and add to it to $PATH")
    if len(error_list)>0:
        print("\n".join(error_list))
        sys.exit()

def get_outdir(out_directory, add_dir=""):
    """generates output directory in case it does not exist."""
    if not isinstance(out_directory, str):
        print(f"\t[!] {out_directory} is NOT a directory! Please specify an output directory")
        sys.exit()
    elif os.path.isfile(out_directory):
        print(f"\t[!] {out_directory} is a File! Please specify an output directory")
        sys.exit()
    elif not os.path.exists(os.path.join(out_directory, add_dir)):
        os.mkdir(os.path.join(out_directory, add_dir))
        return os.path.abspath(os.path.join(out_directory, add_dir))
    else:
        return os.path.abspath(os.path.join(out_directory, add_dir))

def check_indir(input_dir):
    """check if directory exists."""
    if not os.path.exists(input_dir):
        print(f"\t[!] FATAL ERROR: '{input_dir}' directory not found")
        sys.exit()
    elif not os.path.isdir(input_dir):
        print(f"\t[!] FATAL ERROR: '{input_dir}' is not a directory")
        sys.exit()
    else:
        return input_dir

def open_file(fname):
    """Can open gzip files."""
    if fname.endswith('.gz'):
        return gzip.open(fname, 'rt')
    return open(fname, 'r', encoding = "UTF-8")

def convert(filei, formati, fileo, formato):
    """Converts betwenn file formats."""
    records = SeqIO.parse(filei, formati)
    SeqIO.write(records, fileo, formato)

def progress(iteration, steps, max_value, no_limit=False):
    """Progress bar woohoo!"""
    if int(iteration) == max_value:
        if no_limit:
            sys.stdout.write('\r')
            print ("[x] \t%d%%" % (100), end='\r')
        else:
            sys.stdout.write('\r')
            print ("[x] \t%d%%" % (100))
    elif int(iteration) % steps == 0:
        sys.stdout.write('\r')
        print ("[x] \t%d%%" % (float(int(iteration) / int(max_value)) * 100), end='\r')
        sys.stdout.flush()
    else:
        pass
