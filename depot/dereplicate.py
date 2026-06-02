#!/usr/bin/env python

"""
Usage:
  R2V2 dereplicate -f <FILE> -o <DIR> -c <FILE> 

Options:
  -h, --help                show help.
  -f, --files_info <FILE>   A file containing files to parse.
  -o, --output <DIR>        output directory.
  -c, --config_file <FILE>  config file
"""

import os, sys

from docopt import docopt

import timeit
import yaml

from depot.R2V2IO import check_programs, get_outdir, open_file
from depot.R2V2Lib import unique_filter, orf_finder
from depot.R2V2Soft import run_hmmer


def main():
    """Performs global dereplication and orf filtering"""

    check_programs("hmmscan")

    args = docopt(__doc__)

    inputfile = args['--files_info']
    output_dir = get_outdir(args['--output'])
    config_yaml = args['--config_file']

    stream = open(config_yaml, 'r')
    config_opts = yaml.safe_load(stream)
    stream.close()

    num_threads = int(config_opts["num_threads"])
    maxee = float(config_opts["maxee"])
    min_un_size = int(config_opts["min_un_size"])
    min_num_libs = int(config_opts["min_num_libs"])
    orf = config_opts["orf"]
    hmm_file = config_opts["hmm_file"]
    len_ratio = float(config_opts["len_ratio"])

    libs, file_lib_len_list = get_lib_file_names(inputfile)

    print ("[+] Running global dereplication.")
    derep_output = os.path.join(output_dir, "Derep.fasta")
    seq_lib_info, derep_res = unique_filter(file_lib_len_list, derep_output, maxee, min_un_size, min_num_libs)
    print(derep_res)

    abundance_file = os.path.join(output_dir, "Abundance.tsv")

    with open(abundance_file, "w", encoding = "UTF-8") as abundance_file_h:
        abundance_file_h.write("#ID\t" + "\t".join(libs) + "\n")
        for sid, sid_libs in seq_lib_info.items():
            abundances = []
            for lib in libs:
                abundances.append(str(sid_libs[lib]))
            abundance_file_h.write(sid + "\t" + "\t".join(abundances) + "\n")

    if orf:

        print ("[+] Filtering based on orfs.")
        orfs_file = os.path.join(output_dir, "Derep.fasta.orfs")
        seqs_nucl, seqs_aa, orf_res = orf_finder(derep_output, orfs_file)
        print(orf_res)

        print ("[+] Running Hmmer.")
        hmm_scan = os.path.join(output_dir, "Derep.hmm.scan")
        valid_orfs, hmm_res = run_hmmer(num_threads, hmm_scan, hmm_file, orfs_file, len_ratio, otu=1)
        print(hmm_res)

        derep_nucl = os.path.join(output_dir, "Derep.valid.fasta")
        derep_aa = os.path.join(output_dir, "Derep.valid.aa.fasta")

        with open(derep_nucl, "w", encoding = "UTF-8") as derep_nucl_h, open(derep_aa, "w", encoding = "UTF-8") as derep_aa_h:
            for valid_id, valid_aa_id in valid_orfs.items():
                derep_nucl_h.write(">" + valid_id + "\n" + seqs_nucl[valid_id] + "\n")
                derep_aa_h.write(">" + valid_id + "\n" + seqs_aa[valid_aa_id] + "\n")


def get_lib_file_names(inputfile):
    """Get lib names"""
    libs = []
    file_lib_len_list = {}
    with open_file(inputfile) as inputfile_h:
        for line in inputfile_h:
            line = line.strip()
            cols = line.split()
            file_lib_len_list[cols[0]]=cols[1]
            libs.append(cols[0])
    return libs, file_lib_len_list

if __name__ == '__main__':
    main()
