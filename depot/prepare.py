#!/usr/bin/env python

"""
Usage:
  R2V2 prepare -f <FILE> -o <DIR> -c <FILE> [--keep]

Options:
  -h, --help                show help.
  -f, --files_info <FILE>   A file containing files to parse.
  -o, --output <DIR>        output directory.
  -c, --config_file <FILE>  config file
  --keep                    keep intermediate files
"""

import os
import glob
import subprocess
import yaml

from pathlib import Path
from docopt import docopt

from Bio import SeqIO
from depot.R2V2IO import check_programs, get_outdir, open_file
from depot.R2V2Soft import run_pear, run_cutadapt, run_vsearch_filt_ee

def process_file(output_dir, l_type, lib, reads, params):
    
    min_ovlp, num_threads, primer_string, cutadapt_error_rate, min_len, max_len = params
    print ("[+] Now processing " + lib + ".")
    lib_folder = os.path.join(output_dir,lib)

    get_outdir(lib_folder)

    print ("[+] Merging reads with PEAR.")
    pear_output = os.path.join(lib_folder, lib + ".pear")
    pear_res = run_pear(reads[0], reads[1] , pear_output, min_ovlp, num_threads)
    print(pear_res)

    merged_reads = os.path.join(lib_folder, lib + ".pear.assembled.fastq")
    
    print ("[+] Adapter trimming with cutadapt.")
    cutadapt_output = os.path.join(lib_folder, lib + ".pear.trim.fastq")
    cutadapt_res = run_cutadapt(primer_string, cutadapt_output, cutadapt_error_rate, num_threads, merged_reads)
    print(cutadapt_res)
    
    print ("[+] Filtering based on length.")
    
    filter_len_output = os.path.join(lib_folder, lib + ".pear.trim.len.fasta")
    vsearch_filt_res = run_vsearch_filt_ee(cutadapt_output, filter_len_output, min_len, max_len, num_threads)
    print(vsearch_filt_res)

    return(filter_len_output)

def gzip_file(filter_len_output):
    gzip_filter_len_output = filter_len_output + ".gz"
    gzip_command = "gzip " + filter_len_output
    subprocess.run(gzip_command, shell=True, capture_output=True, check=True)
    return(gzip_filter_len_output)


def main():
    """Prepare files"""
    check_programs("pear", "cutadapt", "vsearch", "gzip")

    args = docopt(__doc__)

    inputfile = args['--files_info']
    output_dir = get_outdir(args['--output'])
    config_yaml = args['--config_file']

    stream = open(config_yaml, 'r')
    config_opts = yaml.safe_load(stream)
    stream.close()

    num_threads = int(config_opts["num_threads"])
    min_ovlp = int(config_opts["min_ovlp"])
    cutadapt_error_rate = float(config_opts["cutadapt_error_rate"])
    fprimer = config_opts["fprimer"]
    rprimer = config_opts["rprimer"]
    min_len = int(config_opts["min_len"])
    max_len = int(config_opts["max_len"])

    primer_string = fprimer + "..." + rprimer

    params = [min_ovlp, num_threads, primer_string, cutadapt_error_rate, min_len, max_len]

    libs = {}
    with open_file(inputfile) as inputfile_h:
        for line in inputfile_h:
            line = line.strip()
            if not len(line) == 0:
                cols = line.split()
                l_type = cols[0]
                lib = cols[1]
                libs[lib]={}
                libs[lib]["type"] = l_type
                libs[lib]["files"] = []
                i = 2
                while i<len(cols):
                    libs[lib]["files"].append(cols[i])
                    i += 1


    libs_prepared_loc = os.path.join(output_dir, "Samples.tsv")

    negative_controls = 0

    with open(libs_prepared_loc, "w", encoding = "UTF-8") as libs_prepared_loc_h:

        negative_controls_data = {}
        for lib in libs:
            l_type = libs[lib]["type"]
            reads = libs[lib]["files"]

            if l_type == "N":
                negative_controls = 1
                filter_len_output = process_file(output_dir, l_type, lib, reads, params)
                with open_file(filter_len_output) as handle:
                    negative_controls_seqs = {}
                    for record in SeqIO.parse(handle, "fasta"):
                        seq = str(record.seq)
                        if seq not in negative_controls_seqs:
                            negative_controls_seqs[seq] = 0
                        negative_controls_seqs[seq] += 1
                for seq, ab in negative_controls_seqs.items():
                    if seq not in negative_controls_data:
                        negative_controls_data[seq] = []
                    negative_controls_data[seq].append(ab)
                gzip_filter_len_output = gzip_file(filter_len_output)

                if not args['--keep']:
                    print ("[+] Deleting intermediate files.")
                    lib_folder = os.path.join(output_dir,lib)
                    files_to_delete = os.path.join(lib_folder, lib + ".*.fastq")
                    for file in glob.glob(files_to_delete):
                        os.remove(file)

        for lib in libs:
            negative_controls_max = {}
            for seq, ab in negative_controls_data.items():
                negative_controls_max[seq] = max(ab)

            l_type = libs[lib]["type"]
            reads = libs[lib]["files"]

            if l_type == "S":
                if len(negative_controls_max) > 0:
                    rep_id = "t"
                    lib_name = f"{lib}_{rep_id}"
                    filter_len_output_t = process_file(output_dir, l_type, lib_name, reads, params)
                    dir_to_remove = Path(os.path.join(output_dir,lib_name))
                    files_to_delete = os.path.join(dir_to_remove,lib_name + ".*.fast*")
                    lib_folder = os.path.join(output_dir,lib)
                    get_outdir(lib_folder)
                    filter_len_output  = os.path.join(lib_folder, lib + ".pear.trim.len.fasta")
                    negative_removed = 0
                    with open(filter_len_output , "w", encoding = "UTF-8") as fo, open_file(filter_len_output_t) as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            header = str(record.id)
                            seq = str(record.seq)
                            to_print = 1
                            if seq in negative_controls_max:
                                if negative_controls_max[seq] > 0:
                                    negative_controls_max[seq] -= 1
                                    negative_removed += 1
                                    to_print = 0
                            if to_print:
                                fo.write(f">{header}\n{seq}\n")
                    print (f"[!] Removed {negative_removed} reads based on negative controls.")
                    if not args['--keep']:
                        print ("[+] Deleting intermediate files.")
                        for file in glob.glob(files_to_delete):
                            os.remove(file)                        
                        dir_to_remove.rmdir()               
                else:
                    filter_len_output = process_file(output_dir, l_type, lib, reads, params)

                gzip_filter_len_output = gzip_file(filter_len_output)
                if not args['--keep']:
                    print ("[+] Deleting intermediate files.")
                    files_to_delete = os.path.join(lib_folder, lib + ".*.fastq")
                    for file in glob.glob(files_to_delete):
                        os.remove(file)

            if l_type == "R":
                i = 0
                rep_id = 1
                files_to_check = []
                seq_data = []
                dirs_to_remove = []                
                while i < len(reads):
                    reads_r = [libs[lib]["files"][i], libs[lib]["files"][i+1]]
                    lib_name = f"{lib}_{rep_id}"
                    dir_to_remove = Path(os.path.join(output_dir,lib_name))
                    dirs_to_remove.append(dir_to_remove)
                    filter_len_output = process_file(output_dir, l_type, lib_name, reads_r, params)
                    files_to_check.append(filter_len_output)
                    i += 2
                    rep_id += 1

                for file in files_to_check:
                    seq_data_set = set()
                    with open_file(file) as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            seq = str(record.seq)
                            seq_data_set.add(seq)
                    seq_data.append(seq_data_set)

                print (f"[+] Evaluating common reads between replicates.")
                new_set = set.intersection(*seq_data)

                lib_folder = os.path.join(output_dir,lib)
                get_outdir(lib_folder)

                filter_len_output  = os.path.join(lib_folder, lib + ".pear.trim.len.fasta")
                common_reads = 0
                total_reads = 0
                negative_removed = 0
                with open(filter_len_output , "w", encoding = "UTF-8") as fo:
                    for file in files_to_check:
                        negative_controls_max = {}
                        for seq, ab in negative_controls_data.items():
                            negative_controls_max[seq] = max(ab)
                        with open_file(file) as handle:
                            for record in SeqIO.parse(handle, "fasta"):
                                header = str(record.id)
                                seq = str(record.seq)
                                to_print = 0
                                if seq in new_set:
                                    to_print = 1
                                    if seq in negative_controls_max:
                                        if negative_controls_max[seq] > 0:
                                            negative_controls_max[seq] -= 1
                                            negative_removed += 1
                                            to_print = 0
                                if to_print:
                                    fo.write(f">{header}\n{seq}\n")
                                    common_reads += 1
                                total_reads += 1

                gzip_filter_len_output = gzip_file(filter_len_output)
                common_per = round(100 * common_reads / total_reads, 2)
                if len(negative_controls_data) > 0:
                    print (f"[!] Kept {common_per} % common reads between replicates. Removed {negative_removed} reads based on negative controls.")
                else:
                    print (f"[!] Kept {common_per} % common reads between replicates.")
                if not args['--keep']:
                    print ("[+] Deleting intermediate files.")
                    i = 0
                    while i<rep_id:
                        i += 1
                        lib_folder_name = f"{lib_folder}_{i}"
                        lib_name = f"{lib}_{i}"
                        files_to_delete = os.path.join(lib_folder_name, lib_name + ".*.fast*")
                        for file in glob.glob(files_to_delete):
                            os.remove(file)

                    for dir_to_remove in dirs_to_remove:
                        dir_to_remove.rmdir()


            if l_type != "N":
                libs_prepared_loc_h.write(lib + "\t" + gzip_filter_len_output + "\n")

            

    print ("[+] Done.")

if __name__ == '__main__':
    main()
