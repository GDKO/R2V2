#!/usr/bin/env python

"""
Usage:
  R2V2 cluster -f <FILE> -a <FILE> -o <DIR> -c <FILE> [-p <FILE>]

Options:
  -h, --help                           show help.
  -f, --dereplicated_nucl_file <FILE>  the dereplicate file in fasta
  -a, --abundance_info <FILE>          file containing the abundance values
  -o, --output <DIR>                   output directory.
  -c, --config_file <FILE>             config file
  -p, --dereplicated_prot_file <FILE>  the dereplicate protein file in fasta
"""

import sys
import os
import yaml

from docopt import docopt
from Bio import SeqIO

from depot.R2V2IO import check_programs, get_outdir, open_file
from depot.R2V2Soft import run_swarm, run_vsearch_chim, run_vsearch_clus, run_vsearch_denoise


def main():
    """Performs global dereplication, orf filtering, and OTU clustering"""

    check_programs("vsearch", "hmmscan" "swarm")

    args = docopt(__doc__)
    derep_file_nucl_loc = args['--dereplicated_nucl_file']
    derep_file_prot_loc = args['--dereplicated_prot_file']
    abundance_file = args['--abundance_info']
    output_dir = get_outdir(args['--output'])
    config_yaml = args['--config_file']

    stream = open(config_yaml, 'r')
    config_opts = yaml.safe_load(stream)
    stream.close()

    num_threads = int(config_opts["num_threads"])
    cl_mode = config_opts["cl_mode"]
    pid = float(config_opts["pid"])
    alpha = int(config_opts["alpha"])
    min_den_size = int(config_opts["min_den_size"])


    valid_cl_modes = ["swarm", "vsearch_otu", "vsearch_denoise"]
    if cl_mode not in valid_cl_modes:
        sys.exit("\t[!] FATAL ERROR: Clustering mode should be either swarm or vsearch_otu or vsearch_denoise")


    if derep_file_prot_loc:
        aa_seqs = get_aa_seqs(derep_file_prot_loc)
    else:
        aa_seqs = False

    libs, seq_lib_info = get_seq_lib_abundance(abundance_file)

    if cl_mode == "swarm":
        print("[+] OTU clustering with swarm.")
        centroids = os.path.join(output_dir, "Derep.valid.cl.fasta")
        uc_output = os.path.join(output_dir, "Derep.valid.cl.uc")
        clus_res = run_swarm(num_threads, uc_output, centroids, derep_file_nucl_loc)
    elif cl_mode == "vsearch_otu":
        print("[+] OTU clustering with vsearch.")
        centroids = os.path.join(output_dir, "Derep.valid.cl.fasta")
        uc_output = os.path.join(output_dir, "Derep.valid.cl.uc")
        clus_res = run_vsearch_clus(derep_file_nucl_loc, pid, centroids, num_threads, uc_output)
    elif cl_mode == "vsearch_denoise":
        print("[+] ASV clustering with vsearch.")
        centroids = os.path.join(output_dir, "Derep.valid.cl.fasta")
        uc_output = os.path.join(output_dir, "Derep.valid.cl.uc")
        clus_res = run_vsearch_denoise(derep_file_nucl_loc, alpha, min_den_size, centroids, num_threads, uc_output)
    print(clus_res)

    print("[+] Chimera detection with vsearch.")
    chim_output = os.path.join(output_dir, "Derep.valid.cl.ch.fasta")
    chim_res = run_vsearch_chim(centroids, chim_output, num_threads)
    print(chim_res)


    ## Get lib stats
    cluster_libs = get_stats_from_uc(uc_output)

    # NEED TO CHECK THIS
    ## Rename to OTU and output lib statistics file
    rename_and_ouput_stats_file(output_dir, libs, chim_output, cluster_libs, seq_lib_info, cl_mode, aa_seqs)

def get_seq_lib_abundance(inputfile):
    """Get lib names"""
    libs = []
    seq_lib_info = {}
    with open_file(inputfile) as inputfile_h:
        for line in inputfile_h:
            line = line.strip()
            cols = line.split()
            if cols[0].startswith("#"):
                i = 1
                while i<len(cols):
                    libs.append(cols[i])
                    i += 1
            else:
                seq_lib_info[cols[0]]={}
                i = 1
                while i<len(cols):
                    seq_lib_info[cols[0]][libs[i-1]]=int(cols[i])
                    i += 1
    return libs, seq_lib_info

def get_stats_from_uc(uc_output):
    """Returns a dictionary with clustering info"""
    cluster_libs = {}
    with open_file(uc_output) as uc_output_h:
        for line in uc_output_h:
            cols = line.split()
            if cols[0] == "S":
                c_id = cols[8].split(";")
                cluster_libs[c_id[0]]=[]
                cluster_libs[c_id[0]].append(c_id[0])
            elif cols[0] == "H":
                s_id = cols[8].split(";")
                c_id = cols[9].split(";")
                cluster_libs[c_id[0]].append(s_id[0])

    return cluster_libs

def rename_and_ouput_stats_file(output_dir, libs, chim_output, cluster_libs, seq_lib_info, cl_mode, aa_seqs):
    """Renames sequences to OTUs and provides lib statistics"""

    if cl_mode == "vsearch_denoise":
        otu_file = os.path.join(output_dir, "ASVs.fasta")
        otu_file_aa = os.path.join(output_dir, "ASVs.aa.fasta")
        otu_file_stats = os.path.join(output_dir, "ASVs.lib.stats.tsv")
    else:
        otu_file = os.path.join(output_dir, "OTUs.fasta")
        otu_file_aa = os.path.join(output_dir, "OTUs.aa.fasta")
        otu_file_stats = os.path.join(output_dir, "OTUs.lib.stats.tsv")


    with (open_file(chim_output) as chim_output_h,
      open(otu_file, "w", encoding="UTF-8") as otu_file_h,
      open(otu_file_stats, "w", encoding="UTF-8") as otu_file_stats_h):

        if cl_mode == "vsearch_denoise":
            naming = "ASV"
        else:
            naming = "OTU"

        otu_file_stats_h.write("#" + naming + "\t" + "\t".join(libs) + "\n")
        i = 0
        
        if aa_seqs:
            otu_file_aa_h = open(otu_file_aa, "w", encoding="UTF-8")
        else:
            otu_file_aa_h = None

        try:
            for record in SeqIO.parse(chim_output_h, "fasta"):
                sid_l = record.id.split(";")
                size_l = sid_l[1].split("=")
                seq = str(record.seq)

                # Get the subset of unique IDs belonging to this specific cluster once
                current_cluster_ids = cluster_libs.get(sid_l[0], [])

                # Optimization: Pre-sum the library abundances for these specific unique IDs
                # This avoids nested loop overhead for every library sequence
                lib_totals = {
                    lib: sum(seq_lib_info[un_id].get(lib, 0) for un_id in current_cluster_ids)
                    for lib in libs
                }
                
                # Map the pre-calculated numbers directly into your string array
                abundances = [str(lib_totals[lib]) for lib in libs]

                i += 1
                new_sid = naming + "_" + str(i) + ";size=" + str(size_l[1])

                otu_file_h.write(">" + new_sid + "\n" + seq.upper() + "\n")
                
                if otu_file_aa_h:
                    otu_file_aa_h.write(">" + new_sid + "\n" + aa_seqs[sid_l[0]].upper() + "\n")
                    
                otu_file_stats_h.write(new_sid + "\t" + "\t".join(abundances) + "\n")
                
        finally:
            if otu_file_aa_h:
                otu_file_aa_h.close()

def get_aa_seqs(derep_file_prot_loc):
    """Return dictionary of orfs"""
    aa_seqs = {}
    with open_file(derep_file_prot_loc) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sid_l = record.id.split(";")
            aa_seqs[sid_l[0]]=str(record.seq)
    return aa_seqs

if __name__ == '__main__':
    main()
