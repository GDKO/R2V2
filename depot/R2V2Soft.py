"""Contains external software functions."""

import subprocess
from depot.R2V2Lib import parse_hmmer

def run_pear(file_f, file_r, pear_output, min_ovlp, num_threads):
    """Runs Pear read merger."""
    pear_command = "pear -f " + file_f + " -r " + file_r + " -o " + pear_output + " -v " + str(min_ovlp) + " -c 93 -j " + str(num_threads)
    pear_result = subprocess.run(pear_command, shell=True, capture_output=True, check=True)
    pear_result_list = list(pear_result.stdout.decode('UTF-8').split("\n"))
    return "[!] " + pear_result_list[-8]

def run_cutadapt(primer_string, cutadapt_output, cutadapt_error_rate, num_threads, merged_reads):
    """Runs Cutadapt adapter trimmer."""
    cutadapt_command = "cutadapt --discard-untrimmed --report minimal --rc -g " + primer_string + " -o " + cutadapt_output + " -e " + str(cutadapt_error_rate) + " -j " + str(num_threads) + " " + merged_reads
    cutadapt_result = subprocess.run(cutadapt_command, shell=True, capture_output=True, check=True)
    cutadapt_result_list = list(cutadapt_result.stdout.decode('UTF-8').split("\t"))
    removed = int(cutadapt_result_list[10]) - int(cutadapt_result_list[16])
    return "[!] Removed " + str(removed) + " sequences."

def run_vsearch_filt_ee(cutadapt_output, filter_len_output, min_len, max_len, num_threads):
    """Filters trimmed reads based on length and adds ee value."""
    vsearch_command = "vsearch --fastx_filter " + cutadapt_output + " -fastaout " + filter_len_output + " --fastq_minlen " + str(min_len) + " --fastq_maxlen " + str(max_len) + " --fasta_width 0 --relabel R --threads " + str(num_threads) + " --fastq_qmax 93 --eeout"
    vsearch_result = subprocess.run(vsearch_command, shell=True, capture_output=True, check=True)
    vsearch_result_list = list(vsearch_result.stderr.decode('UTF-8').split("\n"))
    return "[!] " + vsearch_result_list[-2]

def run_vsearch_filt(cutadapt_output, filter_len_output, min_len, max_len, num_threads):
    """Filters trimmed reads based on length"""
    vsearch_command = "vsearch --fastx_filter " + cutadapt_output + " -fastaout " + filter_len_output + " --fastq_minlen " + str(min_len) + " --fastq_maxlen " + str(max_len) + " --fasta_width 0 --threads " + str(num_threads)
    vsearch_result = subprocess.run(vsearch_command, shell=True, capture_output=True, check=True)
    vsearch_result_list = list(vsearch_result.stderr.decode('UTF-8').split("\n"))
    return "[!] " + vsearch_result_list[-2]

def run_hmmer(num_threads, hmm_scan, coi_hmm, orfs_file, len_ratio, otu):
    """Runs hmmer."""
    hmmer_command = "hmmscan --cpu " + str(num_threads) + " -o " + hmm_scan + " " + coi_hmm  + " " + orfs_file
    hmmer_result = subprocess.run(hmmer_command, shell=True, capture_output=True, check=True)
    valid_orfs = parse_hmmer(hmm_scan,len_ratio,otu)
    hmm_res = "[!] " + str(len(valid_orfs)) + " sequences with valid orfs."
    return valid_orfs, hmm_res

def run_swarm(num_threads, uc_output, centroids, derep_file_loc):
    """Runs Swarm clustering."""
    swarm_command = "swarm -f -z -t " + str(num_threads) + " -u " + uc_output + " -w " + centroids + " " + derep_file_loc
    swarm_result = subprocess.run(swarm_command, shell=True, capture_output=True, check=True)
    swarm_result_list = list(swarm_result.stderr.decode('UTF-8').split("\n"))
    clus_res = "[!] " + swarm_result_list[-4]
    return clus_res

def run_vsearch_clus(derep_file_loc, pid, centroids, num_threads, uc_output):
    """Runs Vsearch otu clustering."""
    clus_command = "vsearch --cluster_size " + derep_file_loc + " --id " + str(pid) + " --fasta_width 0 --sizein --sizeout --centroids " + centroids + " --threads " + str(num_threads) + " --uc " + uc_output
    clus_result = subprocess.run(clus_command, shell=True, capture_output=True, check=True)
    clus_result_list = list(clus_result.stderr.decode('UTF-8').split("\n"))
    clus_res = "[!] " + clus_result_list[-3]
    return clus_res

def run_vsearch_denoise(derep_file_loc, alpha, min_den_size, centroids, num_threads, uc_output):
    """Runs Vsearch denoising."""
    clus_command = "vsearch --cluster_unoise " + derep_file_loc + " --unoise_alpha " + str(alpha) + " --minsize " + str(min_den_size) + " --fasta_width 0 --sizein --sizeout --centroids " + centroids + " --threads " + str(num_threads) + " --uc " + uc_output
    clus_result = subprocess.run(clus_command, shell=True, capture_output=True, check=True)
    clus_result_list = list(clus_result.stderr.decode('UTF-8').split("\n"))
    clus_res = "[!] " + clus_result_list[-3]
    return clus_res

def run_vsearch_chim(swarm_output, chim_output, num_threads):
    """Runs Vsearch chimera detection."""
    chim_command = "vsearch --uchime3_denovo " + swarm_output + " --fasta_width 0 --nonchimeras " + chim_output + " --threads " + str(num_threads)
    chim_result = subprocess.run(chim_command, shell=True, capture_output=True, check=True)
    chim_result_list = list(chim_result.stderr.decode('UTF-8').split("\n"))
    chim_res = "[!] " + chim_result_list[-6] + "\n" + "[!] " + chim_result_list[-4] + "\n" + "[!] " + chim_result_list[-3]
    return chim_res
