"""Contains processing functions."""

import warnings
import json #GK
import jsonpickle
from json import JSONEncoder
import timeit

from Bio import SeqIO, SearchIO, BiopythonWarning
from depot.R2V2IO import open_file, progress
warnings.simplefilter('ignore', BiopythonWarning)


def create_lca_dict(lca):
    lca["species"]=[]
    lca["genus"]=[]
    lca["family"]=[]
    lca["order"]=[]
    lca["class"]=[]
    return lca

def unique_filter(filel, fileo, maxee, min_un_size, min_num_libs):
    """Dereplicates reads."""
    ambiguous_dna = ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
    num_seqs = 0
    data = {}

    start = timeit.default_timer()
    jb = 0
    # -1 is ee, -2 is size
    for lib in filel:
        dna_removed = 0
        with open_file(filel[lib]) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                header = record.id.split("=")
                ee = float(header[-1])
                num_seqs += 1
                seq = str(record.seq)
                if any(b in seq for b in ambiguous_dna):
                    dna_removed += 1
                else:               
                    if seq not in data:
                        data[seq]=[0] * (len(filel)+2)
                    data[seq][-1] += ee
                    data[seq][-2] += 1
                    data[seq][jb] += 1
        #print(f"Removed {dna_removed} ambiguous dna seqs. Total dict size is {len(data)}.")
        jb += 1
        progress(jb,1,len(filel))

    stop = timeit.default_timer()
    #print('Time: ', stop - start)

    num_clusters = len(data)

    rem_clus_size = 0
    rem_seqs_size = 0
    rem_clus_lib = 0
    rem_seqs_lib = 0
    rem_clus_ee = 0
    rem_seqs_ee = 0

    for seq in list(data.keys()):
        mean_ee = data[seq][-1]/data[seq][-2]
        num_libs_passed = 0
        k = 0
        for lib in filel:
            if data[seq][k] >= min_un_size:
                num_libs_passed +=1
            k+=1

        if mean_ee > maxee or num_libs_passed < min_num_libs:
            if data[seq][-2] < min_un_size:
                rem_clus_size += 1
                rem_seqs_size += data[seq][-2]
            elif num_libs_passed < min_num_libs:
                rem_clus_lib += 1
                rem_seqs_lib += data[seq][-2]
            else:
                rem_clus_ee += 1
                rem_seqs_ee += data[seq][-2]
            del data[seq]


    data_size_sorted = dict(sorted(data.items(), key=lambda t:t[1][-2],reverse=True))

    stop = timeit.default_timer()
    #print('Time: ', stop - start)
    i = 0
    seq_lib_info = {}
    with open(fileo, "w", encoding = "UTF-8") as fo:
        for seq in data_size_sorted:
            i+=1
            sid = "UN_" + str(i)
            seq_lib_info[sid] = {}
            k = 0
            for lib in filel:
                seq_lib_info[sid][lib] = data_size_sorted[seq][k]
                k+=1
            fo.write(">" + sid + ";size=" + str(data_size_sorted[seq][-2]) + ";\n" + seq + "\n")

    kept_clus = num_clusters - rem_clus_size - rem_clus_ee - rem_clus_lib
    kept_seqs = num_seqs - rem_seqs_size - rem_seqs_ee - rem_seqs_lib

    # print("Removed " + str(rem_clus_size) + " clusters containing " + str(rem_seqs_size)
    # + " sequences (cluster_size<" + str(min_un_size) + ").")
    # print("Removed " + str(rem_clus_ee) + " clusters containing " + str(rem_seqs_ee)
    # + " sequences (ee>" + str(maxee) + ").")
    derep_res = "[!] Kept " + str(kept_clus) + " clusters with " + str(kept_seqs) + " sequences."
    stop = timeit.default_timer()
    #print('Time: ', stop - start)
    return seq_lib_info, derep_res

def orf_finder(filei,fileo):
    """Finds orfs."""
    seqs_nucl = {}
    seqs_aa = {}
    dna_removed = 0

    orf_stats = [0,0,0,0]

    ambiguous_dna = ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
    ambiguous_aa = ["*", "B", "X", "J", "Z"]

    with open_file(filei) as handle, open(fileo, "w", encoding = "UTF-8") as fo:
        for record in SeqIO.parse(handle, "fasta"):
            found = 0
            orf_seqs = []
            frames = []
            seq = record.seq
            seqs_nucl[record.id] = str(seq)
            if any(b in seq for b in ambiguous_dna):
                dna_removed += 1
            else:
                for frame_start in range(3):
                    orf_seq = seq[frame_start:].translate(table=5)
                    if not any(b in orf_seq for b in ambiguous_aa):
                        orf_seqs.append(orf_seq)
                        frames.append(frame_start)
                        found +=1

                orf_stats[found] += 1

            for orf_seq, frame in zip(orf_seqs, frames):
                orf_id = record.id + ";f=" + str(frame)+ ";orfs=" + str(found)
                seqs_aa[orf_id] = str(orf_seq)
                fo.write(">" + orf_id + "\n" + str(orf_seq) + "\n")


    orf_res = "[!] Ambiguous:" + str(dna_removed) + "|No orfs:" + str(orf_stats[0]) + "|One orf:" + str(orf_stats[1]) + "|Two orfs:" + str(orf_stats[2]) + "|Three orfs:" + str(orf_stats[3])
    return seqs_nucl, seqs_aa, orf_res


def parse_hmmer(filei,len_ratio,otu):
    """Parses hmmer output DB - TODO what if 2 orfs are valid."""
    valid_orfs = {}
    for qresult in SearchIO.parse(filei,"hmmer3-text"):
        if len(qresult)>0:
            for hit in qresult:
                for hsp in hit:
                    per = (hsp.query_end - hsp.query_start) / qresult.seq_len
                    if per>=len_ratio:
                        id_list = qresult.id.split(";")
                        if otu:
                            sid=id_list[0]+";"+id_list[1]+";"
                        else:
                            sid=id_list[0]
                        valid_orfs[sid]=qresult.id

    return valid_orfs

def parse_sintax(sintaxfile, sintax_cutoff):
    sintax = {}
    sintax2 = {}
    sorted_otus = {}
    with open_file(sintaxfile) as sintaxfile_h:
        for line in sintaxfile_h:
            line_cols = line.strip().split()
            header, size = line_cols[0].split(";")
            o, k = header.split("_")
            sorted_otus[int(k)] = line_cols[0]
            if len(line_cols) == 3:
                taxs = line_cols[1].split(",")
                taxa = []
                confs = []
                for tax in taxs:
                    tax_i = tax.split(":")
                    tax_a = tax_i[1].split("(")
                    tax_conf = tax_a[1].replace(")","")
                    taxon = tax_a[0]
                    if float(tax_conf) >= sintax_cutoff:
                        taxa.append(taxon)
                        confs.append(tax_conf)

                sintax[line_cols[0]] = "|".join(taxa) + "\t" + "|".join(confs)

            else:
                sintax[line_cols[0]] = "NA\tNA"
    sorted_otus = dict(sorted(sorted_otus.items()))
    for otu in sorted_otus.values():
        sintax2[otu] = sintax[otu]
    return sintax2

def parse_gk_sintax(sintaxfile, sintax_cutoff):
    sintax = {}
    sintax2 = {}
    sorted_otus = {}
    with open_file(sintaxfile) as sintaxfile_h:
        for line in sintaxfile_h:
            line_cols = line.strip().split()
            header, size = line_cols[0].split(";")
            o, k = header.split("_")
            sorted_otus[int(k)] = line_cols[0]
            if len(line_cols) == 3:

                sintax[line_cols[0]] = line_cols[1] + "\t" + line_cols[2]

            else:
                sintax[line_cols[0]] = "NA\tNA"
    sorted_otus = dict(sorted(sorted_otus.items()))
    for otu in sorted_otus.values():
        sintax2[otu] = sintax[otu]
    return sintax2

def nucl_aln_from_aa_aln(file_aa_aln, file_nt_fna):
    aln_seq = {}
    aa_seq = {}
    nucl_aln_seq = {}

    with open(file_aa_aln) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id
            seq = record.seq.replace("-","")
            seq_aln = "-" + record.seq + "-" #adding these gaps for overhang nts
            aln_seq[seq_id]=seq_aln.upper()
            aa_seq[seq_id]=seq.upper()


    with open(file_nt_fna) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id

            seq = record.seq
            correct_frame = 0
            for frame_start in range(3):
                my_seq = seq[frame_start:].translate(table=5)
                if aa_seq[seq_id] == my_seq:
                    correct_frame = frame_start

            begin_nucl = correct_frame

            end_nucl = (len(seq)-begin_nucl)%3


            ff=0

            if begin_nucl == 0:
                ff = 0
            elif begin_nucl == 1:
                first_triplet = "--" + str(seq[0:1])
                ff =1 
            else:
                first_triplet = "-" + str(seq[0:2])
                ff = 1

            ee=0
            if end_nucl == 0:
                ee = 0
            elif end_nucl == 1:
                end_triplet = str(seq[-1:]) + "--"
                ee=1
            else:
                end_triplet = str(seq[-2:]) + "-"
                ee=1

            aln_nn = []
                
            f = 0
            a = 0
            while f<len(aln_seq[seq_id]):
                triplet = str(seq[begin_nucl:begin_nucl+3])
                if aa_seq[seq_id][a:a+1] == aln_seq[seq_id][f:f+1]:
                    a += 1
                    f += 1
                    begin_nucl += 3
                    aln_nn.append(""+triplet)
                else:
                    aln_nn.append("---")
                    f+=1

            #adding overhang nts
            if ff:
                i = 0
                found = 0
                while i<len(aln_nn):
                    if (aln_nn[i] == "---"):
                        fs = i
                        found = 1
                    else:
                        i = len(aln_nn)
                    i+=1
                if found:
                    aln_nn[fs]=first_triplet
            if ee:
                found = 0
                i=-1
                while i<0:
                    if aln_nn[i] == "---":
                        fs = i
                        found = 1
                    else:
                        i=1
                    i-=1
                if found:
                    aln_nn[fs]=end_triplet
            
            nucl_aln_seq[seq_id] =  "".join(aln_nn)
    i = 0
    residue_per_column = {}
    for seq in nucl_aln_seq.values():
        if i == 0:
            num_residues = len(seq)
            for k in range(num_residues):
                residue_per_column[k] = set()
            i += 1
        residues = list(seq)
        for k in range(num_residues):
            residue_per_column[k].add(residues[k])

    nucl_alignment = {}
    for seq_id, seq in nucl_aln_seq.items():
        residues = list(seq)
        final_seq = ""
        for k in range(len(residues)):
            if residue_per_column[k] != set("-"):
                final_seq = final_seq + residues[k]
        nucl_alignment[seq_id] = final_seq

    return nucl_alignment





