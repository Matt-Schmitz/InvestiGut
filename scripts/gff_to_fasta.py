#!/usr/bin/env python3

import re
from Bio import SeqIO
from collections import defaultdict
import gzip
from investigut_utils import zopen, translation
import sys

original_fasta_files = sys.argv[1]
dir_base = sys.argv[2]
fasta_files = sys.argv[3].split()
predictions = sys.argv[4].split()

fasta_gff_list=[]
for fasta_file in fasta_files:
    org_group, codon_table, *_ = fasta_file.split("_")
    gff_list = []
    for gff_file in predictions:
        tool_gff, org_group_gff, codon_table_gff, *_ = gff_file.split("_")
        if org_group == org_group_gff and codon_table == codon_table_gff:
            gff_list.append(gff_file)
    fasta_gff_list.append([fasta_file, gff_list])


database = defaultdict(list)
get_filename_without_extension = lambda filepath: re.search(r'([^/]+?)(?:\.fa|\.fna|\.fasta|\.gff)?(?:\.gz)?$', filepath).group(1)
get_filename_after_second_underscore = lambda filepath: filepath[filepath.index('_', filepath.index('_') + 1) + 1:]
file=get_filename_after_second_underscore(get_filename_without_extension(fasta_files[0]))


def rev_comp(seq):
    return seq.upper()[::-1].translate(str.maketrans("ATCGRYKMVBHD","TAGCYRMKBVDH"))

def parse_gff(file, tool):
    res = []
    with (gzip.open(file, 'rt') if file.endswith('.gz') else open(file, 'r')) as file:
        if tool == "MetaGeneAnnotator":
            index=-10
            name = "ERROR: NO CONTIG NAME"
            reg_exp = re.compile(r"([^\t]+)\t(\d+)\t(\d+)\t([+-])\t")
            for i, line in enumerate(file):
                m=re.match(r"# ([^\t]+)", line)
                if line[0]=="#" and i-index>2:
                    name = re.match(r"\S+", m[1])[0]
                    index=i
                m=reg_exp.match(line)
                if m:
                    start = m[2]
                    stop = m[3]
                    strand = m[4]
                    res.append([name, [(int(start), int(stop))], strand])
        elif tool == "Augustus":
            contig_name = ""
            contig_strand = ""
            one_seq=[]
            has_start_stop = [False, False]
            reg_exp = re.compile(r"(\S+)[^\t]*\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t")
            for line in file:
                m=reg_exp.match(line)
                if m:
                    name = m[1]
                    seq_type = m[3]
                    start = m[4]
                    stop = m[5]
                    strand = m[7]
                    if name != contig_name:
                        contig_name = name
                        contig_strand = strand
                        one_seq = []
                        has_start_stop = [False, False]
                    if seq_type == "gene":
                        contig_name = name
                        contig_strand = strand
                        one_seq=[]
                        has_start_stop = [False, False]
                    elif seq_type == "CDS":
                        if contig_strand != strand:
                            raise Exception(f"strand of CDS different from gene for file: {file}     tool: {tool}     contig: {name}")
                        one_seq.append((int(start), int(stop)))
                    elif seq_type == "start_codon" and strand == contig_strand: has_start_stop[0] = True
                    elif seq_type == "stop_codon" and strand == contig_strand: has_start_stop[1] = True
                    if has_start_stop == [True, True]:
                        res.append([name, one_seq, strand])
                        contig_strand = ""
                        one_seq = []
                        has_start_stop = [False, False]
        elif tool == "Snap": 
            gene_num = ""
            contig_name = ""
            contig_strand = ""
            one_seq=[]
            has_start = False
            reg_exp = re.compile(r"(\S+)[^\t]*\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t([^\t]+)\t([^\t\n]+)\n")
            for line in file:
                m=reg_exp.match(line)
                if m:
                    name = m[1]
                    seq_type = m[3]
                    start = m[4]
                    stop = m[5]
                    strand = m[7]
                    this_gene_num = m[9]
                    if this_gene_num != gene_num:
                        gene_num = this_gene_num
                        contig_strand = strand
                        one_seq = []
                        has_start = False
                    if name != contig_name:
                        contig_name = name
                        contig_strand = strand
                        one_seq = []
                        has_start = False
                    if seq_type == "Esngl":
                        res.append([name, [(int(start), int(stop))], strand])
                        one_seq=[]
                        has_start = False
                    elif seq_type == "Einit":
                        has_start = True
                        one_seq=[(int(start), int(stop))]
                        contig_name = name
                        contig_strand = strand
                    elif seq_type == "Eterm":
                        if has_start and strand == contig_strand:
                            one_seq.append((int(start), int(stop)))
                            res.append([name, one_seq, strand])
                        has_start = False
                        one_seq = []
                    elif seq_type == "Exon": 
                        if has_start and strand == contig_strand:
                            one_seq.append((int(start), int(stop)))
                        else:
                            has_start = False
                            one_seq = []
        elif tool == "GeneMarkS":
            reg_exp = re.compile(r"(\S+)[^\t]*\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t")
            for line in file:
                m=reg_exp.match(line)
                if m:
                    name = m[1]
                    coding_seq = m[3]
                    if coding_seq != "CDS": continue
                    start = m[4]
                    stop = m[5]
                    strand = m[7]
                    res.append([name, [(int(start), int(stop))], strand])
        else:
            reg_exp = re.compile(r"(\S+)[^\t]*\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t")
            for line in file:
                m=reg_exp.match(line)
                if m:
                    name = m[1]
                    start = m[4]
                    stop = m[5]
                    strand = m[7]
                    res.append([name, [(int(start), int(stop))], strand])
    return res


def get_genes(one_metagenome_list):

    fasta_index = 1
    #print(len(human_contigs), subfolders, file)


    with open(file + ".fa", 'w') as new_file:   
        for gen_file, gff_files in one_metagenome_list:
            start_stop_dict = dict() ## example: {((100),((150,200),(250,300)),"+", "name"):[50, ["tool1","tool2"]],   ((50,100),(150),"-"):[200, ["tool1","tool2"]]}
            fasta_dict = dict()
            with (gzip.open(gen_file, 'rt') if gen_file.endswith('.gz') else open(gen_file, 'r')) as handle:
                for seq_record in SeqIO.parse(handle, "fasta"):
                    fasta_dict[seq_record.id] = seq_record.seq
            for gff_file in gff_files:
                tool, org_group, codon_table, *_ = gff_file.split("_")
                gff_data = parse_gff(gff_file, tool)
                for name, start_stop_list, strand in gff_data:
                    start_stop_list = sorted(tuple(sorted([int(start), int(stop)])) for start, stop in sorted(start_stop_list))
                    if strand == "+":
                        key = (start_stop_list[0][1],tuple(start_stop_list[1:]), strand, name)
                        if start_stop_dict.get(key):
                            start_so_far, tool_list = start_stop_dict[key]
                            start_stop_dict[key] = [start_so_far + [start_stop_list[0][0]], tool_list+[tool]] # for troubleshooting : [tool+name+str(start_stop_list)]
                        else:
                            start_stop_dict[key] = [[start_stop_list[0][0]], [tool]] # for troubleshooting : [tool+name+str(start_stop_list)]
                    elif strand == "-":
                        key = (start_stop_list[-1][0], tuple(start_stop_list[:-1]), strand, name)
                        if start_stop_dict.get(key):
                            start_so_far, tool_list = start_stop_dict[key]
                            start_stop_dict[key] = [start_so_far + [start_stop_list[-1][1]], tool_list+[tool]] # for troubleshooting : [tool+name+str(start_stop_list)]
                        else:
                            start_stop_dict[key] = [[start_stop_list[-1][1]], [tool]] # for troubleshooting : [tool+name+str(start_stop_list)]
                    else:
                        raise Exception(f"Invalid strand ({strand}). It should be - or +.")             
            
            correct_gene = re.compile("([TCAG]{3}){2,}$")
            for (second_val, start_stop_vals, strand, name), [first_val_list, tools_list] in start_stop_dict.items():
                seq_name = fasta_dict[name]
                for first_val in sorted(first_val_list, reverse = (strand=="-")):
                    if strand == "+":
                        start_stop_list = [(first_val, second_val), *start_stop_vals]
                    if strand == "-":
                        start_stop_list = [*start_stop_vals, (second_val, first_val)]
                    if len(start_stop_list)==1:
                        start, stop = start_stop_list[0]
                        start, stop = int(start), int(stop)
                        if not fasta_dict.get(name):##############
                            raise Exception(f"{tools_list} {org_group} {codon_table} {file} {name}") ##########
                        seq = str(seq_name[start-1:stop])
                        if len(seq)%3 != 0: continue
                        if strand == "-": seq = rev_comp(seq)
                        #print(seq, f'>{file}|{dataset}|{subdataset}|{org_group}|{codon_table}|{"-".join(sorted(tools_list))}|||{name}')
                        seq_len=len(seq_name)
                        edge_start = int(min(int(first_val), int(second_val)))-1
                        edge_stop = seq_len- int(max(int(first_val), int(second_val)))
                        if strand == "-":
                            edge_start,edge_stop=edge_stop,edge_start
                        
                        if not correct_gene.match(seq):
                            continue
                        aa_seq, wrong_start, wrong_stop = translation(seq, int(codon_table))
                        if len(aa_seq)<5:
                            print(seq, aa_seq)
                        aa_seq= "M" + aa_seq[1:-1]
                        if "*" in aa_seq: continue
                        ############################################### name of contig added
                        if (not wrong_start) and (not wrong_stop) and (not edge_start <= 2):# and (not edge_stop <= 2):
                            new_file.write(f'>{fasta_index}|{file}|{org_group}|{codon_table}|{"-".join(sorted(tools_list))}|{name}\n{aa_seq}\n') # translation(seq, int(codon_table))     {name} is the contig name          |{seq_len}|||{second_val}|{start_stop_vals}|{strand}|{first_val}|{name} |||{edge_start}|{edge_stop}|||{wrong_start}|{wrong_stop}
                            fasta_index +=1
                            break
                    else:
                        whole_seq = []
                        #print(sorted(start_stop_list))
                        for start, stop in sorted(start_stop_list):
                            start, stop = sorted([int(start), int(stop)])
                            seq = str(seq_name[start-1:stop])
                            whole_seq.append(seq)
                        whole_seq = "".join(whole_seq)
                        if len(whole_seq)%3 != 0: continue
                        if strand == "-": whole_seq = rev_comp(whole_seq)
                        #print(seq, f'>{file}|{dataset}|{subdataset}|{org_group}|{codon_table}|{"-".join(sorted(tools_list))}|||multi {name}')
                        seq_len=len(seq_name)
                        edge_start = int(min(int(first_val), int(second_val)))-1
                        edge_stop = seq_len- int(max(int(first_val), int(second_val)))
                        if strand == "-":
                            edge_start,edge_stop=edge_stop,edge_start
                        aa_seq, wrong_start, wrong_stop = translation(whole_seq, int(codon_table))
                        if not correct_gene.match(whole_seq):
                            continue
                        aa_seq= "M" + aa_seq[1:-1]
                        if "*" in aa_seq: continue
                        ############################################### name of contig added
                        if (not wrong_start) and (not wrong_stop) and (not edge_start <= 2):# and (not edge_stop <= 2):
                            new_file.write(f'>{fasta_index}|{file}|{org_group}|{codon_table}|{"-".join(sorted(tools_list))}|{name}\n{aa_seq}\n') # translation(whole_seq, int(codon_table))b       {name} is the contig name             |||{edge_start}|{edge_stop}|||{wrong_start}|{wrong_stop}|multiple_exons\n{whole_seq}
                            fasta_index +=1
                            break


get_genes(fasta_gff_list)
