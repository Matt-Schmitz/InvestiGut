#!/usr/bin/env python3

import re
import os
from collections import defaultdict
import pickle
import sys
from investigut_utils import zopen
from os.path import basename

fasta_file = str(sys.argv[1])
dir_base = sys.argv[2]



def load_data(filename):
    with open(filename, 'rb') as file:
        data = pickle.load(file)
    return data


tax_table_pairs = load_data(os.path.join(dir_base, "tax_table_pairs.pkl"))


def get_gen_code(taxid):
    if int(taxid)==0: return 11 # Default bacterial codon table for unclassified samples
    try: 
        return tax_table_pairs[int(taxid)]
    except:
        raise Exception(f'taxid {taxid} is not in genetic codes table')

mag_taxids = load_data(os.path.join(dir_base, "mag_taxids.pkl"))

def get_mag_taxid(name):
    return mag_taxids[name]





rep_org_group={"Bacteria":[], "Archaea":[]}
with open("/home/horneflablinux/Documents/workflows/investigut_pipeline/pipeline/gtdbtk.bac120.summary.tsv", 'r') as f:
    for line in f:
        m = re.match(r"(Rep_\d+)\t(\w__[^\t]+)\t", line)
        if m: 
            rep_name = m[1]
            rep_org_group["Bacteria"].append(rep_name)
with open("/home/horneflablinux/Documents/workflows/investigut_pipeline/pipeline/gtdbtk.ar53.summary.tsv", 'r') as f:
    for line in f:
        m = re.match(r"(Rep_\d+)\t(\w__[^\t]+)\t", line)
        if m: 
            rep_name = m[1]
            rep_org_group["Archaea"].append(rep_name)



def split_mags(fasta_file):
    d = defaultdict(list)
    with zopen(fasta_file, 'rt') as f:
        pattern = r'^(.+?)(?:\.fa|\.fna|\.fasta)(?:\.gz)?$'
        org_name=re.match(pattern, basename(fasta_file))[1].split("__")[1]
        print(org_name)
        if org_name in rep_org_group["Bacteria"]: org_group = "Bacteria"
        elif org_name in rep_org_group["Archaea"]: org_group = "Archaea"
        else:
            raise Exception('MAG neither in Bacteria or Archaea')

        whole_fasta_file = []
        while (one_line:=f.readline()):
            whole_fasta_file.append(one_line)
        gen_code =str(get_gen_code(get_mag_taxid(org_name)))
        d[(fasta_file, org_group, gen_code)].append("".join(whole_fasta_file))
    #(?:\.fa|\.fna|\.fasta)
    pattern = r'^(.+?)(?:\.gz)?$'
    for (file, org_group, gen_code), contents in d.items():         
        with open(f'{org_group}_{gen_code}_{re.match(pattern, basename(file))[1]}', 'w') as new_file:
            new_file.write("".join(contents))


split_mags(fasta_file)

