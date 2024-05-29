#!/usr/bin/env python3

import re
import os
from ete3 import NCBITaxa
from collections import defaultdict
import pickle
import sys
from investigut_utils import zopen
from os.path import basename

fasta_file = str(sys.argv[1])
kraken_tree = sys.argv[2]
kraken_taxid = sys.argv[3]
dir_base = sys.argv[4]

def load_data(filename):
    with open(filename, 'rb') as file:
        data = pickle.load(file)
    return data

print("!!!!!!!!!!!!!!!!!!!!!!!")
print(dir_base)
print(os.path.join(dir_base, "/scripts/tax_table_pairs.pkl"))
tax_table_pairs = load_data(os.path.join(dir_base, "scripts/tax_table_pairs.pkl"))


def get_gen_code(taxid):
    if int(taxid)==0: return 11 # Default bacterial codon table for unclassified samples
    try: 
        return tax_table_pairs[int(taxid)]
    except:
        raise Exception(f'taxid {taxid} is not in genetic codes table')



ncbi = NCBITaxa()

def split_metagenomes(fasta_file):
    if os.path.getsize(kraken_tree)==0: return
    d = defaultdict(list)
    prev_fasta_desc = ""
    tax_reg = re.compile(r"[CU]\t(.*?)\t.*?\(taxid (\d+)\)") 
    with zopen(fasta_file, 'rt') as f:
        with open(kraken_taxid, 'r') as t:
            for taxid in t:
                taxid = taxid.rstrip()
                taxid_re = tax_reg.match(taxid)
                taxid_name = taxid_re[1]
                taxid_id = taxid_re[2]
                lineage = ncbi.get_lineage(taxid_id) if taxid_id != "0" else [0]
                org_group = "unknown" # Taxid 0 - unclassified     Taxid 131567 - cellular organisms
                if 9696 in lineage: continue # skip human contigs
                if 10239 in lineage: org_group = "Viruses"
                elif 2157 in lineage: org_group = "Archaea"
                elif 2759 in lineage: org_group = "Eukaryota"
                elif 2 in lineage: org_group = "Bacteria"

                fasta_desc = prev_fasta_desc if prev_fasta_desc else f.readline()
                fasta_name = re.match(r">(.*?)\s", fasta_desc)[1]
                if taxid_name != fasta_name: 
                    print("ids don't match!!!!")
                    continue
                whole_fasta_seq = []
                while (fasta_seq:=f.readline()) and fasta_seq[0]!=">":
                    whole_fasta_seq.append(fasta_seq)
                prev_fasta_desc=fasta_seq
                gen_code = str(get_gen_code(taxid_id))
                d[(fasta_file, org_group, gen_code)].append(fasta_desc)
                d[(fasta_file, org_group, gen_code)].append("".join(whole_fasta_seq))
    # (?:\.fa|\.fna|\.fasta)
    pattern = r'^(.+?)(?:\.gz)?$'
    for (file, org_group, gen_code), contents in d.items():         
        with open(f'{org_group}_{gen_code}_{re.match(pattern, basename(file))[1]}', 'w') as new_file:
            new_file.write("".join(contents))



split_metagenomes(fasta_file)
