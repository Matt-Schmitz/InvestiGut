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

rep_org_group={"Bacteria":[], "Archaea":[]}
with open(f"{dir_base}/data/metadata/gtdbtk.bac120.summary.tsv", 'r') as f:
    for line in f:
        m = re.match(r"(Rep_\d+)\t(\w__[^\t]+)\t", line)
        if m: 
            rep_name = m[1]
            rep_org_group["Bacteria"].append(rep_name)
with open(f"{dir_base}/data/metadata/gtdbtk.ar53.summary.tsv", 'r') as f:
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
        if org_name[:4]!="Rep_": raise Exception("wrong org name")
        gen_code ="4" if org_name in ["Rep_1090", "Rep_1284", "Rep_1667", "Rep_1668"] else "11"
        d[(fasta_file, org_group, gen_code)].append("".join(whole_fasta_file))
    #(?:\.fa|\.fna|\.fasta)
    pattern = r'^(.+?)(?:\.gz)?$'
    for (file, org_group, gen_code), contents in d.items():         
        with open(f'{org_group}_{gen_code}_{re.match(pattern, basename(file))[1]}', 'w') as new_file:
            new_file.write("".join(contents))


split_mags(fasta_file)

