#!/usr/bin/env python3

import re
from ete3 import NCBITaxa
import pickle
import multiprocessing
import sys

dir_base = sys.argv[1]
kraken_files = sys.argv[2:]

all_taxids = set()
lock = multiprocessing.Lock()

# adds taxids from Kraken2 metagenome predictions to all_taxids
def get_taxids(kraken_file):
    taxid_pattern = re.compile(r"C\t(.*?)\t.*?\(taxid (\d+)\)") 
    with open(kraken_file, 'r') as t:
        return [int(taxid_re.group(2)) for taxid in t if (taxid_re := taxid_pattern.match(taxid))]

with multiprocessing.Pool() as pool:
    results = pool.map(get_taxids, kraken_files)
    for taxids in results:
        with lock:
            all_taxids.update(taxids)

mag_taxids=dict()
# adds taxids from Leviatan MAGs to all_taxids
ncbi = NCBITaxa()
metadata_files = [f"{dir_base}/gtdbtk.bac120.summary.tsv", f"{dir_base}/gtdbtk.ar53.summary.tsv"]
for metadata_file in metadata_files:
    with open(metadata_file, 'r') as f:
        for line in f:
            m = re.match(r"(Rep_\d+)\t(\w__[^\t]+)\t", line)
            if m: 
                tax_names = m[2].split(";")
                tax_names = [re.sub(r"_\w","", re.sub(r"\w__", "", x)) for x in tax_names[::-1]]
                tax_dict = ncbi.get_name_translator(tax_names)
                for name in tax_names:
                    if name in tax_dict:
                        all_taxids.add(tax_dict[name][0])
                        mag_taxids[m[1]] =tax_dict[name][0]
                        break

with open("mag_taxids.pkl", 'wb') as f:
    pickle.dump(mag_taxids, f)


# makes conversion dictionary between taxids and codon tables
tax_gc_convert = dict()
with open(f"{dir_base}/nodes.dmp", 'r') as f:
    for line in f:
        m = re.match(r"([^\t\|]*?)\t\|\t[^\t\|]*?\t\|\t[^\t\|]*?\t\|\t[^\t\|]*?\t\|\t[^\t\|]*?\t\|\t[^\t\|]*?\t\|\t([^\t\|]*?)\t\|\t", line)
        node_taxid = int(m[1])
        node_gc = int(m[2])
        tax_gc_convert[node_taxid] = node_gc

tax_table_pairs = dict()

for tax in all_taxids:
    if tax == 0: continue
    tax_table_pairs[tax]=tax_gc_convert[tax]


with open('tax_table_pairs.pkl', 'wb') as file:
    pickle.dump(tax_table_pairs, file)