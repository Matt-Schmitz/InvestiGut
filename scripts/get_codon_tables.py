#!/usr/bin/env python3

import re
import pickle
import multiprocessing
import sys
from os.path import dirname, join

dir_base = sys.argv[1]
nodes_dmp = sys.argv[2]
kraken_files = sys.argv[3:]

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


# makes conversion dictionary between taxids and codon tables
tax_gc_convert = dict()
with open(nodes_dmp, 'r') as f:
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