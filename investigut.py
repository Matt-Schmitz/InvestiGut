import glob
import scipy.stats
import numpy as np
from scipy import stats
from numpy import median
from statsmodels.stats.multitest import multipletests
import subprocess
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from itertools import combinations
from countrycode import country_convert
from ete3 import NCBITaxa
import re
from os.path import exists
from scipy.stats import mannwhitneyu
import os
import multiprocessing
import pickle
from investigut_utils import disk_cache
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
import seaborn as sns
import numpy as np
from statistics import mean, median
from scipy.stats import gmean
import argparse
from datetime import datetime


def main():
    parser = argparse.ArgumentParser(description='Investigut (Version 1.0.0)')
    parser.add_argument('-i', metavar='{INPUT}', required=True, help='input fasta file')
    parser.add_argument('-o', metavar='{OUTPUT}', help='output file')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', action='store_true', help='single mode (default)')
    group.add_argument('-m', action='store_true', help='multi mode')
    parser.add_argument('-d', metavar='{DMND_TSV_OUTPUT}', help='optional DIAMOND .tsv output file to skip alignment')


    args = parser.parse_args()

    combine_qseqids = False
    save_folder="."
    qseqid_file = args.i
    save_folder = datetime.now().strftime(f"%y_%m_%d_%H_%M_%S_investigut_output")
    if args.o:
        save_folder=args.o
    if args.s:
        pass
    if args.m:
        combine_qseqids = True
    dmnd_tsv=""
    if args.d:
        dmnd_tsv = args.d

    # Get Metagenome Lengths
    @disk_cache(include_func_in_hash=False)
    def metagenome_files_set():
        genome_base = "//home/matt/DATA/Proteins_small4/"
        fasta_files = glob.glob(genome_base + '**/*.fa', recursive=True)
        metagenome_files = set()
        for fasta_file in fasta_files:
            if "Leviatan" in fasta_file: continue
            study_name, sample_id = fasta_file.split("/")[-1].rsplit(".",1)[0].split("__")
            command = f"tail -n 2 {fasta_file}"
            process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell = False)
            stdout, stderr = process.communicate()
            try:
                file_length = int(stdout.split("|")[0][1:])
                if file_length >= 1000: metagenome_files.add((study_name, sample_id))
            except:
                print(fasta_file)
        return metagenome_files
        # len(metagenome_files) == 9419
    metagenome_files=metagenome_files_set()
    print("Finished obtaining metagenome lengths")

    # Get MAG Lengths
    len_lst=[]

    @disk_cache
    def mag_files_set():
        unread=0
        genome_base = "//home/matt/DATA/Proteins_small4/Isolate_collections/"
        fasta_files = glob.glob(genome_base + '**/*.fa', recursive=True)
        metagenome_files = set()
        for fasta_file in fasta_files:
            study_name, sample_id = fasta_file.split("/")[-1].rsplit(".",1)[0].split("__")
            command = f"tail -n 2 {fasta_file}"
            process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell = False)
            stdout, stderr = process.communicate()
            try:
                file_length = int(stdout.split("|")[0][1:])
                metagenome_files.add((study_name, sample_id))
                len_lst.append((file_length,sample_id))
            except:
                unread+=1
        if unread: print(f'{unread} MAGs could not be read.')
        return metagenome_files
        #len(metagenome_files) == 9419
    mag_files=mag_files_set()
    print("Finished obtaining MAG lengths")

    # Get Metadata
    genome_base = "/DATA/Matt/Genomes/Special_datasets/Human_passolli_assemblies/Meta-data/"
    meta_files = list(glob.iglob(genome_base + '**/*.tsv*', recursive=True))

    col_lengths=[]
    all_cols=[]
    for meta_file in meta_files:
        with open(meta_file, 'r') as f:
            first = f.readline().rstrip()
            cols = first.split("\t")
            if "ALT" in cols:
                cols[cols.index("ALT")] = "alt"
            if "non_westernized" in cols:###
                cols[cols.index("non_westernized")] = "westernized"###
            col_lengths.append(len(cols))
            all_cols.extend(cols)
    columns = Counter(all_cols).keys()
    col_name_to_num = {column:i for i, column in enumerate(columns)}
    all_patient_rows=[]
    for meta_file in meta_files:
        with open(meta_file, 'r') as f:
            first = f.readline().rstrip().split("\t")
            if "non_westernized" in first: first[first.index("non_westernized")] = "westernized"
            for line in f:
                patient = line.rstrip().split("\t")
                patient_row = [None]*len(columns)
                for col_val, pat_val in zip(first, patient):
                    if (not pat_val) or (pat_val =="NA"): pat_val = None
                    if col_val == "ALT": col_val = "alt"
                    if col_val == "westernized":
                        if pat_val=="yes": pat_val="no"
                        elif pat_val=="no": pat_val="yes"
                    patient_row[col_name_to_num[col_val]]=pat_val
                all_patient_rows.append(patient_row)

    def is_integer(some_string):
        try:
            int(some_string)
            return True
        except ValueError: 
            return False
        
    def is_float(some_string):
        try: float(some_string); return True
        except ValueError: return False

    for i, col in enumerate(list(zip(*all_patient_rows))):
        if all((x is None) or (is_integer(x)) for x in col): 
            for x in all_patient_rows:
                if x[i]: x[i] = int(x[i])
        elif all((x is None) or (is_float(x)) for x in col):
            for x in all_patient_rows:
                if x[i]: x[i] = float(x[i])

    print("Finished obtaining metadata")

    # Metadata Extractor
    class meta:

        def __init__(self, toindex: list):
            self.index = defaultdict(list)
            toindex_cols = [col_name_to_num[name] for name in toindex]
            for i, patient_row in enumerate(all_patient_rows):
                self.index[tuple(patient_row[c] for c in toindex_cols)].append(i)

        def keys(self, filter_none=False):
            return [k for k in self.index.keys() if all(e!=None for e in k)] if filter_none else self.index.keys() 
        def count(self, filter_none=False):
            return [(k, len(v)) for k,v in self.index.items() if (all(e!=None for e in k) if filter_none else True)]
        
        def get_metadata(self, wanted_index: tuple, columns = ()):
            if columns:
                res = []
                col_inds = [col_name_to_num[col] for col in columns]
                for x in self.index[wanted_index]:
                    one_patient_row = all_patient_rows[x]
                    res.append(tuple(one_patient_row[i] for i in col_inds))
                return res
            else:
                return list(tuple(all_patient_rows[x]) for x in self.index[wanted_index])

    print("Finished creating metadata extractor")

    # Name and ID Metadata
    meta_index = meta(["study_name", "sample_id"])
    print("Finished indexing metadata name/ID")

    # Get DMND Results
    #dmnd_tsv="/DATA/Matt/Diamond_databases/seaweed_50.tsv"
    #qseqid_file="/DATA/Matt/Diamond_databases/seaweed.tsv"
    all_qseqids=[]
    with open(qseqid_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                all_qseqids.append(line[1:].split()[0])

    def get_metadata_count(columns: tuple, all_not_none=False, keep_only=dict(), multi=combine_qseqids, sample_id_set=False):
        all_res = defaultdict(dict)
        keep_only_list=list(keep_only.values())
        for study_name, sample_id in metagenome_files:
            data = meta_index.get_metadata((study_name, sample_id), tuple(list(columns) + list(keep_only.keys())))
            if data and data[0][0] and (all(val!=None for val in data[0]) if all_not_none else True) and (all(x in keep_only_list[i] for i,x in enumerate(data[0][-len(keep_only):])) if keep_only else True):
                if sample_id_set:
                    if all_res[study_name].get(tuple(data[0][:-len(keep_only)] if keep_only else data[0])): all_res[study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])].add(sample_id)
                    else: all_res[study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])]={sample_id}          
                else:
                    if all_res[study_name].get(tuple(data[0][:-len(keep_only)] if keep_only else data[0])): all_res[study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])]+=1
                    else: all_res[study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])]=1
        def defaultdict2():
            return defaultdict(dict)
        dmnd_res=defaultdict(defaultdict2)
        dmnd_res2=defaultdict(defaultdict2)
        #print(sum(v for d in all_res.values() for v in d.values()))        # total count of specified columns in all metadata
        with open(dmnd_tsv, 'r') as file:
            seen = set()
            for line in file:
                qseqid, sseqid, _, _, _, _, _, _, _, _, evalue, bitscore, *_, aa_seq = line.rstrip().split("\t")
                evalue = float(evalue)
                bitscore = float(bitscore)
                fasta_index, file, org_group, codon_table, tools = sseqid.split("|")
                dataset, individual = file.split("__")
                if (dataset != "Leviatan") or ("Yachida" in file): ##################################################################
                    study_name, sample_id = file.split("__")
                    if (study_name, sample_id, qseqid) not in seen:
                        
                        if evalue < 1e-3 and (study_name, sample_id) in metagenome_files: 
                            seen.add((study_name, sample_id, qseqid))
                            if multi and (not all((study_name, sample_id, one_qseqid) in seen for one_qseqid in all_qseqids)): continue
                            if multi: qseqid = " + ".join(all_qseqids)
                            data = meta_index.get_metadata((study_name, sample_id), tuple(list(columns) + list(keep_only.keys())))
                            if data and data[0][0] and (all(val!=None for val in data[0]) if all_not_none else True) and (all(x in keep_only_list[i] for i,x in enumerate(data[0][-len(keep_only):])) if keep_only else True):
                                #if sample_id_set:
                                if dmnd_res2[qseqid][study_name].get(tuple(data[0][:-len(keep_only)] if keep_only else data[0])): dmnd_res2[qseqid][study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])].add(sample_id)
                                else: dmnd_res2[qseqid][study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])]={sample_id}
                                #else:
                                if dmnd_res[qseqid][study_name].get(tuple(data[0][:-len(keep_only)] if keep_only else data[0])): dmnd_res[qseqid][study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])]+=1
                                else: dmnd_res[qseqid][study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])]=1
                                if dmnd_res[qseqid][study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])] !=len(dmnd_res2[qseqid][study_name][tuple(data[0][:-len(keep_only)] if keep_only else data[0])]):
                                    raise Exception("something went wrong!")
            #print(sum(v for d in one_res.values() for v in d.values()))   # total count of specified columns in dmnd result metadata
        
        for study, count_dict in all_res.items():
            for k in count_dict.keys():
                for qseqid in all_qseqids:
                    if multi: qseqid = " + ".join(all_qseqids)
                    if not dmnd_res[qseqid][study].get(k):
                        dmnd_res[qseqid][study][k]= 0
                    if not dmnd_res2[qseqid][study].get(k):
                        dmnd_res2[qseqid][study][k]= set()
        return (all_res, dmnd_res2 if sample_id_set else dmnd_res)
    #get_metadata_count(("disease", "disease_subtype",))
    print("Finished obtaining DMND results")

    # Global Prevalence
    #meta_index.get_metadata(("SID103_12M", "stool"), ("antibiotics_current_use",))[0]
    all_res, dmnd_res = get_metadata_count(("study_name", "sample_id",), multi=combine_qseqids)
    for qseqid, one_res in dmnd_res.items():      
        found_count = sum(1 for study_name_dict in one_res.values() for name_id_count in study_name_dict.values() if name_id_count)
        global_count = sum(len(x) for x in all_res.values())
        #print(qseqid)
        #print(f"Global prevalence: {found_count}/{global_count} ({found_count/global_count*100 :.02f} %)")
        os.makedirs(f"./{save_folder}/{qseqid}", exist_ok=True)
        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            f.write(qseqid +"\n")
            f.write(f"Global prevalence: {found_count}/{global_count} ({found_count/global_count*100 :.02f} %)\n")
            f.write("\nStudy\tDisease\tPrevalence_Healthy\tPrevalence_Healthy_%\tPrevalence_Disease\tPrevalence_Disease_%\tp-value\tBH-adj._p-value\n")

    print("Finished writing global prevalences")

    # Helper Classes
    class MyTuple:
        def __init__(self, value1=0, value2=0):
            self.value1 = value1
            self.value2 = value2
        def add(self, x, y):
            self.value1 += x
            self.value2 += y
        def __getitem__(self, index):
            if index == "dmnd":
                return self.value1
            elif index == "all":
                return self.value2
            else:
                raise IndexError("MyTuple index not in ['dmnd','all']")
    class MyFisher:
        def __init__(self):
            self.value1 = np.array([[0,0],[0,0]])
        def add(self, x):
            self.value1 += x
        def get(self):
            return self.value1

    def add_sig(ylims, sig_combs):
        bottom, top = ylims
        y_range = top - bottom
        for i, significant_combination in enumerate(sig_combs):
            # Columns corresponding to the datasets of interest
            x1 = significant_combination[0][0]
            x2 = significant_combination[0][1]
            # What level is this bar among the bars above the plot?
            level = len(sig_combs) - i
            # Plot the bar
            bar_height = (y_range * 0.07 * level) + top
            bar_tips = bar_height - (y_range * 0.02)
            plt.plot(
                [x1, x1, x2, x2],
                [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k', scalex=False
            )
            # Significance level
            p = significant_combination[1]
            sig_symbol = 'NS'
            if p < 0.001:  sig_symbol = '***'
            elif p < 0.01: sig_symbol = '**'
            elif p < 0.05: sig_symbol = '*'
            text_height = bar_height + (y_range * 0)
            plt.text((x1 + x2) * 0.5, text_height, sig_symbol, ha='center', va='bottom', c='k')

    def pairs(lst):
        ls = list(range(1, len(lst) + 1))
        return ((lst[ls[x]-1], lst[ls[x + y-1]]) for y in reversed(ls) for x in range((len(ls) - y)))
    print("Finished loading helper classes/functions")


    # Disease Data
    ######################## Edit ################################
    keep_japan = {}#{"country":["JPN"]}

    ##############################################################

    all_res, dmnd_res = get_metadata_count(("disease", "disease_subtype",), multi=combine_qseqids, keep_only = keep_japan)

    def multi_get_disease(data, disease_cond, exclude_other_diseases=False, check_other_disease=[], subtype_cond=None, exclude_other_subtypes=False, check_other_subtype=[]):
        """
        When check_other_diseases is None, all additional diseases are accepted. When not None, then only these other diseases are accepted.
        """
        lst=[]#####
        res = 0
        for (disease, disease_subtype), value in data.items():
            if disease_cond in disease.split(";") and (disease_cond==disease if exclude_other_diseases else True):      #       (all(dst in disease_subtype.split(";") for dst in disease_)):
                if (check_other_disease==[]) or (set(check_other_disease + [disease_cond]) >= set(disease.split(";"))):
                    if (subtype_cond is None) or ((subtype_cond in disease_subtype.split(";")) and (subtype_cond==disease_subtype if exclude_other_subtypes else True)):
                        if (check_other_subtype==[]) or (set(check_other_subtype + [subtype_cond]) >= (set(disease_subtype.split(";")) if disease_subtype else set())):
                            lst.append((disease, disease_subtype, value))####
                            res += value
        return res,lst
    qseqid_list = []
    stored_studies = []
    stored_stats = []
    groups_list = []
    prev_list = []

    qseqid_list_multi = []
    stored_stats_multi = []
    groups_list_multi = []
    prev_list_multi = []
    # disease_cond, exclude_other_diseases=False, check_other_disease=[], subtype_cond=None, exclude_other_subtypes=False, check_other_subtype=[]
    run_list = {"CD":   ["IBD",         False,  ["perianal_fistula"],   "CD",   False,  []],
                "UC":   ["IBD",         False,  ["perianal_fistula"],   "UC",   False,  []],
            "CRC":   ["CRC",         True,   [],                     None,   False,  []],
            "CDI":   ["CDI",         True,   [],                     None,   False,  []],
            "T2D":   ["T2D",         True,   [],                     None,   False,  []],
                "RA":   ["RA",          True,   [],                     None,   False,  []],
    #"fatty_liver":   ["fatty_liver", True,   [],                     None,   False,  []]
    }

    def get_prevalence(array): # -> (healthy prev %, disease prev %)
        return {'healthy':(f'{array[0][1]}/{(array[1][1]+array[0][1])}', array[0][1]/(array[1][1]+array[0][1])*100), "disease":(f'{array[0][0]}/{(array[1][0]+array[0][0])}', array[0][0]/(array[1][0]+array[0][0])*100)}
    for qseqid, one_res in dmnd_res.items():
        #print(qseqid)
        multi_study_cd = np.array([[0,0],[0,0]])
        multi_study_uc = np.array([[0,0],[0,0]])
        disease_dict = defaultdict(MyFisher)
        healthy_dict = defaultdict(MyTuple)
        for study, data in all_res.items():
            if not data.get(("healthy", None)): continue
            dmnd_healthy = (one_res[study].get(("healthy",None)) or 0)
            all_healthy = data[("healthy", None)]
            healthy_perc = (dmnd_healthy/all_healthy)*100
            if multi_get_disease(data, *run_list["CD"])[0] or multi_get_disease(data, *run_list["UC"])[0]:
                healthy_dict["IBD"].add(dmnd_healthy,all_healthy)
            for key in run_list.keys():
                if multi_get_disease(data, *run_list[key])[0]:
                    healthy_dict[key].add(dmnd_healthy,all_healthy)

            for key, run_commands in run_list.items():
                if multi_get_disease(data, *run_commands)[0]: #, disease_cond=run_commands[0], exclude_other_diseases=run_commands[1], check_other_disease=run_commands[2], subtype_cond=run_commands[3], exclude_other_subtypes=run_commands[4], check_other_subtype=run_commands[5]
                    dmnd_disease,dmnd_disease_stats =  multi_get_disease(one_res[study], *run_commands)                #(one_res[study].get(("IBD","CD")) or 0)
                    all_disease,all_disease_stats = multi_get_disease(data, *run_commands)           #data[("IBD", "CD")]
                    disease_perc = (dmnd_disease/all_disease)*100
                    observed = [dmnd_disease,dmnd_healthy] # [Number of times specific protein appears in disease patients in study, Number of healthy patients in study]
                    Negative_samples = [all_disease-dmnd_disease,all_healthy-dmnd_healthy] 
                    disease_dict[key].add(np.array([observed,Negative_samples]))
                    if key=="CD":
                        multi_study_cd += np.array([observed,Negative_samples])
                    elif key=="UC":
                        multi_study_uc += np.array([observed,Negative_samples])
                    disease_stats = scipy.stats.fisher_exact(np.array([observed,Negative_samples]))[1]
                    qseqid_list.append(qseqid)
                    stored_studies.append(study)
                    stored_stats.append(disease_stats)
                    groups_list.append(key)
                    #prev_list.append((healthy_perc, CD_perc))
                    prev_list.append(get_prevalence(np.array([observed,Negative_samples])))
                    if disease_stats < 5:
                        pass
                        #print("   "+study)
                        #print (f'\tEnriched in {["healthy",key][disease_perc > healthy_perc] if disease_perc!=healthy_perc else "neither"}: ',disease_stats)
                        #print ('\tHealthy prevalence: ', healthy_perc)
                        #print (f'\t{key} prevalence: ', disease_perc)
                        #print(f'\t{get_prevalence(np.array([observed,Negative_samples]))}')
                        #print("\tdmnd", dmnd_disease_stats)
                        #print("\tall", all_disease_stats)
                        #print([observed,Negative_samples])
                

        multi_study_ibd = disease_dict["CD"].get()+disease_dict["UC"].get()
        multi_study_ibd[0][1]=healthy_dict["IBD"]["dmnd"]
        multi_study_ibd[1][1]=healthy_dict["IBD"]["all"]-healthy_dict["IBD"]["dmnd"]
        multi_study_ibd_stats = scipy.stats.fisher_exact(multi_study_ibd)[1]
        qseqid_list_multi.append(qseqid)
        stored_stats_multi.append(multi_study_ibd_stats)
        groups_list_multi.append("IBD")
        prev_list_multi.append(get_prevalence(multi_study_ibd))
        #print("multi IBD", get_prevalence(multi_study_ibd), multi_study_ibd_stats)

        for key in run_list.keys():
            multi_study_disease_stats = scipy.stats.fisher_exact(disease_dict[key].get())[1]
            qseqid_list_multi.append(qseqid)
            stored_stats_multi.append(multi_study_disease_stats)
            groups_list_multi.append(key)
            prev_list_multi.append(get_prevalence(disease_dict[key].get()))    
            #print(f"multi {key}", get_prevalence(disease_dict[key].get()), multi_study_disease_stats)
    print("Finished obtaining disease data")

    # Write Disease Data
    for qseqid, study, group, prev, orig, corr in zip(qseqid_list, stored_studies, groups_list, prev_list, stored_stats, list(multipletests(stored_stats, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1])):
        if corr < 5:   
            x_vals, y_vals, bar_labels = [], [], []
            x_vals.append("Healthy")
            x_vals.append(group)
            y_vals.append(prev["healthy"][1])
            y_vals.append(prev["disease"][1])
            bar_labels.append(prev["healthy"][0])
            bar_labels.append(prev["disease"][0])
            plt.figure(figsize=(2, 5))
            fig = sns.barplot(x=x_vals, y=y_vals, color='#D04356', width=0.8) ### leaving out bw_adjust leads to a very smooth curve but also negative gene lengths
            fig.bar_label(fig.containers[0], labels = bar_labels)
            fig.set(xlabel="Disease Status", ylabel="Prevalence (%)", title=f"Prevalence of {qseqid} in {group}\n({study})") #Prevalence of {qseqid} by Disease ({study})
            add_sig(fig.get_ylim(),[[(0,1),corr]])
            #plt.xlim(-100,5100)
            #plt.xticks(rotation=90)
            os.makedirs(f"./{save_folder}/{qseqid}/disease/one_study/", exist_ok=True)
            plt.savefig(f"./{save_folder}/{qseqid}/disease/one_study/Prevalence_of_{qseqid}_in_{group}_({study}).svg")
            plt.close()
            #plt.show() 
            #print (qseqid, study, group, prev, orig, corr)
            with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
                f.write(f'{study}\t{group}\t{prev["healthy"][0]}\t{prev["healthy"][1]}\t{prev["disease"][0]}\t{prev["disease"][1]}\t{orig}\t{corr}\n')
    for qseqid, group, prev, orig, corr in zip(qseqid_list_multi, groups_list_multi, prev_list_multi, stored_stats_multi, list(multipletests(stored_stats_multi, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1])):
        if corr < 5: 
            
            x_vals, y_vals, bar_labels = [], [], []
            x_vals.append("Healthy")
            x_vals.append(group)
            y_vals.append(prev["healthy"][1])
            y_vals.append(prev["disease"][1])
            bar_labels.append(prev["healthy"][0])
            bar_labels.append(prev["disease"][0])
            plt.figure(figsize=(2, 5))
            fig = sns.barplot(x=x_vals, y=y_vals, color='#D04356', width=0.8) ### leaving out bw_adjust leads to a very smooth curve but also negative gene lengths
            fig.bar_label(fig.containers[0], labels = bar_labels)
            fig.set(xlabel="Disease Status", ylabel="Prevalence (%)", title=f"Prevalence of {qseqid} in {group}") #Prevalence of {qseqid} in {group}         Bai Protein Set 
            add_sig(fig.get_ylim(),[[(0,1),corr]])
            #plt.xlim(-100,5100)
            #plt.xticks(rotation=90)
            os.makedirs(f"./{save_folder}/{qseqid}/disease/combined_studies/", exist_ok=True)
            plt.savefig(f"./{save_folder}/{qseqid}/disease/combined_studies/Prevalence_of_{qseqid}_in_{group}.svg")
            plt.close()
            #plt.show() 
            #print ("multi", qseqid, group, prev, orig, corr)
            with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
                f.write(f'multi_study\t{group}\t{prev["healthy"][0]}\t{prev["healthy"][1]}\t{prev["disease"][0]}\t{prev["disease"][1]}\t{orig}\t{corr}\n')

    print("Finished writing disease data")

    # Write Country Data
    all_res, dmnd_res = get_metadata_count(("country",), multi=combine_qseqids)
    stored_stats = []
    groups_list = []
    qseqid_list = []
    prev_list = []
    figure_data = dict()
    x_pos_dict = dict()
    for qseqid, one_res in dmnd_res.items():
        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            f.write("\nCountry_1\tCountry_2\tPrevalence_Country_1\tPrevalence_Country_1_%\tPrevalence_Country_2\tPrevalence_Country_2_%\tp-value\tBH-adj._p-value\n")
        combined_data_all = defaultdict(int)
        for study, data in all_res.items():
            for columns, count in data.items():
                combined_data_all[columns] += count
        combined_data_dmnd = defaultdict(int)
        for study, data in one_res.items():
            for columns, count in data.items():
                combined_data_dmnd[columns] += count

        fig_data = []
        for columns, all_count in combined_data_all.items():
            if all_count < 10: continue
            dmnd_count = combined_data_dmnd[columns]
            col_perc = (dmnd_count / all_count) * 100
            fig_data.append((columns, col_perc, (dmnd_count, all_count)))
        fig_data=sorted(fig_data, key=lambda x: (-x[1],x[0]))

        x_vals = [country_convert(" +\n".join(x)) for x,y,z in fig_data]
        x_pos_dict[qseqid]={i:x_val for x_val,i in enumerate(x_vals)}
        y_vals = [y for x,y,z in fig_data]
        fischer_vals = [[z[0], z[1]-z[0]] for x,y,z in fig_data]
        bar_labels = [f'{z[0]}/{z[1]}' for x,y,z in fig_data]
        figure_data[qseqid]=(x_vals, y_vals, bar_labels)
        for (cols1, fisch1), (cols2, fisch2) in pairs(list(zip(x_vals, fischer_vals))):
            stat = scipy.stats.fisher_exact(np.array([fisch1,fisch2]))[1]
            prev_list.append([fisch1,fisch2])
            groups_list.append((cols1, cols2))
            stored_stats.append(stat)
            qseqid_list.append(qseqid)

    corr_p_val_dict=defaultdict(list)
    for qseqid, group, orig, corr in zip(qseqid_list, groups_list, stored_stats, list(multipletests(stored_stats, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1])):
        if corr < 0.05:
            corr_p_val_dict[qseqid].append([(x_pos_dict[qseqid][group[0]],x_pos_dict[qseqid][group[1]]),corr])

    previous_qseqid=""

    sigletters=defaultdict(list)
    whichletter=dict()
    allletters="abcdefghijklmnopqrstuvwxyz"
    for qseqid, group, orig, corr in zip(qseqid_list, groups_list, stored_stats, list(multipletests(stored_stats, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1])):
        if corr >= 0.05: 
            first, second = group
            if second not in sigletters[first]: sigletters[first].append(second)
            if first not in sigletters[second]: sigletters[second].append(first)
    #print(sigletters)
    for qseqid, prev, group, orig, corr in zip(qseqid_list, prev_list, groups_list, stored_stats, list(multipletests(stored_stats, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1])):
        if qseqid != previous_qseqid:
            #print(qseqid)
            previous_qseqid = qseqid
            x_vals, y_vals, bar_labels = figure_data[qseqid]
            bar_labels_new=[]
            i=0
            for x_val, bar_label in zip(x_vals, bar_labels):
                next_letters=[]
                for same_group_x_val in sigletters[x_val]:
                    let = whichletter.get(same_group_x_val)
                    if let: 
                        i+=1
                        next_letters.append(allletters)
                if x_val == "Denmark": 
                    pass#print("test",next_letters)
                if not next_letters: 
                    next_letters.append(allletters[i])
                    whichletter[x_val]=allletters[i]
                    i+=1
                bar_labels_new.append("".join(next_letters)+"\n"+bar_label)

            plt.figure(figsize=(17, 8))
            fig = sns.barplot(x=x_vals, y=y_vals, color='#E29953', width=0.8) ### leaving out bw_adjust leads to a very smooth curve but also negative gene lengths
            fig.bar_label(fig.containers[0], labels = bar_labels) ########################### bar_labels_new
            fig.set(xlabel="Country", ylabel="Prevalence (%)", title=f"Prevalence of {qseqid} by Country")
            #plt.xlim(-100,5100)
            plt.xticks(rotation=90)
            #add_sig(fig.get_ylim(),corr_p_val_dict[qseqid])
            os.makedirs(f"./{save_folder}/{qseqid}/country/", exist_ok=True)
            plt.savefig(f"./{save_folder}/{qseqid}/country/Prevalence_of_{qseqid}_by_Country.svg")
            #plt.show() 
            plt.close()
        [[c1_yes, c1_no],[c2_yes,c2_no]] = prev
        
        #if corr < 0.05: print ("  ", f'{group[0]}\t{group[1]}\t{c1_yes}/{c1_yes+c1_no}\t{c1_yes/(c1_yes+c1_no)}\t{c2_yes}/{c2_yes+c2_no}\t{c2_yes/(c2_yes+c2_no)}\t{orig}\t{corr}')   

        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            f.write(f'{group[0]}\t{group[1]}\t{c1_yes}/{c1_yes+c1_no}\t{100*c1_yes/(c1_yes+c1_no)}\t{c2_yes}/{c2_yes+c2_no}\t{100*c2_yes/(c2_yes+c2_no)}\t{orig}\t{corr}\n')


    # Write Other Demographic Factors
    global_features = [("smoker","ever_smoker",),("BMI",),("gender",),("westernized",),("age",),("antibiotics_current_use",),("ajcc",)]##age_category   ,("ajcc",),("ctp",)    ("smoker","ever_smoker",),("BMI",),("gender",),("westernized",),("age",),("antibiotics_current_use",),("ajcc",)
    color_map= {("ctp",):"#494E6B",("smoker","ever_smoker",):"#494E6B",("BMI",):"#494E6B",("gender",):"#494E6B",("westernized",):"#494E6B",("age",):"#494E6B",("antibiotics_current_use",):"#494E6B",("ajcc",):"#494E6B",("ctp",):"#494E6B"}
    stored_stats = []
    groups_list = []
    feature_list = []
    qseqid_list = []
    prev_list = []
    figure_data = dict()
    x_pos_dict=dict()
    ######################## Edit ################################
    keep_japan = {}#{"country":["JPN"]}
    ##############################################################
    for qseqid, one_res in dmnd_res.items():
        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            f.write("\nDemographic_Factor\tValue_1\tValue_2\tPrevalence_Value_1\tPrevalence_Value_1_%\tPrevalence_Value_2\tPrevalence_Value_2_%\tp-value\tBH-adj._p-value\n")
    for feature in global_features:
        all_res, dmnd_res = get_metadata_count(feature, all_not_none=True, keep_only= keep_japan, multi=combine_qseqids)
        #print(dmnd_res)
        for qseqid, one_res in dmnd_res.items():
            combined_data_all = defaultdict(int)
            combined_data_dmnd = defaultdict(int)
            for res, combined_data in [(all_res, combined_data_all), (one_res, combined_data_dmnd)]:
                for study, data in res.items():
                    for columns, count in data.items():
                        if feature == ("BMI",):
                            bmi = columns[0]
                            if bmi < 18.5: combined_data[("Underweight\n(<18.5)",)] += count
                            elif bmi < 25: combined_data[("Normal\n(18.5 - <25)",)] += count
                            elif bmi < 30: combined_data[("Overweight\n(25 - <30)",)] += count
                            else: combined_data[("Obese\n(≥30)",)] += count
                        elif feature == ("ctp",):
                            score = columns[0]
                            if score <= 6: combined_data[("A (5-6 Points)",)] += count
                            elif score <= 9: combined_data[("B (7-9 Points)",)] += count
                            elif score <= 15: combined_data[("C (10-15 Points)",)] += count
                        elif feature == ("age",):
                            age = columns[0]
                            decade, _ = divmod(age, 10)
                            combined_data[(str(int(decade)),)] += count
                        else:
                            combined_data[columns] += count

            fig_data = []
            for columns, all_count in combined_data_all.items():
                if all_count < 10: continue   #Minimum number of counts per column
                dmnd_count = combined_data_dmnd[columns]
                col_perc = (dmnd_count / all_count) * 100
                fig_data.append((columns, col_perc, (dmnd_count, all_count)))
            if feature == ("age",): fig_data=sorted(fig_data, key=lambda x: (x[0]))
            elif feature == ("BMI",): fig_data=sorted(fig_data, key=lambda x: ["Underweight\n(<18.5)","Normal\n(18.5 - <25)","Overweight\n(25 - <30)","Obese\n(≥30)"].index(x[0][0]))
            elif feature == ("ajcc",): fig_data=sorted(fig_data, key=lambda x: ["0","i","ii","iii","iv"].index(x[0][0]))
            elif feature == ("smoker","ever_smoker",): fig_data=sorted(fig_data, key=lambda x: [("no","no"),("no","yes"),("yes","yes")].index((x[0][0],x[0][1])))
            else: fig_data=sorted(fig_data, key=lambda x: (-x[1],x[0]))
            #print(fig_data)
            x_vals = ["\n+ ".join(tuple(str(obj) for obj in x)) for x,y,z in fig_data]
            x_pos_dict[(qseqid, feature)]={i:x_val for x_val,i in enumerate(x_vals)}
            y_vals = [y for x,y,z in fig_data]
            fischer_vals = [[z[0], z[1]-z[0]] for x,y,z in fig_data]
            bar_labels = [f'{z[0]}/{z[1]}' for x,y,z in fig_data]
            figure_data[(qseqid, feature)]=(x_vals, y_vals, bar_labels)
            for (cols1, fisch1), (cols2, fisch2) in pairs(list(zip(x_vals, fischer_vals))):
                stat = scipy.stats.fisher_exact(np.array([fisch1,fisch2]))[1]
                feature_list.append(feature)
                groups_list.append((cols1, cols2))
                prev_list.append([fisch1,fisch2])
                stored_stats.append(stat)
                qseqid_list.append(qseqid)

    corr_p_val_dict=defaultdict(list)
    for qseqid, feature, group, orig, corr in zip(qseqid_list, feature_list, groups_list, stored_stats, list(multipletests(stored_stats, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1])):
        if corr < 0.05:
            corr_p_val_dict[(qseqid, feature)].append([(x_pos_dict[(qseqid, feature)][group[0]],x_pos_dict[(qseqid, feature)][group[1]]),corr])

    previous_qseqid=""
    previous_feature=""

    for qseqid, prev, feature, group, orig, corr in zip(qseqid_list, prev_list, feature_list, groups_list, stored_stats, list(multipletests(stored_stats, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1])):

        if qseqid != previous_qseqid or feature != previous_feature:
            #print(qseqid)
            previous_qseqid = qseqid
            previous_feature = feature
            x_vals, y_vals, bar_labels = figure_data[(qseqid, feature)]
            plt.figure(figsize=(len(x_vals), 6))
            fig = sns.barplot(x=x_vals, y=y_vals, color=color_map[feature], width=0.8) ### leaving out bw_adjust leads to a very smooth curve but also negative gene lengths
            fig.bar_label(fig.containers[0], labels = bar_labels)
            join_cols = " + ".join(" ".join(y[0].upper()+y[1:] for y in x.split("_")) for x in feature)
            if join_cols=="Age": join_cols="Decade of Life"
            if join_cols=="Ctp": join_cols="CTP Class"
            if join_cols=="Antibiotics Current Use": join_cols="Current Antibiotic Use"
            fig.set(xlabel=f"{join_cols}", ylabel="Prevalence (%)", title=f"Prevalence of {qseqid} by {join_cols}") #Prevalence of Methyl-Coenzyme M Reductase (Alpha Subunit) MCR (αβ Subunit)   Japanese Prevalence of {qseqid} by {join_cols}         MCR (α, β, and γ  Subunits)
            #plt.xlim(-100,5100)
            if any(len(x)>8 for x in x_vals):
                plt.xticks(rotation=90)
            if len(x_vals)<6:
                #print(corr_p_val_dict[(qseqid,feature)])
                add_sig(fig.get_ylim(),corr_p_val_dict[(qseqid,feature)])
            #plt.savefig(f"./figures/archaea/Prevalence of MCR (α, β, and γ  Subunits) by {join_cols}.svg")
            os.makedirs(f"./{save_folder}/{qseqid}/other_demographic_factors/", exist_ok=True)
            plt.savefig(f"./{save_folder}/{qseqid}/other_demographic_factors/Prevalence_of_{qseqid}_by_{join_cols}.svg")
            #plt.show()
            plt.close()  
        #if corr < 1.05: print('{:>30} {:>30} {:>30} {:>30}'.format(*(str(feature), str(group),orig,corr)))
        [[c1_yes, c1_no],[c2_yes,c2_no]] = prev
        new_line="\n"
        #if corr < 0.05: print ("  ", f'{group[0]}\t{group[1]}\t{c1_yes}/{c1_yes+c1_no}\t{c1_yes/(c1_yes+c1_no)}\t{c2_yes}/{c2_yes+c2_no}\t{c2_yes/(c2_yes+c2_no)}\t{orig}\t{corr}')   
        #print(f'{"+".join(feature)}\t{group[0].replace(new_line,"")}\t{group[1].replace(new_line,"")}\t{c1_yes}/{c1_yes+c1_no}\t{100*c1_yes/(c1_yes+c1_no)}\t{c2_yes}/{c2_yes+c2_no}\t{100*c2_yes/(c2_yes+c2_no)}\t{orig}\t{corr}\n')
        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            f.write(f'{"+".join(feature)}\t{group[0].replace(new_line,"")}\t{group[1].replace(new_line,"")}\t{c1_yes}/{c1_yes+c1_no}\t{100*c1_yes/(c1_yes+c1_no)}\t{c2_yes}/{c2_yes+c2_no}\t{100*c2_yes/(c2_yes+c2_no)}\t{orig}\t{corr}\n')
    #print('{:>30} {:>30} {:>30} {:>30}'.format(*("Metadata Feature", "Groups Compared","Original p-Value","BH corrected p-Value"))+"\n")


    # Read MAG Abundance
    @disk_cache
    def read_leviatan_abundance():
        return pd.read_excel('/DATA/Matt/mastersproject/Leviatan_species_abundance.xlsx', index_col=0) 

    df = read_leviatan_abundance()

    # Get DMND MAGs
    #print(sum(v for d in all_res.values() for v in d.values()))        # total count of specified columns in all metadata

    #lev_files = []
    lev_files_dict = defaultdict(list)
    multi_lev=combine_qseqids
    with open(dmnd_tsv, 'r') as file:
        seen_lev = set()
        one_res = defaultdict(dict)
        for line in file:
            qseqid, sseqid, _, _, _, _, _, _, _, _, evalue, bitscore, *_ = line.rstrip().split("\t")
            evalue = float(evalue)
            bitscore = float(bitscore)
            fasta_index, file, org_group, codon_table, tools = sseqid.split("|")
            dataset, individual = file.split("__")
            if dataset == "LeviatanS_2022":
                if (dataset, individual, qseqid) not in seen_lev:
                    if evalue < 1e-3: 
                        seen_lev.add((dataset, individual, qseqid))
                        if multi_lev and (not all((dataset, individual, one_qseqid) in seen_lev for one_qseqid in all_qseqids)): continue
                        if multi_lev: qseqid = " + ".join(all_qseqids)
                        lev_files_dict[qseqid].append(individual)
        #print(sum(v for d in one_res.values() for v in d.values()))   # total count of specified columns in dmnd result metadata
    for k,v in lev_files_dict.items():
        lev_files_dict[k] = list(set(v))
    #print(lev_files_dict)

    # Get Taxa
    

    tax_table_pairs = []
    i=0
    ncbi = NCBITaxa()
    # gtdbtk.bac120.summary.tsv      gtdbtk.ar53.summary.tsv
    lev_metadata = dict()
    to_rank = {"d__":"domain","k__":"kingdom","p__":"phylum","c__":"class","o__":"order","f__":"family","g__":"genus","s__":"species"}
    lev_metadata_files = ["/home/matt/DATA/metadata/gtdbtk.bac120.summary.tsv", "/home/matt/DATA/metadata/gtdbtk.ar53.summary.tsv"]
    for lev_metadata_file in lev_metadata_files:
        with open(lev_metadata_file, 'r') as f:
            for line in f:
                m = re.match(r"(Rep_\d+)\t([^\t]+)\t", line)
                if m: 
                    lev_metadata[m[1]]=dict()
                    tax_names = m[2].split(";")
                    #tax_names = [re.sub(r"_\w","", re.sub(r"\w__", "", x)) for x in tax_names[::-1]]
                    filtered_tax_names_rank = []
                    filtered_tax_names = []
                    for tax_name in tax_names:
                        mat = re.match("(\w__)(.*)", tax_name)
                        if mat:
                            filtered_tax_names_rank.append(to_rank[mat[1]])
                            filtered_tax_names.append(re.sub(r"_\w+","",mat[2]))


                    for tax_rank,tax_name in zip(filtered_tax_names_rank,filtered_tax_names):

                        taxid = ncbi.get_name_translator([tax_name]).get(tax_name)
                        lev_metadata[m[1]][tax_rank]=(tax_name, taxid)

    # Write MAG Species

    stats_all, stats_ind = get_metadata_count(("study_name","sample_id"))
    total_mags = sum(sum(d.values()) for d in stats_all.values() )
    def total_qseqid(id):
        if multi_lev:
            all_sets=[]
            for qseqid,study_dict in stats_ind.items():
                all_sets.append(set())
                for k,v in study_dict.items():
                    for k2,v2 in v.items():
                        if v2>1: raise Exception("Error!")
                        if v2: all_sets[-1].add(k2)
            intersection_set = all_sets[0]
            for s in all_sets[1:]:
                intersection_set &= s
            return len(intersection_set) ############################### CHECK IF CORRECT
        return sum(sum(d.values()) for d in stats_ind[id].values() )

    c=[]#########
    map_rep = dict()


    for qseqid, _ in dmnd_res.items():#for qseqid in all_qseqids:
        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            if lev_files_dict.get(qseqid):
                lev_files=lev_files_dict.get(qseqid)
            
                #print(qseqid)
                f.write(f"\nLeviatan species genomes:\n")        ##{total_qseqid(qseqid)}/{total_mags}     1629/9634
                for rep in lev_files:
                    res = []
                    for i, tax_rank in enumerate(["species","genus","family","order","class","phylum","domain"]):
                        tax_name, taxid = lev_metadata[rep][tax_rank]
                        if " " in tax_name: c.append(tax_name) #######
                        if tax_name and not map_rep.get(rep): map_rep[rep]=tax_name 
                        if taxid:
                            res.append(" "*(2*i+2)+f'{tax_rank}: {tax_name} ({taxid[0]})')  ##also include leviatan Rep number with rep?
                            break
                        else:
                            res.append(" "*(2*i+2)+f'{tax_rank} not in database: {tax_name if tax_name else "N/A"}')
                        
                    f.write("\n".join(res)+"\n")
            else:
                f.write(f"\nLeviatan species genomes:\n\tN/A\n")


    # Write MAG Taxa
    for qseqid, _ in dmnd_res.items():#for qseqid in all_qseqids:
        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            if lev_files_dict.get(qseqid):
                lev_files=lev_files_dict.get(qseqid)
                f.write("\nTaxonomic rank abundance:\n")
                for rank in ["phylum","family","genus"]:
                    pos_rank_dict = defaultdict(int)
                    all_rank_dict = defaultdict(int)

                    for rep in lev_files:
                        tax_name, taxid = lev_metadata[rep][rank]
                        if tax_name:
                            pos_rank_dict[tax_name]+=1
                    for rep in df.columns.to_list():
                        tax = lev_metadata[rep][rank][0]
                        if tax in pos_rank_dict:
                            all_rank_dict[tax]+=1
                    f.write("\t"+rank+"\n")
                    for rep, count in pos_rank_dict.items():    
                        f.write(f'\t\t{rep}: {count}/{all_rank_dict[rep]} ({count/all_rank_dict[rep]*100 :.2f} %)\n')
            else:
                f.write(f"\nTaxonomic rank abundance:\n\tN/A\n")


    lev_violin = dict()

    for qseqid, lev_files in lev_files_dict.items():
        #lev_files.remove('Rep_1238')
        #print(lev_files)######################################
        fpf = df[lev_files].sum(axis=1)
        fpf_nonzero = fpf[fpf != 0]*100

        fpf_il = [value for label, value in fpf.items() if label.startswith('IL')]
        fpf_nonzero_il = [value for label, value in fpf_nonzero.items() if label.startswith('IL')]

        fpf_nld = [value for label, value in fpf.items() if label.startswith('NLD')]
        fpf_nonzero_nld = [value for label, value in fpf_nonzero.items() if label.startswith('NLD')]

        lev_violin[qseqid] = [fpf, fpf_nonzero, fpf_il, fpf_nonzero_il, fpf_nld, fpf_nonzero_nld]

    # Write MAG inclusion
    total_il = df.index.str.startswith('IL').sum()
    total_nld = df.index.str.startswith('NLD').sum()
    for qseqid, lev_files in lev_files_dict.items():
        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            f.write("Total MAG inclusion\n")
            for rep_num,rep_non_zero_count in df[lev_files].count().items():
                f.write(f'\t{map_rep[rep_num]}\t({rep_num})\t{rep_non_zero_count}/{total_il+total_nld}\n')

            f.write("Israel MAG inclusion\n")
            for rep_num,rep_non_zero_count in df[lev_files][df[lev_files].index.str.startswith('IL')].count().items():
                f.write(f'\t{map_rep[rep_num]}\t({rep_num})\t{rep_non_zero_count}/{total_il}\n')

            f.write("Netherlands MAG inclusion\n")
            for rep_num,rep_non_zero_count in df[lev_files][df[lev_files].index.str.startswith('NLD')].count().items():
                f.write(f'\t{map_rep[rep_num]}\t({rep_num})\t{rep_non_zero_count}/{total_nld}\n')

    # Write FPF


    for qseqid, [fpf, fpf_nonzero, fpf_il, fpf_nonzero_il, fpf_nld, fpf_nonzero_nld] in lev_violin.items():
        with open(f"./{save_folder}/{qseqid}/overview.txt", 'a+') as f:
            u_statistic, p_value = mannwhitneyu(fpf_nonzero_il, fpf_nonzero_nld)
            f.write(f'Function Positive Fraction:\n')
            f.write(f'\tMean: Israel:{mean(fpf_il)}\tIsrael-Nonzero: {mean(fpf_nonzero_il)}\tNeatherlands: {mean(fpf_nld)}\t:Neatherlands-Nonzero: {mean(fpf_nonzero_nld)}\n')
            f.write(f'\tGeometric Mean: Israel:{gmean(fpf_il)}\tIsrael-Nonzero: {gmean(fpf_nonzero_il)}\tNeatherlands: {gmean(fpf_nld)}\t:Neatherlands-Nonzero: {gmean(fpf_nonzero_nld)}\n')
            f.write(f'\tMedian: Israel:{median(fpf_il)}\tIsrael-Nonzero: {median(fpf_nonzero_il)}\tNeatherlands: {median(fpf_nld)}\t:Neatherlands-Nonzero: {median(fpf_nonzero_nld)}\n')
            f.write(f'Israel vs Netherlands:\n\tU-statistic: {u_statistic}\n\tp-value: {p_value}\n')
            data = [fpf_nonzero, fpf_nonzero_il, fpf_nonzero_nld]

            log_data = [[np.log10(d) for d in row] for row in data]
            fig, ax = plt.subplots(figsize=(14, 8))#9,6

            
            
            #this is geometric mean
            sns.boxplot(data=log_data, width=0.25, 
                    showfliers=False, showmeans=True, 
                    meanprops=dict(marker='o', markerfacecolor='white', markeredgecolor="white",
                                    markersize=5, zorder=3),
                    boxprops=dict(facecolor=(0,0,0,0), 
                                    linewidth=0, zorder=3),
                    whiskerprops=dict(linewidth=0, zorder=3),
                    capprops=dict(linewidth=0, zorder=3),
                    medianprops=dict(linewidth=2, zorder=3, color="white"))
            sns.violinplot(data=log_data, zorder=1, ax=ax, inner=None, color='lightgray', density_norm='count')

            #sns.swarmplot(data=log_data, zorder=2, s=3, ax=ax, palette='dark:black')

            
            sns.swarmplot(data=log_data, zorder=2, s=1.5, ax=ax, palette='dark:black')


            ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
            ymin, ymax = ax.get_ylim()
            tick_range = np.arange(np.floor(ymin), ymax)
            ax.yaxis.set_ticks(tick_range)
            ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)

            # Set x-axis labels
            labels = [f"All Samples\n{len(fpf_nonzero)}/{len(fpf)}", f"Israel\n{len(fpf_nonzero_il)}/{len(fpf_il)}", f"Netherlands\n{len(fpf_nonzero_nld)}/{len(fpf_nld)}"]
            ax.set_xticks([0,1,2])
            ax.set_xticklabels(labels)

            plt.xlabel('Samples')
            plt.ylabel('Function Positive Fraction (%)')
            plt.title(f'Function Positive Fraction of {qseqid}')#of Methyl-Coenzyme M Reductase (Alpha Subunit)    Function Positive Fraction of MCR (α, β, and γ  Subunits) Combined Violin Plot and Swarm Plot for {qseqid
            plt.show()
            #plt.savefig(f"./figures/Function Positive Fraction of {qseqid}.svg")
            os.makedirs(f"./{save_folder}/{qseqid}/FPF/", exist_ok=True)
            plt.savefig(f"./{save_folder}/{qseqid}/FPF/Function Positive Fraction of {qseqid}.svg")
            plt.tight_layout()

if __name__ == "__main__":
    main()