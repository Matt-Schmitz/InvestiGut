#!/usr/bin/env python3

import re
import subprocess
import sys
import json
from os.path import exists, dirname, getsize

fasta_file = sys.argv[1]
dir_base = sys.argv[2]
pathToGeneMarkS = sys.argv[3]
pathToGeneMarkS2 = sys.argv[4]
pathToGeneMarkST = sys.argv[5]
pathToMetaGeneMark = sys.argv[6]
pathToMetaGeneMark2 = sys.argv[7]
sorted_mag_metagenomes = sys.argv[8].split()
tools_for_org_group = json.loads(sys.argv[9])

mgm_mod = dirname(pathToMetaGeneMark) + "/MetaGeneMark_v1.mod"

print("!!!!!!!!!!!!!!!!!!")
print(fasta_file)
print(dir_base)
print(sorted_mag_metagenomes)
print(type(sorted_mag_metagenomes))
print("!!!!!!!!!!!!!!!!!!")



tools_commands = {
                  "GlimmerHMM": "glimmerhmm {gen_file} "+f"{dir_base}/gene_prediction/TrainGlimmM_S_cerevisiae"+" -o {pred_file} -g",
                  "Phanotate": "phanotate.py {gen_file}",
                  "Augustus": "augustus --gff3=on --species=saccharomyces_cerevisiae_S288C --stopCodonExcludedFromCDS=false {gen_file}",
                  "Snap": "snap "+f"{dir_base}/gene_prediction/SNAP/Saccharomyces_cerevisiae.hmm"+" {gen_file} -gff",
                  "FragGeneScan": "FragGeneScan -s {gen_file} -o {pred_file_no_ext} -w 0 -t complete",
                  "Pyrodigal": "pyrodigal -i {gen_file} -o {pred_file} -g {gen_code} -p {mode} -f gff --min-gene 21 --min-edge-gene 21 --max-overlap 20",
                  "GeneMarkS": pathToGeneMarkS + " --format GFF3 --output {pred_file} --gcode {gen_code} {tool_org} {gen_file}",
                  "GeneMarkS2": pathToGeneMarkS2 + " --genome-type auto --format gff --output {pred_file} --gcode {gen_code} --seq {gen_file}",
                  "GeneMarkST": pathToGeneMarkST + " --format GFF --output {pred_file} --gcode {gen_code} --prok {gen_file}", # not specifying prok doesn't work???
                  "MetaGeneAnnotator": "mga {gen_file} -m",
                  "MetaGeneMark": pathToMetaGeneMark + " -m "+mgm_mod+" -f G -o {pred_file} {gen_file}",
                  "MetaGeneMark2": pathToMetaGeneMark2 + " --seq {gen_file} --out {pred_file} --format gff3",
                  "Prodigal": "prodigal -i {gen_file} -f gff -o {pred_file} -g {gen_code}"}  # -p meta
available_tables = {"GeneMarkS":[1,4,11,25], "GeneMarkST":[1,4,11], "GeneMarkS2":[4,11,15,25], "Prodigal":[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31,33], "Pyrodigal":[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31,33]}
#tools_for_org_group = {"Archaea":["FragGeneScan","GeneMarkST","MetaGeneMark"], "Bacteria":["MetaGeneAnnotator","MetaGeneMark","Pyrodigal"], "Eukaryota":["Augustus","Pyrodigal","Snap"], "Viruses":["MetaGeneAnnotator","GeneMarkS","GeneMarkS2"]}

#tools_for_org_group = {"Archaea":["Prodigal", "MetaGeneAnnotator"], "Bacteria":["Prodigal", "MetaGeneAnnotator"], "Eukaryota":["Prodigal", "MetaGeneAnnotator"], "Viruses":["Prodigal", "MetaGeneAnnotator", "GeneMarkS2"]}

stdout_output = ["MetaGeneAnnotator", "Phanotate", "Augustus", "Snap"]

def extract_filename(file_path):
    pattern = r'(.+?)(?:\.fna|\.fa|\.fasta)?(?:\.gz)?$'
    return re.match(pattern, file_path)[1]

for sorted_mag_metagenome in sorted_mag_metagenomes:
    org_group, codon_table, *_ = sorted_mag_metagenome.split("_")
    file = extract_filename(sorted_mag_metagenome)

    paramlist = [[tool, tools_commands[tool], sorted_mag_metagenome] for tool in tools_for_org_group["Bacteria" if org_group == "unknown" else org_group]]

    result = {}
    for sublist in paramlist:
        key = sublist[0]
        result.setdefault(key, []).append(sublist)

    def create_predictions(params):
        tool, pre_command, gen_file = params
        if tool == "FragGeneScan":
            frag_command = 'ln -s "$(dirname "$(which FragGeneScan)")/train" train'
            process_frag1 = subprocess.Popen(frag_command, shell=True)
            process_frag1.communicate()
            if process_frag1.returncode != 0:
                raise Exception(f"Command failed with return code {process_frag1.returncode}")
        pred_folder = "./"
        pred_file_no_ext = pred_folder + tool+"_"+ file
        pred_file = pred_folder + tool+"_"+ file + ".gff"
        gen_code = codon_table if tool in available_tables and int(codon_table) in available_tables[tool] else "11"
        mode = "meta" if tool == "Pyrodigal" and "metagenome" in fasta_file.lower() else  "single"
        org = {"GeneMarkS": {"Eukaryota":"--euk", "Bacteria":"--prok", "Viruses":"--virus", "Phages":"--phage"},
                "GeneMarkST":{"Bacteria":"--prok"}}
        o=org.get(tool)
        tool_org = (o.get(org_group) if o else "")
        d = {"{gen_file}":gen_file, "{pred_file}": pred_file, "{pred_file_no_ext}":pred_file_no_ext, "{gen_code}": gen_code, "{tool_org}": tool_org, "{mode}": mode}
        command = re.sub("{.*?}",lambda m: (d.get(m[0]) or m[0]), pre_command)

        #command = [(o.get(orga) if o else "") if c=="tool_org" else c for c in command]
        command = re.sub(" +"," ", command)
        
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell = False)
        stdout, stderr = process.communicate()
        if tool in stdout_output:
            with open(pred_file, 'w+') as f:
                if stdout:
                    f.write(stdout)
        if tool == "FragGeneScan":
            frag_command2 = [
                'awk',
                'BEGIN{print "##gff-version 3";} {s=substr($1,1,1); if (s==">") {seqid=substr($1,2);} else {s=split($0, t, "\\t"); id="ID=" seqid "_" t[1] "_" t[2] "_" t[3] ";product=predicted protein"; print seqid "\\tFGS\\tCDS\\t" t[1] "\\t" t[2] "\\t.\\t" t[3] "\\t" int(t[4]-1) "\\t" id; }}',
                f"{pred_file_no_ext}.out"
            ]
            output_file = f"{pred_file_no_ext}.gff"

            with open(output_file, "w") as output:
                process_frag2 = subprocess.Popen(frag_command2, stdout=output)
                stdout2, stderr2 = process_frag2.communicate()
            if process_frag2.returncode != 0:
                raise Exception(f"Command failed with return code {process_frag2.returncode}\n{' '.join(frag_command2)}\n{stderr2}")

        if not exists(pred_file)or (getsize(pred_file) == 0 and tool!="Snap"):# or getsize(pred_file) == 0: <- this causes SNAP to fail on files without proteins
            raise Exception(f"Error: {tool} failed to create a prediction file ({pred_file}) for {gen_file}.\n{command}\n{stderr}")        



    for p in paramlist:
        create_predictions(p)