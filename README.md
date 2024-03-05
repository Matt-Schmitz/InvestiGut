# InvestiGut

## Instructions
1. Clone the repository.  
`git clone https://github.com/Matt-Schmitz/InvestiGut.git`
2. Enter the InvestiGut folder.  
`cd InvestiGut`
3. Create the conda environment.  
`conda env create -f environment.yml -n investigut`
4. Activate the environment.  
`conda activate investigut`
5. Download proteins and create DIAMOND database. The compressed data take up 79GB. The DIAMOND database will be 311GB. After the DIAMOND database has been created, the compressed fasta data can be deleted.  
`./download_fasta.sh`
6. Add the InvestiGut folder to your `~/.bashrc` and apply the changes.  
`export PATH="/PATH/TO/InvestiGut:$PATH"`  
`source ~/.bashrc`
8. Run InvestiGut on chosen proteins.  
Basic usage: `investigut.py -i /PATH/TO/FASTA/FILE`
  
Options list:
```
-h, --help            show this help message and exit
-i {INPUT}            input fasta file
-o {OUTPUT}           output folder
-s                    single mode (default)
-m                    multi mode
-d {DMND_TSV_OUTPUT}  optional DIAMOND .tsv output file to skip alignment
--query-cover {NUM}   DIAMOND query-coverage (default 90.0)
--subject-cover {NUM}
                      DIAMOND subject-coverage (default 90.0)
--id {NUM}            DIAMOND %id (default 90)
--low                 overrides default DIAMOND query/sujbect coverage and %id giving all a
                      value of 50
```
