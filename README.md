# InvestiGut

## Instructions
1. Clone the repository.  
```bash
git clone https://github.com/Matt-Schmitz/InvestiGut.git
```
3. Enter the InvestiGut folder.  
```bash
cd InvestiGut
```
4. Create the conda environment.  
```bash
conda env create -f environment.yml -n investigut
```
5. Activate the environment.  
```bash
conda activate investigut
```
6. Download proteins and create DIAMOND database. The compressed data take up 79GB. The DIAMOND database will be 311GB. After the DIAMOND database has been created, the compressed fasta data can be deleted.  
```bash
./download_fasta.sh
```
8. Add the InvestiGut folder to your `~/.bashrc` and apply the changes.  
- Add to .bashrc: `export PATH="/PATH/TO/InvestiGut:$PATH"`
- Enter in terminal:
```bash
source ~/.bashrc
```
- Reactivate the conda environment:
```bash
conda activate investigut
```

9. Run InvestiGut on chosen proteins.  
Basic usage:

```bash
investigut.py -i /PATH/TO/FASTA/FILE
```
  
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
