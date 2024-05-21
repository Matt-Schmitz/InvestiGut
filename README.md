<img src="https://github.com/Matt-Schmitz/InvestiGut/assets/34464190/ce6ea7aa-7ac1-4ee2-9e3a-c2de169f531e" alt="investigut" width="555" height="128"/>

## Instructions
1. Clone the repository.  
```bash
git clone https://github.com/Matt-Schmitz/InvestiGut.git
```
3. Enter the InvestiGut folder.  
```bash
cd InvestiGut
```
4. Create the environment with [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).  
```bash
mamba create --no-channel-priority -n investigut -c bioconda -c conda-forge "python=3.11" "numpy=1.24.3" "scipy=1.10.1" "matplotlib=3.7.1" "seaborn=0.13.0" "pandas=1.5.3" "statsmodels=0.13.5" "ete3=3.1.2" "openpyxl=3.0.10" "bioconda::diamond=2.1.8"
```
5. Activate the environment.  
```bash
conda activate investigut
```
6. Download proteins and create DIAMOND database.  After the DIAMOND database has been created, the compressed fasta data can be deleted.
$${\color{red}WARNING: \space The \space compressed \space data \space take \space up \space 79GB\space (two \space files). The \space DIAMOND \space database \space will \space be \space 311GB.}$$
```bash
./download_fasta.sh
```
7. Add the InvestiGut folder to your `~/.bashrc` and apply the changes.  
- Add to .bashrc: `export PATH="/PATH/TO/InvestiGut:$PATH"`
- Enter in terminal:
```bash
source ~/.bashrc
```
- Reactivate the mamba environment:
```bash
mamba activate investigut
```

8. Run InvestiGut on chosen proteins.  
Basic usage:

```bash
investigut.py -i /PATH/TO/FASTA/FILE
```
  
Options list:

Any arguments not recognized by InvesitGut will be passed on to DIAMOND (e.g., `--threads`, `--faster`, `--ultra-sensitive`).
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
