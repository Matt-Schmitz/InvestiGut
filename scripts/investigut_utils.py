# Caching
import pickle
import os
import hashlib
import inspect

class DiskCache:
    def __init__(self, cache_message:str="Results saved", load_message:str="Results loaded", cache_folder:str="", filename_override:str="", include_func_in_hash:bool=False):
        self.cache_message = cache_message
        self.load_message = load_message
        self.cache_folder = cache_folder if cache_folder else os.getcwd()
        self.filename_override = filename_override
        self.include_func_in_hash = include_func_in_hash

        os.makedirs(self.cache_folder, exist_ok=True)

        if self.filename_override and ("." not in self.filename_override):
            self.filename_override += ".pkl"  

    def __call__(self, *args, cache_message:str="", load_message:str="", cache_folder:str="", filename_override:str="", include_func_in_hash:bool=""):
        os.makedirs(cache_folder or self.cache_folder, exist_ok=True)
        
        if filename_override and ("." not in filename_override):
            filename_override += ".pkl"
        
        def decorator(func):
            def wrapper(*args, **kwargs):
                # Generate a unique cache file name based on the function name and its arguments
                to_hash = (args, kwargs) if ((include_func_in_hash==False) or ((include_func_in_hash != True) and (self.include_func_in_hash==False))) else (args, kwargs, inspect.getsource(func).split("\n")[1:]) 
                arg_hash = hashlib.sha256(pickle.dumps(to_hash)).hexdigest()[:32]
                cache_key = filename_override or self.filename_override or f"{func.__name__}_{arg_hash}.pkl"
                cache_file = os.path.join((cache_folder or self.cache_folder), cache_key)

                if os.path.exists(cache_file):
                    # Load cached results from file
                    with open(cache_file, "rb") as f:
                        cached_results = pickle.load(f)
                    if load_message or self.load_message: 
                        print(load_message or self.load_message)
                    return cached_results
                else:
                    # Execute the function and store the results
                    results = func(*args, **kwargs)
                    with open(cache_file, "wb") as f:
                        pickle.dump(results, f)
                    if cache_message or self.cache_message:
                        print(cache_message or self.cache_message)
                    return results

            return wrapper
        if len(args) == 1 and callable(args[0]):
            return decorator(args[0])
        else:
            return decorator

    def set_load_message(self, message="Results loaded"):
        self.load_message = message

    def set_cache_message(self, message="Results saved"):
        self.cache_message = message

    def set_cache_folder(self, cache_folder):
        self.cache_folder = cache_folder if cache_folder else os.getcwd()
        
        if not os.path.exists(self.cache_folder):
            os.makedirs(self.cache_folder)

disk_cache = DiskCache()

# File I/O
import gzip
import bz2
import lzma
import warnings
from os.path import exists

def zopen(file_path, mode='', write_format=''):

    if 'r' in mode or ('a' in mode and exists(file_path)):
        if 'r' in mode and write_format:
            warnings.warn(f"write_format specified in read mode.", RuntimeWarning)

        with open(file_path, 'rb') as file:
            signature = file.read(6)

            if signature[:2] == b'\x1f\x8b':  # gz
                return gzip.open(file_path, mode)
            elif signature[:3] == b'BZh':  # bz2
                return bz2.open(file_path, mode)
            elif signature[:6] == b'\xfd7zXZ\x00':  # xz
                return lzma.open(file_path)
            else:
                return open(file_path, mode)
    elif 'w' in mode or 'x' in mode or ('a' in mode and not exists(file_path)):
        if (ending := file_path[-len(write_format):]) != write_format:
            warnings.warn(f"write_format = '{write_format}', but the path ends with '{ending}'.", RuntimeWarning)

        if write_format == 'gz':
            return gzip.open(file_path, mode)
        elif write_format == 'bz2':
            return bz2.open(file_path, mode)
        elif write_format == 'xz':
            return lzma.open(file_path, mode)
        elif write_format == '':
            return open(file_path, mode)
        else:
            raise ValueError("Unsupported write_format. Use 'gz', 'bz2', 'xz', or leave empty.")
    else:
        raise ValueError("Unsupported mode. Use 'r' for reading, 'w' for writing, 'x' for exclusive creation, or 'a' for appending.")
    


#Translation
from itertools import product
import re

codon_table ={
1:  ["Standard", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "---M------**--*----M---------------M---------------------------- "],
2:  ["Vertebrate Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG ", "----------**--------------------MMMM----------**---M------------ "],
3:  ["Yeast Mitochondrial", "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "----------**----------------------MM---------------M------------ "],
4:  ["Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "--MM------**-------M------------MMMM---------------M------------ "],
5:  ["Invertebrate Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG ", "---M------**--------------------MMMM---------------M------------ "],
6:  ["Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear", "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "--------------*--------------------M---------------------------- "],
9:  ["Echinoderm Mitochondrial; Flatworm Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG ", "----------**-----------------------M---------------M------------ "],
10: ["Euplotid Nuclear", "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "----------**-----------------------M---------------------------- "],
11: ["Bacterial, Archaeal and Plant Plastid", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "---M------**--*----M------------MMMM---------------M------------ "],
12: ["Alternative Yeast Nuclear", "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "----------**--*----M---------------M---------------------------- "],
13: ["Ascidian Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG ", "---M------**----------------------MM---------------M------------ "],
14: ["Alternative Flatworm Mitochondrial", "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG ", "-----------*-----------------------M---------------------------- "],
15: ["Blepharisma Macronuclear", "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "----------*---*--------------------M---------------------------- "],
16: ["Chlorophycean Mitochondrial", "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "----------*---*--------------------M---------------------------- "],
21: ["Trematode Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG ", "----------**-----------------------M---------------M------------ "],
22: ["Scenedesmus obliquus mitochondrial", "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "------*---*---*--------------------M---------------------------- "],
23: ["Thraustochytrium mitochondrial code", "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "--*-------**--*-----------------M--M---------------M------------ "],
24: ["Rhabdopleuridae Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M------**-------M---------------M---------------M------------"],
25: ["Candidate Division SR1 and Gracilibacteria", "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "---M------**-----------------------M---------------M------------ "],
26: ["Pachysolen tannophilus Nuclear", "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "----------**--*----M---------------M---------------------------- "],
27: ["Karyorelict Nuclear", "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "--------------*--------------------M---------------------------- "],
28: ["Condylostoma Nuclear", "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "----------**--*--------------------M---------------------------- "],
29: ["Mesodinium Nuclear", "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "--------------*--------------------M---------------------------- "],
30: ["Peritrich Nuclear", "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "--------------*--------------------M---------------------------- "],
31: ["Blastocrithidia Nuclear", "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "----------**-----------------------M---------------------------- "],
32: ["Balanophoraceae Plastid", "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ", "---M------*---*----M------------MMMM---------------M------------ "],
33: ["Cephalodiscidae Mitochondrial", "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG ", "---M-------*-------M---------------M---------------M------------ "]}

codons = ["".join(x) for x in product("TCAG", repeat=3)]

translation_dicts = dict()
start_stop_dicts = dict()

for table, [organisms, AAs, starts_stops] in codon_table.items():
    translation_dicts[table] = {codon: aa for codon, aa in zip(codons, AAs)}
    start_stop_dicts[table] = {codon: start_stop for codon, start_stop in zip(codons, starts_stops)}

def translation(seq, table):
    if not re.match("([TCAG]{3})+$", seq):
        raise Exception(f"Error processing sequence ({seq}). The input sequence should have a length that is a multiple of 3 and only contain T, C, A, and G.")
    translation_dict = translation_dicts[table]
    start_stop_dict = start_stop_dicts[table]
    wrong_start = False 
    wrong_stop = False 
    if start_stop_dict[seq[:3]] != "M":
        wrong_start = True
    if start_stop_dict[seq[-3:]] != "*":
        wrong_stop = True
    return (f'M{re.sub("...", lambda m: translation_dict[m.group()], seq[3:-3])}', wrong_start, wrong_stop)