# Caching
import pickle
import os
import hashlib
import inspect

class DiskCache:
    def __init__(self, cache_message:str="Results cached", load_message:str="Results loaded from cache", cache_folder:str="", filename_override:str="", include_func_in_hash:bool=False):
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

    def set_load_message(self, message="Results loaded from cache"):
        self.load_message = message

    def set_cache_message(self, message="Results cached"):
        self.cache_message = message

    def set_cache_folder(self, cache_folder):
        self.cache_folder = cache_folder if cache_folder else os.getcwd()
        
        if not os.path.exists(self.cache_folder):
            os.makedirs(self.cache_folder)

disk_cache = DiskCache()


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