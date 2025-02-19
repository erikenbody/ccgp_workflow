#!/bin/python3
import gzip
import shutil
import argparse
from pathlib import Path
from typing import TextIO
import sys
# User provides list of sample names, 1 per line
# User provides path to where fastq files are. Assume paired end and that file name has sample name in it
# User provides path to reference genome. Need to copy this to proper path

def read_sample_list(sample_fh: TextIO) -> list:
    return sample_fh.read().splitlines()

def find_sample_fastqs(samples: list, fastq_dir: Path) -> dict:
    """Searches fastq_dir for sample names and associates in a dict"""
    sample_fastq_paths = {}
    cant_find = []
    for samp_name in samples:
        fqs = sorted(list(fastq_dir.glob(f"*{samp_name}*")))  # Hoping that sorting will make fq1 first. 
        if len(fqs) != 2:
            cant_find.append(samp_name)
        else:
            sample_fastq_paths[samp_name] = fqs
    return sample_fastq_paths, cant_find

def copy_reference(ref: Path) -> str:
    exts = ['.fna', '.fa', '.fasta']
    for ext in exts:
        if ext in ref.name:
            ref_name = ref.name.split(ext)[0]
    if Path('..', 'data', 'genome', ref_name + ".fna").exists():
        return ref_name
    if not Path("../data/genome").exists():
        try:
            Path("../data/genome").mkdir(parents=True)
        except PermissionError as e:
            print(e)
            sys.exit("Must run script from workflow directory.")
    if ref.suffix == ".gz":
        with gzip.open(ref, 'rb') as f_in:
            with open(Path('..', 'data', 'genome', ref_name + ".fna"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copyfile(ref, Path('..', 'data', 'genome', ref_name + ".fna"))
    return ref_name

def write_sample_sheet(sample_dict: dict, ref_name: str, organism: str, ) -> None:
    """Writes the sample sheet"""
    org = organism.replace(" ", "_")
    with open(Path("../sample_sheets", f"{org}_samples.csv"), "w") as out:
        out.write("BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2\n")
        for i, (k, v) in enumerate(sample_dict.items()):
            out.write(f"{k},lib_{k},{ref_name},{i},{organism},NaN,{v[0]},{v[1]}\n")
            

def main() -> None:
    
    parser = argparse.ArgumentParser(description='Write sample files.')
    parser.add_argument('-s', '--sample_list', dest='samp', required=True, help="Specify path to sample list")
    parser.add_argument('-f', '--fastq_dir', dest='fastq', required=True, help="Specify path to fastq dir")
    parser.add_argument('-r', '--ref', dest='ref', required=True, help="Specify path to reference genome")
    parser.add_argument('-o', '--org', dest='org', required=True, help="Specify organism name")
    args = parser.parse_args()
    
    sample_list = args.samp
    fastq_dir = Path(args.fastq)
    ref = Path(args.ref)
    organism = args.org
    
    with open(sample_list, "r") as f:
        samples = read_sample_list(f)
        
    sample_dict, cant_find = find_sample_fastqs(samples, fastq_dir)
    ref_name = copy_reference(ref)
    write_sample_sheet(sample_dict, ref_name, organism)
    if cant_find:
        print("Couldnt' find fastqs for these files:")
        for name in cant_find:
            print(name)

main()