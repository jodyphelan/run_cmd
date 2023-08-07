import pysam
import os
import joblib
from joblib import Parallel, delayed
import subprocess as sp
# from tqdm.rich import trange, tqdm
from tqdm import tqdm

__version__ = "0.0.1"

def sanitize_region(region):
    return region.replace(":","_").replace("-","_")
def genome_job(cmd,region):
    region_safe = sanitize_region(region)
    out = sp.run(cmd.format(region=region,region_safe=region_safe),shell=True,stderr=sp.PIPE,stdout=sp.PIPE)
    return out


def get_genome_chunks(fasta,nchunks):
    """Split genome into n chunks"""
    genome = pysam.FastaFile(fasta)
    total_length = sum(genome.lengths)
    chunk_length = int(total_length/nchunks)
    chunks = []
    chunk_start = 0
    chunk_end = 0
    for chrom in genome.references:
        while chunk_end < genome.get_reference_length(chrom):
            chunk_end += chunk_length
            if chunk_end > genome.get_reference_length(chrom):
                chunk_end = genome.get_reference_length(chrom)
            chunks.append([chrom,chunk_start,chunk_end])
            chunk_start = chunk_end
    regions = [f"{r[0]}:{r[1]}-{r[2]}" for r in chunks]
    return regions


def load_bed_regions(bed_file):
    regions = []
    with open(bed_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom,start,end = line.strip().split("\t")[:3]
            regions.append(f"{chrom}:{start}-{end}")
    return regions

def run_cmd_parallel_on_genome(cmd,genome,threads=2,task=None,bed_file=None):
    if bed_file:
        regions = load_bed_regions(bed_file)
    else:
        regions = get_genome_chunks(genome,nchunks=threads)
    
    parallel = Parallel(n_jobs=threads, return_as="generator")
    task = task if task else "Running command in parallel..."
    results = [r for r in tqdm(parallel(delayed(genome_job)(cmd,r) for r in regions),total=len(regions),desc=task)]
    print(results)


def run_cmd(cmd,log=None):
    print(cmd)
    output = open(log,"w") if log else sp.PIPE
    sp.run(cmd,shell=True,check=True,stderr=output,stdout=output)
