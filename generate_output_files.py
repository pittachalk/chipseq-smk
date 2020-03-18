# this file generates the output files for the rule all

from snakemake.io import expand

def bam_indexes(bamdir, fastq_info):
    bamindices = expand(bamdir + "{name}-rep{replic}.merged.sorted.bam.bai", 
                        name=fastq_info["name"].unique(), replic=fastq_info["rep"].unique())
    return bamindices

def align_stats(bamdir, fastq_info):
    alignstats = expand(bamdir + "{name}-rep{replic}.merged.alignstat.txt", 
                        name=fastq_info["name"].unique(), replic=fastq_info["rep"].unique())
    return alignstats

def tdf_files(bamdir, fastq_info):
    tdf_files = expand(bamdir + "{name}-rep{replic}.merged.sorted.tdf", 
                       name=fastq_info["name"].unique(), replic=fastq_info["rep"].unique())
    return tdf_files