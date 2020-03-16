# this file generates the output files for the rule all

from snakemake.io import expand

def bam_indexes(bamdir, fastq_info):
    bamindices = expand(bamdir+"{name}-rep{replic}.merged.sorted.bam.bai", 
                        name=fastq_info["name"].unique(), replic=fastq_info["rep"].unique())
    return bamindices
