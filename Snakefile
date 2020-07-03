configfile: "config.yaml"

import pandas as pd
import itertools
from os.path import join
import generate_output_files as gof

# obtain directories from the config file
indir       = config["sampledir"]
outdir      = "output/"
logdir      = "log/"

# subdirectories of outdir
trimdir     = join(outdir, "trim/", "")
qcdir       = join(outdir, "qc/", "")
bamdir      = join(outdir, "bam/", "")
bwdir       = join(outdir, "bw/", "")
peaksdir    = join(outdir, "peaks/", "") #join(outdir, "macs2_individual/", "")

# read which fastq files belong to each sample
fastq_info = pd.read_csv("info-fastq.csv")
fastq_info["rep"] = fastq_info["rep"].astype(str)
fastq_info = fastq_info.set_index(["name", "rep", "lane"], drop=False)

# read which control belong to each sample
control_info = pd.read_csv("info-control.csv")
control_info["idrep"]      = control_info["idrep"].astype(str)
control_info["treatrep"]   = control_info["treatrep"].astype(str)
control_info["controlrep"] = control_info["controlrep"].astype(str)
control_info = control_info.set_index(["id", "idrep"], drop=False)

# output files
set1 = set(expand("{name}-rep{replic}", zip, name=fastq_info["name"], replic=fastq_info["rep"]))
set2 = set(expand("{id}-rep{idrep}", zip, id=control_info["id"], idrep=control_info["idrep"]))
rule all:
    input:
        expand(bamdir + "{name_rep}.merged.alignstat.txt", name_rep=set1),
        expand(bamdir + "{name_rep}.merged.sorted.bam", name_rep=set1),
        expand(bamdir + "{name_rep}.merged.sorted.bam.bai", name_rep=set1),
        expand(bamdir + "{name_rep}.merged.sorted.tdf", name_rep=set1),
        expand(bamdir + "{name_rep}.merged.sorted.tdf", name_rep=set1),
        expand(peaksdir + "{id_idrep}.relaxed.bed", id_idrep=set2)

# rules to include
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/callpeak.smk"
include: "rules/igv.smk"
