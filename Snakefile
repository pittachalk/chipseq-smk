configfile: "config.yaml"  # put quotes

# End users should not change anything below this line
# Parameters for the run should be modified in the CONFIG file

######################################################################
######################################################################
#     Setting up
######################################################################

import pandas as pd
import itertools
from os.path import join

# obtain directories from the CONFIG file
indir       = config["sampledir"]
outdir      = "output/"
logdir      = "log/"

# subdirectories of outdir
trimdir     = join(outdir, "trim/", "")
qcdir       = join(outdir, "qc/", "")
bamdir      = join(outdir, "bam/", "")
bwdir       = join(outdir, "bw/", "")
peaksdir    = join(outdir, "macs2_individual/", "")
pairsdir    = join(outdir, "macs2_sample/", "")
summarydir  = join(outdir, "macs2_summary/", "")
overalldir  = join(outdir, "overall/", "")

# reading CSVs about fastq file and treated-control information
fastq_info = pd.read_csv("info-fastq.csv")
fastq_info["rep"] = fastq_info["rep"].astype(str)
fastq_info = fastq_info.set_index(["name", "rep", "lane"], drop=False)

control_info = pd.read_csv("info-control.csv")
control_info["idrep"]      = control_info["idrep"].astype(str)
control_info["treatrep"]   = control_info["treatrep"].astype(str)
control_info["controlrep"] = control_info["controlrep"].astype(str)
control_info = control_info.set_index(["id", "idrep"], drop=False)

######################################################################
######################################################################
#     The rule all
######################################################################

"""
Final output files
Note that the user can comment out what is not required
"""

rule all:
    input:
        expand(summarydir + "{id}.commonpeaks.bed", id = control_info["id"].unique()),
        expand(summarydir + "{id}.intersect.bed", id = control_info["id"].unique()),
        expand(summarydir + "{id}-peaks-matrix-heatmap.png", id = control_info["id"].unique()),
        overalldir + "overall.unionpeaks.bed",
        overalldir + "overallbw-peaks-matrix-heatmap.png"

include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/callpeak.smk"
include: "rules/igv.smk"
include: "rules/comparepairs.smk"
include: "rules/summarize.smk"
#include: "rules/testrandomstuff.smk"


