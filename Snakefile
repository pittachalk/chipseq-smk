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
indir       = config["input"]
outdir      = config["out"]
logdir      = config["log"]

# subdirectories of outdir
trimdir     = join(outdir, "trim/", "")
qcdir       = join(outdir, config["qc"], "")
bamdir      = join(outdir, config["bam"], "")
peaksdir    = join(outdir, config["peaks"], "")
pairsdir    = join(outdir, config["pairs"], "")
summarydir  = join(outdir, config["summary"], "")

# reading CSVs about fastq file and treated-control information
fastq_info = pd.read_csv("info-fastq.csv")
fastq_info["rep"] = fastq_info["rep"].astype(str)
fastq_info["lane"] = fastq_info["lane"].astype(str)
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
        #bamdir + "SING-rep1-1.sai",
        bamdir + "PAIR_A-rep1-1.sai",

        #trimdir + "PAIR_A-rep1-L1.pair1.fq.gz"
        #bamdir + "PAIR_A-rep1-L1.sai"
        #

        #bamdir + "{name}-rep{replic}-{lane}.sai"

        # # alignment files for individuals:
        # # sorted BAM files and quality control reports
        # expand(["{bamdir}{sample}_sorted.bam", "{qcdir}{sample}_fastqc.html"], 
        #    sample=config["samples"], bamdir=bamdir, qcdir=qcdir),

        # # MACS2 output for individuals: BED files for peaks, binary TDF for IGV
        # expand(["{peaksdir}{id}_peaks.narrowPeak","{peaksdir}{id}_linearFE_sorted.tdf", 
        #         "{peaksdir}{id}_logLR_sorted.tdf"], id=config["ids"], peaksdir=peaksdir),

        # # combined output from pairs of replicates:
        # # common peaks BED, IDR values, heatmap of peak profiles
        # expand(["{pairsdir}{combined}_commonpeaks.bed",
        #        "{pairsdir}{combined}_idrValues.txt",
        #        "{pairsdir}{combined}-peaks-matrix-heatmap.png"], 
        #    pairsdir=pairsdir, combined=config["combined"] ),

        # # summary files for everything:
        # # union of all peaks BED, PCA, correlation and heatmap of peak profiles
        # expand(["{summarydir}summary-unionpeaks.bed",
        #        "{summarydir}summarybw-peak-pca.png",
        #        "{summarydir}summarybw-peak-corr-heatmap.png",
        #        "{summarydir}summarybw-peaks-matrix-heatmap.png"],
        #    summarydir=summarydir)

        expand("test/{id}.commonpeaks.bed", id = control_info["id"].unique())
        # #"testRoutput.pdf"

#include: "rules/preprocessing.smk"
#include: "rules/align.smk"
#include: "rules/callpeak.smk"
#include: "rules/igv.smk"
#include: "rules/comparepairs.smk"
#include: "rules/summarize.smk"
#include: "rules/testrandomstuff.smk"
include: "rules/prototype.smk"


