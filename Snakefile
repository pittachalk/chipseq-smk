configfile: "config.yaml"  # put quotes

# End users should not change anything below this line
# Parameters for the run should be modified in the CONFIG file

######################################################################
######################################################################
#     Setting up
######################################################################

import pandas as pd
import itertools

# obtain directories from the CONFIG file
sampledir   = config["sampledir"]
outputdir   = config["outputdir"]
tempdir     = outputdir + config["subdir"]["tmp"]
logdir      = outputdir + config["subdir"]["log"]
qcdir       = outputdir + config["subdir"]["qc"]
bamdir      = outputdir + config["subdir"]["bam"]
macs2dir    = outputdir + config["subdir"]["macs2"]
summarydir  = macs2dir + config["subdir"]["summary"]  # subdir for all samples
combineddir = macs2dir + config["subdir"]["combined"] # subdir for combined replicates

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
        expand("test/{id}.commonpeaks.bed", id = control_info["id"].unique())

        #"testRoutput.pdf"

        ## alignment files for individuals:
        ## sorted BAM files and quality control reports
        #expand(["{outputdir}{bamdir}{sample}_sorted.bam",
        #        "{outputdir}{qc}{sample}_fastqc.html"], 
        #    sample=config["samples"], outputdir=config["outputdir"], 
        #    bamdir=config["subdir"]["bam"], qc = config["subdir"]["qc"]),
#
        ## MACS2 output for individuals: BED files for peaks, binary TDF for IGV
        #expand(["{outputdir}{macs2}{id}_peaks.narrowPeak",
        #        "{outputdir}{macs2}{id}_linearFE_sorted.tdf", 
        #        "{outputdir}{macs2}{id}_logLR_sorted.tdf"],
        #    outputdir=config["outputdir"], id=config["ids"], macs2=config["subdir"]["macs2"]),
#
        ## combined output from pairs of replicates:
        ## common peaks BED, IDR values, heatmap of peak profiles
        #expand(["{outputdir}{macs2}{combineddir}{combined}_commonpeaks.bed",
        #        "{outputdir}{macs2}{combineddir}{combined}_idrValues.txt",
        #        "{outputdir}{macs2}{combineddir}{combined}-peaks-matrix-heatmap.png"],
        #    outputdir=config["outputdir"], macs2=config["subdir"]["macs2"], 
        #    combineddir=config["subdir"]["combined"], combined=config["combined"] ),
#
        ## summary files for everything:
        ## union of all peaks BED, PCA, correlation and heatmap of peak profiles
        #expand(["{outputdir}{macs2}{summarydir}summary-unionpeaks.bed",
        #        "{outputdir}{macs2}{summarydir}summarybw-peak-pca.png",
        #        "{outputdir}{macs2}{summarydir}summarybw-peak-corr-heatmap.png",
        #        "{outputdir}{macs2}{summarydir}summarybw-peaks-matrix-heatmap.png"],
        #    outputdir=config["outputdir"],
        #    macs2=config["subdir"]["macs2"], summarydir=config["subdir"]["summary"] )

#include: "rules/preprocessing.smk"
#include: "rules/align.smk"
#include: "rules/callpeak.smk"
#include: "rules/igv.smk"
#include: "rules/comparepairs.smk"
#include: "rules/summarize.smk"
#include: "rules/testrandomstuff.smk"
include: "rules/prototype.smk"


