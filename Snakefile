configfile: "config.yaml"  # put quotes

# End users should not change anything below this line
# Parameters for the run should be modified in the CONFIG file

######################################################################
######################################################################
#     Setting up
######################################################################
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

import pandas as pd

fastq_info = pd.read_csv("info-fastq.csv")
fastq_info["rep"] = fastq_info["rep"].astype(str)
fastq_info = fastq_info.set_index(["name", "rep", "lane"], drop=False)

control_info = pd.read_csv("info-control.csv")
control_info["idrep"]      = control_info["idrep"].astype(str)
control_info["treatrep"]   = control_info["treatrep"].astype(str)
control_info["controlrep"] = control_info["controlrep"].astype(str)

control_info = control_info.set_index(["id", "idrep"], drop=False)
#control_info = control_info.set_index("species", drop=False)

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


def is_single_end(name, replic, lane):
    """
    Check if there is a second fastq file for the given sample.
    """
    return pd.isnull(fastq_info.loc[(name, replic, lane), "fq2"])

def get_fq_se(wildcards):
    if(is_single_end(wildcards.name, wildcards.replic, wildcards.lane)):
        return "data/samples/" + fastq_info.loc[(wildcards.name, wildcards.replic, wildcards.lane), "fq1"]
    else:
        raise ValueError("This is a paired-end sample")

def get_fq_pe(wildcards):
    if(is_single_end(wildcards.name, wildcards.replic, wildcards.lane)):
        raise ValueError("This is a single-end sample")
    else:
        return "data/samples/" + fastq_info.loc[(wildcards.name, wildcards.replic, wildcards.lane), ["fq1", "fq2"]]

rule proto_trimSE:
    # output: {name}-{replicate}-{lane}.txt   
    input:
        get_fq_se
    output:
        fq = "test/{name}-rep{replic}-{lane}.fq.gz"
	#log:
	#	logdir + "test/trimPE/{name}-rep{replic}-{lane}.log"
	#params:
	#	settings=config["trimmomatic"]["settings"]
    shell:
        "cat {input} > {output}"
		#"trimmomatic PE {input} {output} {params.settings} 2>{log}"

rule proto_trimPE:
    # output: {name}-{replicate}-{lane}.1.txt, {name}-{replicate}-{lane}.2.txt
    input:
        get_fq_pe
    output:
        fq1 = "test/{name}-rep{replic}-{lane}.1.fq.gz",
        fq2 = "test/{name}-rep{replic}-{lane}.2.fq.gz"
	#log:
    #	logdir + "test/trimPE/{name}-{replic}-{lane}.log"
	#params:
	#	settings=config["trimmomatic"]["settings"]
    shell:
        "cat {input} > {output}"
		#"trimmomatic PE {input} {output} {params.settings} 2>{log}"
    

#rule proto_catlanes:
    # output: {name}-{replicate}.1.txt, {name}-{replicate}.2.txt
    # output: {name}-{replicate}.txt'
    # NOT SURE if necessary (maybe needed for QC)
    #input: 

def get_trimmed_fq(wildcards):
    if is_single_end(wildcards.name, wildcards.replic, wildcards.lane):
        return "test/{name}-rep{replic}-{lane}.fq.gz".format(name = wildcards.name, replic = wildcards.replic, lane = wildcards.lane)
    else:
        return expand("test/{name}-rep{replic}-{lane}.{pair}.fq.gz", 
        name = wildcards.name, replic = wildcards.replic, lane = wildcards.lane, pair = [1, 2])

rule proto_bwa:
    # bwa is run on each lane separately
    # this is the step where paired end reads are combined
    input:
        fq = get_trimmed_fq
    output:
        "test/{name}-rep{replic}-{lane}.bam"
    shell:
        "cat {input} > {output}"


# def get_bam_units(wildcards):
#     # return "mapped/{name}-*.bam".format(name = wildcards.name)
#     units = unit_info[unit_info["name"] == wildcards.name]["unit"]
#     return expand("mapped/{name}-{unit}.bam", name = wildcards.name, unit = units)


def get_bam_lanes(wildcards):
    lanes   = fastq_info[(fastq_info["name"] == wildcards.name) &  (fastq_info["rep"] == wildcards.replic)]["lane"]
    return expand("test/{name}-rep{replic}-{lane}.bam", name = wildcards.name, replic = wildcards.replic, lane = lanes)

rule proto_cat_bam:
   # now we cat the bam files from separate lanes
   # this will be replaced with conversion to sam, so the cat is just placeholder
   # output: {name}.{replicate}-aligned.txt
   input:
       fq = get_bam_lanes
   output:
       "test/{name}-rep{replic}.merged.bam"
   shell:
       "cat {input} > {output}"

def get_treatvscontrol_bam(wildcards):
    # this should give a single row describing the treatment vs control pairing
    x = control_info[(control_info["id"] == wildcards.id) & (control_info["idrep"] == wildcards.idrep)]

    treat = expand("test/{treatname}-rep{treatrep}.merged.bam", treatname = x["treatname"], treatrep = x["treatrep"])
    control = expand("test/{controlname}-rep{controlrep}.merged.bam", controlname = x["controlname"], controlrep = x["controlrep"])
    
    assert (len(treat) == 1), "treated list not length 1, are id-idrep pairs unique?"
    assert (len(control) == 1), "control list not length 1, are id-idrep pairs unique?"

    return (treat[0], control[0])

rule proto_vsinput:
    # compare vs input
    # output: {name}-controlled.txt
    input:
        bam = get_treatvscontrol_bam
    output:
        "test/{id}-{idrep}.narrowPeak"
    shell:
        "cat t {input.bam[0]} c {input.bam[1]} > {output}"



def are_there_replicates(id):
    # check if there are replicates for a given id
    return (len(control_info.loc[id]) > 1)

#def get_replicate_pairs(wildcards):
    # returns all pairwise combinations of replicates for a given id, if those exist

def get_replicate_list(wildcards):
    # returns a list of of replicates for a given id, if those exist
    if(are_there_replicates(wildcards.id)):
        return expand("test/{id}-{idrep}.narrowPeak", id = wildcards.id, idrep = control_info.loc[wildcards.id]["idrep"])
    else:
        raise ValueError("There are no replicates for this ID.")


rule proto_combinereplicates:
    # equivalent to calling common peaks
    # output: {name}-commonpeaks.txt
    input:
        replicates = get_replicate_list
    output:
        "test/{id}.commonpeaks.bed"
    shell:
        "cat a {input.replicates[0]} b {input.replicates[1]} > {output}"