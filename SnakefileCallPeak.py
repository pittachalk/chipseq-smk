CONFIG = "config.yaml" # change this, put quotes

# End users should not change anything below this line
# Parameters for the run should be modified in the CONFIG file

######################################################################
######################################################################
#     Setting up
######################################################################
configfile: CONFIG # do NOT change this
print(config)

# obtain directories from the CONFIG file
sampledir = config["sampledir"]
outputdir = config["outputdir"]
tempdir   = outputdir + config["subdir"]["tmp"]
logdir    = outputdir + config["subdir"]["log"]
qcdir     = outputdir + config["subdir"]["qc"]
macs2dir  = outputdir + config["subdir"]["macs2"]

######################################################################
######################################################################
#     The rule all
######################################################################

"""
Final output files (in outputdir):
- sorted BAM files (and their index)
- log files for bwa, bwasamse, fastqc, trimmomatic
- qc files for FASTQC and flagstat
- a copy of the config.yaml file
"""

rule all:
	input:
	    #expand("{outputdir}config.yaml", outputdir = config["outputdir"]),
	    expand(["{outputdir}{macs2}{id}_linearFE_sorted.tdf", "{outputdir}{macs2}{id}_logLR_sorted.tdf"], 
	    	outputdir = config["outputdir"], id=config["ids"], macs2 = config["subdir"]["macs2"]),

"""

######################################################################
######################################################################
#     Preprocessing of fastq.gz files and quality control
######################################################################

rule cat:
# concatenate fastq.gz files of each sample
	input:
		lambda x: map(lambda y: sampledir + y, config["samples"][x.sample])
	output:
		temp(tempdir + "{sample}_untrimmed.fastq.gz")
	shell:
		"cat {input} > {output}"

rule trim:
# trim Illumina adapters from fastq.gz
	input:
		tempdir + "{sample}_untrimmed.fastq.gz"
	output:
		temp(tempdir + "{sample}.fastq.gz")
	log:
		logdir + "trimmomatic/{sample}.log"
	params:
		bin=config["trimmomatic"]["bin"],
		settings=config["trimmomatic"]["settings"]
	shell:
		"java -jar {params.bin} SE "
		"{input} {output} {params.settings} 2>{log}"

rule decompress:
# unzip fastq.gz for BWA
	input:
		tempdir + "{sample}.fastq.gz"
	output:
		temp(tempdir + "{sample}.fastq")
	shell:
		"gunzip --keep {input}"

rule qctrim:
# run FastQC on trimmed fastq.gz files
	input:
		tempdir + "{sample}.fastq.gz"
	output:
		qcdir + "{sample}_fastqc.html",
		qcdir + "{sample}_fastqc.zip"
	log:
		logdir + "fastqc/{sample}.log"
	shell:
		"fastqc -o {qcdir} {input} 2>{log}"


######################################################################
######################################################################
#     BWA alignment and samtools
######################################################################

rule bwa_map:
# run bwa aln to find the SA coordinates of the input reads (.sai file)
# for paired-end, this part needs to be modified
	input:
		tempdir + "{sample}.fastq"
	output:
	    temp(tempdir + "{sample}.sai")
	log:
		logdir + "bwa/{sample}.log"
	threads: 8
	params:
		ref=config["refgenome"]
	shell:
		"bwa aln -t {threads} {params.ref} {input} 2>{log} >{output} "

rule bwa_samse:
# generate alignments from .sai file in the .sam format
# for paired-end, this part needs to be modified to use bwa_sampe
	input:
		tempdir + "{sample}.sai", 
		tempdir + "{sample}.fastq"
	output:
		temp(tempdir + "{sample}.sam")
	log:
		logdir + "bwa_samse/{sample}.log"
	shell:
		"bwa samse -n 50 data/genome.fa {input} 2>{log} >{output}"

rule samtools:
# get alignment stat with flagstat, sort SAM file, convert to BAM, index BAM file
	input:
		tempdir + "{sample}.sam"
	output:
		flagstat = qcdir + "{sample}_alignstat.txt",
		sortedbam = protected(outputdir + "{sample}_sorted.bam"),
		bai = protected(outputdir + "{sample}_sorted.bam.bai")
	shell:
		"samtools flagstat {input} > {output.flagstat}; "
		"samtools view -bS {input} | "
		"samtools sort - -o {output.sortedbam}; "
		"samtools index {output.sortedbam}"

"""

######################################################################
######################################################################
#     ChIP-seq: MACS2
######################################################################

rule macs2:
# call peaks with MACS2 
	input:
		sample = outputdir + "{id}_sorted.bam",
		#control = lambda x: map(lambda y: outputdir + y + "_sorted.bam", config["ids"][x.id])
		control = lambda x: outputdir + config["ids"][x.id] + "_sorted.bam"
	output:
		#map(lambda macs2outfile: macs2dir + "{id}-" + config["ids"][x.id] + macs2outfile, 
		#	["_peaks.narrowPeak"]),
		map(lambda macs2outfile: macs2dir + "{id}_" + macs2outfile, 
			["peaks.narrowPeak", "treat_pileup.bdg", "control_lambda.bdg"])
	log:
		logdir + "macs2/{id}.log"
	shell:
		"macs2 callpeak -t {input.sample} -c {input.control} "
		"--name {wildcards.id} --outdir " + macs2dir + " "
		"--gsize 8.8e8 --extsize 147 --nomodel -q 0.01 -B --cutoff-analysis"

rule bedgraph:
# prepare bedgraph of linear fold enrichment and log10 likelihood
	input: 
		sample = macs2dir + "{id}_treat_pileup.bdg",
		control = macs2dir + "{id}_control_lambda.bdg"
	output:
		FE = macs2dir + "{id}_linearFE.bdg",
		logLR = macs2dir + "{id}_logLR.bdg"
	shell:
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output.FE} -m FE; "
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output.logLR} -m logLR --pseudocount 0.00001"

rule igvsort:
# sort bedgraph
	input: 
		FE = macs2dir + "{id}_linearFE.bdg",
		logLR = macs2dir + "{id}_logLR.bdg"
	output: 
		FE = temp(tempdir + "{id}_linearFE_sorted.bdg"),
		logLR = temp(tempdir + "{id}_logLR_sorted.bdg")
	params:
		bin=config["igvtools"]["bin"]
	shell:
		"{params.bin} sort {input.FE} {output.FE}; "
		"{params.bin} sort {input.logLR} {output.logLR}"

rule igvtotdf:
# convert bedgraph to TDF for IGV
	input: 
		FE = tempdir + "{id}_linearFE_sorted.bdg",
		logLR = tempdir + "{id}_logLR_sorted.bdg"
	output: 
		FE = macs2dir + "{id}_linearFE_sorted.tdf",
		logLR = macs2dir + "{id}_logLR_sorted.tdf"
	params:
		bin=config["igvtools"]["bin"],
		ref=config["refgenome"]
	shell:
		"{params.bin} toTDF {input.FE} {output.FE} {params.ref}; "
		"{params.bin} toTDF {input.logLR} {output.logLR} {params.ref}"

"""
######################################################################
######################################################################
#     Miscellaneous
######################################################################

rule copyconfig:
# copy config.yaml into the outputdir
	input:
		CONFIG
	output:
		protected(outputdir + "config.yaml")
	shell:
		"cp {input} {output}"


rule chipseq:
# echo a chipseq command
	input:
		sample = outputdir + "{id}_sorted.bam.bai",
		control = lambda x: map(lambda y: outputdir + y + "_sorted.bam.bai", config["ids"][x.id])
	output:
		outputdir + "{id}_chipseq.txt"
	shell:
		"echo -t {input.sample} -c {input.control} > {output}"
"""
