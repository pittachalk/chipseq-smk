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
	    expand("{outputdir}{sample}_sorted.bam", sample=config["samples"],
	    	outputdir = config["outputdir"]),
	    expand("{outputdir}{qc}{sample}_fastqc.html", sample=config["samples"], 
	    	outputdir = config["outputdir"], qc = config["subdir"]["qc"]),
	    expand("{outputdir}config.yaml", outputdir = config["outputdir"]),
	    expand("{outputdir}{id}_chipseq.txt", outputdir = config["outputdir"], id=config["ids"])

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


######################################################################
######################################################################
#     ChIP-seq: MACS2
######################################################################

# coming soon!


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
