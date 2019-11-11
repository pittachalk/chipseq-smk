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
	    expand("{outputdir}{sample}_T_sorted.bam", sample=config["samples"],
	    	outputdir = config["outputdir"]),
	    expand("{outputdir}{qc}{sample}_T_fastqc.html", sample=config["samples"], 
	    	outputdir = config["outputdir"], qc = config["subdir"]["qc"]),
	    expand("{outputdir}config.yaml", outputdir = config["outputdir"])


######################################################################
######################################################################
#     Preprocessing of fastq.gz files and quality control
######################################################################

rule cat:
# concatenate fastq.gz files of each sample
	input:
		lambda x: map(lambda y: sampledir + y, config["samples"][x.sample])
	output:
		temp(tempdir + "{sample}.fastq.gz")
	shell:
		"cat {input} > {output}"

rule trim:
# trim Illumina adapters from fastq.gz
	input:
		tempdir + "{sample}.fastq.gz"
	output:
		temp(tempdir + "{sample}_T.fastq.gz")
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
		tempdir + "{sample}_T.fastq.gz"
	output:
		temp(tempdir + "{sample}_T.fastq")
	shell:
		"gunzip --keep {input}"

rule qctrim:
# run FastQC on trimmed fastq.gz files
	input:
		tempdir + "{sample}_T.fastq.gz"
	output:
		qcdir + "{sample}_T_fastqc.html",
		qcdir + "{sample}_T_fastqc.zip"
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
		tempdir + "{sample}_T.fastq"
	output:
	    temp(tempdir + "{sample}_T.sai")
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
		tempdir + "{sample}_T.sai", 
		tempdir + "{sample}_T.fastq"
	output:
		temp(tempdir + "{sample}_T.sam")
	log:
		logdir + "bwa_samse/{sample}.log"
	shell:
		"bwa samse -n 50 data/genome.fa {input} 2>{log} >{output}"

rule samtools:
# get alignment stat with flagstat, sort SAM file, convert to BAM, index BAM file
	input:
		tempdir + "{sample}_T.sam"
	output:
		flagstat = qcdir + "{sample}_T_alignstat.txt",
		sortedbam = protected(outputdir + "{sample}_T_sorted.bam"),
		bai = protected(outputdir + "{sample}_T_sorted.bam.bai")
	shell:
		"samtools flagstat {input} > {output.flagstat}; "
		"samtools view -bS {input} | "
		"samtools sort - -o {output.sortedbam}; "
		"samtools index {output.sortedbam}"


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
