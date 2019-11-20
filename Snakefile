configfile: "config.yaml"  # put quotes

# End users should not change anything below this line
# Parameters for the run should be modified in the CONFIG file

######################################################################
######################################################################
#     Setting up
######################################################################
# obtain directories from the CONFIG file
sampledir = config["sampledir"]
outputdir = config["outputdir"]
tempdir   = outputdir + config["subdir"]["tmp"]
logdir    = outputdir + config["subdir"]["log"]
qcdir     = outputdir + config["subdir"]["qc"]
bamdir    = outputdir + config["subdir"]["bam"]
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

the user can comment out what is not required
"""

rule all:
	input:
	    expand("{outputdir}{bamdir}{sample}_sorted.bam", sample=config["samples"],
	    	outputdir = config["outputdir"], bamdir = config["subdir"]["bam"]),
	    expand("{outputdir}{qc}{sample}_fastqc.html", sample=config["samples"], 
	    	outputdir = config["outputdir"], qc = config["subdir"]["qc"]),
	    expand("{outputdir}{macs2}{id}_peaks.narrowPeak",
	    	outputdir = config["outputdir"], id=config["ids"], macs2 = config["subdir"]["macs2"]),
	    expand(["{outputdir}{macs2}{id}_linearFE_sorted.tdf", "{outputdir}{macs2}{id}_logLR_sorted.tdf"], 
	    	outputdir = config["outputdir"], id=config["ids"], macs2 = config["subdir"]["macs2"])


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
	params:
		ref=config["refgenome"]
	shell:
		"bwa samse -n 50 {params.ref} {input} 2>{log} >{output}"

rule samtools:
# get alignment stat with flagstat, sort SAM file, convert to BAM, index BAM file
	input:
		tempdir + "{sample}.sam"
	output:
		flagstat = qcdir + "{sample}_alignstat.txt",
		sortedbam = protected(bamdir + "{sample}_sorted.bam"),
		bai = protected(bamdir + "{sample}_sorted.bam.bai")
	shell:
		"samtools flagstat {input} > {output.flagstat}; "
		"samtools view -bS {input} | "
		"samtools sort - -o {output.sortedbam}; "
		"samtools index {output.sortedbam}"


######################################################################
######################################################################
#     ChIP-seq: MACS2
######################################################################

rule macs2:
# call peaks with MACS2 
	input:
		sample = bamdir + "{id}_sorted.bam",
		control = lambda x: bamdir + config["ids"][x.id] + "_sorted.bam"
	output:
		map(lambda macs2outfile: macs2dir + "{id}_" + macs2outfile, 
			["peaks.narrowPeak", "treat_pileup.bdg", "control_lambda.bdg"])
	log:
		logdir + "macs2/{id}.log"
	params:
		settings=config["macs2"]["settings"]
	shell:
		"macs2 callpeak -t {input.sample} -c {input.control} "
		"--name {wildcards.id} --outdir " + macs2dir + " "
		"{params.settings} 2>{log}"

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


######################################################################
######################################################################
#     Convert bedgraphs to the TDF and BigWig binary formats
######################################################################

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

# coming soon: bedgraphtobigwig, mergebigwig


######################################################################
######################################################################
#     Get common peaks between the replicates
######################################################################

# coming soon: getcommmonpeaks, extendcommonpeaks


######################################################################
######################################################################
#     Summary files
######################################################################

# coming soon: idr, multibigwigsummary, pcacorr, heatmap (maybe?)



