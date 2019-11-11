# specify the config file
CONFIG = "config.yaml" # change this, put quotes
configfile: CONFIG # do NOT change this
print(config)

sampledir = config["sampledir"]
outputdir = config["outputdir"]
tempdir   = outputdir + config["subdir"]["tmp"]
logdir    = outputdir + config["subdir"]["log"]
qcdir     = outputdir + config["subdir"]["qc"]

rule all:
	input:
		# expand takes place in initialisation, we do not know wildcard values expanded upon above
		# we have to take what is present in the config file
	    expand("{outputdir}{sample}_T_sorted.bam.bai", sample=config["samples"],
	    	outputdir = config["outputdir"]),
	    #expand("{outputdir}{qc}{sample}_fastqc.html", sample=config["samples"],
	    #	outputdir = config["outputdir"], qc = config["subdir"]["qc"]),
	    expand("{outputdir}{qc}{sample}_T_fastqc.html", sample=config["samples"], 
	    	outputdir = config["outputdir"], qc = config["subdir"]["qc"]),
	    expand("{outputdir}config.yaml", outputdir = config["outputdir"])

rule copyconfig:
	input:
		CONFIG
	output:
		protected(outputdir + "config.yaml")
	shell:
		"cp {input} {output}"
 
rule cat:
	# a nested list of lists
	input:
		lambda x: map(lambda y: sampledir + y, config["samples"][x.sample])
	output:
		temp(tempdir + "{sample}.fastq.gz")
	shell:
		"cat {input} > {output}"

rule trim:
	input:
		tempdir + "{sample}.fastq.gz"
	output:
		temp(tempdir + "{sample}_T.fastq.gz")
	log:
		logdir + "trimmomatic/{sample}.log"
	shell:
		"java -jar /mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/trimmomatic-0.39.jar SE "
		"{input} {output} "
		"ILLUMINACLIP:/mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 "
		"ILLUMINACLIP:/mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 "
		"SLIDINGWINDOW:4:28 MINLEN:20  2>{log}"

rule decompress:
	input:
		tempdir + "{sample}_T.fastq.gz"
	output:
		temp(tempdir + "{sample}_T.fastq")
	shell:
		"gunzip --keep {input}"

"""
rule qc:
	input:
		tempdir + "{sample}.fastq.gz"
	output:
		qcdir + "{sample}_fastqc.html",
		qcdir + "{sample}_fastqc.zip"
	shell:
		"fastqc -o qc/ {input}"
"""

rule qctrim:
	input:
		tempdir + "{sample}_T.fastq.gz"
	output:
		qcdir + "{sample}_T_fastqc.html",
		qcdir + "{sample}_T_fastqc.zip"
	log:
		logdir + "fastqc/{sample}.log"
	shell:
		"fastqc -o {qcdir} {input} 2>{log}"

rule bwa_map:
	input:
		"data/genome.fa",
		tempdir + "{sample}_T.fastq"
	output:
	    temp(tempdir + "{sample}_T.sai")
	log:
		logdir + "bwa/{sample}.log"
	threads: 8
	shell:
		"bwa aln -t {threads} {input} 2>{log} >{output} "

rule bwa_samse: # for paired-end, may need to do expand for both
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

#for i in $(ls *_R1_* | sed 's/_L00[1|2]_R1_001.fastq.gz//' | sort -u ); 
#do cat ${i}*L001_R1* ${i}*L002_R1* > ${i}_R1.fastq.gz   & done

"""
rule bwa_map:
	input:
		"data/genome.fa",
		lambda wildcards: config["samples"][wildcards.sample]
	output:
	    temp("mapped_reads/{sample}.bam")
	log:
		"logs/bwa_mem/{sample}.log"
	threads: 8
	shell:
		"bwa mem -t {threads} {input} | "
		"samtools view  -Sb - > {output} 2> {log}"

rule samtools_sort:
	input:
		"mapped_reads/{sample}.bam"
	output:
		protected("sorted_reads/{sample}.bam")
	shell:
		"samtools sort -T sorted_reads/{wildcards.sample} "
		"-O bam {input} > {output}"

rule samtools_index:
	input:
		"sorted_reads/{sample}.bam"
	output:
		"sorted_reads/{sample}.bam.bai"
	shell:
		"samtools index {input}"
"""
