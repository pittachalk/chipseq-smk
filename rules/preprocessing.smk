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
		settings=config["trimmomatic"]["settings"]
	shell:
		"trimmomatic SE {input} {output} {params.settings} 2>{log}"

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
