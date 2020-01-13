######################################################################
######################################################################
#     Preprocessing of fastq.gz files and quality control
######################################################################

# rule cat:
# # concatenate fastq.gz files of each sample
# 	input:
# 		lambda x: map(lambda y: indir + y, config["samples"][x.sample])
# 	output:
# 		temp(trimdir + "{sample}_untrimmed.fastq.gz")
# 	conda:
# 		"../envs/py3.yml"
# 	shell:
# 		"cat {input} > {output}"

# rule trim:
# # trim Illumina adapters from fastq.gz
# 	input:
# 		trimdir + "{sample}_untrimmed.fastq.gz"
# 	output:
# 		temp(trimdir + "{sample}.fastq.gz")
# 	log:
# 		logdir + "trimmomatic/{sample}.log"
# 	params:
# 		settings=config["trimmomatic"]["settings"]
# 	conda:
# 		"../envs/py3.yml"
# 	shell:
# 		"trimmomatic SE {input} {output} {params.settings} 2>{log}"

# rule decompress:
# # unzip fastq.gz for BWA
# 	input:
# 		trimdir + "{sample}.fastq.gz"
# 	output:
# 		temp(trimdir + "{sample}.fastq")
# 	conda:
# 		"../envs/py3.yml"
# 	shell:
# 		"gunzip --keep {input}"

rule qctrim:
# run FastQC on trimmed fastq.gz files
	input:
		trimdir + "{sample}.fastq.gz"
	output:
		qcdir + "{sample}_fastqc.html",
		qcdir + "{sample}_fastqc.zip"
	log:
		logdir + "fastqc/{sample}.log"
	conda:
		"../envs/py3.yml"
	shell:
		"fastqc -o {qcdir} {input} 2>{log}"
