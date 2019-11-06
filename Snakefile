configfile: "config.yaml"

rule all:
	input:
	    "haha.txt"

rule cat:
	input:
		lambda wildcards: config["samples"][wildcards.sample]
	output:
		"0_raw/{sample}.fastq.gz" # temporary
	shell:
		"cat {input} > {output}"

rule trim:
	input:
		"0_raw/{sample}.fastq.gz"
	output:
		"0_raw/{sample}_T.fastq.gz"
	log:
		"log/trimmomatic/{sample}.log"
	shell:
		"java -jar /mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/trimmomatic-0.39.jar SE "
		"{input} {output} "
		"ILLUMINACLIP:/mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 "
		"ILLUMINACLIP:/mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 "
		"SLIDINGWINDOW:4:28 MINLEN:20  2> {log}"

rule decompress:
	input:
		"0_raw/{sample}_T.fastq.gz"
	output:
		"0_raw/{sample}_T.fastq" # temporary
	shell:
		"gunzip {input}"

rule qc:
	input:
		"0_raw/{sample}.fastq.gz"
	output:
		"qc/{sample}_fastqc.html"
		"qc/{sample}_fastqc.zip"
	shell:
		"fastqc -o qc/ {input}"

rule qctrim:
	input:
		"0_raw/{sample}_T.fastq.gz"
	output:
		"qc/{sample}_T_fastqc.html"
		"qc/{sample}_T_fastqc.zip"
	shell:
		"fastqc -o  {input}"

rule bwa_map:
	input:
		"data/genome.fa",
		"0_raw/{sample}_T.fastq"
	output:
	    temp("2_mapped/{sample}_T.sai")
	log:
		"logs/bwa/{sample}.log"
	threads: 8
	shell:
		"bwa aln -t {threads} {input} > {output}"

rule bwa_samse: # for paired-end, may need to do expand for both
	input:
		"2_mapped/{sample}_T.sai", 
		"0_raw/{sample}_T.fastq"
	output:
		"2_mapped/{sample}_T.sam"
	shell:
		"bwa samse -n 50 data/genome.fa {input} >  {output}"

rule samtools_sort:
	input:
		"2_mapped/{sample}.bam"
	output:
		bam=protected("2_mapped/{sample}.bam"),
		sortedbam=protected("2_mapped/{sample}_sorted.bam"),
		bai=protected("2_mapped/{sample}_sorted.bam.bai")
	shell:
		"samtools_mt view -bS -o {output.bam} {input};"
		"samtools_mt sort {output.bam} -o {output.sortedbam};"
		"samtools_mt index {output.sortedbam}"

rule countfile:
	input:
		gz=expand("2_mapped/{sample}_sorted.bam.bai", sample=config["samples"])
	output:
		"haha.txt"
	shell:
		"ls {input} | wc > {output}"



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
