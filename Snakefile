configfile: "config.yaml"

sampledir = config["sampledir"]
outputdir = config["outputdir"] # config["outputdir"]
tempdir   = outputdir + "temp/"

rule all:
	input:
		# expand takes place in initialisation, we do not know wildcard values expanded upon above
		# we have to take what is present in the config file
	    expand("{outputdir}{sample}_T_sorted.bam.bai", sample=config["samples"],
	    	outputdir = config["outputdir"])
	    #expand("qc/{sample}_fastqc.html", sample=config["samples"]),
	    #expand("qc/{sample}_T_fastqc.html", sample=config["samples"])
	    
rule cat:
	# a nested list of lists
	input:
		lambda x: map(lambda y: sampledir + y, config["samples"][x.sample])
	output:
		temp("0_raw/{sample}.fastq.gz")
	shell:
		"cat {input} > {output}"

rule trim:
	input:
		"0_raw/{sample}.fastq.gz"
	output:
		temp("0_raw/{sample}_T.fastq.gz")
	log:
		"logs/trimmomatic/{sample}.log"
	shell:
		"java -jar /mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/trimmomatic-0.39.jar SE "
		"{input} {output} "
		"ILLUMINACLIP:/mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 "
		"ILLUMINACLIP:/mnt/home1/miska/jlp76/programs/Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 "
		"SLIDINGWINDOW:4:28 MINLEN:20  2>{log}"

rule decompress:
	input:
		"0_raw/{sample}_T.fastq.gz"
	output:
		temp("0_raw/{sample}_T.fastq")
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
		"fastqc -o qc/ {input}"

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
		"bwa aln -t {threads} {input} 2>{log} >{output} "

rule bwa_samse: # for paired-end, may need to do expand for both
	input:
		"2_mapped/{sample}_T.sai", 
		"0_raw/{sample}_T.fastq"
	output:
		temp("2_mapped/{sample}_T.sam")
	log:
		"logs/bwa_samse/{sample}.log"
	shell:
		"bwa samse -n 50 data/genome.fa {input} 2>{log} >{output}"

rule samtools_sort:
	input:
		"2_mapped/{sample}_T.sam"
	output:
		sortedbam=protected(outputdir + "{sample}_T_sorted.bam"),
		bai=protected(outputdir + "{sample}_T_sorted.bam.bai")
	shell:
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
