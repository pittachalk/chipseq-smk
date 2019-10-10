configfile: "config.yaml"

rule all:
	input:
	    "plots/quals.svg"

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

rule bcftools_call:
	input:
		fa="data/genome.fa",
		bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
		bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
	output:
		"calls/all.vcf"
	log:
		"logs/bcftools_call/bcftools.log"
	params:
		u=config["mutationrate"]
	shell:
		"samtools mpileup -g -f {input.fa} {input.bam} | "
		"bcftools call -P {params.u} -mv - > {output} 2> {log}"

rule plot_quals:
	input:
		"calls/all.vcf"
	output:
		"plots/quals.svg"
	script:
		"scripts/plot-quals.py"

