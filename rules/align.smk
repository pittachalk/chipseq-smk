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
		sortedbam = bamdir + "{sample}_sorted.bam",
		bai = bamdir + "{sample}_sorted.bam.bai"
	shell:
		"samtools flagstat {input} > {output.flagstat}; "
		"samtools view -bS {input} | "
		"samtools sort - -o {output.sortedbam}; "
		"samtools index {output.sortedbam}"