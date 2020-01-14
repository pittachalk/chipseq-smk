rule cat:
# concatenate fastq.gz files of each sample
	input:
		lambda x: map(lambda y: indir + y, config["samples"][x.sample])
	output:
		temp(trimdir + "{sample}_untrimmed.fastq.gz")
	conda:
		"../envs/py3.yml"
	shell:
		"cat {input} > {output}"

rule trim:
# trim Illumina adapters from fastq.gz
	input:
		trimdir + "{sample}_untrimmed.fastq.gz"
	output:
		temp(trimdir + "{sample}.fastq.gz")
	log:
		logdir + "trimmomatic/{sample}.log"
	params:
		settings=config["trimmomatic"]["settings"]
	conda:
		"../envs/py3.yml"
	shell:
		"trimmomatic SE {input} {output} {params.settings} 2>{log}"

rule decompress:
# unzip fastq.gz for BWA
	input:
		trimdir + "{sample}.fastq.gz"
	output:
		temp(trimdir + "{sample}.fastq")
	conda:
		"../envs/py3.yml"
	shell:
		"gunzip --keep {input}"

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

rule bwa_map:
# run bwa aln to find the SA coordinates of the input reads (.sai file)
# for paired-end, this part needs to be modified
	input:
		trimdir + "{sample}.fastq"
	output:
	    temp(bamdir + "{sample}.sai") #bam
	log:
		logdir + "bwa/{sample}.log"
	threads: 8
	params:
		ref=config["refgenome"]
	conda:
		"../envs/py3.yml"
	shell:
		"bwa aln -t {threads} {params.ref} {input} 2>{log} >{output} "

rule bwa_samse:
# generate alignments from .sai file in the .sam format
# for paired-end, this part needs to be modified to use bwa_sampe
	input:
		bamdir + "{sample}.sai", 
		bamdir + "{sample}.fastq"
	output:
		temp(bamdir + "{sample}.sam")
	log:
		logdir + "bwa_samse/{sample}.log"
	params:
		ref=config["refgenome"]
	conda:
		"../envs/py3.yml"
	shell:
		"bwa samse -n 50 {params.ref} {input} 2>{log} >{output}"

def get_replicate_pairs(wildcards):
    # returns all pairwise combinations of replicates for a given id, if those exist
    # returns all pairwise combinations of replicates for a given id, if those exist
    # THIS NEEDS TO BE IN A SEPARATE FUNCTION
    # GENERATING PAIRS AS INPUT DOESNT WORK BECAUSE THOSE HAVE TO BE THE OUTPUT OF SOMETHING
    if(are_there_replicates(wildcards.id)):
        replicate_list = control_info.loc[wildcards.id]["idrep"]
        pairs = itertools.combinations(replicate_list, 2)
        return expand("test/{id}-{x[0]}.narrowPeak test/{id}-{x[1]}.narrowPeak", id = wildcards.id, x = pairs)
    else:
        raise ValueError("There are no replicates for this ID.")

rule proto_comparereplicatecombinations:
    # equivalent to calling common peaks
    # probably need a separate script for this
    # output: {name}-commonpeaks.txt
    # THIS NEEDS TO BE IN A SEPARATE FUNCTION
    # GENERATING PAIRS AS INPUT DOESNT WORK BECAUSE THOSE HAVE TO BE THE OUTPUT OF SOMETHING
    input:
        pairs = get_replicate_pairs
    output:
        "test/{id}.pairs.bed"
    shell:
        "cat {input.pairs} > {output}"

rule bigwigmerge:
# merge bigwig files into bedgraph
# NOTE: mergedbigwig for visualisation only, not recommended as analysis
	input:
		lambda x: map(lambda y: peaksdir + y + "_linearFE.bw", config["combined"][x.combined])
	output:
		temp(pairsdir + "{combined}_combined.bedGraph")
	conda:
		"../envs/py3.yml"
	shell:
		"bigWigMerge {input} {output}"

rule converttomultibigwig:
# convert merged bedgraph into bigwig file
# NOTE: mergedbigwig for visualisation only, not recommended as analysis
	input:
		pairsdir + "{combined}_combined.bedGraph"
	output:
		pairsdir + "{combined}_combined.bw"
	params:		
		ref=config["refchromsizes"]
	conda:
		"../envs/py3.yml"
	shell:
		"bedGraphToBigWig {input} {params.ref} {output}"

rule idr:
# calculate IDR statistic for the two replicates
	input:
		peakfiles=lambda x: map(lambda y: peaksdir + y + "_peaks.narrowPeak", config["combined"][x.combined]),
		overlap=pairsdir + "{combined}_commonpeaks.bed"
	output:
		pairsdir + "{combined}_idrValues.txt"
	log:
		logdir + "idr/{combined}.log"
	conda:
		"../envs/py3.yml"
	shell:
		"idr --samples {input.peakfiles} "
		"--output-file {output} "
		"--peak-list {input.overlap} --plot 2>{log}"

def get_replicate_list(wildcards):
    # returns a list of of replicates for a given id, if those exist
    # this might crash if there are no replicates (to debug)
    if(are_there_replicates(wildcards.id)):
        return expand(peaksdir + "{id}-{idrep}.narrowPeak", id = wildcards.id, idrep = control_info.loc[wildcards.id]["idrep"])
    else:
        raise ValueError("There are no replicates for this ID.")

rule getcommonpeaks:
	input:
		lambda x: map(lambda y: peaksdir + y + "_peaks.narrowPeak", config["combined"][x.combined])
	output:
		a=temp(pairsdir + "{combined}_commonpeaks1.bed"),
		b=temp(pairsdir + "{combined}_commonpeaks2.bed")

	shell:
		"bedtools intersect -a {input[0]} -b {input[1]} -wa > {output.a}; "
		"bedtools intersect -a {input[1]} -b {input[0]} -wa > {output.b}"


rule extendcommonpeaks:
# extend common peaks between two replicates
# note: used 'count' for summit (column 10), because later steps need this to be an integer
	input:
		a=pairsdir + "{combined}_commonpeaks1.bed",
		b=pairsdir + "{combined}_commonpeaks2.bed"
	output:
		pairsdir + "{combined}_commonpeaks.bed"
	conda:
		"../envs/py3.yml"
	shell:
		"cat {input} | bedtools sort | "
		"bedtools merge -c 4,5,6,7,8,9,10 -o collapse,mean,collapse,mean,mean,mean,count | "
		"""awk '$6="."' FS="\t" OFS="\t" """
		"> {output}"

