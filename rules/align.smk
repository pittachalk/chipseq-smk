def get_trimmed_fq(wildcards):
    if is_single_end(wildcards.name, wildcards.replic, wildcards.lane):
        return trimdir + "{name}-rep{replic}-{lane}.fq.gz".format(name = wildcards.name, replic = wildcards.replic, lane = wildcards.lane)
    else:
        return expand(trimdir + "{name}-rep{replic}-{lane}.{pair}.fq.gz", 
            name = wildcards.name, replic = wildcards.replic, lane = wildcards.lane, pair = [1, 2])

def get_bam_lanes(wildcards):
    lanes   = fastq_info[(fastq_info["name"] == wildcards.name) &  (fastq_info["rep"] == wildcards.replic)]["lane"]
    return expand(bamdir + "{name}-rep{replic}-{lane}.sorted.bam", name = wildcards.name, replic = wildcards.replic, lane = lanes)

rule bwamem:
    # bwa is run on each lane separately
    # this is the step where paired end reads are combined
    input:
        get_trimmed_fq
    output:
        bamdir + "{name}-rep{replic}-{lane}.sam"
    log:
        logdir + "bwa/{name}-rep{replic}-{lane}.log"
    threads: 8
    params:
        ref=config["refgenome"]
    conda:
        "../envs/align.yml"
    shell:
        "bwa mem -t {threads} {params.ref} {input} 2>{log} >{output}"

rule sort_bam_lanes:
    # get alignment stat with flagstat, sort SAM file, convert to BAM, index BAM file
	input:
		bamdir + "{name}-rep{replic}-{lane}.sam"
	output:
		sortedbam = bamdir + "{name}-rep{replic}-{lane}.sorted.bam"
	conda:
		"../envs/align.yml"
	shell:
		"samtools view -bS {input} | samtools sort - -o {output.sortedbam}"

rule merge_bam_lanes:
    # now we cat the bam files from separate lanes
    # this will be replaced with conversion to sam, so the cat is just placeholder
    # output: {name}.{replicate}-aligned.txt
    input:
        fq = get_bam_lanes
    output:
        bamdir + "{name}-rep{replic}.merged.bam"
    conda:
        "../envs/align.yml"
    shell:
        "samtools merge {output} {input}"

rule flagstat:
    # get alignment stat with flagstat
    input:
        bamdir + "{name}-rep{replic}.merged.bam"
    output:
        bamdir + "{name}-rep{replic}.merged.alignstat.txt"
    conda:
        "../envs/align.yml"
    shell:
        "samtools flagstat {input} > {output}"

rule sort_mergedbam:
    # sort SAM file, convert to BAM
    input:
        bamdir + "{name}-rep{replic}.merged.bam"
    output:
        bamdir + "{name}-rep{replic}.merged.sorted.bam"
    conda:
        "../envs/align.yml"
    shell:
        "samtools view -bS {input} | samtools sort - -o {output}"

rule index_mergedbam:
    # index sorted bam files
    # alignstat file (as input) is to "trick" Snakemake
    # so I do not need to specify them in rule all
    input:
        bam = bamdir + "{name}-rep{replic}.merged.sorted.bam",
        alignstat = bamdir + "{name}-rep{replic}.merged.alignstat.txt"
    output:
        bai = bamdir + "{name}-rep{replic}.merged.sorted.bam.bai"
    conda:
        "../envs/align.yml"
    shell:
        "samtools index {input.bam}"