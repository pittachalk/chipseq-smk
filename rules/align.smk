def get_trimmed_fq(wildcards):
    if is_single_end(wildcards.name, wildcards.replic, wildcards.lane):
        return trimdir + "{name}-rep{replic}-{lane}.fq.gz".format(name = wildcards.name, replic = wildcards.replic, lane = wildcards.lane)
    else:
        return expand(trimdir + "{name}-rep{replic}-{lane}.{pair}.fq.gz", 
            name = wildcards.name, replic = wildcards.replic, lane = wildcards.lane, pair = [1, 2])

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

def get_bam_lanes(wildcards):
    lanes   = fastq_info[(fastq_info["name"] == wildcards.name) &  (fastq_info["rep"] == wildcards.replic)]["lane"]
    return expand(bamdir + "{name}-rep{replic}-{lane}.sorted.bam", name = wildcards.name, replic = wildcards.replic, lane = lanes)

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

rule convert_bam_to_bedgraph:
    # convert merged bam files into a bedgraph (needed for SEACR)
    # also, bedgraph format is needed by igv to convert to TDF
    # need to add scaling factors later on
    input:
        bamdir + "{name}-rep{replic}.merged.sorted.bam"
    output:
        bamdir + "{name}-rep{replic}.merged.sorted.bedgraph"
    conda:
        "../envs/analyzepeak.yml"
    shell:
        "bamCoverage -b {input} -o {output} -of bedgraph --scaleFactor 1"

rule igv_to_tdf:
    # convert bedgraph to TDF for IGV
    input: 
        bamdir + "{name}-rep{replic}.merged.sorted.bedgraph"
    output:
        bamdir + "{name}-rep{replic}.merged.sorted.tdf"
    params:
        ref=config["refgenome"]
    conda:
        "../envs/igv.yml"
    shell:
        "igvtools toTDF {input} {output} {params.ref}"


def get_treatvscontrol_bedgraph(wildcards):
    # this should give a single row describing the treatment vs control pairing
    x = control_info[(control_info["id"] == wildcards.id) & (control_info["idrep"] == wildcards.idrep)]

    treat = expand(bamdir + "{treatname}-rep{treatrep}.merged.sorted.bedgraph",
        treatname = x["treatname"], treatrep = x["treatrep"])
    control = expand(bamdir + "{controlname}-rep{controlrep}.merged.sorted.bedgraph", 
        controlname = x["controlname"], controlrep = x["controlrep"])
    
    assert (len(treat) == 1), "treated list not length 1, are id-idrep pairs unique?"
    assert (len(control) == 1), "control list not length 1, are id-idrep pairs unique?"
    return (treat[0], control[0])

rule seacr:
    # call peaks using SEACR
    input:
        bedgraph = get_treatvscontrol_bedgraph
    output:
        relaxed = peaksdir + "{id}-rep{idrep}.relaxed.bed"
    log:
        logdir + "seacr/{id}-{idrep}-relaxed.log"
    shell:
        "bash script/SEACR_1.3.sh {input.bedgraph[0]} {input.bedgraph[1]} 'norm' 'relaxed' {output.relaxed} 2> {log}; "
        "mv {output.relaxed}.relaxed.bed {output.relaxed}"
# should add a rename option
# also, the log is not redicreting for some reason
