def is_single_end(name, replic, lane):
    """
    Check if there is a second fastq file for the given sample.
    """
    return pd.isnull(fastq_info.loc[(name, replic, lane), "fq2"])

def get_fastq(wildcards):
    return indir + fastq_info.loc[(wildcards.name, wildcards.replic, wildcards.lane), ["fq1", "fq2"]].dropna()

rule trimSE:
    input:
        get_fastq
    output:
        fq = trimdir + "{name}-rep{replic}-{lane}.fq.gz"
    log:
        logdir + "trimmomatic/{name}-rep{replic}-{lane}.log"
    params:
        settings = config["trimmomatic"]["settings"]
    threads: 8
    shell:
        "trimmomatic SE -threads {threads} {input} {output} {params.settings} 2>{log}"

rule trimPE:
    input:
        get_fastq
    output:
        fq1 = trimdir + "{name}-rep{replic}-{lane}.1.fq.gz",
        fq2 = trimdir + "{name}-rep{replic}-{lane}.2.fq.gz",
        fq1_unpaired = trimdir + "{name}-rep{replic}-{lane}.1.unpaired.fq.gz",
        fq2_unpaired = trimdir + "{name}-rep{replic}-{lane}.2.unpaired.fq.gz"
    log:
        logdir + "trimmomatic/{name}-rep{replic}-{lane}.log"
    params:
        settings=config["trimmomatic"]["settings"]
    threads: 8
    shell:
        "trimmomatic PE -threads {threads} {input} {output.fq1} {output.fq1_unpaired} {output.fq2} {output.fq2_unpaired} {params.settings} 2>{log}"

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
        "../envs/py3.yml"
    shell:
        "bwa mem -t {threads} {params.ref} {input} 2>{log} >{output}"

rule sort_bam_lanes:
    # get alignment stat with flagstat, sort SAM file, convert to BAM, index BAM file
	input:
		bamdir + "{name}-rep{replic}-{lane}.sam"
	output:
		sortedbam = bamdir + "{name}-rep{replic}-{lane}.sorted.bam"
	conda:
		"../envs/py3.yml"
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
    shell:
        "samtools merge {output} {input}"

rule sortindex_mergedbam:
    # get alignment stat with flagstat, sort SAM file, convert to BAM, index BAM file
    input:
        bamdir + "{name}-rep{replic}.merged.bam"
    output:
        flagstat = bamdir + "{name}-rep{replic}.merged.alignstat.txt",
        sortedbam = bamdir + "{name}-rep{replic}.merged.sorted.bam",
        bai = bamdir + "{name}-rep{replic}.merged.sorted.bam.bai"
    conda:
        "../envs/py3.yml"
    shell:
        "samtools flagstat {input} > {output.flagstat}; "
        "samtools view -bS {input} | "
        "samtools sort - -o {output.sortedbam}"

def get_treatvscontrol(wildcards):
    # this should give a single row describing the treatment vs control pairing
    x = control_info[(control_info["id"] == wildcards.id) & (control_info["idrep"] == wildcards.idrep)]

    treat = expand(bamdir + "{treatname}-rep{treatrep}.merged.bam",
        treatname = x["treatname"], treatrep = x["treatrep"])
    control = expand(bamdir + "{controlname}-rep{controlrep}.merged.bam", 
        controlname = x["controlname"], controlrep = x["controlrep"])
    
    assert (len(treat) == 1), "treated list not length 1, are id-idrep pairs unique?"
    assert (len(control) == 1), "control list not length 1, are id-idrep pairs unique?"

    return (treat[0], control[0])


rule macs2:
    # call peaks with MACS2 2.1.2, run in a Python 2 environment
    input:
        bam = get_treatvscontrol
    output:
        peaksdir + "{id}-{idrep}.narrowPeak", #change to include bdg too
        peaksdir + "{id}-{idrep}_treat_pileup.bdg",
        peaksdir + "{id}-{idrep}_control_lambda.bdg"
    log:
        logdir + "macs2/{id}-{idrep}.log"
    params:
        settings=config["macs2"]["settings"]
    conda:
        "../envs/macs.yml"
    shell:
        "macs2 callpeak -t {input.bam[0]} -c {input.bam[1]} "
        "--name {wildcards.id}-{wildcards.idrep} --outdir " + peaksdir + " "
        "{params.settings} 2>{log}"

rule bedgraph:
    # prepare bedgraph of fold enrichment, run in a Python 2 environment
	input: 
		sample = peaksdir + "{id}-{idrep}_treat_pileup.bdg",
		control = peaksdir + "{id}-{idrep}_control_lambda.bdg"
	output:
		FE = peaksdir + "{id}-{idrep}_linearFE.bdg",
		logLR = peaksdir + "{id}-{idrep}_logLR.bdg"
	conda:
		"../envs/macs.yml"
	shell:
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output.FE} -m FE; "
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output.logLR} -m logLR --pseudocount 0.00001"



# def are_there_replicates(id):
#     # check if there are replicates for a given id
#     return (len(control_info.loc[id]) > 1)

# def get_replicate_list(wildcards):
#     # returns a list of of replicates for a given id, if those exist
#     # CURRENTLY THIS DOESNT WORK IF THERE ARE INDEED NO REPLICATES
#     # I THINK PE AND SE WORKS BECAUSE THEY COMPLEMENT EACH OTHER
#     # this function might not even be necessary anymore, can merge with proto_combineallreplicates
#     if(are_there_replicates(wildcards.id)):
#         return expand(peaksdir + "{id}-{idrep}.narrowPeak", id = wildcards.id, idrep = control_info.loc[wildcards.id]["idrep"])
#     else:
#         raise ValueError("There are no replicates for this ID.")

# rule proto_combineallreplicates:
#     # equivalent to calling common peaks
#     # output: {name}-commonpeaks.txt
#     input:
#         replicates = get_replicate_list
#     output:
#         summarydir + "{id}.commonpeaks.bed"
#     shell:
#         "cat a {input.replicates[0]} b {input.replicates[1]} > {output}"



# def get_replicate_pairs(wildcards):
#     # returns all pairwise combinations of replicates for a given id, if those exist
#     # returns all pairwise combinations of replicates for a given id, if those exist
#     # THIS NEEDS TO BE IN A SEPARATE FUNCTION
#     # GENERATING PAIRS AS INPUT DOESNT WORK BECAUSE THOSE HAVE TO BE THE OUTPUT OF SOMETHING
#     if(are_there_replicates(wildcards.id)):
#         replicate_list = control_info.loc[wildcards.id]["idrep"]
#         pairs = itertools.combinations(replicate_list, 2)
#         return expand("test/{id}-{x[0]}.narrowPeak test/{id}-{x[1]}.narrowPeak", id = wildcards.id, x = pairs)
#     else:
#         raise ValueError("There are no replicates for this ID.")

# rule proto_comparereplicatecombinations:
#     # equivalent to calling common peaks
#     # probably need a separate script for this
#     # output: {name}-commonpeaks.txt
#     # THIS NEEDS TO BE IN A SEPARATE FUNCTION
#     # GENERATING PAIRS AS INPUT DOESNT WORK BECAUSE THOSE HAVE TO BE THE OUTPUT OF SOMETHING
#     input:
#         pairs = get_replicate_pairs
#     output:
#         "test/{id}.pairs.bed"
#     shell:
#         "cat {input.pairs} > {output}"