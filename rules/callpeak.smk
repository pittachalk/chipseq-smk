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

def get_treatvscontrol_mergedbam(wildcards):
    # this should give a single row describing the treatment vs control pairing
    x = control_info[(control_info["id"] == wildcards.id) & (control_info["idrep"] == wildcards.idrep)]

    treat = expand(bamdir + "{treatname}-rep{treatrep}.merged.sorted.bam",
        treatname = x["treatname"], treatrep = x["treatrep"])
    control = expand(bamdir + "{controlname}-rep{controlrep}.merged.sorted.bam", 
        controlname = x["controlname"], controlrep = x["controlrep"])

    assert (len(treat) == 1), "treated list not length 1, are id-idrep pairs unique?"
    assert (len(control) == 1), "control list not length 1, are id-idrep pairs unique?"
    return (treat[0], control[0])

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

rule macs2:
    # call peaks with MACS2 2.1.2, run in a Python 2 environment
    input:
        bam = get_treatvscontrol_mergedbam
    output:
        peaksdir + "{id}-rep{idrep}.narrowPeak",
        peaksdir + "{id}-rep{idrep}_treat_pileup.bdg",
        peaksdir + "{id}-rep{idrep}_control_lambda.bdg"
    log:
        logdir + "macs2/{id}-rep{idrep}.log"
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
		sample = peaksdir + "{id}-rep{idrep}_treat_pileup.bdg",
		control = peaksdir + "{id}-rep{idrep}_control_lambda.bdg"
	output:
		peaksdir + "{id}-rep{idrep}_linearFE.bdg"
	conda:
		"../envs/macs.yml"
	shell:
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output} -m FE"