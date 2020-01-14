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
