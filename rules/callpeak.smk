######################################################################
######################################################################
#     ChIP-seq: MACS2
######################################################################

rule macs2:
# call peaks with MACS2 2.1.2
# this needs to be run in a Python 2 environment
	input:
		sample = bamdir + "{id}_sorted.bam",
		control = lambda x: bamdir + config["ids"][x.id] + "_sorted.bam"
	output:
		map(lambda macs2outfile: macs2dir + "{id}_" + macs2outfile, 
			["peaks.narrowPeak", "treat_pileup.bdg", "control_lambda.bdg"])
	log:
		logdir + "macs2/{id}.log"
	params:
		settings=config["macs2"]["settings"]
	conda:
		"envs/py2.yml"
	shell:
		"macs2 callpeak -t {input.sample} -c {input.control} "
		"--name {wildcards.id} --outdir " + macs2dir + " "
		"{params.settings} 2>{log}"

rule bedgraph:
# prepare bedgraph of linear fold enrichment and log10 likelihood
# this needs to be run in a Python 2 environment
	input: 
		sample = macs2dir + "{id}_treat_pileup.bdg",
		control = macs2dir + "{id}_control_lambda.bdg"
	output:
		FE = macs2dir + "{id}_linearFE.bdg",
		logLR = macs2dir + "{id}_logLR.bdg"
	conda:
		"envs/py2.yml"
	shell:
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output.FE} -m FE; "
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output.logLR} -m logLR --pseudocount 0.00001"


