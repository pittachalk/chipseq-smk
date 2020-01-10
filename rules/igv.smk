######################################################################
######################################################################
#     Convert bedgraphs to the TDF and BigWig binary formats
######################################################################

rule igvsort:
# sort bedgraph
	input: 
		FE = peaksdir + "{id}_linearFE.bdg",
		logLR = peaksdir + "{id}_logLR.bdg"
	output: 
		FE = temp(peaksdir + "{id}_linearFE_sorted.bdg"),  # temp this later on!
		logLR = temp(peaksdir + "{id}_logLR_sorted.bdg")
	conda:
		"../envs/py3.yml"
	shell:
		"igvtools sort {input.FE} {output.FE}; "
		"igvtools sort {input.logLR} {output.logLR}"

rule igvtotdf:
# convert bedgraph to TDF for IGV
	input: 
		FE = peaksdir + "{id}_linearFE_sorted.bdg",
		logLR = peaksdir + "{id}_logLR_sorted.bdg"
	output: 
		FE = peaksdir + "{id}_linearFE_sorted.tdf",
		logLR = peaksdir + "{id}_logLR_sorted.tdf"
	params:
		ref=config["refgenome"]
	conda:
		"../envs/py3.yml"
	shell:
		"igvtools toTDF {input.FE} {output.FE} {params.ref}; "
		"igvtools toTDF {input.logLR} {output.logLR} {params.ref}"

rule bedgraphtobigwig:
# convert bedgraph to bigwig for deeptools
# this is only done for the linear fold enrichment bedgraph
	input: 
		peaksdir + "{id}_linearFE_sorted.bdg"
	output: 
		peaksdir + "{id}_linearFE.bw"
	params:
		ref=config["refchromsizes"]
	conda:
		"../envs/py3.yml"
	shell:
		"bedGraphToBigWig {input} {params.ref} {output}"

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
