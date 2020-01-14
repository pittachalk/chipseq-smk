######################################################################
######################################################################
#     Convert bedgraphs to the TDF and BigWig binary formats
######################################################################

rule igvsort:
# sort bedgraph
	input: 
		FE = peaksdir + "{id}-{idrep}_linearFE.bdg",
		logLR = peaksdir + "{id}-{idrep}_logLR.bdg"
	output: 
		FE = temp(bwdir + "{id}-{idrep}_linearFE_sorted.bdg"), 
		logLR = temp(bwdir + "{id}-{idrep}_logLR_sorted.bdg")
	conda:
		"../envs/py3.yml"
	shell:
		"igvtools sort {input.FE} {output.FE}; "
		"igvtools sort {input.logLR} {output.logLR}"

rule igvtotdf:
# convert bedgraph to TDF for IGV
	input: 
		FE = bwdir + "{id}-{idrep}_linearFE_sorted.bdg",
		logLR = bwdir + "{id}-{idrep}_logLR_sorted.bdg"
	output: 
		FE = bwdir + "{id}-{idrep}_linearFE_sorted.tdf",
		logLR = bwdir + "{id}-{idrep}_logLR_sorted.tdf"
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
		bwdir + "{id}-{idrep}_linearFE_sorted.bdg"
	output: 
		bwdir + "{id}-{idrep}_linearFE.bw"
	params:
		ref=config["refchromsizes"]
	conda:
		"../envs/py3.yml"
	shell:
		"bedGraphToBigWig {input} {params.ref} {output}"


