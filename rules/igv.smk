######################################################################
######################################################################
#     Convert bedgraphs to the TDF and BigWig binary formats
######################################################################

rule igvsort:
# sort bedgraph
	input: 
		FE = macs2dir + "{id}_linearFE.bdg",
		logLR = macs2dir + "{id}_logLR.bdg"
	output: 
		FE = temp(tempdir + "{id}_linearFE_sorted.bdg"),  # temp this later on!
		logLR = temp(tempdir + "{id}_logLR_sorted.bdg")
	conda:
		"../envs/py3.yml"
	shell:
		"igvtools sort {input.FE} {output.FE}; "
		"igvtools sort {input.logLR} {output.logLR}"

rule igvtotdf:
# convert bedgraph to TDF for IGV
	input: 
		FE = tempdir + "{id}_linearFE_sorted.bdg",
		logLR = tempdir + "{id}_logLR_sorted.bdg"
	output: 
		FE = macs2dir + "{id}_linearFE_sorted.tdf",
		logLR = macs2dir + "{id}_logLR_sorted.tdf"
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
		tempdir + "{id}_linearFE_sorted.bdg"
	output: 
		macs2dir + "{id}_linearFE.bw"
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
		lambda x: map(lambda y: macs2dir + y + "_linearFE.bw", config["combined"][x.combined])
	output:
		temp(tempdir + "{combined}_combined.bedGraph")
	conda:
		"../envs/py3.yml"
	shell:
		"bigWigMerge {input} {output}"

rule converttomultibigwig:
# convert merged bedgraph into bigwig file
# NOTE: mergedbigwig for visualisation only, not recommended as analysis
	input:
		tempdir + "{combined}_combined.bedGraph"
	output:
		combineddir + "{combined}_combined.bw"
	params:		
		ref=config["refchromsizes"]
	conda:
		"../envs/py3.yml"
	shell:
		"bedGraphToBigWig {input} {params.ref} {output}"
