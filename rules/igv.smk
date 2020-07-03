rule igv_convert_bedgraph_to_tdf:
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

rule igv_sort_bdgcmp:
# sort bedgraph
	input: 
		FE = peaksdir + "{id}-{idrep}_linearFE.bdg",
		logLR = peaksdir + "{id}-{idrep}_logLR.bdg"
	output: 
		FE = temp(bwdir + "{id}-{idrep}_linearFE_sorted.bdg"), 
		logLR = temp(bwdir + "{id}-{idrep}_logLR_sorted.bdg")
	conda:
		"../envs/igv.yml"
	shell:
		"igvtools sort {input.FE} {output.FE}; "
		"igvtools sort {input.logLR} {output.logLR}"

rule igv_convert_bdgcmp_to_tdf:
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
		"../envs/igv.yml"
	shell:
		"igvtools toTDF {input.FE} {output.FE} {params.ref}; "
		"igvtools toTDF {input.logLR} {output.logLR} {params.ref}"


rule bedgraphtobigwig:
# UNNEDDED FOR NOT, REMEMBER TO ADD "REP" BEFORE IDREP
# convert bedgraph to bigwig for deeptools
# this is only done for the linear fold enrichment bedgraph
	input: 
		bwdir + "{id}-{idrep}_linearFE_sorted.bdg"
	output: 
		bwdir + "{id}-{idrep}_linearFE.bw"
	params:
		ref=config["refchromsizes"]
	conda:
		"../envs/analyzepeak.yml"
	shell:
		"bedGraphToBigWig {input} {params.ref} {output}"


