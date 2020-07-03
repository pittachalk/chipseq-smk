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
		peaksdir + "{id}-rep{idrep}_linearFE.bdg"
	output: 
		temp(peaksdir + "{id}-rep{idrep}_linearFE_sorted.bdg")
	conda:
		"../envs/igv.yml"
	shell:
		"igvtools sort {input} {output}"

rule igv_convert_bdgcmp_to_tdf:
# convert bedgraph to TDF for IGV
	input: 
		peaksdir + "{id}-rep{idrep}_linearFE_sorted.bdg"
	output: 
		peaksdir + "{id}-rep{idrep}_linearFE_sorted.tdf"
	params:
		ref=config["refgenome"]
	conda:
		"../envs/igv.yml"
	shell:
		"igvtools toTDF {input} {output} {params.ref}"
