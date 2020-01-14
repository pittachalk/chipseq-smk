######################################################################
######################################################################
#     Compare individual pairs of replicates
######################################################################

def are_there_replicates(id):
    # check if there are replicates for a given id/species
    return (len(control_info.loc[id]) > 1)

def get_first_sample(wildcards):
    # returns the first sample from a list of of replicates for a given id, if those exist
    # this might crash if there are no replicates (to debug)
    if(are_there_replicates(wildcards.id)):
        replist = expand(peaksdir + "{id}-{idrep}.narrowPeak", id = wildcards.id, idrep = control_info.loc[wildcards.id]["idrep"])
        return replist[0]
    else:
        raise ValueError("There are no replicates for this ID.")

def get_other_replicates(wildcards):
    # returns a list of replicates excluding the first for a given id, if those exist
    # this might crash if there are no replicates (to debug)
    if(are_there_replicates(wildcards.id)):
        replist = expand(peaksdir + "{id}-{idrep}.narrowPeak", id = wildcards.id, idrep = control_info.loc[wildcards.id]["idrep"])
        replist.pop(0)
        return replist
    else:
        raise ValueError("There are no replicates for this ID.")

def get_replicate_list(wildcards):
    # returns a list of of replicates for a given id, if those exist
    # this might crash if there are no replicates (to debug)
    if(are_there_replicates(wildcards.id)):
        return expand(peaksdir + "{id}-{idrep}.narrowPeak", id = wildcards.id, idrep = control_info.loc[wildcards.id]["idrep"])
    else:
        raise ValueError("There are no replicates for this ID.")

rule getsampleintersect:
    # find common peaks between all replicates for a given id
    input:
        first = get_first_sample,
        others = get_other_replicates
    output:
        summarydir + "{id}.intersect.bed"
    conda:
        "../envs/py3.yml"
    shell:
        "bedtools intersect -a {input.first} -b {input.others} -wa > {output}"

rule getsampleunion:
    # extend the common peaks
    input:
        get_replicate_list
    output:
        summarydir + "{id}.commonpeaks.bed"
    conda:
        "../envs/py3.yml"
    shell:
        "cat {input} | bedtools sort | "
        "bedtools merge -c 4,5,6,7,8,9,10 -o mean,mean,mean,mean,mean,mean,count | "
		"""awk '$6="."' FS="\t" OFS="\t" """
		"> {output}"

def get_replicate_bws(wildcards):
    # returns a list of bw of replicates for a given id, if those exist
    # this might crash if there are no replicates (to debug)
    if(are_there_replicates(wildcards.id)):
        return expand(bwdir + "{id}-{idrep}_linearFE.bw", id = wildcards.id, idrep = control_info.loc[wildcards.id]["idrep"])
    else:
        raise ValueError("There are no replicates for this ID.")

rule computematrixbysample:
# calculate scores per genome regions and prepares an intermediate file for plotHeatmap and plotProfiles
    input:
        bwfiles = get_replicate_bws,
        overlap = summarydir + "{id}.commonpeaks.bed"
    output:
        gzipped = temp(summarydir + "{id}-peaks-matrix.mat.gz"),
        tab = summarydir + "{id}-peaks-matrix.tab"
    conda:
        "../envs/py3.yml"
    shell:
        "computeMatrix scale-regions "
        "-S {input.bwfiles} -R {input.overlap} "
        "--beforeRegionStartLength 3000 --afterRegionStartLength 3000 "
        "--regionBodyLength 5000 --skipZeros "
        "-o {output.gzipped} --outFileNameMatrix {output.tab}"

rule plotheatmapbysample:
# create heatmap for scores associated with common peaks (by sample)
# for quality check --- see consistency of peak profiles
	input:
		summarydir + "{id}-peaks-matrix.mat.gz"
	output:
		summarydir + "{id}-peaks-matrix-heatmap.png"
	conda:
		"../envs/py3.yml"
	shell:
		"plotHeatmap -m {input} -out {output} "
		"--heatmapHeight 12 --heatmapWidth 8 --colorMap RdBu "
		"--startLabel 'peak start' --endLabel 'peak end' "
		"--yAxisLabel 'average signal' --xAxisLabel ' ' "
		"--legendLocation none  --kmeans 3"
