######################################################################
######################################################################
#     Summary files
######################################################################

def get_replicate_bws(wildcards):
    # returns a list of bw of replicates for a given id, if those exist
    # this might crash if there are no replicates (to debug)
    if(are_there_replicates(wildcards.id)):
        return expand(bwdir + "{id}-{idrep}_linearFE.bw", id = wildcards.id, idrep = control_info.loc[wildcards.id]["idrep"])
    else:
        raise ValueError("There are no replicates for this ID.")


rule getoverallunion:
# get the union of peaks between all samples
# note: used 'count' for summit (column 10), because later steps need this to be an integer
	input:
		expand(summarydir + "{id}.commonpeaks.bed", id = control_info["id"].unique())
	output:
		overalldir + "overall.unionpeaks.bed"
	conda:
		"../envs/py3.yml"
	shell:
		"cat {input} | bedtools sort | bedtools merge -c 4,5,6,7,8,9,10 -o mean,mean,mean,mean,mean,mean,count | "
		"""awk '$6="."' FS="\t" OFS="\t" """
		"> {output}"

rule multibigwigsummary:
	# computes the average scores for every bw in every genomic region
	# for the entire genome (bins mode), and for consensus peaks
	input:
		bwfiles = expand(bwdir + "{id}-{idrep}_linearFE.bw", zip, id=control_info["id"], idrep=control_info["idrep"]),
		overlap = overalldir + "overall.unionpeaks.bed"
	output:
		npzbins = overalldir + "overallbw-bins.npz",
		tabbins = overalldir + "overallbw-bins.tab",
		npzpeak = overalldir + "overallbw-peak.npz", 
		tabpeak = overalldir + "overallbw-peak.tab"
	conda:
		"../envs/py3.yml"
	shell:
		"multiBigwigSummary bins -b {input.bwfiles} -o {output.npzbins} "
		"--outRawCounts {output.tabbins}; "
		"multiBigwigSummary BED-file -b {input.bwfiles} -o {output.npzpeak} "
		"--outRawCounts {output.tabpeak} --BED {input.overlap}"

rule pcacorr:
# plot PCA and correlation plots for replicates of ALL samples
	input:
		npzbins = overalldir + "overallbw-bins.npz",
		npzpeak = overalldir + "overallbw-peak.npz"
	output:
		pcabins  = overalldir + "overallbw-bins-pca.png",
		corrbins = [overalldir + "overallbw-bins-corr-heatmap" + f for f in [".png", ".tab"]],
		pcapeak  = overalldir + "overallbw-peak-pca.png",
		corrpeak = [overalldir + "overallbw-peak-corr-heatmap" + f for f in [".png", ".tab"]]
	conda:
		"../envs/py3.yml"
	shell:
		"plotPCA -in {input.npzbins} -o {output.pcabins} -T 'PCA: whole genome bins'; "
		
		"plotCorrelation -in {input.npzbins} --corMethod spearman --skipZeros "
		"--plotTitle 'Spearman Correlation: whole genome bins' "
		"--whatToPlot heatmap --colorMap RdYlBu --plotNumbers "
		"-o {output.corrbins[0]} --outFileCorMatrix {output.corrbins[1]}; "

		"plotPCA -in {input.npzpeak} -o {output.pcapeak} -T 'PCA: union of sample peaks'; "
		
		"plotCorrelation -in {input.npzpeak} --corMethod spearman --skipZeros "
		"--plotTitle 'Spearman Correlation: union of sample peaks' "
		"--whatToPlot heatmap --colorMap RdYlBu --plotNumbers "
		"-o {output.corrpeak[0]} --outFileCorMatrix {output.corrpeak[1]}"

rule computematrixall:
# calculate scores per genome regions across all samples
	input:
		bwfiles = expand(bwdir + "{id}-{idrep}_linearFE.bw", zip, id=control_info["id"], idrep=control_info["idrep"]),
		overlap = overalldir + "overall.unionpeaks.bed"
	output:
		gzipped = overalldir + "overallbw-peaks-matrix.mat.gz",
		tab = overalldir + "overallbw-peaks-matrix.tab"
	conda:
		"../envs/py3.yml"
	shell:
		"computeMatrix scale-regions "
		"-S {input.bwfiles} -R {input.overlap} "
		"--beforeRegionStartLength 3000 --afterRegionStartLength 3000 "
		"--regionBodyLength 5000 --skipZeros "
		"-o {output.gzipped} --outFileNameMatrix {output.tab}"

rule plotheatmapall:
# create heatmap for scores associated with the union of peak regions
	input:
		overalldir + "overallbw-peaks-matrix.mat.gz"
	output:
		overalldir + "overallbw-peaks-matrix-heatmap.png"
	conda:
		"../envs/py3.yml"
	shell:
		"plotHeatmap -m {input} -out {output} "
		"--heatmapHeight 12 --heatmapWidth 8 --colorMap RdBu "
		"--startLabel 'peak start' --endLabel 'peak end' "
		"--yAxisLabel 'average signal' --xAxisLabel ' ' "
		"--legendLocation none  --kmeans 3"

