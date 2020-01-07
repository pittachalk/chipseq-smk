######################################################################
######################################################################
#     Summary files
######################################################################

rule compilepeakunion:
# get the union of peaks between all samples
# note: used 'count' for summit (column 10), because later steps need this to be an integer
	input:
		map(lambda x: combineddir + x + "_commonpeaks.bed", config["combined"])
	output:
		summarydir + "summary-unionpeaks.bed"
	shell:
		"cat {input} | bedtools sort | bedtools merge -c 4,5,6,7,8,9,10 -o mean,mean,mean,mean,mean,mean,count | "
		"""awk '$6="."' FS="\t" OFS="\t" """
		"> {output}"

rule multibigwigsummary:
# computes the average scores for each replicate in every genomic region
# for the entire genome (bins mode), and for consensus peaks
	input:
		bwfiles=map(lambda x: macs2dir + x + "_linearFE.bw", config["ids"]),
		overlap=summarydir + "summary-unionpeaks.bed"
	output:
		npzbins=summarydir + "summarybw-bins.npz",
		tabbins=summarydir + "summarybw-bins.tab",
		npzpeak=summarydir + "summarybw-peak.npz", 
		tabpeak=summarydir + "summarybw-peak.tab"
	shell:
		"multiBigwigSummary bins -b {input.bwfiles} -o {output.npzbins} "
		"--outRawCounts {output.tabbins}; "
		"multiBigwigSummary BED-file -b {input.bwfiles} -o {output.npzpeak} "
		"--outRawCounts {output.tabpeak} --BED {input.overlap}"

rule pcacorr:
# plot PCA and correlation plots for replicates of ALL samples
	input:
		npzbins=summarydir + "summarybw-bins.npz",
		npzpeak=summarydir + "summarybw-peak.npz"
	output:
		pcabins=summarydir + "summarybw-bins-pca.png",
		corrbins=[summarydir + "summarybw-bins-corr-heatmap" + f for f in [".png", ".tab"]],
		pcapeak=summarydir + "summarybw-peak-pca.png",
		corrpeak=[summarydir + "summarybw-peak-corr-heatmap" + f for f in [".png", ".tab"]],
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
		bwfiles=map(lambda x: macs2dir + x + "_linearFE.bw", config["ids"]),
		overlap=summarydir + "summary-unionpeaks.bed"
	output:
		gzipped=summarydir + "summarybw-peaks-matrix.mat.gz",
		tab=summarydir + "summarybw-peaks-matrix.tab"
	shell:
		"computeMatrix scale-regions "
		"-S {input.bwfiles} -R {input.overlap} "
		"--beforeRegionStartLength 3000 --afterRegionStartLength 3000 "
		"--regionBodyLength 5000 --skipZeros "
		"-o {output.gzipped} --outFileNameMatrix {output.tab}"

rule plotheatmapall:
# create heatmap for scores associated with the union of peak regions
	input:
		summarydir + "summarybw-peaks-matrix.mat.gz"
	output:
		summarydir + "summarybw-peaks-matrix-heatmap.png"
	shell:
		"plotHeatmap -m {input} -out {output} "
		"--heatmapHeight 12 --heatmapWidth 8 --colorMap RdBu "
		"--startLabel 'peak start' --endLabel 'peak end' "
		"--yAxisLabel 'average signal' --xAxisLabel ' ' "
		"--legendLocation none  --kmeans 3"

