######################################################################
######################################################################
#     Compare individual pairs of replicates
######################################################################

rule getcommonpeaks:
# find common peaks between two replicates (it doesn't work for >2 atm)
	input:
		lambda x: map(lambda y: macs2dir + y + "_peaks.narrowPeak", config["combined"][x.combined])
	output:
		a=temp(combineddir + "{combined}_commonpeaks1.bed"),
		b=temp(combineddir + "{combined}_commonpeaks2.bed")
	conda:
		"../envs/py3.yml"
	shell:
		"bedtools intersect -a {input[0]} -b {input[1]} -wa > {output.a}; "
		"bedtools intersect -a {input[1]} -b {input[0]} -wa > {output.b}"

rule extendcommonpeaks:
# extend common peaks between two replicates
# note: used 'count' for summit (column 10), because later steps need this to be an integer
	input:
		a=combineddir + "{combined}_commonpeaks1.bed",
		b=combineddir + "{combined}_commonpeaks2.bed"
	output:
		combineddir + "{combined}_commonpeaks.bed"
	conda:
		"../envs/py3.yml"
	shell:
		"cat {input} | bedtools sort | "
		"bedtools merge -c 4,5,6,7,8,9,10 -o collapse,mean,collapse,mean,mean,mean,count | "
		"""awk '$6="."' FS="\t" OFS="\t" """
		"> {output}"

rule idr:
# calculate IDR statistic for the two replicates
	input:
		peakfiles=lambda x: map(lambda y: macs2dir + y + "_peaks.narrowPeak", config["combined"][x.combined]),
		overlap=combineddir + "{combined}_commonpeaks.bed"
	output:
		combineddir + "{combined}_idrValues.txt"
	log:
		logdir + "idr/{combined}.log"
	conda:
		"../envs/py3.yml"
	shell:
		"idr --samples {input.peakfiles} "
		"--output-file {output} "
		"--peak-list {input.overlap} --plot 2>{log}"

rule computematrixbysample:
# calculate scores per genome regions and prepares an intermediate file for plotHeatmap and plotProfiles
	input:
		bwfiles=lambda x: map(lambda y: macs2dir + y + "_linearFE.bw", config["combined"][x.combined]),
		overlap=combineddir + "{combined}_commonpeaks.bed"
	output:
		gzipped=temp(combineddir + "{combined}-peaks-matrix.mat.gz")
	conda:
		"../envs/py3.yml"
	shell:
		"computeMatrix scale-regions "
		"-S {input.bwfiles} -R {input.overlap} "
		"--beforeRegionStartLength 3000 --afterRegionStartLength 3000 "
		"--regionBodyLength 5000 --skipZeros "
		"-o {output.gzipped}"

rule plotheatmapbysample:
# create heatmap for scores associated with common peaks (by sample)
# for quality check --- see consistency of peak profiles
	input:
		combineddir + "{combined}-peaks-matrix.mat.gz"
	output:
		combineddir + "{combined}-peaks-matrix-heatmap.png"
	conda:
		"../envs/py3.yml"
	shell:
		"plotHeatmap -m {input} -out {output} "
		"--heatmapHeight 12 --heatmapWidth 8 --colorMap RdBu "
		"--startLabel 'peak start' --endLabel 'peak end' "
		"--yAxisLabel 'average signal' --xAxisLabel ' ' "
		"--legendLocation none  --kmeans 3"


