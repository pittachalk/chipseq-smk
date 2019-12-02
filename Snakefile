configfile: "config.yaml"  # put quotes

# End users should not change anything below this line
# Parameters for the run should be modified in the CONFIG file

######################################################################
######################################################################
#     Setting up
######################################################################
# obtain directories from the CONFIG file
sampledir   = config["sampledir"]
outputdir   = config["outputdir"]
tempdir     = outputdir + config["subdir"]["tmp"]
logdir      = outputdir + config["subdir"]["log"]
qcdir       = outputdir + config["subdir"]["qc"]
bamdir      = outputdir + config["subdir"]["bam"]
macs2dir    = outputdir + config["subdir"]["macs2"]
summarydir  = macs2dir + config["subdir"]["summary"]  # subdir for all samples
combineddir = macs2dir + config["subdir"]["combined"] # subdir for combined replicates

######################################################################
######################################################################
#     The rule all
######################################################################

"""
Final output files
Note that the user can comment out what is not required
"""

rule all:
	input:
		# alignment files for individuals:
		# sorted BAM files and quality control reports
	    expand(["{outputdir}{bamdir}{sample}_sorted.bam",
	    		"{outputdir}{qc}{sample}_fastqc.html"], 
	    	sample=config["samples"], outputdir=config["outputdir"], 
	    	bamdir=config["subdir"]["bam"], qc = config["subdir"]["qc"]),

	    # MACS2 output for individuals: BED files for peaks, binary TDF for IGV
	    expand(["{outputdir}{macs2}{id}_peaks.narrowPeak",
				"{outputdir}{macs2}{id}_linearFE_sorted.tdf", 
				"{outputdir}{macs2}{id}_logLR_sorted.tdf"],
	    	outputdir=config["outputdir"], id=config["ids"], macs2=config["subdir"]["macs2"]),

	    # combined output from pairs of replicates:
	    # common peaks BED, IDR values, heatmap of peak profiles
	    expand(["{outputdir}{macs2}{combineddir}{combined}_commonpeaks.bed",
	    	    "{outputdir}{macs2}{combineddir}{combined}_idrValues.txt",
	    	    "{outputdir}{macs2}{combineddir}{combined}-peaks-matrix-heatmap.png"],
			outputdir=config["outputdir"], macs2=config["subdir"]["macs2"], 
			combineddir=config["subdir"]["combined"], combined=config["combined"] ),

	    # summary files for everything:
	    # union of all peaks BED, PCA, correlation and heatmap of peak profiles
	    expand(["{outputdir}{macs2}{summarydir}summary-unionpeaks.bed",
	    	    "{outputdir}{macs2}{summarydir}summarybw-peak-pca.png",
	    	    "{outputdir}{macs2}{summarydir}summarybw-peak-corr-heatmap.png",
	    	    "{outputdir}{macs2}{summarydir}summarybw-peaks-matrix-heatmap.png"],
			outputdir=config["outputdir"],
			macs2=config["subdir"]["macs2"], summarydir=config["subdir"]["summary"] )


######################################################################
######################################################################
#     Preprocessing of fastq.gz files and quality control
######################################################################

rule cat:
# concatenate fastq.gz files of each sample
	input:
		lambda x: map(lambda y: sampledir + y, config["samples"][x.sample])
	output:
		temp(tempdir + "{sample}_untrimmed.fastq.gz")
	shell:
		"cat {input} > {output}"

rule trim:
# trim Illumina adapters from fastq.gz
	input:
		tempdir + "{sample}_untrimmed.fastq.gz"
	output:
		temp(tempdir + "{sample}.fastq.gz")
	log:
		logdir + "trimmomatic/{sample}.log"
	params:
		bin=config["trimmomatic"]["bin"],
		settings=config["trimmomatic"]["settings"]
	shell:
		"trimmomatic SE {input} {output} {params.settings} 2>{log}"

rule decompress:
# unzip fastq.gz for BWA
	input:
		tempdir + "{sample}.fastq.gz"
	output:
		temp(tempdir + "{sample}.fastq")
	shell:
		"gunzip --keep {input}"

rule qctrim:
# run FastQC on trimmed fastq.gz files
	input:
		tempdir + "{sample}.fastq.gz"
	output:
		qcdir + "{sample}_fastqc.html",
		qcdir + "{sample}_fastqc.zip"
	log:
		logdir + "fastqc/{sample}.log"
	shell:
		"fastqc -o {qcdir} {input} 2>{log}"


######################################################################
######################################################################
#     BWA alignment and samtools
######################################################################

rule bwa_map:
# run bwa aln to find the SA coordinates of the input reads (.sai file)
# for paired-end, this part needs to be modified
	input:
		tempdir + "{sample}.fastq"
	output:
	    temp(tempdir + "{sample}.sai")
	log:
		logdir + "bwa/{sample}.log"
	threads: 8
	params:
		ref=config["refgenome"]
	shell:
		"bwa aln -t {threads} {params.ref} {input} 2>{log} >{output} "

rule bwa_samse:
# generate alignments from .sai file in the .sam format
# for paired-end, this part needs to be modified to use bwa_sampe
	input:
		tempdir + "{sample}.sai", 
		tempdir + "{sample}.fastq"
	output:
		temp(tempdir + "{sample}.sam")
	log:
		logdir + "bwa_samse/{sample}.log"
	params:
		ref=config["refgenome"]
	shell:
		"bwa samse -n 50 {params.ref} {input} 2>{log} >{output}"

rule samtools:
# get alignment stat with flagstat, sort SAM file, convert to BAM, index BAM file
	input:
		tempdir + "{sample}.sam"
	output:
		flagstat = qcdir + "{sample}_alignstat.txt",
		sortedbam = bamdir + "{sample}_sorted.bam",
		bai = bamdir + "{sample}_sorted.bam.bai"
	shell:
		"samtools flagstat {input} > {output.flagstat}; "
		"samtools view -bS {input} | "
		"samtools sort - -o {output.sortedbam}; "
		"samtools index {output.sortedbam}"


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
        "envs/py2.yaml"
	shell:
		"macs2 callpeak -t {input.sample} -c {input.control} "
		"--name {wildcards.id} --outdir " + macs2dir + " "
		"{params.settings} 2>{log}"

rule bedgraph:
# prepare bedgraph of linear fold enrichment and log10 likelihood
	input: 
		sample = macs2dir + "{id}_treat_pileup.bdg",
		control = macs2dir + "{id}_control_lambda.bdg"
	output:
		FE = macs2dir + "{id}_linearFE.bdg",
		logLR = macs2dir + "{id}_logLR.bdg"
	shell:
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output.FE} -m FE; "
		"macs2 bdgcmp -t {input.sample} -c {input.control} -o {output.logLR} -m logLR --pseudocount 0.00001"


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
		FE = tempdir + "{id}_linearFE_sorted.bdg",  # temp this later on!
		logLR = temp(tempdir + "{id}_logLR_sorted.bdg")
	params:
		bin=config["igvtools"]["bin"]
	shell:
		"{params.bin} sort {input.FE} {output.FE}; "
		"{params.bin} sort {input.logLR} {output.logLR}"

rule igvtotdf:
# convert bedgraph to TDF for IGV
	input: 
		FE = tempdir + "{id}_linearFE_sorted.bdg",
		logLR = tempdir + "{id}_logLR_sorted.bdg"
	output: 
		FE = macs2dir + "{id}_linearFE_sorted.tdf",
		logLR = macs2dir + "{id}_logLR_sorted.tdf"
	params:
		bin=config["igvtools"]["bin"],
		ref=config["refgenome"]
	shell:
		"{params.bin} toTDF {input.FE} {output.FE} {params.ref}; "
		"{params.bin} toTDF {input.logLR} {output.logLR} {params.ref}"

rule bedgraphtobigwig:
# convert bedgraph to bigwig for deeptools
# this is only done for the linear fold enrichment bedgraph
	input: 
		tempdir + "{id}_linearFE_sorted.bdg"
	output: 
		macs2dir + "{id}_linearFE.bw"
	params:
		ref=config["refchromsizes"]
	shell:
		"bedGraphToBigWig {input} {params.ref} {output}"

rule bigwigmerge:
# merge bigwig files into bedgraph
# NOTE: mergedbigwig for visualisation only, not recommended as analysis
	input:
		lambda x: map(lambda y: macs2dir + y + "_linearFE.bw", config["combined"][x.combined])
	output:
		temp(tempdir + "{combined}_combined.bedGraph")
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
	shell:
		"bedGraphToBigWig {input} {params.ref} {output}"


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
	shell:
		"plotHeatmap -m {input} -out {output} "
		"--heatmapHeight 12 --heatmapWidth 8 --colorMap RdBu "
		"--startLabel 'peak start' --endLabel 'peak end' "
		"--yAxisLabel 'average signal' --xAxisLabel ' ' "
		"--legendLocation none  --kmeans 3"


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
