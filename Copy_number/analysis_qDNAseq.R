#################################################################
# Authors: Friday
# Modified Date: 30-05-2023
# Version: 1.0

# This script is designed to run copy number analysis using QDNAseq.  It requires all the analysis-ready BAM files to be stored within the same folder of the working diretory and enter the path at "path=" at the readCounts step. For details, please refer to https://github.com/ccagc/QDNAseq

##Installation
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# install QDNAseq bin pre-build annotations adopted from QDNAseq.hg19
#devtools::install_github("asntech/QDNAseq.hg38@main")
#BiocManager::install(c("QDNAseq","future"))
#################################################################

library(QDNAseq)
library(QDNAseq.hg38)

# call wrapper package future version 1.32.0 for sequential and parallel processing
library(future)
require(parallelly)

# instantize parallelization of estimateCorrection(), segmentBins(), createBins()and calculateBlacklist(). binReadCounts(); default is to use all cores available; ask for more nodes from schedulers and turns it into a multinodal job during qsub to high-performance compute (HPC) clusters with e.g. #PBS -l nodes=16:ppn=1
future::plan("multisession", workers=parallelly::availableCores()/2)

# Set working directory
setwd("/lustre/alice3/scratch/spectre/y/yyl23/IRP")
#setwd("/lustre/alice3/scratch/spectre/y/yyl23/IRP_merged")
#setwd("/lustre/alice3/scratch/spectre/y/yyl23/IRP/dataDir/filtered_BAM")
#setwd("/home/yyl23/Desktop/BS7130_Independent_Research_Project/QDNAseq/PTCL/1000bin")

# set bin annotation variable
bins <- getBinAnnotations(binSize=1000, genome="hg38")
#bins <- getBinAnnotations(binSize=500, genome="hg38")
#bins <- getBinAnnotations(binSize=100, genome="hg38")
#bins <- getBinAnnotations(binSize=50, genome="hg38")

# Load sequence data from BAM files ending in .bam from the subdirectory "finalBam" and return with object of class QDNAseqReadCounts
readCounts <- binReadCounts(bins, path="90-160")

# plot raw copy number profile (read counts across the genome) and highlight bins that will be removed with default filtering
plot(readCounts, logTransform=FALSE, ylim=c(-50, 20000))
highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
save(readCounts, file="../qdnaseq/readCounts.RData")

# apply filters and plot median read counts as a function of GC content and mappability
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
isobarPlot(readCountsFiltered)

# estimate the correction for GC content and mappability, and make a plot for the relationship between the observed standard deviation in the data and its read depth
# The theoretical expectation is a linear relationship, which is shown in the plot with a black line. Samples with low-quality DNA will be noisier than expected and appear further above the line than good-quality samples. 
readCountsFiltered <- estimateCorrection(readCountsFiltered)
noisePlot(readCountsFiltered)
save(readCountsFiltered, file="../qdnaseq/readCountsFiltered.RData")

# --> 1st Copy number profile after correcting for GC content and mappability

# apply the correction for GC content and mappability; it returns  a QDNAseq-CopyNumbers object, which we then normalize, smooth outliers, and plot the copy number profile
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers

copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
plot(copyNumbersSmooth)

# export dataset for downstream analysis
saveRDS(copyNumbersSmooth, file="../qdnaseq/PTCL_copyNumbersSmooth.rds")

## Downstream analysis
# Segmentation with the CBS algorithm from DNAcopy , and calling copy number aberrations with CGHcall or cutoffs have been implemented for convenience.
# By default, segmentation uses a log2-transformation, but a sqrt(x + 3/8) can also be used as it stabilizes the variance of a Poisson distribution (Anscombe transform):
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)
saveRDS(copyNumbersSegmented, file="../qdnaseq/PTCL_copyNumbersSegmented.rds")
# Figure 5: Copy number profile after segmenting.


# call aberrations and plot
# use N/2 with availableCores() to keep overhead in-check and return at least one core
# alternatively use # ncores <- max(1L, detectCores(), na.rm = TRUE) as detectCores() may return a missing value/1 and it does not give the number of “allowed” cores
copyNumbersCalled <- callBins(copyNumbersSegmented, ncpus= parallelly::availableCores()/2)
plot(copyNumbersCalled)

# Figure 6: Copy number profile after calling gains and losses.
# save abberations for further manipulation
saveRDS(copyNumbersCalled, file="../qdnaseq/PTCL_copyNumbersCalled.rds")

# callBins() was developed based on a mixture model, and when there are not enough aberrations present in the data, model fitting can fail. if fails, it can also perform simple cutoff-based calling by setting parameter method="cutoff". based on the assumption of uniform cell populations, and in case of cancer samples will most likely need calibration by adjusting parameter cutoffs.

# convert to a cghCall object for other downstream analyses
cgh <- makeCgh(copyNumbersCalled)
save(cgh, file="../qdnaseq/cgh.RData")

save(copyNumbersSegmented, copyNumbersCalled, copyNumbersSmooth, file="../qdnaseq/PTCL_copyNumbersCombined.RData")

# save output to pdf
pdf(file="../qdnaseq/Fig4_copyNumbersSmooth.pdf")
par(cex.main=0.9,cex.lab=0.9, cex.axis=0.8)
plot(copyNumbersSmooth)
dev.off()

pdf(file="../qdnaseq/Fig5_copyNumbersSegmented.pdf")
par(cex.main=0.9,cex.lab=0.9, cex.axis=0.8)
plot(copyNumbersSegmented)
dev.off()

pdf(file="../qdnaseq/Fig6_copyNumbersCalledpdf")
par(cex.main=0.9,cex.lab=0.9, cex.axis=0.8)
plot(copyNumbersCalled)
dev.off()

##################################
setwd("/home/yyl23/Desktop/BS7130_Independent_Research_Project/QDNAseq/PBMC_90-160")

if (file.exists("/home/yyl23/Desktop/BS7130_Independent_Research_Project/QDNAseq/PBMC_90-160/500bin/copyNumbersSegmented.RData"))
{
  load("/home/yyl23/Desktop/BS7130_Independent_Research_Project/QDNAseq/PBMC_90-160/500bin/copyNumbersSegmented.RData")
} else{}

# save to .rds
saveRDS(copyNumbersSegmented, file = "PBMC_90-160_500bin_copyNumbersSegmented.rds")

# Restore the object
readRDS(file = "PTCL_1000bin_copyNumbersSegmented.rds")
