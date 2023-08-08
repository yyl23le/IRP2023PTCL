#!/bin/bash

#***************************************************************************************************************
# Modified by:Friday
# Modified Date:30-05-2023 BST 00:00
# Version: 1

# This analysis script is designed to analyse the fragmentation profiles using ichorCNA.  It uses the analysis-ready bam files, and produce WIG files as an input for ichorCNA. Subsequently, we re-run the analysis on different subsetted BAMs.

# Prior to running the script, please ensure u have hmmcopy_utils and readCounter hmmcopy_utils installed and compiled using C++ and cmake.

##installation of ichorCNA
#install.packages("devtools")
#library(devtools)
#install_github("broadinstitute/ichorCNA")

#Install the HMMcopy suite from https://github.com/shahcompbio/hmmcopy_utils
#cmake .
#make
#***************************************************************************************************************
# Load R
module load R
module load samtools/1.9

# Get the system time
now="$(date +"%c")"

IRPDIR="/lustre/alice3/scratch/spectre/y/yyl23/IRP"
dataDir="/lustre/alice3/scratch/spectre/y/yyl23/IRP/dataDir"

echo ----
echo "Job started at $now"
echo "current wd = $IRPDIR"
echo ----

while read prefix; do
# Create ouput directories
# awk parses the filename - it uses the underscore character to separate fields. A few shells can also do this, like bash. Field $1 is the one you want, so awk prints the the $1 field. $[field number] is the awk syntax for each field in a file - the default field separator is tab/space. The -F '_' tells awk to use underscore.
tmp=$(echo "${prefix}" | awk -F '_' '{print $1}' )
mkdir -p $IRPDIR/ichorCNA/PTCL_${tmp}

#Generate Read Count File; use --list .bam to list all chrom name
echo "generating read count file for ${prefix}..."
/lustre/alice3/scratch/spectre/y/yyl23/IRP/hmmcopy_utils-master/bin/readCounter --window 1000000 --quality 20 -c "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" $IRPDIR/finalBam/${prefix}_dD-RG-BQSR_FINAL.bam > $IRPDIR/ichorCNA/PTCL_${tmp}.wig

echo "running ichorCNA on ${prefix}..."
Rscript /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/scripts/runIchorCNA.R --id PTCL_${tmp} \
  --WIG $IRPDIR/ichorCNA/PTCL_${tmp}.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/gc_hg38_1000kb.wig \
  --mapWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/map_hg38_1000kb.wig \
  --centromere /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --normalPanel /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $IRPDIR/ichorCNA/PTCL_${tmp}/

done <$dataDir/sample_list.txt

##50-100
while read prefix; do
# Create ouput directories
tmp=$(echo "${prefix}" | awk -F '_' '{print $1}' )
mkdir -p $IRPDIR/ichorCNA/50-100/PTCL_${tmp}

# Generate index files
echo "indexing ${prefix}_is-50-100..."
samtools index $dataDir/filtered_BAM/50-100/${prefix}_is-50-100.bam

#Generate Read Count File; use --list .bam to list all chrom name
echo "generating read count file of ${prefix}..."
/lustre/alice3/scratch/spectre/y/yyl23/IRP/hmmcopy_utils-master/bin/readCounter --window 1000000 --quality 20 -c "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" $dataDir/filtered_BAM/50-100/${prefix}_is-50-100.bam > $IRPDIR/ichorCNA/50-100/PTCL_${tmp}_50-100.wig;

echo "running ichorCNA on ${prefix}..."
Rscript /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/scripts/runIchorCNA.R --id PTCL_50-100_${tmp} \
  --WIG $IRPDIR/ichorCNA/50-100/PTCL_${tmp}_50-100.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/gc_hg38_1000kb.wig \
  --mapWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/map_hg38_1000kb.wig \
  --centromere /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --normalPanel /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $IRPDIR/ichorCNA/50-100/PTCL_${tmp}/
done <$dataDir/sample_list.txt

##100-150
while read prefix; do
# Create ouput directories
tmp=$(echo "${prefix}" | awk -F '_' '{print $1}' )
mkdir -p $IRPDIR/ichorCNA/100-150/PTCL_${tmp}

# Generate index files
echo "indexing ${prefix}_is-100-150..."
samtools index $dataDir/filtered_BAM/100-150/${prefix}_is-100-150.bam

#Generate Read Count File; use --list .bam to list all chrom name
echo "generating read count file of ${prefix}..."
/lustre/alice3/scratch/spectre/y/yyl23/IRP/hmmcopy_utils-master/bin/readCounter --window 1000000 --quality 20 -c "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" $dataDir/filtered_BAM/100-150/${prefix}_is-100-150.bam > $IRPDIR/ichorCNA/100-150/PTCL_${tmp}_100-150.wig;

echo "running ichorCNA on ${prefix}..."
Rscript /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/scripts/runIchorCNA.R --id PTCL_100-150_${tmp} \
  --WIG $IRPDIR/ichorCNA/100-150/PTCL_${tmp}_100-150.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/gc_hg38_1000kb.wig \
  --mapWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/map_hg38_1000kb.wig \
  --centromere /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --normalPanel /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $IRPDIR/ichorCNA/100-150/PTCL_${tmp}/
done <$dataDir/sample_list.txt


##50-150
while read prefix; do
# Create ouput directories
tmp=$(echo "${prefix}" | awk -F '_' '{print $1}' )
mkdir -p $IRPDIR/ichorCNA/50-150/PTCL_${tmp}

# Generate index files
echo "indexing ${prefix}_is-50-150..."
samtools index $dataDir/filtered_BAM/50-150/${prefix}_is-50-150.bam

#Generate Read Count File; use --list .bam to list all chrom name
echo "generating read count file of ${prefix}..."
/lustre/alice3/scratch/spectre/y/yyl23/IRP/hmmcopy_utils-master/bin/readCounter --window 1000000 --quality 20 -c "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" $dataDir/filtered_BAM/50-150/${prefix}_is-50-150.bam > $IRPDIR/ichorCNA/50-150/PTCL_${tmp}_50-150.wig;

echo "running ichorCNA on ${prefix}..."
Rscript /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/scripts/runIchorCNA.R --id PTCL_50-150_${tmp} \
  --WIG $IRPDIR/ichorCNA/50-150/PTCL_${tmp}_50-150.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/gc_hg38_1000kb.wig \
  --mapWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/map_hg38_1000kb.wig \
  --centromere /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --normalPanel /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $IRPDIR/ichorCNA/50-150/PTCL_${tmp}/
done <$dataDir/sample_list.txt


##90-160
while read prefix; do
# Create ouput directories
tmp=$(echo "${prefix}" | awk -F '_' '{print $1}' )
mkdir -p $IRPDIR/ichorCNA/90-160/PTCL_${tmp}

# Generate index files
echo "indexing ${prefix}_is-90-160..."
samtools index $dataDir/filtered_BAM/90-160/${prefix}_is-90-160.bam

#Generate Read Count File; use --list .bam to list all chrom name
echo "generating read count file of ${prefix}..."
/lustre/alice3/scratch/spectre/y/yyl23/IRP/hmmcopy_utils-master/bin/readCounter --window 1000000 --quality 20 -c "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" $dataDir/filtered_BAM/90-160/${prefix}_is-90-160.bam > $IRPDIR/ichorCNA/90-160/PTCL_${tmp}_90-160.wig;

echo "running ichorCNA on ${prefix}..."
Rscript /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/scripts/runIchorCNA.R --id PTCL_90-160_${tmp} \
  --WIG $IRPDIR/ichorCNA/90-160/PTCL_${tmp}_90-160.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/gc_hg38_1000kb.wig \
  --mapWig /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/map_hg38_1000kb.wig \
  --centromere /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --normalPanel /home/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2/ichorCNA/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $IRPDIR/ichorCNA/90-160/PTCL_${tmp}/
done <$dataDir/sample_list.txt
