library(chromVAR)
library(GenomicRanges)
library(tidyverse)

# Define sample metadata and select histone mark
metadata <- read.table("~/metadata_CnT.txt", header = TRUE, sep = "\t")
metadata <- subset(metadata, USE !="No" & QC != "No")

metadata <- subset(metadata, histone=="H3K4me3")

histone <- unique(metadata$histone)
condition <- unique(metadata$condition)
replicate <- unique(metadata$replicate)
tissue <- unique(metadata$tissue)
GMO <- unique(metadata$GMO)
mergedPeaks <- unique(metadata$peakFile)
samples <- metadata$ID

#Load blacklist region filtered peaks (these are all peaks --> we will select the peaks for a specifi hPTM)
filtPeaks   <- readRDS("~/filtPeaks.rds")
peakSummary <- filtPeaks$peakSummary

# List peak files to be used
condPeaks <- unique(metadata$ID)

# Merge all called peaks for all samples
mPeak <- GRanges()

for(cp in condPeaks)
{
  peakRes <- subset(peakSummary, ID == cp)
  mPeak   <- GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}

# Reduce genomic ranges (only for genome wide quantification)
masterPeak <- GenomicRanges::reduce(mPeak)

# List samples
samples <- metadata$ID

# Initialize count matrix
countMat <- matrix(NA, length(masterPeak), length(samples))
dim(countMat)

# Overlap with bam files to get count table
for(i in (1:nrow(metadata)))
{
  bam <- metadata[i, "bamFile"]
  fragment_counts <- chromVAR::getCounts(bam, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
  countMat[, i]   <- counts(fragment_counts)[,1]
}

colnames(countMat) <- samples

df <- data.frame(masterPeak)
rownames(countMat) <- paste0(df$seqnames, ":", df$start, "-", df$end)

saveRDS(countMat, "~/H3K4me3_filtPeakCounts.rds")
### these counts can be used for differential analysis or be subsetted for specific regulatory elements (eg. promoters)
