library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(dplyr)
library(edgeR)
library(ggfortify)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(tibble)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

## Load CnT sample metadata
metadata <- read.table("~/metadata_CnT.txt", header = TRUE, sep = "\t")

# Select only info for the specific histone 
metadata <- subset(metadata, histone=="H3K4me3")

# Load raw CnT peak counts of filtered peaks
rawCnTCounts <- readRDS("~/H3K4me3_filtPeakCounts.rds")

rawCnTCounts <- rawCnTCounts[,metadata$ID]

# List unique factors
CnTsamples <- colnames(rawCnTCounts)

# Annotate peaks
rawCnTCounts <- data.frame(rawCnTCounts)
rawCnTCounts$loci <- rownames(rawCnTCounts)

CnTpeakRawCounts <- separate(rawCnTCounts, "loci", c("chr", "pos"), sep=":")
CnTpeakRawCounts <- separate(CnTpeakRawCounts, "pos", c("start", "end"), sep="-")


CnTpeakRawCounts <- makeGRangesFromDataFrame(CnTpeakRawCounts,
                                             keep.extra.columns=TRUE,
                                             ignore.strand=FALSE,
                                             seqinfo=NULL,
                                             seqnames.field=c("chromosome", "chrom",
                                                              "chr", "chromosome_name"),
                                             start.field="start",
                                             end.field=c("end", "stop"),
                                             strand.field="strand",
                                             starts.in.df.are.0based=FALSE)

CnTpeakRawCounts <- renameSeqlevels(CnTpeakRawCounts, mapSeqlevels(seqlevels(CnTpeakRawCounts), "UCSC"))

CnTpeakAnno <- annotatePeak(CnTpeakRawCounts, tssRegion=c(-2000, 2000), 
                            TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")



# Convert the peak annotations into data frame
CnTpeakAnno_df <- data.frame(CnTpeakAnno)

# Summarize peak read counts in promoter
CnTpeakPr <- subset(CnTpeakAnno_df, abs(distanceToTSS) <= 2000)
CnTpeakPr <- CnTpeakPr[-grep("predicted gene|microRNA|RIKEN|Riken|pseudogene", CnTpeakPr$GENENAME),]

CnTpeakPr <- CnTpeakPr[,c("geneId", CnTsamples)]

CnTpeakPr_rawCounts <- CnTpeakPr %>% 
  group_by(geneId) %>% 
  summarise(across(everything(), sum))

CnTpeakPr_rawCounts <- column_to_rownames(CnTpeakPr_rawCounts, "geneId")
saveRDS(CnTpeakPr_rawCounts, "~/H3K4me3_peakPr_rawCounts.rds")




