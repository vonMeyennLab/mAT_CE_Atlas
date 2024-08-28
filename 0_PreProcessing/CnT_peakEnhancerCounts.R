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
library(tidyr)

##generate enhancer bed files per GMO
#We ran ChromHMM to identify the enhancers

###Ucp1ERCre
# enhancer states are 7 and 8
# Load bed files from chroHMM output and select correct states to use as enhancers
CE <- read.table("~/CE_segments.bed")
CERT <- read.table("~/CERT_segments.bed")
CETN <- read.table("~/CETN_segments.bed")

enhancers_CE <- subset(CE, CE$V4 == c("E7", "E8"))
enhancers_CERT <- subset(CERT, CERT$V4 == c("E7", "E8"))
enhancers_CETN <- subset(CETN, CETN$V4 == c("E7", "E8"))

Ucp1ERCre_enhancers <- rbind(enhancers_CE,enhancers_CERT, enhancers_CETN)
colnames(Ucp1ERCre_enhancers) <- c('seqnames', 'start', 'end', 'name')
### save as one bed file 
write.table(Ucp1ERCre_enhancers, "~/Enhancers_Ucp1ERCre.bed", sep = "\t")

### AdipoCre 
#enhancer state 9
CE_DT <- read.table("~/CE_DT_9_segments.bed")
RT <- read.table("~/RT_9_segments.bed")
TN <- read.table("~/TN_9_segments.bed")

enhancers_CE_DT <- subset(CE_DT, CE_DT$V4 == c("E9"))
enhancers_RT <- subset(RT, RT$V4 == c("E9"))
enhancers_TN <- subset(TN, TN$V4 == c("E9"))
AdipoCre_enhancers <- rbind(enhancers_CE_DT,enhancers_RT, enhancers_TN)
colnames(AdipoCre_enhancers) <- c('seqnames', 'start', 'end', 'name')
write.table(AdipoCre_enhancers, "~/Enhancers_AdipoCre.bed", sep = "\t")


## Load enhancers for Ucp1ERCre

enhancers <- read.table("~/Enhancers_Ucp1ERCre.bed", header = TRUE, sep = "\t") 

# Make GRanges from df
EN_GR <- makeGRangesFromDataFrame(enhancers,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="seqnames",
                                  start.field="start",
                                  end.field="end")

EN_GR <- renameSeqlevels(EN_GR, mapSeqlevels(seqlevels(EN_GR), "UCSC"))


## Load CnT sample data
metadata <- read.table("~/metadata_CnT.txt", header = TRUE, sep = "\t")

# Select only info for the specific histone 
metadata <- subset(metadata, USE !="No"  & QC != "No")

metadata <- subset(metadata, histone=="H3K4me1")

# Load raw CnT peak counts
rawCnTCounts <- readRDS("~/H3K4me1_filtPeakCounts.rds")

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


# Find peak regions overlapping with enhancers
table(!is.na(findOverlaps(CnTpeakRawCounts, EN_GR, select="arbitrary")))

EN_overlapsGR <- GenomicRanges::findOverlaps(CnTpeakRawCounts, EN_GR)
EN_ol <- data.frame(CnTpeakRawCounts[unique(EN_overlapsGR@from)])  
EN_ol$loci <- paste0(EN_ol$seqnames, ":", EN_ol$start, "-", EN_ol$end)

EN_ol_GR <- GRanges(seqnames = EN_ol$seqnames, IRanges(start = EN_ol$start, end = EN_ol$end), strand = "*", keep.extra.columns=TRUE)

# Annotate CnT peaks
CnTpeakAnno <- annotatePeak(EN_ol_GR, tssRegion=c(-2000, 2000), 
                            TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

# Convert the peak annotations into data frame
CnTpeakAnno_df <- data.frame(CnTpeakAnno)

# Select peaks outside promoter region
CnTpeakEn <- subset(CnTpeakAnno_df, abs(distanceToTSS) > 2000)
CnTpeakEn <- CnTpeakEn[-grep("predicted gene|microRNA|RIKEN|Riken|pseudogene", CnTpeakEn$GENENAME),]

CnTpeakEn$loci <- paste0(CnTpeakEn$seqnames, ":", CnTpeakEn$start, "-", CnTpeakEn$end)

EN_ol_peakCounts <- subset(EN_ol, loci %in% CnTpeakEn$loci)

CnTpeakEn <- EN_ol_peakCounts[,c("loci", CnTsamples)]
rownames(CnTpeakEn) <- NULL
CnTpeakEn <- column_to_rownames(CnTpeakEn, "loci")

saveRDS(CnTpeakEn, "~/H3K4me1_filtPeakUcp1ERCreEnCounts.rds")


###
## Load enhancers for AdipoCre

enhancers <- read.table("~/Enhancers_AdipoCre.bed", header = TRUE, sep = "\t") 

# Make GRanges from df
EN_GR <- makeGRangesFromDataFrame(enhancers,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="seqnames",
                                  start.field="start",
                                  end.field="end")

EN_GR <- renameSeqlevels(EN_GR, mapSeqlevels(seqlevels(EN_GR), "UCSC"))


## Load CnT sample data
metadata <- read.table("~/metadata_CnT.txt", header = TRUE, sep = "\t")

# Select only info for the specific histone 
metadata <- subset(metadata, USE !="No" &  QC != "No")

metadata <- subset(metadata, histone=="H3K27ac")

# Load raw CnT peak counts
rawCnTCounts <- readRDS("~/H3K27ac_filtPeakCounts.rds")

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


# Find peak regions overlapping with enhancers
table(!is.na(findOverlaps(CnTpeakRawCounts, EN_GR, select="arbitrary")))

EN_overlapsGR <- GenomicRanges::findOverlaps(CnTpeakRawCounts, EN_GR)
EN_ol <- data.frame(CnTpeakRawCounts[unique(EN_overlapsGR@from)])  
EN_ol$loci <- paste0(EN_ol$seqnames, ":", EN_ol$start, "-", EN_ol$end)

EN_ol_GR <- GRanges(seqnames = EN_ol$seqnames, IRanges(start = EN_ol$start, end = EN_ol$end), strand = "*", keep.extra.columns=TRUE)

# Annotate CnT peaks
CnTpeakAnno <- annotatePeak(EN_ol_GR, tssRegion=c(-2000, 2000), 
                            TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

# Convert the peak annotations into data frame
CnTpeakAnno_df <- data.frame(CnTpeakAnno)

# Select peaks outside promoter region
CnTpeakEn <- subset(CnTpeakAnno_df, abs(distanceToTSS) > 2000)
CnTpeakEn <- CnTpeakEn[-grep("predicted gene|microRNA|RIKEN|Riken|pseudogene", CnTpeakEn$GENENAME),]

CnTpeakEn$loci <- paste0(CnTpeakEn$seqnames, ":", CnTpeakEn$start, "-", CnTpeakEn$end)

EN_ol_peakCounts <- subset(EN_ol, loci %in% CnTpeakEn$loci)

CnTpeakEn <- EN_ol_peakCounts[,c("loci", CnTsamples)]
rownames(CnTpeakEn) <- NULL
CnTpeakEn <- column_to_rownames(CnTpeakEn, "loci")

saveRDS(CnTpeakEn, "~/H3K27ac_filtPeakAdipoCreEnCounts.rds")
