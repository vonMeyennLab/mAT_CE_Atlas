library(tidyverse)
library(biomaRt)
library(dplyr)
library(edgeR)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(tibble)
library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

## Diff analysis 

  ### CountQC for H3K4me1 quantified in enhancer peaks

## Load CnT sample data
metadata <- read.table("~metadata_CnT.txt", header = TRUE, sep = "\t")

# Select only info for the specific histone and conditions to be analysed
metadata <- subset(metadata, USE !="No" & QC != "No")

metadata <- subset(metadata, histone=="H3K4me1" &tissue == "ingAT" & tissue = "BAT")
# Load filtered peak counts (check for which enhancers)
rawCounts <- readRDS("~/public/LauraHinte/BRITE/Enhancers/Counts/H3K4me1_filtPeakUcp1ERCreEnCounts.rds")

rawCounts <- rawCounts[,metadata$ID]
rawCounts <- rawCounts[rowSums(rawCounts) > 0, ]

# Extract sample info
samples <- colnames(rawCounts)
histone <- as.factor(metadata$histone)
condition <- as.factor(metadata$condition)
replicate <- as.factor(metadata$replicate)
tissue <- as.factor(metadata$tissue)
GMO <- as.factor(metadata$GMO)

# Add sample information
group <- data.frame(histone, condition, tissue, GMO, row.names = samples)

##### For different tissues of origin we use FDR and for from same tissue we use adj p-val

# Create DGEList object
y <- edgeR::DGEList(counts=rawCounts, group=group$condition)

# Remove genes/peaks with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts <- edgeR::cpm(y, log=TRUE)

  
  #### Sample clustering


# Run PCA
pcDat  <- prcomp(t(normCounts), scale. = TRUE)

# Tissue and GMO
p1<-autoplot(pcDat, data=group,
             colour="condition", shape = "GMO", size=4, frame = F, label = F) + theme_classic()

# GMO and condition

p2<-autoplot(pcDat, data=group,
             colour="condition", shape="tissue", size=4, label = F) + theme_classic()

p1 
p2


#### Clustering of specific comparison

## Load sample data
metadata_sub <- subset(metadata, condition %in% c("CE", "CE"))
metadata_sub <- subset(metadata_sub, tissue %in% c("ingAT", "BAT"))

# Load raw peak counts
rawCounts_sub <- rawCounts[,metadata_sub$ID]
samples_sub <- colnames(rawCounts_sub)

# Add sample information
group_sub <- data.frame(tissue=as.factor(metadata_sub$tissue), row.names = samples_sub)

# Create DGEList object
y <- DGEList(counts=rawCounts_sub, group=group_sub$tissue)

# Remove genes with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts_sub <- edgeR::cpm(y, log=TRUE)


##### PCA

# Run PCA
pcDat  <- prcomp(t(normCounts_sub), scale. = TRUE)

# PC1 vs PC2
autoplot(pcDat, data=group_sub,
         colour="tissue", size=4) + theme_classic()

  
  #### Volcano plot


# Construct design matrix
design <- model.matrix(~ 0 + group_sub$tissue)
colnames(design) <- levels(group_sub$tissue)

# Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)                  

# Estimate quasi-likelihood (QL) dispersions
fitGlm <- glmQLFit(y, design, robust = TRUE)

# Define contrast of interest
contr <- makeContrasts(ingAT - epiAT, levels=design)

# Apply QL F-test for differential expression
qlf <- glmQLFTest(fitGlm, contrast=contr)

# Extract DE results
res <- topTags(qlf, n = Inf, sort.by = "PValue")$table                         

# Annotate genes based on log2FC / adjusted p-value thresholds
res$regulation <- ""
res[which(res$logFC >= 1 & res$FDR <= 0.05),"regulation"] <- "up"
res[which(res$logFC <= -1 & res$FDR <= 0.05),"regulation"] <- "down"
res[which(abs(res$logFC) < 1 | res$FDR > 0.05),"regulation"] <- "non-significant"

table(res$regulation)


# Select significant DEGs based on FDR corrected p values
sig_peaks <- subset(res, FDR < 0.05 & abs(logFC) > 1)

# Add annotation
sig_peaks <- rownames_to_column(sig_peaks, "loci")

sig_peaks <- separate(sig_peaks, "loci", c("chr", "pos"), sep=":")
sig_peaks <- separate(sig_peaks, "pos", c("start", "end"), sep="-")
sig_peaks$start <- as.integer(sig_peaks$start)
sig_peaks$end <- as.integer(sig_peaks$end)

# Save significant diff peaks data.frame as csv to check peaks in Seqmonk
write.csv(sig_peaks, "~/H3K4me1_sigDP_CE_ingAT_vs_CE_BAT_Ucp1ERCre.csv")
### rember that for GSEA the distance has to be in there

# Update for peak fold enrichmnet
sig_peaks$chr <- gsub("chr", "", sig_peaks$chr)

# Split significant diff peaks into up and down
upPeaks <- subset(sig_peaks, regulation == "up")
downPeaks <- subset(sig_peaks, regulation == "down")

histDiffPeaks <- list("up" = upPeaks, 
                      "down" = downPeaks)


```
#### Peak annotation
<br/>
  Significant diff peaks: log2FC > 1 ; FDR < 0.05

```{r echo=FALSE, message=FALSE, warning=FALSE}

histDiffPeaks <- lapply(histDiffPeaks, function(x) {makeGRangesFromDataFrame(x,
                                                                             keep.extra.columns=TRUE,
                                                                             ignore.strand=FALSE,
                                                                             seqinfo=NULL,
                                                                             seqnames.field=c("chromosome", "chrom",
                                                                                              "chr", "chromosome_name"),
                                                                             start.field="start",
                                                                             end.field=c("end", "stop"),
                                                                             strand.field="strand",
                                                                             starts.in.df.are.0based=FALSE)})

histDiffPeaks <- lapply(histDiffPeaks, function(x) {renameSeqlevels(x, mapSeqlevels(seqlevels(x), "UCSC"))})

# Annotate the peaks and save (for GSEA for example)
histDiffPeaks_anno <- lapply(histDiffPeaks, function(x) {annotatePeak(x, tssRegion=c(-2000, 2000), 
                                                                      TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")})

write.csv(data.frame(histDiffPeaks_anno$up), "~/H3K4me1_CE_ingAT_vs_CE_BAT_Ucp1ERCre_upPeaks.csv")
write.csv(data.frame(histDiffPeaks_anno$down), "~/H3K4me1_CE_ingAT_vs_CE_BAT_Ucp1ERCre_downPeaks.csv")

