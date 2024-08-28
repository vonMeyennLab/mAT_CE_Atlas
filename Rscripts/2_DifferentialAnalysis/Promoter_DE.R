

library(dbplyr)
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


## Diff analysis 

 ### CountQC

## Load CnT sample data
metadata <- read.table("~/metadata_CnT.txt", header = TRUE, sep = "\t")

# Select only info for the specific histone and maybe GMOs 
metadata <- subset(metadata, USE !="No" & QC != "No")

metadata <- subset(metadata, histone=="H3K27me3"& GMO == "Ucp1ERCre")
# Load promoter counts
rawCounts <- readRDS("~/H3K27me3_peakPr_rawCounts.rds")

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

### for gene annotation later

ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")
### initial inspection
# Create DGEList object
y <- edgeR::DGEList(counts=rawCounts, group=group$condition)

# Remove genes with low expression in > 50% samples
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
             colour="condition", shape = "tissue", size=4, frame = F, label = F) + theme_classic()

# GMO and condition

p2<-autoplot(pcDat, data=group,
             colour="condition", shape="GMO", size=4, label = F) + theme_classic()

p1 
p2

### Differential analysis between two conditions/sample groups
  ### Example for comparing brown and beige cold exposed adipocytes
  ### Ucp1ERCre brown CE vs Ucp1ERCre beige CE
  
  #### Clustering

## Load sample data
metadata_sub <- subset(metadata, condition %in% "CE")
metadata_sub <- subset(metadata_sub, tissue %in% c("BAT", "ingAT"))

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

# Construct design matrix
design <- model.matrix(~ 0 + group_sub$tissue)
colnames(design) <- levels(group_sub$tissue)

# Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)                  

# Estimate quasi-likelihood (QL) dispersions
fitGlm <- glmQLFit(y, design, robust = TRUE)

# Define contrast of interest
contr <- makeContrasts(BAT - ingAT, levels=design)

# Apply QL F-test for differential expression
qlf <- glmQLFTest(fitGlm, contrast=contr)

# Extract DE results
res <- topTags(qlf, n = Inf, sort.by = "PValue")$table                         

# Annotate genes based on log2FC / adjusted p-value thresholds
res$regulation <- ""
res[which(res$logFC >= 1 & res$FDR <= 0.05),"regulation"] <- "up"
res[which(res$logFC <= -1 & res$FDR <= 0.05),"regulation"] <- "down"
res[which(abs(res$logFC) < 1 | res$FDR > 0.05),"regulation"] <- "non-significant"

table(res$regulation) # check numbers and decide if you need to filter more or less (expecially logFC could be altered)

# Volcano plot
ggplot(res, aes(x = logFC, y = -log10(FDR), color = regulation)) +
  geom_point() + theme_minimal() +
  scale_color_manual(values = c("#00AFBB","#DCDCDC","#FC4E07")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        legend.title = element_text(size=10, face="italic")) +
  labs(x="log fold change", y="-log(FDR)") +
  geom_vline(xintercept=c(-1, 1), lty=2, colour = "#323232") +
  geom_hline(yintercept=-log10(0.05), lty=2, colour = "#323232" ) 

# Select significant DEGs based on FDR corrected p values; cut offs for FDR or p.adjust and logFC can be changed if needed (e.g. too many hits)
sigGenes <- subset(res, FDR < 0.05 & abs(logFC) > 1)

# Add gene names and save DF

geneAnnotation <- getBM(attributes=c("entrezgene_id", "external_gene_name"), 
                        filters="entrezgene_id", values=rownames(sigGenes), mart=ensembl)

sigGenesPr <- merge(geneAnnotation, sigGenes, by.x="entrezgene_id", by.y="row.names")

write.csv(sigGenesPr, "~/H3K27me3_sigDEGPr_CE_BAT_vs_CE_ingAT.csv")

# Combine with normalized counts 
sigGenesPr <- merge(sigGenesPr, normCounts, by.x="entrezgene_id", by.y="row.names")

# Split significant DEGs into up and down
upGenes <- subset(sigGenesPr, regulation == "up")
downGenes <- subset(sigGenesPr, regulation == "down")
### can be used to generate gene/promoter specific heatmaps (though not in this MS)
```
