
library(ChIPseeker)
library(ChIPpeakAnno)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

## Load sample data
## metadata contains information on conditions, sequencing QC, sample ID, hPTM, GMO, biological replicate etc.
metadata <- read.table("~/metadata_CnT.txt", header = TRUE, sep = "\t")
metadata <- subset(metadata, USE !="No" & QC != "No") ## subset what you want to annotate

# Load blacklist regions (with BAC or not BAC) The bed file that includes the ~300 kb region (BAC) spanning Ucp1 used to create the Ucp1ErCre mice was used for H3K9me3
blReg <- read.table("~/mm10_blacklist_nochr.bed", sep = "\t", header = FALSE)
blRegGR <- GRanges(blReg$V1, IRanges(start = blReg$V2, end = blReg$V3), strand = "*")

# List unique factors
histone <- unique(metadata$histone)
condition <- unique(metadata$condition)
replicate <- unique(metadata$replicate)
tissue <- unique(metadata$tissue)
GMO <- unique(metadata$GMO)

# Initialize data frames
peakN       <- c()    # peak number
peakSummary <- c()    # peak summary

# Calculate peak number and peak width for each sample
for(i in (1:nrow(metadata)))
{
  sampleInfo <- strsplit(metadata[i,"ID"], "_")[[1]]
  peakInfo   <- read.table(metadata[i, "peakFile"], header = FALSE, fill = TRUE)
  peakInfo$ID <- paste0(peakInfo$V1, ":", peakInfo$V2, "-", peakInfo$V3)
  
  # Calculate peak width
  peakInfo$ID           <- metadata[i,"ID"]
  peakInfo$condition    <- sampleInfo[1]
  peakInfo$tissue          <- sampleInfo[2]
  peakInfo$GMO      <- sampleInfo[3]
  peakInfo$histone    <- sampleInfo[4]
  peakInfo$replicate    <- sampleInfo[5]
  peakInfo$width        <- abs(peakInfo$V3 - peakInfo$V2) + 1
  
  # Calculate peak % overlapping with blacklist region
  histGR <- GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
  peakInfo$numBLoverlap <- GenomicRanges::countOverlaps(histGR, blRegGR)
  
  countBLoverlap <- length(which(peakInfo$numBLoverlap>0))
  
  # Summarize peak statistics
  peakSummary <- rbind(peakSummary, peakInfo)
  peakN       <- data.frame(sampleID = metadata[i,"ID"],
                            condition = sampleInfo[1], tissue = sampleInfo[2], GMO = sampleInfo[3], histone = sampleInfo[4], replicate = sampleInfo[5], 
                            peakN = nrow(peakInfo) - countBLoverlap, countBLoverlap = countBLoverlap) %>% rbind(peakN, .)
}

# Remove the peaks overlapping the blacklist region
peakSummary <- subset(peakSummary, numBLoverlap == 0)

saveRDS(list(peakSummary = peakSummary, peakN = peakN), "CnT_wBAC_filtPeaks.rds")

################################################################################################################################

## Peak annotation

# Initialize empty list to store the peak annotations
# needs to be done per mark to finish
# to get annotated peaks per mark
### or can be done for all marks but it can be too many plos to look at


sampleL <- list()
metadata2 <- subset(metadata, histone == "H3K4me3")
samples <- metadata2$ID


# Calculate overlapping peaks for each mark
for(sample in samples)
{

  hist <- subset(peakSummary, ID == sample)

  histGR <- GRanges(hist$V1, IRanges(start = hist$V2, end = hist$V3), strand = "*")

  histGR <- renameSeqlevels(histGR, mapSeqlevels(seqlevels(histGR), "UCSC"))
  sampleL[[sample]] <- histGR
}

# Annotate the peaks
peakAnno <- lapply(sampleL, function(x){annotatePeak(x, tssRegion=c(-2000, 2000),
                                                      TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")})

saveRDS(peakAnno, "CnT_filtPeakAnno_H3K4me3.rds")

#promoter <- getPromoters(TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, upstream=3000, downstream=3000)
#tagMatrix <- lapply(sampleL, function(x){getTagMatrix(x, windows=promoter)})
 
#saveRDS(tagMatrix, "CnT_filtPeakTagMatrix.rds")

################################################################################################################################

# Load mm10 genomic annotation
genome_annotation <- readRDS("~/genomeAnnotation_mm10.rds")

# Adding CGI and Non-CGI promoters to the annotation list
genome_annotation$`CGI Promoters`     <- genome_annotation$Promoters %>% .[.$CGI == TRUE]
genome_annotation$`Non-CGI Promoters` <- genome_annotation$Promoters %>% .[.$CGI == FALSE]

genome_annotation_subset <-
  genome_annotation %>% .[c("CGI Promoters", "Non-CGI Promoters", "3' UTRs", "5' UTRs",
                            "Exons", "Introns", "Intergenic regions", 
                            "PLS", "pELS", "dELS", "CTCF", "DNase_H3K4me3" )] %>% GRangesList()
                            
# Function calculating genome size based on genic and intergenic regions
genome_size <- 
  sum(
    GenomicRanges::reduce(genome_annotation$`Gene bodies`, ignore.strand = TRUE) %>% 
      GenomicRanges::width() %>% sum,
    GenomicRanges::reduce(genome_annotation$`Intergenic regions`, ignore.strand = TRUE) %>%
      GenomicRanges::width() %>% sum
  )

calc_fold_enrichment <- 
  function(gr1, gr2, genome_size){
    
    A <- 
      GenomicRanges::intersect(
        GenomicRanges::reduce(gr2, ignore.strand = TRUE),
        GenomicRanges::reduce(gr1, ignore.strand = TRUE)
      ) %>% GenomicRanges::width() %>% sum
    
    B <- 
      GenomicRanges::reduce(gr1, ignore.strand = TRUE) %>%
      GenomicRanges::width() %>% sum
    
    C <- 
      GenomicRanges::reduce(gr2, ignore.strand = TRUE) %>% 
      GenomicRanges::width() %>% sum
    
    log2((A/B)/(C/genome_size))
  }

## Calculate peak fold enrichment
sampleFE <- list()

for (sample in samples)
{
  hist <- subset(peakSummary, ID == sample)
  hist$chr <- gsub("chr", "", hist$V1)
  histGR <- GRanges(hist$V1, IRanges(start = hist$V2, end = hist$V3), strand = "*")

  histFE <- lapply(genome_annotation_subset, function(x){calc_fold_enrichment(histGR, x, genome_size = genome_size)})
  sampleFE[[sample]] <- histFE
}

sampleFE <- lapply(sampleFE, function(x){x <- do.call("rbind", x)})

# Peak fold enrichment values of all samples in table format
comb_FE <- do.call("cbind", sampleFE)
colnames(comb_FE) <- names(sampleFE)

saveRDS(comb_FE, "CnT_filtPeakFE.rds")
### can be used to plot the peak fold enrichment in a heatmap