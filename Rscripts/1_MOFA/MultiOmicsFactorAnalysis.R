
# This R script has been used to generate the following figures:
### Figure 1-5 and associated Extended Data

library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(MOFAcellulaR)

# Load sample metadata for H3K4me3, H3K27me3, H3K27ac, H3K4me1 and H3K9me3
# Metadata table has the following columns:
# ID | omics | replicate | tissue | GMO | condition | bamFilePath | peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Define samples and conditions
condition <- metadata$condition
samples <- metadata$ID

################################################################################

# Following code chunk has been used for each histone mark
# to identify the most variable 3,000 peak regions

# Load union peak based raw counts (for each histone mark)
# Peak regions are represented in (chr):(start)-(end)_(annotated_closest_gene) format
peakCounts <- readRDS("../filtPeakCounts.rds")

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(peakCounts), 
                              colData = metadata,
                              design = ~condition)

# Remove peak regions with < 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Apply variance stabilizing transformation
vsd <- vst(dds)

# Extract normalized peak counts
normCounts <- assay(vsd)

# Identify most variable 3,000 peak regions
varPeaks <- rownames(normCounts)[head(order(rowSds(as.matrix(normCounts)),
                                            decreasing=TRUE),3000)]

# Select normalized counts of most variable peak regions
varPeaks_normCounts <- normCounts[varPeaks,]

# Transform normalized count matrix into data frame
h3k27ac <- data.frame(varGenes_normCounts)

# Add NA to missing replicate(s)
setdiff(sampleL, colnames(h3k27ac))
h3k27ac$CETN_ingAT_Ucp1ERCre_4 <- NA

# Sort data frame in sample order
h3k27ac <- h3k27ac[,sampleL]

################################################################################

# Following code chunk has been used to prepare the input
# and run Multi-Omics Factor Analysis (MOFA)

# Combine most variable features from different omics assays into a list
# Ensure for each assay, samples are in the same order
matList <- list("H3K9me3" = h3k9me3,
                "H3K4me3"=h3k4me3,
                "H3K27me3"=h3k27me3,
                "H3K27ac"=h3k27ac,
                "H3K4me1"=h3k4me1)

# Convert into list of matrices
matList <- lapply(matList, function(x) { x <- as.matrix(x)})

# Create MOFA object
MOFAobject <- create_mofa(matList)

# Check summary of different omics assays
plot_data_overview(MOFAobject)

# Define MOFA options
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

# Build model with 5 latent factors
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5

# Train model
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$weight_views <- TRUE

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts)

# Run MOFA
MOFAobject <- run_mofa(MOFAobject)

# Add sample metadata to the model
metadata <- rownames_to_column(metadata, "sample")
samples_metadata(MOFAobject) <- metadata

# Define group info
samples <- colnames(h3k4me1)
group <- metadata[,c("cond", "tissue", "condition", "GMO")]
group <- unique(group)

rownames(group) <- NULL
group <- column_to_rownames(group, "cond")
group <- group[sampleL,]

################################################################################

# Following code chunk has been used to explore and visualize MOFA results

# Check if the factors are largely not correlated
# Indicator of good model fit
plot_factor_cor(MOFAobject)

# Visualize % of variance explained by each factor across each data modality 
p <- plot_variance_explained(MOFAobject, plot_total = FALSE) 
p + scale_color_brewer(palette="Reds")

# Total variance explained per assay
plot_variance_explained(MOFAobject, plot_total = T)[[2]]

# Cumulative variance explained by factors
r2 <- model@cache$variance_explained$r2_per_factor[[1]]

r2.dt <- r2 %>%
  as.data.table %>% .[,factor:=as.factor(1:model@dimensions$K)] %>%
  melt(id.vars=c("factor"), variable.name="view", value.name = "r2") %>%
  .[,cum_r2:=cumsum(r2), by="view"]

# Visualize cumulative variance as explained by each hPTM
ggplot(r2.dt, aes(x=view, y=cum_r2, fill=view)) + 
  geom_bar(stat="identity",  width=0.9) +
  xlab("") + ylab("Variance explained (%)") +
  scale_fill_manual(values = colPalette) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color="black"),
    axis.text.y = element_text(color="black"),
    axis.title.y = element_text(color="black"),
    axis.line = element_line(color="black"),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text()
  )

# Statistical test to determine the association between categorical variable and latent factors
# Tissue: F1
get_associations(model = MOFAobject, metadata = metadata,
                 sample_id_column = "sample", test_variable = "tissue", 
                 test_type = "categorical", group = FALSE)

# Condition: F2
get_associations(model = MOFAobject, metadata = metadata,
                 sample_id_column = "sample", test_variable = "condition", 
                 test_type = "categorical", group = FALSE)
            
assoc_list <- list(Tissue = tissue_assoc,
                   Condition = cond_assoc)

# Visualize model summary as a heatmap
plot_MOFA_hmap(model = MOFAobject,
               group = FALSE,
               metadata = meta,
               sample_id_column = "sample",
               sample_anns = c("tissue", "condition"),
               assoc_list = assoc_list)

# Characterizing  latent factors
# Plot Factor1 vs Factor2
p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "condition",
                  shape_by = "tissue",
                  dot_size = 6,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

# Visualize latent factor scores 
plot_factor(MOFAobject, factors = 1, color_by = "condition", shape_by = "tissue", 
            dodge = TRUE, add_violin = TRUE)
