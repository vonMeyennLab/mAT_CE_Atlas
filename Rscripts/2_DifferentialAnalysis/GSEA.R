#### GSEA 
library(SCpubr)
library(enrichR)
library(Seurat)
library(readr)
library(openxlsx)
library(readxl)
library(janitor)
require(openxlsx)
library(ggpubr)
library(data.table)

#### load functions

suppressMessages({
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  options(modEnrichR.use = TRUE)
  options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))
  
  # Set the search to Human genes.
  enrichR::setEnrichrSite(site = "Enrichr")
  
  websiteLive <- TRUE
  dbs <- enrichR::listEnrichrDbs()
  # Get all the possible databases to query.
  dbs <- sort(dbs$libraryName)
})
### choose which ones you would like to query; we often query several
dbs_use <- c("Reactome_2022",
             "GO_Biological_Process_2021", 
             "GO_Cellular_Component_2021", 
             "GO_Molecular_Function_2021",
             "WikiPathways_2019_Mouse")

### function to read output later
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#### load DE promoter and enhancer files from comparisons for speciic hPTMs

H3K27me3 <- read_csv("H3K27me3_sigDEGPr_ingAT_vs_epiAT.csv")
enhancers_down <- read_csv("H3K4me1_ingATvsepiAT_downPeaks.csv")
## if needed you can also subset these dfs again for FDR or p-value if not done before 

enriched_terms <- enrichR::enrichr(unlist(subset(H3K27me3, H3K27me3$regulation == "down")$external_gene_name), dbs_use)
write.xlsx(enriched_terms,file = "~/ingATvepiAT/H3K27me3_down_ingAT.xlsx")
SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms, plot.title = NULL, 
                              plot.subtitle = "H3K27me3 enriched in ingAT adipocytes over epiAT adipocytes",  
                              text_labels_size = 4, colors.use = c("red", "#FFFACD"), nchar_wrap = 15, nterms = 8)

#### enhancers
enriched_terms <- enrichR::enrichr(unlist(subset(enhancers_down, abs(enhancers_down$distanceToTSS)<20000)$SYMBOL), dbs_use)
write.xlsx(enriched_terms,file = "~/ingATvepiAT/H3K4me1_enhancers_down_ingAT.xlsx")
## usually we save all outputs first and then proceed to plotting

## enhancers brown and beige (final lists are derived from intersects of DE analysis)
## load GSEA files (were done with 20 000 bp cutoff)

X<- read_excel_allsheets("~/public/LauraHinte/BRITE/Enhancers/DiffAnalysis/Ucp1ERCre_Enhancers/H3K4me1/BATvBeige/Brown_Enhancers_GO.xlsx")
SCpubr::do_TermEnrichmentPlot(X, plot.title = NULL, 
                              #plot.subtitle = "brown specific enhancers",  
                              text_labels_size = 3, colors.use = c("red", "#FFFACD"), nchar_wrap = 15, nterms = 5)

Y<- read_excel_allsheets("~/public/LauraHinte/BRITE/Enhancers/DiffAnalysis/Ucp1ERCre_Enhancers/H3K4me1/BATvBeige/Beige_Enhancers_GO.xlsx")
SCpubr::do_TermEnrichmentPlot(Y, plot.title = NULL, 
                              #plot.subtitle = "beige specific enhancers",  
                              text_labels_size = 3, colors.use = c("red", "#FFFACD"), nchar_wrap = 15, nterms = 5)



