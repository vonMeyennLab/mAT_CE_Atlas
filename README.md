# An Epigenome Atlas of Mouse Adipocytes

This respository contains code and files related to the manuscript [An Epigenome Atlas of Mouse Adipocytes](https://doi.org/10.1016/j.molmet.2025.102197).

## Abstract
Epigenetic modifications can influence gene expression in adipocytes, potentially contributing to metabolic dysfunctions, obesity, and insulin resistance. Despite recent advances in the characterization of the mouse adipocyte epigenome, epigenetic characterization of adipocytes in vivo has been challenging, particularly across different adipose depots and several epigenetic modifications. Here we generate paired single mouse datasets and conduct an integrative multimodal analysis of five histone marks – H3K4me3, H3K27me3, H3K4me1, H3K27ac and H3K9me3 – in beige, brown, and white adipocytes from three distinct mouse adipose tissue depots. Our analysis reveals that enhancers distinguish adipocytes by their tissue of origin, with H3K4me1 deposition differentiating between beige and brown adipocytes. Diphtheria toxin-mediated ablation of beige adipocytes shows that non-beigeing white adipocytes in inguinal adipose tissue and beige adipocytes are not inherently epigenetically different suggesting that they share a common developmental progenitor. These paired multimodal data comprise an extensive resource for the further exploration of the mouse adipocyte epigenome which will enable discovery of regulatory elements governing adipocyte identity and gene regulation. 


## Contents of this Repository
#### 1. :file_folder: ```Rscripts```</p>
&emsp;&emsp;:file_folder: ```Epigenetics ```</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```0_PreProcessing```&ensp;*Peak QC, annotation, filtering and quantification*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```1_MOFA```&ensp;*Cut&Tag based multi-omics factor analysis*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```2_DifferentialAnalysis```&ensp;*Differential analysis of promoters, enhancers, and GSEA*</p>

#### 2. :file_folder: ```ChromHMM```</p>
&emsp;&emsp;:file_folder: ```ChromHMM```&ensp;*Commands and outputs of ChroHMM analysis (enhancers)*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```tracks```&ensp;*Enhancer bed files for adipocytes*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```scripts```&ensp;*Commands for Cut&Tag based ChromHMM analysis*</p>

#### 3. :file_folder: ```Enhancer_Tracks```</p>
Contains bed files for enhancers from ChromHMM analysis for white (AdipoCre) and beige/brown (Ucp1ERCre) adipocytes. 

## Associated Repositories 
1. [NextFlow Pipeline for CUT&Tag](https://github.com/vonMeyennLab/nf_cutntag)


## DEGs
[App to explore differential enrichment results for promoters and enhancers](https://nme.ethz.ch/mAT_CEAtlas.html)

## Accession Codes
GEO: GSE262213
