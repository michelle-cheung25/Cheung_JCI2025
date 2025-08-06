# Cheung_JCI2025
Code used to analyze single-cell RNA-sequencing data (GEX and BCR) associated with the following manuscript:

Cheung et al. Reduced vaccine-induced germinal center outputs in inflammatory bowel disease patients treated with anti-TNF biologics. Journal of Clinical Investigation 2025. https://doi.org/10.1172/JCI192589.

Raw GEX and BCR data can be downloaded from NCBI GEO under accession number: GSE290006.

## GEX data analysis (proceed in the order listed below)
### Load and cluster GEX data 
load_and_QR.R
### Clustering of cells and annotation 
clustering_and_annotation.R
### DotPlots and Violin Plots
dotplots_violinplots.R
### Proportion Test
proportion_test.R
### Comparing the transcriptome of cell clusters between study groups
DEGs_between_studygroups.R

## BCR VDJ data analysis (proceed in the order listed below)
### Loading BCR data and combining with GEX data 
BCR_load_and_mergeGEX.R
### Assigning clonal groups and reconstructing germline sequences 
BCR_assign_clonal_groups.R
### BCR somatic hypermutation and gene usage analyses 
BCR_SHM_gene_usage_analyses.R


Questions? Contact Michelle Cheung (michellew.cheung@mail.utoronto.ca).
