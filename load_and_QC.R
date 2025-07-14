## Loading and QC of GEX Data ##
## Cheung et al. JCI 2025 ##

### Set working directory, load libraries ####
setwd("path/to/projects")
library(Seurat)
library(clustree)
library(Seurat)
library(clustree)
library(ggrepel)
library(scProportionTest)
library(harmony)
###

#### Importing and quality Control ####

# Read in CITE-seq data (10X Genomics cellranger pipeline outputs: filtered feature matrices) 
# Files can be found under NCBI GEO accession GSE290006
# T4 refers to Timepoint 4: 3-4 months post dose 2
# T6 referes to Timepoint 6: 3-4 months post dose 3

T4_MBC_CITEseq <- Read10X("20240913_LH00244_0159_A22NKN7LT3_Watts_Michelle_CITE_Seq/Watts_Michelle_24_08_01__Timepoint4_MBC_HT_Citeseq/outs/filtered_feature_bc_matrix")
T6_MBC_CITEseq <- Read10X("20240503_LH00244_0092_A22K3FHLT3_Watts_Michelle_CITE-Seq/Watts_Michelle__Timepoint6-MBC_HTCiteseq/outs/filtered_feature_bc_matrix")

# Creating Seurat object using GEX data 
T4_MBC_CITEseq_obj <- CreateSeuratObject(counts = T4_MBC_CITEseq$`Gene Expression`, project = "MBC CITE-seq")
T6_MBC_CITEseq_obj <- CreateSeuratObject(counts = T6_MBC_CITEseq$`Gene Expression`, project = "MBC CITE-seq")

# Adding HTO data (sample hashtagging) to Seurat Object
T4_MBC_CITEseq_obj[['HTO']] <- CreateAssayObject(counts = T4_MBC_CITEseq$`Antibody Capture`)
T6_MBC_CITEseq_obj[['HTO']] <- CreateAssayObject(counts = T6_MBC_CITEseq$`Antibody Capture`)

# Calculating the percentage of mitochondrial genes for each donor
T4_MBC_CITEseq_obj[["percent.mt"]] <- PercentageFeatureSet(T4_MBC_CITEseq_obj, pattern="^MT-", assay="RNA")
T6_MBC_CITEseq_obj[["percent.mt"]] <- PercentageFeatureSet(T6_MBC_CITEseq_obj, pattern="^MT-", assay="RNA")

# QC: Keep only cells with 200 < features < 6500 and <10% mitochondrial genes
T4_MBC_CITEseq_obj <- subset(T4_MBC_CITEseq_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 10)
T6_MBC_CITEseq_obj <- subset(T6_MBC_CITEseq_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 10)

#### Data normalization #### 

# Set the seed to enable statistical learning methods to produce the same solution.
set.seed(123)

# Normalize GEX data with LogNormalize method
T4_MBC_CITEseq_obj <- NormalizeData(T4_MBC_CITEseq_obj, assay = "RNA")
T6_MBC_CITEseq_obj <- NormalizeData(T6_MBC_CITEseq_obj, assay = "RNA")

# Normalize HTO data with CLR method, margin = 2 normalize across cells (vs normalizing across tags margin = 1)
T4_MBC_CITEseq_obj <- NormalizeData(T4_MBC_CITEseq_obj, assay = "HTO", normalization.method = "CLR", margin = 2)
T6_MBC_CITEseq_obj <- NormalizeData(T6_MBC_CITEseq_obj, assay = "HTO", normalization.method = "CLR", margin = 2)


#### De-multiplexing samples based on TotalSeq-C Hashtags ####
T4_MBC_CITEseq_obj <- HTODemux(T4_MBC_CITEseq_obj, assay = "HTO", positive.quantile = 0.999, seed=25)
T6_MBC_CITEseq_obj <- HTODemux(T6_MBC_CITEseq_obj, assay = "HTO", positive.quantile = 0.999, seed=25)

#Group cells based on the max HTO signal 
Idents(T4_MBC_CITEseq_obj) <- "HTO_maxID"
Idents(T6_MBC_CITEseq_obj) <- "HTO_maxID"

# Keep only cells with a minimum HTO_margin to be 2.0 between HTO maxID and secondary HTO maxID to differentiate doublets 
T4_MBC_CITEseq_obj_filtered <- subset(T4_MBC_CITEseq_obj, subset = HTO_margin >= 2.0)
T6_MBC_CITEseq_obj_filtered <- subset(T6_MBC_CITEseq_obj, subset = HTO_margin >= 2.0)

#### Addition of metadata to objects  ####

# Timepoint 4 
# Hashtags 1-3 - Healthy controls
# Hashtags 4-7 - IL12/23 IBD
# Hashtags 8-12 - TNF IBD
# Adding "studygroup" annotation to each cell

T4_MBC_CITEseq_obj_filtered[['studygroup']] <- "Healthy Controls"

for (index in 1:length(T4_MBC_CITEseq_obj_filtered$HTO_maxID)) {
  
  if ((T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0254-anti-human-Hashtag-4") |
      (T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0255-anti-human-Hashtag-5") | 
      (T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0256-anti-human-Hashtag-6") |
      (T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0257-anti-human-Hashtag-7")) {
    T4_MBC_CITEseq_obj_filtered$studygroup[index] <- "Anti-IL12-23 IBD" }
  else if ((T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0258-anti-human-Hashtag-8") |
           (T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0259-anti-human-Hashtag-9") | 
           (T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0260-anti-human-Hashtag-10") |
           (T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0296-anti-human-Hashtag-11") |
           (T4_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0262-anti-human-Hashtag-12")) {
    T4_MBC_CITEseq_obj_filtered$studygroup[index] <- "Anti-TNF IBD" }
}

# Timepoint 6
# Hashtags 1-3 - Healthy controls
# Hashtags 4-7 - IL12/23 IBD
# Hashtags 8-11 - TNF IBD
# Adding "studygroup" annotation to each cell


T6_MBC_CITEseq_obj_filtered[['studygroup']] <- "Healthy Controls"

for (index in 1:length(T6_MBC_CITEseq_obj_filtered$HTO_maxID)) {
  
  if ((T6_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0254-anti-human-Hashtag-4") |
      (T6_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0255-anti-human-Hashtag-5") | 
      (T6_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0256-anti-human-Hashtag-6") |
      (T6_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0257-anti-human-Hashtag-7")) {
    T6_MBC_CITEseq_obj_filtered$studygroup[index] <- "Anti-IL12-23 IBD" }
  else if ((T6_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0258-anti-human-Hashtag-8") |
           (T6_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0259-anti-human-Hashtag-9") | 
           (T6_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0260-anti-human-Hashtag-10") |
           (T6_MBC_CITEseq_obj_filtered$HTO_maxID[index] == "TotalSeq-C-0296-anti-human-Hashtag-11")) {
    T6_MBC_CITEseq_obj_filtered$studygroup[index] <- "Anti-TNF IBD" }
}

# Add "Timepoint" annotation to each cell
T4_MBC_CITEseq_obj_filtered$Timepoint <- 4
T6_MBC_CITEseq_obj_filtered$Timepoint <- 6

MBC_CITEseq_obj_merged_filtered$cluster_ident <- Idents(MBC_CITEseq_obj_merged_filtered)

MBC_CITEseq_obj_merged_filtered$donor_ID <- ""
for (index in 1:length(MBC_CITEseq_obj_merged_filtered$HTO_maxID)) {
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0251-anti-human-Hashtag-1") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "HC-01"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0252-anti-human-Hashtag-2") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "HC-02"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0253-anti-human-Hashtag-3") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "HC-03"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0254-anti-human-Hashtag-4") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "IL12_23-01"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0255-anti-human-Hashtag-5") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "IL12_23-02"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0256-anti-human-Hashtag-6") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "IL12_23-03"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0257-anti-human-Hashtag-7") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "IL12_23-04"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0258-anti-human-Hashtag-8") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-01"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0259-anti-human-Hashtag-9") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-02"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0260-anti-human-Hashtag-10") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-03"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0296-anti-human-Hashtag-11") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-04"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0262-anti-human-Hashtag-12") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-05"
  }
}


#### Merge T4 and T6 MBC objects ####
MBC_CITEseq_obj_merged <- merge(T4_MBC_CITEseq_obj_filtered, y = T6_MBC_CITEseq_obj_filtered,
                                merge.data = TRUE, add.cell.ids = c("T4", "T6"))

#### Add Donor_ID to metadata ####

MBC_CITEseq_obj_merged_filtered$donor_ID <- ""
for (index in 1:length(MBC_CITEseq_obj_merged_filtered$HTO_maxID)) {
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0251-anti-human-Hashtag-1") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "HC-01"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0252-anti-human-Hashtag-2") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "HC-02"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0253-anti-human-Hashtag-3") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "HC-03"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0254-anti-human-Hashtag-4") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "IL12_23-01"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0255-anti-human-Hashtag-5") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "IL12_23-02"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0256-anti-human-Hashtag-6") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "IL12_23-03"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0257-anti-human-Hashtag-7") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "IL12_23-04"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0258-anti-human-Hashtag-8") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-01"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0259-anti-human-Hashtag-9") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-02"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0260-anti-human-Hashtag-10") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-03"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0296-anti-human-Hashtag-11") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-04"
  }
  if (MBC_CITEseq_obj_merged_filtered$HTO_maxID[index] == "TotalSeq-C-0262-anti-human-Hashtag-12") {
    MBC_CITEseq_obj_merged_filtered$donor_ID[index] <- "TNF-05"
  }
}
