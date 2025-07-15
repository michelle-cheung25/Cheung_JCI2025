# Loading BCR data and combining with GEX data #
# Cheung et al. JCI 2025

# BCR VDJ data can be downloaded from NCBI GEO accession GSE290006
# filtered_contig.fasta and filtered_contig_annotations.csv outputted by cellranger were used to assign V,D, and J genes using IgBlast as described here: https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html#assign-v-d-and-j-genes-using-igblast.
# In brief, V(D)J contigs were processed using the Immcantation framework (v4.5.0). Contigs were aligned to the IMGT reference using IgBLAST (v4.5.0, installed in the Docker container “immcantation/suite”).
# Outputs of alingment are noted below: 
# Timepoint 4 (post dose 2) file: GSM8803450_Timepoint4_BCR_igblast_db-pass.tsv.gz
# Timepoint 6 (post dose 3) file: GSM8803453_Timepoint6_BCR_igblast_db-pass.tsv.gz

#### Load libraries ####
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dowser))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scoper))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scRepertoire))
suppressPackageStartupMessages(library(ggtree))
####

#### Read in BCR Data & QC #### 

#Read in the data. Specify data types of non AIRR-C standard fields 
#We assign integer type to the *_length fields

T4_bcr_data <- read_rearrangement("path/to/Timepoint4 tsv/file",
                                  aux_types=c("v_germline_length"="i",
                                              "d_germline_length"="i",
                                              "j_germline_length"="i",
                                              "day"="i"))

T6_bcr_data <- read_rearrangement("path/to/Timepoint6 tsv/file", 
                                  aux_types=c("v_germline_length"="i",
                                              "d_germline_length"="i",
                                              "j_germline_length"="i",
                                              "day"="i"))

# Assign timepoint to metadata
T4_bcr_data$Timepoint <- 4
T6_bcr_data$Timepoint <- 6

# Add "T4_" or "T6_" to the front of each cell ID
T4_bcr_data$cell_id <- paste("T4_", T4_bcr_data$cell_id, sep="")
T6_bcr_data$cell_id <- paste("T6_", T6_bcr_data$cell_id, sep="")

# Append T4 and T6 BCR data into one merged tibble 
bcr_data_merged <- T4_bcr_data %>% bind_rows(T6_bcr_data)

#Remove cells with non-productive sequences (keep only productive).
bcr_data_merged <- bcr_data_merged %>% dplyr::filter(productive)

#Remove cells with multiple heavy chains
multi_heavy <- table(dplyr::filter(bcr_data_merged, locus == "IGH")$cell_id)
multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]
bcr_data_merged <- dplyr::filter(bcr_data_merged, !cell_id %in% multi_heavy_cells)

#Note which cells have multiple light chains (for paired H and L chain gene usage)
kappa <- dplyr::filter(bcr_data_merged, locus == "IGK")$cell_id
lambda <- dplyr::filter(bcr_data_merged, locus == "IGL")$cell_id
total <- append(kappa, lambda)
total_table <- table(total)
multi_light_cells <- names(total_table)[total_table > 1]
bcr_data_merged$multi_light <- 0
bcr_data_merged$multi_light[which(bcr_data_merged$cell_id %in% multi_light_cells)] <- 1
#1 indicates has multiple light chains 

#Remove cells without heavy chains (only have light chains)
#Split cells by heavy and light chains
heavy_cells <- dplyr::filter(bcr_data_merged, locus == "IGH")$cell_id
light_cells <- dplyr::filter(bcr_data_merged, locus == "IGK" | locus == "IGL")$cell_id
no_heavy_cells <- light_cells[which(!light_cells %in% heavy_cells)]

bcr_data <- dplyr::filter(bcr_data_merged, !cell_id %in% no_heavy_cells)

#### Load GEX data from Seurat object ####
gex_data <- readRDS("path/to/RDS/file/containing/Seurat/object/GEXdata")

#Match GEX and BCR cell ids - match step to identify the GEX data in the BCR object - some BCRs don't have GEX information 
#Due to BCRs that are covered but didn't pass GEX processing and QC thresholds 

#Match indices to find position of the BCR cells in the GEX data 
match.index <- match(bcr_data_merged$cell_id, Cells(gex_data))

#Proportion of BCRs that do not have GEX information 
mean(is.na(match.index))

#Transfer cell type annotations into the BCR data 
#Add annotations to BCR data 
cell.annotation <- as.character(Idents(gex_data))
bcr_data_merged$gex_annotation <-
  unlist(lapply(match.index, function(x) {
    ifelse(is.na(x), NA, cell.annotation[x])
  }))
bcr_data_merged$gex_annotation[1:10]

#Transfer hashtag annotations into the BCR data 
hashtag.annotation <- as.character(gex_data@meta.data$HTO_maxID)

bcr_data_merged$hashtag_annotation <- 
  unlist(lapply(match.index, function(x) {
    ifelse(is.na(x), NA, hashtag.annotation[x])
  }))
bcr_data_merged$hashtag_annotation[1:10]

#Transfer studygroup annotations into the BCR data 
studygroup.annotation <- as.character(gex_data@meta.data$studygroup)

bcr_data_merged$studygroup_annotation <- 
  unlist(lapply(match.index, function(x) {
    ifelse(is.na(x), NA, studygroup.annotation[x])
  }))
bcr_data_merged$studygroup_annotation[1:10]

#Remove cells from BCR data that lack corresponding partners in GEX data 
bcr_data_merged <- dplyr::filter(bcr_data_merged, !is.na(gex_annotation))
