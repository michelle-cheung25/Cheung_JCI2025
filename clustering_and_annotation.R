## Clustering and Annotation of GEX Data ##
## Cheung et al. JCI 2025 ##

#### Set working directory, load libraries ####
setwd("path/to/projects")
library(Seurat)
library(clustree)
library(Seurat)
library(clustree)
library(ggrepel)
library(scProportionTest)
library(harmony)
####


#### Finding Variable Features ####
MBC_CITEseq_obj_merged <- FindVariableFeatures(MBC_CITEseq_obj_merged, 
                                               selection.method="vst", nfeatures=2000, assay = "RNA")
#### Scaling Data and PCA Analyses #### 

# Scale Data (gene expression)
MBC_CITEseq_obj_merged <- ScaleData(MBC_CITEseq_obj_merged, assay="RNA")

clonotype_genes <- grep(pattern = "^IG[HKL][VDJ]|^IGH[MDE]|^IGHA[1-2]|^IGHG[1-4]|^IGKC|^IGLC[1-7]",
                        rownames(MBC_CITEseq_obj_merged),
                        value=TRUE, invert=FALSE) #get list of genes not related to clonotype specific gene usage 

# Run PCA, regressing out clonotype genes so cells cluster by transcriptome and not isotype
MBC_CITEseq_obj_merged <- RunPCA(MBC_CITEseq_obj_merged, features = VariableFeatures(object = MBC_CITEseq_obj_merged), 
                                 vars.to.regress = clonotype_genes, nfeatures.print = 10)

# Graphs the output of PCA reduction, cells are colored by identity class
DimPlot(MBC_CITEseq_obj_merged, reduction = "pca")

# Display genes heavily associated with the first two dimensions
VizDimLoadings(MBC_CITEseq_obj_merged, dims = 1:2, reduction = "pca")

# Plots the variance of principle components for identification of the elbow in the graph
ElbowPlot(MBC_CITEseq_obj_merged, ndims=40)

#### Harmony Integration post PCA Analyses #### 
MBC_CITEseq_obj_merged <- RunHarmony(MBC_CITEseq_obj_merged, "Timepoint")

#### Finding Neighbors and Clustering ####

#Based on the Elbow plot, use the first 30 PCs to compute the shared nearest neighbors  
MBC_CITEseq_obj_merged <- FindNeighbors(MBC_CITEseq_obj_merged, dims = 1:30, reduction="harmony") 

#Decide the appropriate clustering resolution 
MBC_CITEseq_obj_merged <- FindClusters(MBC_CITEseq_obj_merged, resolution = c(0, 0.2, 0.4, 0.6, 0.8, 1.2))

clustree(MBC_CITEseq_obj_merged, prefix = "RNA_snn_res.")

#Find clusters with the selected resolution of 0.6
MBC_CITEseq_obj_merged <- FindClusters(MBC_CITEseq_obj_merged, resolution = 0.6)

#### UMAP Analyses #### 

#Run UMAP dimensional reduction
MBC_CITEseq_obj_merged <- RunUMAP(MBC_CITEseq_obj_merged, dims = 1:30, reduction="harmony")

#Graph output of UMAP on a scatter plot, samples combined or split by Timepoint/studygroup
DimPlot(MBC_CITEseq_obj_merged, reduction="umap")
DimPlot(MBC_CITEseq_obj_merged, reduction="umap", split.by="Timepoint")
DimPlot(MBC_CITEseq_obj_merged, reduction="umap", split.by="studygroup")
DimPlot(MBC_CITEseq_obj_merged, reduction = "umap", split.by="Timepoint", group.by="studygroup")

#### Cluster Annotation ####

#Cluster (resolution 0.6) annotation
new.cluster.ids <- c("C0 SwIg MBCs", "C1 SwIg MBCs", "C2 IgM+ MBCs", "C3 SwIg MBCs",
                     "C4 Naive B cells", "C5 Atypical MBCs", "C6 NK Cells")
names(new.cluster.ids) <- levels(MBC_CITEseq_obj_merged)
MBC_CITEseq_obj_merged <- RenameIdents(MBC_CITEseq_obj_merged, new.cluster.ids)

#Subcluster on memory B cells: remove Naive B cells and NK cells from Seurat object
MBC_CITEseq_obj_merged_filtered <- subset(MBC_CITEseq_obj_merged, idents = c("C0 SwIg MBCs", "C1 SwIg MBCs", "C2 IgM+ MBCs", "C3 SwIg MBCs", "C5 Atypical MBCs"))

#Repeat: Scaling, PCA, harmony integration, clustering, and UMAP without Naive B cells and NK cells
MBC_CITEseq_obj_merged_filtered <- FindVariableFeatures(MBC_CITEseq_obj_merged_filtered, selection.method="vst", nfeatures=2000, assay = "RNA")
all.genes = rownames(MBC_CITEseq_obj_merged_filtered)
MBC_CITEseq_obj_merged_filtered <- ScaleData(MBC_CITEseq_obj_merged_filtered, assay="RNA", features = all.genes)
MBC_CITEseq_obj_merged_filtered <- RunPCA(MBC_CITEseq_obj_merged_filtered, features = VariableFeatures(object = MBC_CITEseq_obj_merged_filtered),
                                          vars.to.regress = clonotype_genes, nfeatures.print = 10)
MBC_CITEseq_obj_merged_filtered <- RunHarmony(MBC_CITEseq_obj_merged_filtered, "Timepoint")
MBC_CITEseq_obj_merged_filtered <- FindNeighbors(MBC_CITEseq_obj_merged_filtered, dims = 1:30) 
MBC_CITEseq_obj_merged_filtered <- FindClusters(MBC_CITEseq_obj_merged_filtered, resolution = c(0, 0.2, 0.4, 0.6, 0.8, 1.2))
clustree(MBC_CITEseq_obj_merged_filtered, prefix = "RNA_snn_res.")
MBC_CITEseq_obj_merged_filtered <- FindClusters(MBC_CITEseq_obj_merged_filtered, resolution = 0.8)
MBC_CITEseq_obj_merged_filtered <- RunUMAP(MBC_CITEseq_obj_merged_filtered, vars.to.regress = clonotype_genes, dims = 1:30, reduction="harmony")
new.cluster.ids2 <- c("C0 SwIg MBCs", "C1 Marginal Zone MBCs", "C2 Atypical MBCs")
names(new.cluster.ids2) <- levels(MBC_CITEseq_obj_merged_filtered)
MBC_CITEseq_obj_merged_filtered <- RenameIdents(MBC_CITEseq_obj_merged_filtered, new.cluster.ids2)

#Add "cluster_ident" to metadata 
MBC_CITEseq_obj_merged_filtered$cluster_ident <- Idents(MBC_CITEseq_obj_merged_filtered)
