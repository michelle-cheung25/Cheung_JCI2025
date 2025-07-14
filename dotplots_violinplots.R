## DotPlots and Violin Plots ##
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
library(extrafont)
####

# Dotplots
loadfonts(device = "win")

genes = c("CD19", "MS4A1", "CD27", "CD38", "CD40", "CD69", "CD71", "CD80", "CD86", "CD180", "CR2", 
                                        "IGHD", "IGHM", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4",
                                        "BACH2", "TCL1A", "SPRY1", "FCER2", "IL4R", "COCH", "HOPX", "TOX", "SELL",
                                        "TNFRSF13B", "MZB1", "GPR183", "CD1C", "CD24", "ZBTB20", "STAT1", "NOTCH2", "MYC", 
                                        "PLD4","PLEK", "FGR", "CIB1","SOX5", "TBX21", "ITGAX", "FCRLA", "FCRL1", "FCRL2", "FCRL3", "FCRL5")

DotPlot(object = MBC_CITEseq_obj_merged_filtered, 
        features = genes, cols=c("lightgrey", "#280137"), dot.scale = 2.8) + RotatedAxis() + coord_flip() +
  xlab("") + ylab("") + theme(text = element_text(family="Arial", size=8), axis.text.x = element_blank(), axis.text.y = element_text(size=8)) 


# Violin Plots
VlnPlot(MBC_CITEseq_obj_merged_filtered, features = c("CD19", "MS4A1", "CD38", "CD27", "IGHD", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1"), stack=TRUE) + scale_fill_brewer(palette="Paired")


