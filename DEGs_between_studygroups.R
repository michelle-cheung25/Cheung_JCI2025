## Differentially Expressed Genes ##
## Cheung et al. JCI 2025 ##

### Set working directory, load libraries ####
setwd("path/to/projects")
library(Seurat)
library(clustree)
library(ggrepel)
####

### Comparing the transcriptome of the same cell cluster (C0 SwIg MBCs , C1 Marginal Zone, C2 Atypical MBCs) between studygroups (Anti-TNF IBD, Anti-IL-12/23 IBD, Healthy Controls) ####

# Join layers
MBC_CITEseq_obj_merged_filtered <- JoinLayers(MBC_CITEseq_obj_merged_filtered)

#C0 SwIg MBCs - TNF vs Healthy
cluster0.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Healthy Controls", ident.2 = "Anti-TNF IBD", group.by="studygroup", subset.ident="C0 SwIg MBCs")
head(cluster0.markers)
cluster0.markers$significance <- "Not significant"
cluster0.markers$significance[cluster0.markers$avg_log2FC > 0.6 & cluster0.markers$p_val_adj < 0.05] <- "Upregulated_HC"
cluster0.markers$significance[cluster0.markers$avg_log2FC < -0.6 & cluster0.markers$p_val_adj < 0.05] <- "Upregulated_TNF"
rownames <- row.names(cluster0.markers)
rownames <- as.data.frame(rownames)

cluster0.markers$gene <- ""
for (index in 1:nrow(cluster0.markers)) {
  cluster0.markers$gene[index] <- rownames$rownames[index]
}
cluster0.markers$delabel <- ifelse((cluster0.markers$gene %in% head(cluster0.markers$gene, 50)) & (cluster0.markers$significance != "Not significant"), cluster0.markers$gene, NA)

ggplot(cluster0.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel, family="sans")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "black", "orange"), 
                     labels = c("Not significant", "Upregulated_HC", "Upregulated_TNF")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point(size=0.75) + theme_classic() + theme(text=element_text(family="sans"), 
                                                  axis.text.x = element_blank(), 
                                                  axis.text.y = element_blank(), 
                                                  axis.title.x = element_blank(), 
                                                  axis.title.y = element_blank()) +
  geom_text_repel(max.overlaps = Inf, color="black", size=1.5) + NoLegend()


#C0 SwIg MBCs - IL12-23 vs Healthy 
cluster0.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Healthy Controls", ident.2 = "Anti-IL12-23 IBD", group.by="studygroup", subset.ident="C0 SwIg MBCs")
head(cluster0.markers)
cluster0.markers$significance <- "Not significant"
cluster0.markers$significance[cluster0.markers$avg_log2FC > 0.6 & cluster0.markers$p_val_adj < 0.05] <- "Upregulated_HC"
cluster0.markers$significance[cluster0.markers$avg_log2FC < -0.6 & cluster0.markers$p_val_adj < 0.05] <- "Upregulated_IL12_23"
rownames <- row.names(cluster0.markers)
rownames <- as.data.frame(rownames)

cluster0.markers$gene <- ""
for (index in 1:nrow(cluster0.markers)) {
  cluster0.markers$gene[index] <- rownames$rownames[index]
}
cluster0.markers$delabel <- ifelse((cluster0.markers$gene %in% head(cluster0.markers$gene, 50)) & (cluster0.markers$significance != "Not significant"), cluster0.markers$gene, NA)


ggplot(cluster0.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel, family="sans")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "black", "blue"), 
                     labels = c("Not significant", "Upregulated_HC", "Upregulated_IL12_23")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point(size=0.75) + theme_classic() + theme(text=element_text(family="sans"), 
                                                  axis.text.x = element_blank(), 
                                                  axis.text.y = element_blank(), 
                                                  axis.title.x = element_blank(), 
                                                  axis.title.y = element_blank()) +
  geom_text_repel(max.overlaps = Inf, color="black", size=1.5) + NoLegend()

ggsave("C:/Users/miche/OneDrive - University of Toronto/IMPACT/MBC scRNA-seq/Affinity Maturation B Cell Paper/Figures/Supplementary Figure 3 GEX comparisons/Volcano_SwIg_HCvsIL12.tiff", 
       width = 3, height = 2, device='tiff', dpi=1200)



#C0 SwIg MBCs - TNF vs IL12-23 
cluster0.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Anti-TNF IBD", ident.2 = "Anti-IL12-23 IBD", group.by="studygroup", subset.ident="C0 SwIg MBCs")
head(cluster0.markers)
cluster0.markers$significance <- "Not significant"
cluster0.markers$significance[cluster0.markers$avg_log2FC > 0.6 & cluster0.markers$p_val_adj < 0.05] <- "Upregulated_TNF"
cluster0.markers$significance[cluster0.markers$avg_log2FC < -0.6 & cluster0.markers$p_val_adj < 0.05] <- "Upregulated_IL12+23"
rownames <- row.names(cluster0.markers)
rownames <- as.data.frame(rownames)

cluster0.markers$gene <- ""
for (index in 1:nrow(cluster0.markers)) {
  cluster0.markers$gene[index] <- rownames$rownames[index]
}
cluster0.markers$delabel <- ifelse((cluster0.markers$gene %in% head(cluster0.markers$gene, 50)) & (cluster0.markers$significance != "Not significant"), cluster0.markers$gene, NA)

ggplot(cluster0.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel, family="sans")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "blue", "orange"), 
                     labels = c("Not significant", "Upregulated_IL12_23", "Upregulated_TNF")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point(size=0.75) + theme_classic() + theme(text=element_text(family="sans"), 
                                                  axis.text.x = element_blank(), 
                                                  axis.text.y = element_blank(), 
                                                  axis.title.x = element_blank(), 
                                                  axis.title.y = element_blank()) +
  geom_text_repel(max.overlaps = Inf, color="black", size=1.5) + NoLegend()


#C1 Marginal Zone MBCs - TNF vs Healthy
MBC_CITEseq_obj_merged_filtered <- JoinLayers(MBC_CITEseq_obj_merged_filtered)

cluster1.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Healthy Controls", ident.2 = "Anti-TNF IBD", group.by="studygroup", subset.ident="C1 Marginal Zone MBCs")
head(cluster1.markers)
cluster1.markers$significance <- "Not significant"
cluster1.markers$significance[cluster1.markers$avg_log2FC > 0.6 & cluster1.markers$p_val_adj < 0.05] <- "Upregulated_HC"
cluster1.markers$significance[cluster1.markers$avg_log2FC < -0.6 & cluster1.markers$p_val_adj < 0.05] <- "Upregulated_TNF"
rownames <- row.names(cluster1.markers)
rownames <- as.data.frame(rownames)

cluster1.markers$gene <- ""
for (index in 1:nrow(cluster1.markers)) {
  cluster1.markers$gene[index] <- rownames$rownames[index]
}
cluster1.markers$delabel <- ifelse((cluster1.markers$gene %in% head(cluster1.markers$gene, 50)) & (cluster1.markers$significance != "Not significant"), cluster1.markers$gene, NA)

ggplot(cluster1.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel, family="sans")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "black", "orange"), 
                     labels = c("Not significant", "Upregulated_HC", "Upregulated_TNF")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point(size=0.75) + theme_classic() + theme(text=element_text(family="sans"), 
                                                  axis.text.x = element_blank(), 
                                                  axis.text.y = element_blank(), 
                                                  axis.title.x = element_blank(), 
                                                  axis.title.y = element_blank()) +
  geom_text_repel(max.overlaps = Inf, color="black", size=1.5) + NoLegend()


#C1 Marginal Zone MBCs - IL12-23 vs Healthy
cluster1.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Healthy Controls", ident.2 = "Anti-IL12-23 IBD", group.by="studygroup", subset.ident="C1 Marginal Zone MBCs")
head(cluster1.markers)
cluster1.markers$significance <- "Not significant"
cluster1.markers$significance[cluster1.markers$avg_log2FC > 0.6 & cluster1.markers$p_val_adj < 0.05] <- "Upregulated_HC"
cluster1.markers$significance[cluster1.markers$avg_log2FC < -0.6 & cluster1.markers$p_val_adj < 0.05] <- "Upregulated_IL12_23"
rownames <- row.names(cluster1.markers)
rownames <- as.data.frame(rownames)

cluster1.markers$gene <- ""
for (index in 1:nrow(cluster1.markers)) {
  cluster1.markers$gene[index] <- rownames$rownames[index]
}
cluster1.markers$delabel <- ifelse((cluster1.markers$gene %in% head(cluster1.markers$gene, 50)) & (cluster1.markers$significance != "Not significant"), cluster1.markers$gene, NA)

ggplot(cluster1.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel, family="sans")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "black", "blue"), 
                     labels = c("Not significant", "Upregulated_HC", "Upregulated_IL12_23")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point(size=0.75) + theme_classic() + theme(text=element_text(family="sans"), 
                                                  axis.text.x = element_blank(), 
                                                  axis.text.y = element_blank(), 
                                                  axis.title.x = element_blank(), 
                                                  axis.title.y = element_blank()) +
  geom_text_repel(max.overlaps = Inf, color="black", size=1.5) + NoLegend()


#C1 Marginal Zone MBCs - IL12-23 vs Anti-TNF
cluster1.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Anti-TNF IBD", ident.2 = "Anti-IL12-23 IBD", group.by="studygroup", subset.ident="C1 Marginal Zone MBCs")
head(cluster1.markers)
cluster1.markers$significance <- "Not significant"
cluster1.markers$significance[cluster1.markers$avg_log2FC > 0.6 & cluster1.markers$p_val_adj < 0.05] <- "Upregulated_TNF"
cluster1.markers$significance[cluster1.markers$avg_log2FC < -0.6 & cluster1.markers$p_val_adj < 0.05] <- "Upregulated_IL12_23"
rownames <- row.names(cluster1.markers)
rownames <- as.data.frame(rownames)

cluster1.markers$gene <- ""
for (index in 1:nrow(cluster1.markers)) {
  cluster1.markers$gene[index] <- rownames$rownames[index]
}
cluster1.markers$delabel <- ifelse((cluster1.markers$gene %in% head(cluster1.markers$gene, 50)) & (cluster1.markers$significance != "Not significant"), cluster1.markers$gene, NA)

ggplot(cluster1.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel, family="sans")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "blue", "orange"), 
                     labels = c("Not significant", "Upregulated_IL12_23", "Upregulated_TNF")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point(size=0.75) + theme_classic() + theme(text=element_text(family="sans"), 
                                                  axis.text.x = element_blank(), 
                                                  axis.text.y = element_blank(), 
                                                  axis.title.x = element_blank(), 
                                                  axis.title.y = element_blank()) +
  geom_text_repel(max.overlaps = Inf, color="black", size=1.5) + NoLegend()

#C2 Atypical MBCs - TNF vs Healthy
cluster2.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Healthy Controls", ident.2 = "Anti-TNF IBD", group.by="studygroup", subset.ident="C2 Atypical MBCs")
head(cluster2.markers)
cluster2.markers$significance <- "Not significant"
cluster2.markers$significance[cluster2.markers$avg_log2FC > 0.6 & cluster2.markers$p_val_adj < 0.05] <- "Upregulated_HC"
cluster2.markers$significance[cluster2.markers$avg_log2FC < -0.6 & cluster2.markers$p_val_adj < 0.05] <- "Upregulated_TNF"
rownames <- row.names(cluster2.markers)
rownames <- as.data.frame(rownames)

cluster2.markers$gene <- ""
for (index in 1:nrow(cluster2.markers)) {
  cluster2.markers$gene[index] <- rownames$rownames[index]
}
cluster2.markers$delabel <- ifelse((cluster2.markers$gene %in% head(cluster2.markers$gene, 50)) & (cluster2.markers$significance != "Not significant"), cluster2.markers$gene, NA)

ggplot(cluster2.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "black", "orange"), 
                     labels = c("Not significant", "Upregulated_HC", "Upregulated_TNF")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point() + theme_classic() + geom_text_repel(max.overlaps = Inf, color="black", size=3) + NoLegend()

#C2 Atypical MBCs - IL12/23 vs Healthy 
cluster2.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Healthy Controls", ident.2 = "Anti-IL12-23 IBD", group.by="studygroup", subset.ident="C2 Atypical MBCs")
head(cluster2.markers)
cluster2.markers$significance <- "Not significant"
cluster2.markers$significance[cluster2.markers$avg_log2FC > 0.6 & cluster2.markers$p_val_adj < 0.05] <- "Upregulated_HC"
cluster2.markers$significance[cluster2.markers$avg_log2FC < -0.6 & cluster2.markers$p_val_adj < 0.05] <- "Upregulated_IL12_23"
rownames <- row.names(cluster2.markers)
rownames <- as.data.frame(rownames)

cluster2.markers$gene <- ""
for (index in 1:nrow(cluster2.markers)) {
  cluster2.markers$gene[index] <- rownames$rownames[index]
}
cluster2.markers$delabel <- ifelse((cluster2.markers$gene %in% head(cluster2.markers$gene, 50)) & (cluster2.markers$significance != "Not significant"), cluster2.markers$gene, NA)

ggplot(cluster2.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "black", "blue"), 
                     labels = c("Not significant", "Upregulated_HC", "Upregulated_IL12_23")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point() + theme_classic() + geom_text_repel(max.overlaps = Inf, color="black", size=3) + NoLegend()

#C2 Atypical MBCs - IL12/23 vs TNF 
cluster2.markers <- FindMarkers(MBC_CITEseq_obj_merged_filtered, ident.1 = "Anti-IL12-23 IBD", ident.2 = "Anti-TNF IBD", group.by="studygroup", subset.ident="C2 Atypical MBCs")
head(cluster2.markers)
cluster2.markers$significance <- "Not significant"
cluster2.markers$significance[cluster2.markers$avg_log2FC > 0.6 & cluster2.markers$p_val_adj < 0.05] <- "Upregulated_IL12_23"
cluster2.markers$significance[cluster2.markers$avg_log2FC < -0.6 & cluster2.markers$p_val_adj < 0.05] <- "Upregulated_TNF"
rownames <- row.names(cluster2.markers)
rownames <- as.data.frame(rownames)

cluster2.markers$gene <- ""
for (index in 1:nrow(cluster2.markers)) {
  cluster2.markers$gene[index] <- rownames$rownames[index]
}
cluster2.markers$delabel <- ifelse((cluster2.markers$gene %in% head(cluster2.markers$gene, 50)) & (cluster2.markers$significance != "Not significant"), cluster2.markers$gene, NA)

ggplot(cluster2.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col=significance, label=delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "blue", "orange"), 
                     labels = c("Not significant", "Upregulated_IL23_23", "Upregulated_TNF")) +
  labs(x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(P-value)"), col=" ") +
  geom_point() + theme_classic() + geom_text_repel(max.overlaps = Inf, color="black", size=3) + NoLegend()


