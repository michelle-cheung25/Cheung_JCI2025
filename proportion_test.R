## Analysis of cluster proportions using permutation tests ##
## Cheung et al. JCI 2025 ##

### Set working directory, load libraries ####
setwd("path/to/projects")
library(Seurat)
library(clustree)
library(Seurat)
library(clustree)
library(ggrepel)
library(scProportionTest)
####


### Proportion Test of Clusters####

# Remove sample with Hashtag 12 (not matched across Timepoints 4 and 6)
MBC_CITEseq_obj_merged_filtered_matched <- subset(MBC_CITEseq_obj_merged_filtered, subset = HTO_maxID != "TotalSeq-C-0262-anti-human-Hashtag-12")

# For Timepoint 6, no atypical MBCs were detected in healthy controls so create object without atypical MBCs for proportion comparisons 
MBC_CITEseq_obj_merged_filtered$cluster_ident <- as.character(MBC_CITEseq_obj_merged_filtered$cluster_ident)
MBC_CITEseq_obj_merged_filtered_no_atMBCs <- subset(MBC_CITEseq_obj_merged_filtered, idents = c("C0 SwIg MBCs", "C1 Marginal Zone MBCs"))

#Proportion test of clusters 
prop_test_T4 <- sc_utils(subset(MBC_CITEseq_obj_merged_filtered_matched, subset= Timepoint==4))
prop_test_T6 <- sc_utils(subset(MBC_CITEseq_obj_merged_filtered_matched, subset= Timepoint==6))
prop_test_T6_no_atMBC <- sc_utils(subset(MBC_CITEseq_obj_merged_filtered_no_atMBCs, subset = Timepoint==6))

#Timepoint 4: Anti-TNF vs HC 
prop_test_TNF_HC <- permutation_test(
  prop_test_T4, cluster_identity = "cluster_ident",
  sample_1 = "Healthy Controls", sample_2 = "Anti-TNF IBD",
  sample_identity = "studygroup"
)

permutation_plot(prop_test_TNF_HC, order_clusters=FALSE) + scale_color_manual(values=c("black", "grey")) + labs(x = "", y = "Observed Log2(Fold Difference)") + NoLegend() +
  theme(text = element_text(family="Arial")) + ylim(-4, 4)

#Timepoint 4: Anti-IL-12/23 vs HC 
prop_test_IL12_23_HC <- permutation_test(
  prop_test_T4, cluster_identity = "cluster_ident",
  sample_1 = "Healthy Controls", sample_2 = "Anti-IL12-23 IBD",
  sample_identity = "studygroup"
)

permutation_plot(prop_test_IL12_23_HC, order_clusters=FALSE) + scale_color_manual(values=c("grey")) + labs(x = "", y = "Observed Log2(Fold Difference)") + NoLegend() +
  theme(text = element_text(family="Arial")) + ylim(-4, 4)


#Timepoint 4: Anti-TNF vs Anti-IL-12/23
prop_test_TNF_IL12_23 <- permutation_test(
  prop_test_T4, cluster_identity = "cluster_ident",
  sample_1 = "Anti-IL12-23 IBD", sample_2 = "Anti-TNF IBD",
  sample_identity = "studygroup"
)

permutation_plot(prop_test_TNF_IL12_23, order_clusters=FALSE) + scale_color_manual(values=c("black", "grey")) + labs(x = "", y = "Observed Log2(Fold Difference)") + NoLegend() +
  theme(text = element_text(family="serif")) + ylim(-4, 4)


#Timepoint 6: Anti-TNF vs HC 
prop_test_TNF_HC_T6 <- permutation_test(
  prop_test_T6_no_atMBC, cluster_identity = "cluster_ident",
  sample_1 = "Healthy Controls", sample_2 = "Anti-TNF IBD",
  sample_identity = "studygroup"
)

permutation_plot(prop_test_TNF_HC_T6, order_clusters=FALSE) + scale_color_manual(values=c("black", "grey")) + labs(x = "", y = "Observed Log2(Fold Difference)") + NoLegend() +
  theme(text = element_text(family="serif")) + ylim(-4, 4)

#Timepoint 6: Anti-IL12/23 vs HC 
prop_test_IL12_23_HC_T6 <- permutation_test(
  prop_test_T6_no_atMBC, cluster_identity = "cluster_ident",
  sample_1 = "Healthy Controls", sample_2 = "Anti-IL12-23 IBD",
  sample_identity = "studygroup"
)

permutation_plot(prop_test_IL12_23_HC_T6, order_clusters=FALSE) + scale_color_manual(values=c("grey")) + labs(x = "", y = "Observed Log2(Fold Difference)") + NoLegend() +
  theme(text = element_text(family="serif")) + ylim(-4, 4)


#Timepoint 6: Anti-TNF vs Anti-IL-12/23 
prop_test_TNF_IL12_23_T6 <- permutation_test(
  prop_test_T6, cluster_identity = "cluster_ident",
  sample_1 = "Anti-IL12-23 IBD", sample_2 = "Anti-TNF IBD",
  sample_identity = "studygroup"
)

permutation_plot(prop_test_TNF_IL12_23_T6, order_clusters=FALSE) + scale_color_manual(values=c("black")) + labs(x = "", y = "Observed Log2(Fold Difference)") + NoLegend() +
  theme(text = element_text(family="serif")) + ylim(-4, 5)
