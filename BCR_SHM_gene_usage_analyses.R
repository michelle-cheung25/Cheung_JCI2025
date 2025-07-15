# BCR Somatic hypermutation and gene usage analyses
# Cheung et al. JCI2025

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
suppressPackageStartupMessages(library(tidyverse))
supppressPackageStartupMuessage(library(reshape2))
####


#### BCR SMH Analyses ####

#Calculate mutational FREQUENCY in the V gene CDR/FRW region of heavy and light chains - includes CDR3 or FWR4
results_heavy <- dplyr::filter(clonal_results, locus == "IGH")
data_mut <- shazam::observedMutations(results_heavy,
                                      sequenceColumn = "sequence_alignment",
                                      germlineColumn = "germline_alignment_d_mask",
                                      regionDefinition = IMGT_VDJ,
                                      frequency = TRUE,
                                      combine = TRUE,
                                      nproc = 1)
data_mut <- data_mut %>% mutate(Timepoint = case_when(
  Timepoint == 4  ~ 'Timepoint 4', 
  Timepoint == 6  ~ 'Timepoint 6')) 


results_light <- dplyr::filter(clonal_results, locus != "IGH")
data_mut_light <- shazam::observedMutations(results_light,
                                            sequenceColumn = "sequence_alignment",
                                            germlineColumn = "germline_alignment_d_mask",
                                            regionDefinition = IMGT_VDJ,
                                            frequency = TRUE,
                                            combine = TRUE,
                                            nproc = 1)
data_mut_light <- data_mut_light %>% mutate(Timepoint = case_when(
  Timepoint == 4  ~ 'Timepoint 4', 
  Timepoint == 6  ~ 'Timepoint 6')) 


#Calculate replace and silent mutation FREQEUNCY in V gene of heavy and light chains, separated by REGION (CDR1/2/3 or FRW)

#mu_freq_cdr1_r, mu_freq_cdr1_s
#mu_freq_cdr2_r, mu_freq_cdr2_s
#mu_freq_cdr3_r, mu_freq_cdr3_s
#mu_freq_fwr1_r, mu_freq_fwr1_s
#mu_freq_fwr2_r, mu_freq_fwr2_s
#mu_freq_fwr3_r, mu_freq_fwr3_s
#mu_freq_fwr4_r, mu_freq_fwr4_s

results_heavy <- dplyr::filter(clonal_results, locus == "IGH")
data_mut_region <- shazam::observedMutations(results_heavy,
                                             sequenceColumn = "sequence_alignment",
                                             germlineColumn = "germline_alignment_d_mask",
                                             regionDefinition = IMGT_VDJ_BY_REGIONS,
                                             frequency = TRUE,
                                             combine = FALSE,
                                             nproc = 1)

results_light <- dplyr::filter(clonal_results, locus != "IGH")
data_mut_region_light <- shazam::observedMutations(results_light,
                                                   sequenceColumn = "sequence_alignment",
                                                   germlineColumn = "germline_alignment_d_mask",
                                                   regionDefinition = IMGT_VDJ_BY_REGIONS,
                                                   frequency = TRUE,
                                                   combine = FALSE,
                                                   nproc = 1)



#### Gene Usage Analyses ####

colors <- c("#008080", "orange", "#746e69", "black")

# Quantify heavy and light chain V gene usage by sample 
genes_heavy <- clonal_results %>% filter(str_detect(v_call, 'IGHV'))
genes_kappa <- clonal_results %>% filter(str_detect(v_call, 'IGK'))
genes_lambda <- clonal_results %>% filter(str_detect(v_call, 'IGL'))

gene_heavy_count <- countGenes(genes_heavy, gene="v_call", groups="studygroup_annotation", mode="gene")
gene_kappa_count <- countGenes(genes_kappa, gene="v_call", groups="studygroup_annotation", mode="gene")
gene_lambda_count <- countGenes(genes_lambda, gene="v_call", groups="studygroup_annotation", mode="gene")

gene_heavy_count$seq_percent <- gene_heavy_count$seq_freq * 100
gene_kappa_count$seq_percent <- gene_kappa_count$seq_freq * 100
gene_lambda_count$seq_percent <- gene_lambda_count$seq_freq * 100


# Heavy chain gene usage 

#Heavy chain V gene usage in CoV-AbDab Database- filter on human and SARS-CoV-2 specific antibodies and genes with a frequency of at least 3%
#Download CoV-AbDAb database as a csv file 
CoV_database <- read.csv("path/to/CoV-AbDab/csv/file", header=TRUE)
CoV_database_human <- CoV_database %>% filter(str_detect(Heavy.V.Gene, '(Human)'))
CoV_database_human <- CoV_database %>% filter(str_detect(Binds.to, 'SARS-CoV2'))
CoV_database_human_count <- countGenes(CoV_database_human, gene="Heavy.V.Gene", mode="gene")
CoV_database_human_count$seq_percent <- CoV_database_human_count$seq_freq * 100
CoV_database_human_count <- CoV_database_human_count %>% filter(seq_percent >= 3)
CoV_database_human_count$studygroup_annotation <- 'CoV-AbDab Database'

#Filter for sample genes with a frequency of at least 3%
gene_heavy_count_IL12_23 <- gene_heavy_count %>% filter(studygroup_annotation=="Anti-IL12-23 IBD" & seq_percent >= 3)
gene_heavy_count_TNF <- gene_heavy_count %>% filter(studygroup_annotation=="Anti-TNF IBD" & seq_percent >= 3)
gene_heavy_count_HC <- gene_heavy_count %>% filter(studygroup_annotation=="Healthy Controls" & seq_percent >= 3)

gene_heavy_count_IL12_23_gene <- gene_heavy_count_IL12_23$gene
gene_heavy_count_HC_gene <- gene_heavy_count_HC$gene
gene_heavy_count_TNF_gene <- gene_heavy_count_TNF$gene
CoV_database_human_count_gene <- CoV_database_human_count$gene

gene_heavy_count_filtered <- gene_heavy_count %>% filter((gene %in% gene_heavy_count_HC_gene) | (gene %in% gene_heavy_count_TNF_gene) 
                                                         | (gene %in% gene_heavy_count_TNF_gene) | (gene %in% CoV_database_human_count_gene))

#Merge database genes and sample genes 
gene_heavy_count_filtered <- merge(gene_heavy_count_filtered, CoV_database_human_count, all=TRUE)

gene_heavy_count_filtered$label_freq <- round(gene_heavy_count_filtered$seq_percent, digits = 1)

#Plot gene usage by studygroup
ggplot(gene_heavy_count_filtered, aes(x=gene, y=seq_percent, fill=studygroup_annotation, label = label_freq)) +
  theme_bw() + theme(text=element_text(family = "Arial"),
                     axis.text.x=element_text(angle = 25, size=6, vjust=1, hjust=1, face="italic"),
                     axis.text.y=element_text(size=6),
                     strip.text.y = element_blank()) + 
  geom_text(position = position_dodge(width = .9),    # move to center of bars
            vjust = -0.5,
            size = 2.3) +
  ylab("") + xlab("") + scale_fill_manual(values=colors) + ylim(c(0, 25)) +
  geom_bar(stat = "identity", position = "dodge") + facet_grid(factor(studygroup_annotation, 
                                                                      levels = c("Anti-IL12-23 IBD", "Anti-TNF IBD", "Healthy Controls", "CoV-AbDab Database"))~.) + NoLegend()


# Light chain kappa gene usage 

#Light chain kappa V gene usage in CoV-AbDab Database- filter on human and SARS-CoV-2 specific antibodies and genes with a frequency of at least 3%
#Download CoV-AbDAb database as a csv file
CoV_database <- read.csv("path/to/CoV-AbDab/csv/file", header=TRUE)
CoV_database_kappa <- CoV_database %>% filter(str_detect(Light.V.Gene, 'IGKV'))
CoV_database_kappa_human <- CoV_database_kappa %>% filter(str_detect(Light.V.Gene, '(Human)'))
CoV_database_kappa_human <- CoV_database_kappa_human %>% filter(str_detect(Binds.to, 'SARS-CoV2'))
CoV_database_kappa_human_count <- countGenes(CoV_database_kappa_human, gene="Light.V.Gene", mode="gene")
CoV_database_kappa_human_count$seq_percent <- CoV_database_kappa_human_count$seq_freq * 100
CoV_database_kappa_human_count <- CoV_database_kappa_human_count %>% filter(seq_percent >= 3)
CoV_database_kappa_human_count$studygroup_annotation <- 'CoV-AbDab Database'


#Filter for sample genes with a frequency of at least 3%
gene_kappa_count_IL12_23 <- gene_kappa_count %>% filter(studygroup_annotation=="Anti-IL12-23 IBD" & seq_percent >= 3)
gene_kappa_count_TNF <- gene_kappa_count %>% filter(studygroup_annotation=="Anti-TNF IBD" & seq_percent >= 3)
gene_kappa_count_HC <- gene_kappa_count %>% filter(studygroup_annotation=="Healthy Controls" & seq_percent >= 3)

gene_kappa_count_IL12_23_gene <- gene_kappa_count_IL12_23$gene
gene_kappa_count_HC_gene <- gene_kappa_count_HC$gene
gene_kappa_count_TNF_gene <- gene_kappa_count_TNF$gene
CoV_database_kappa_human_count_gene <- CoV_database_kappa_human_count$gene

gene_kappa_count_filtered <- gene_kappa_count %>% filter((gene %in% gene_kappa_count_HC_gene) | (gene %in% gene_kappa_count_TNF_gene) 
                                                         | (gene %in% gene_kappa_count_TNF_gene) | (gene %in% CoV_database_kappa_human_count_gene))

#Merge database genes and sample genes 
gene_kappa_count_filtered <- merge(gene_kappa_count_filtered, CoV_database_kappa_human_count, all=TRUE)

gene_kappa_count_filtered$label_freq <- round(gene_kappa_count_filtered$seq_percent, digits = 1)

#Plot gene usage by studygroup
ggplot(gene_kappa_count_filtered, aes(x=gene, y=seq_percent, fill=studygroup_annotation, label = label_freq)) +
  theme_bw() + theme(text=element_text(family = "Arial"),
                     axis.text.x=element_text(angle = 25, size=6, vjust = 1, hjust= 1, face="italic"),
                     axis.text.y=element_text(size=6),
                     strip.text.y = element_blank()) + 
  geom_text(position = position_dodge(width = .9),    # move to center of bars
            vjust = -0.5,
            size = 2.3) +
  ylab("") + xlab("") + scale_fill_manual(values=colors) + ylim(c(0, 25)) +
  geom_bar(stat = "identity", position = "dodge") + facet_grid(factor(studygroup_annotation, 
                                                                      levels = c("Anti-IL12-23 IBD", "Anti-TNF IBD", "Healthy Controls", "CoV-AbDab Database"))~.) + NoLegend()


# Light chain lambda gene usage
#Light chain kappa V gene usage in CoV-AbDab Database- filter on human and SARS-CoV-2 specific antibodies and genes with a frequency of at least 3%
#Download CoV-AbDAb database as a csv file

CoV_database <- read.csv("path/to/CoV-AbDab/csv/file", header=TRUE)
CoV_database_lambda <- CoV_database %>% filter(str_detect(Light.V.Gene, 'IGL'))
CoV_database_lambda_human <- CoV_database_lambda %>% filter(str_detect(Light.V.Gene, '(Human)'))
CoV_database_lambda_human <- CoV_database_lambda_human %>% filter(str_detect(Binds.to, 'SARS-CoV2'))
CoV_database_lambda_human_count <- countGenes(CoV_database_lambda_human, gene="Light.V.Gene", mode="gene")
CoV_database_lambda_human_count$seq_percent <- CoV_database_lambda_human_count$seq_freq * 100
CoV_database_lambda_human_count <- CoV_database_lambda_human_count %>% filter(seq_percent >= 3)
CoV_database_lambda_human_count$studygroup_annotation <- 'CoV-AbDab Database'

#Filter for sample genes with a frequency of at least 3%
gene_lambda_count_IL12_23 <- gene_lambda_count %>% filter(studygroup_annotation=="Anti-IL12-23 IBD" & seq_percent >= 3)
gene_lambda_count_TNF <- gene_lambda_count %>% filter(studygroup_annotation=="Anti-TNF IBD" & seq_percent >= 3)
gene_lambda_count_HC <- gene_lambda_count %>% filter(studygroup_annotation=="Healthy Controls" & seq_percent >= 3)

gene_lambda_count_IL12_23_gene <- gene_lambda_count_IL12_23$gene
gene_lambda_count_HC_gene <- gene_lambda_count_HC$gene
gene_lambda_count_TNF_gene <- gene_lambda_count_TNF$gene
CoV_database_lambda_human_count_gene <- CoV_database_lambda_human_count$gene

gene_lambda_count_filtered <- gene_lambda_count %>% filter((gene %in% gene_lambda_count_HC_gene) | (gene %in% gene_lambda_count_TNF_gene) 
                                                           | (gene %in% gene_lambda_count_TNF_gene) | (gene %in% CoV_database_lambda_human_count_gene))

#Merge database genes and sample genes 
gene_lambda_count_filtered <- merge(gene_lambda_count_filtered, CoV_database_lambda_human_count, all=TRUE)

gene_lambda_count_filtered$label_freq <- round(gene_lambda_count_filtered$seq_percent, digits = 1)

#Plot gene usage by studygroup
ggplot(gene_lambda_count_filtered, aes(x=gene, y=seq_percent, fill=studygroup_annotation, label = label_freq)) +
  theme_bw() + theme(text=element_text(family = "Arial"),
                     axis.text.x=element_text(angle = 25, size=6,  vjust = 1, hjust= 1, face="italic"),
                     axis.text.y=element_text(size=6),
                     strip.text.y = element_blank()) + 
  geom_text(position = position_dodge(width = .9),    # move to center of bars
            vjust = -0.5,
            size = 2.3) +
  ylab("") + xlab("") + scale_fill_manual(values=colors) + ylim(c(0, 25)) +
  geom_bar(stat = "identity", position = "dodge") + facet_grid(factor(studygroup_annotation, 
                                                                      levels = c("Anti-IL12-23 IBD", "Anti-TNF IBD", "Healthy Controls", "CoV-AbDab Database"))~.) + NoLegend()

ggsave(file="C:/Users/miche/OneDrive - University of Toronto/IMPACT/MBC scRNA-seq/Affinity Maturation B Cell Paper/Figures/Figure 5 BCR Gene Usage/Figure5_VL_Lambda_Gene_Usage.svg",
       plot=image, width=3.5, height=3.2)

#Paired VH and VL Gene Usage 

#Heatmap of paired VH and VL
#Re-shape from long to wide
#Keep only cells with one heavy and one light chain (multi_light==0)
clonal_results_no_multi_light <- clonal_results[clonal_results$multi_light==0, ]
clonal_results_no_multi_light <- transform(clonal_results_no_multi_light, chain= ifelse(locus=="IGH", "heavy", "light"))
data_wide <- dcast(clonal_results_no_multi_light, cell_id + studygroup_annotation + hashtag_annotation + clone_id + Timepoint ~ chain, value.var="v_call_10x")
data_wide_matrix <- data.frame(table(data_wide[, c("studygroup_annotation", "heavy", "light")]))
data_wide_hashtag_matrix <- data.frame(table(data_wide[, c("hashtag_annotation", "heavy", "light")]))

data_wide_matrix_wide <- dcast(data_wide_matrix, heavy + light ~ studygroup_annotation, value.var="Freq")
shared_color_scale <- c("white", "black", "orange", "#28a4ac", "#2c37e8", "#7c501d", "darkgreen","purple")
data_wide_matrix_wide_colored <- data_wide_matrix_wide %>%
  mutate(shared=case_when(
    `Anti-IL12-23 IBD`==0 & `Anti-TNF IBD`==0 & `Healthy Controls`==0 ~ 0,
    `Anti-IL12-23 IBD`==0 & `Anti-TNF IBD`==0 & `Healthy Controls`>=1 ~ 1,
    `Anti-IL12-23 IBD`==0 & `Anti-TNF IBD`>=1 & `Healthy Controls`==0 ~ 2,
    `Anti-IL12-23 IBD`>=1 & `Anti-TNF IBD`==0 & `Healthy Controls`==0 ~ 3,
    `Anti-IL12-23 IBD`>=1 & `Anti-TNF IBD`==0 & `Healthy Controls`>=1 ~ 4,
    `Anti-IL12-23 IBD`==0 & `Anti-TNF IBD`>=1 & `Healthy Controls`>=1 ~ 5,
    `Anti-IL12-23 IBD`>=1 & `Anti-TNF IBD`>=1 & `Healthy Controls`==0 ~ 6,
    `Anti-IL12-23 IBD`>=1 & `Anti-TNF IBD`>=1 & `Healthy Controls`>=1 ~ 7,
  ))

data_wide_matrix_wide_ALLshared <- data_wide_matrix_wide_colored[data_wide_matrix_wide_colored$shared==7, ]

#Tile plot showing SHARED VH and VL (yes/no dichotomous based on "shared" status)
image_shared = ggplot(data_wide_matrix_wide_colored, aes(heavy, light)) + geom_tile(aes(fill=factor(shared)), color = "gray") + theme_bw() +
  theme(text = element_text(family = "sans", size = 6), 
        axis.text.x = element_text(family="sans", size=6, angle=90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=6), 
        axis.title.x = element_blank(), axis.title.y=element_blank(),
        strip.text.y = element_text(size=6)) + labs(fill='Frequency (n)') + scale_fill_manual(values=shared_color_scale) + NoLegend()
image_shared

#Paired VH and VL tileplots, separated by studygroup, showing frequncy of the paired VH and VL genes found in ALL THREE GROUPS 
image_purple_1223 = ggplot(data_wide_matrix_wide_ALLshared, aes(heavy, light)) + geom_tile(aes(fill=`Anti-IL12-23 IBD`)) + theme_bw() +
  theme(text = element_text(family = "sans", size = 6), 
        panel.grid.major = element_line(size = 0.15),
        axis.text.x = element_text(family="sans", size=6, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=6), 
        axis.title.x = element_blank(), axis.title.y=element_blank(), 
        strip.text.x = element_text(size=8)) + labs(fill='Frequency (n)') + scale_fill_gradient(low="#e8e8e8", high="darkblue", limits=c(0, 15)) + NoLegend()
image_purple_1223

image_purple_TNF = ggplot(data_wide_matrix_wide_ALLshared, aes(heavy, light)) + geom_tile(aes(fill=`Anti-TNF IBD`)) + theme_bw() +
  theme(text = element_text(family = "sans", size = 6), 
        panel.grid.major = element_line(size = 0.15),
        axis.text.x = element_text(family="sans", size=6, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=6), 
        axis.title.x = element_blank(), axis.title.y=element_blank(), 
        strip.text.x = element_text(size=8)) + labs(fill='Frequency (n)') + scale_fill_gradient(low="#e8e8e8", high="darkblue", limits=c(0, 15)) + NoLegend()
image_purple_TN


image_purple_HC = ggplot(data_wide_matrix_wide_ALLshared, aes(heavy, light)) + geom_tile(aes(fill=`Healthy Controls`)) + theme_bw() +
  theme(text = element_text(family = "sans", size = 6), 
        panel.grid.major = element_line(size = 0.15),
        axis.text.x = element_text(family="sans", size=6, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=6), 
        axis.title.x = element_blank(), axis.title.y=element_blank(), 
        strip.text.x = element_text(size=8)) + labs(fill='Frequency (n)') + scale_fill_gradient(low="#e8e8e8", high="darkblue", limits=c(0, 15)) + NoLegend()
image_purple_HC
