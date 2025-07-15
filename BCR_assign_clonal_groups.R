# Assigning clonal groups and reconstructing germline sequences #
# Cheung et al. JCI 2025

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

#### Clonal Analyses - Immcantation Step 1 Clustering Threshold #### 

#Step 1: Determine clonal clustering threshold; sequences under this cut-off are clonally related 
#Determine an appropriate threshold for trimming hierarchial clustering into B cell clones 
#Ideal threshold for separating clonal groups is the value that separates the two modes of the nearest-neighbor distance distribution
#Chosen threshold: 0.3 

# Check by Timepoint
dist_sc_4 <- distToNearest(subset(bcr_data_merged, subset = Timepoint==4), model=c("ham"), cellIdColumn="cell_id", locusColumn="locus",  
                           VJthenLen=FALSE, onlyHeavy=FALSE, normalize="len")

dist_sc_6 <- distToNearest(subset(bcr_data_merged, subset = Timepoint==6), model=c("ham"), cellIdColumn="cell_id", locusColumn="locus",  
                           VJthenLen=FALSE, onlyHeavy=FALSE, normalize="len")

# Find threshold using density or gamma method
output <- findThreshold(dist_sc_4$dist_nearest, method="gmm", model="gamma-gamma")
output <- findThreshold(dist_sc_6$dist_nearest, method="gmm", model="gamma-gamma")

output_4 <- findThreshold(dist_sc_4$dist_nearest, method="density")
output_6 <- findThreshold(dist_sc_6$dist_nearest, method="density")
plot(output_4, title="Density Method")
plot(output_6, title="Density Method")

# Generate Hamming distance histograms
p1 <- ggplot(subset(dist_sc_4, !is.na(dist_nearest)), 
             aes(x=dist_nearest)) + 
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.3, color="firebrick", linetype=2) +
  labs(x = "Hamming distance", y = "Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  theme_bw()
plot(p1)

p2 <- ggplot(subset(dist_sc_6, !is.na(dist_nearest)), 
             aes(x=dist_nearest)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.3, color="firebrick", linetype=2) +
  labs(x = "Grouped Hamming distance", y = "Count") + 
  facet_grid(hashtag_annotation ~ ., scales="free_y") + 
  theme_bw()
plot(p2)

#### Clonal Analyses - Immcantation Step 2 Assign Clonal Groups #### 

#Step 2: Assign clonal groups; add annotation (clone_id) that can be used to identify a group of seqs that came from the same naive cell 
#Provides a hierarchial agglomerative clustering approach to infer clonal relationships in high-throughput AIRR-seq data
#Clusters B/TCR sequences based on junction region sequence similarity within partitions that share the same V, J gene and junction length
#Define clonal groups 
#Call clones using hierarchicalClones

clonal_results <- hierarchicalClones(bcr_data_merged,
                                     cell_id = "cell_id", 
                                     threshold = 0.3,
                                     locus = "locus",
                                     only_heavy = TRUE, split_light = TRUE,
                                     summarize_clones = FALSE, cdr3 = TRUE,
                                     method="nt")

#### Clonal Analyses - Immcantation Step 3 Reconstruct Germline Sequences #### 

#Step 3: reconstruct germline sequences; figure out germline seq of common ancestor, before mutations are introduced during clonal expansion and SHM 
#Goal to reconstruct the sequence of unmutated ancestor of each clone using a reference database of known alleles (IMGT)

#Read in IMGT data downloaded off of docker and update `dir` to use the path to human\vdj folder
references <- readIMGT(dir = "path/to/IMGT/reference")
clonal_results <- clonal_results[!is.na(clonal_results$clone_id), ]
clonal_results <- createGermlines(clonal_results, references, fields = "Timepoint", nproc = 1)




