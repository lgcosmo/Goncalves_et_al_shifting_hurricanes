#--------------------------------------------------------------#
#-PCA of the topology descriptors------------------------------#
#--------------------------------------------------------------#

# PCA of the five topology descriptors of the 11 island networks: leading
# eigenvalue of the multilayer network, of the pollination layer and of the
# seed-dispersal layer, proportion of dual-role species and proportion of
# interactions involving dual-role species. The first two axes are kept and
# oriented. PC1 is global network integration and PC2 is pollination
# topological dominance.
#
# Inputs:  data/main_summary/macro_structure.csv
# Outputs: data/main_summary/pca_scores.csv
#          data/main_summary/pca_loadings.csv

library(here)
library(tidyverse)
library(FactoMineR)

#--------------------------------------------------------------#
#-Topology descriptors of the 11 islands-----------------------#
#--------------------------------------------------------------#

macro <- read_csv(here("data", "main_summary", "macro_structure.csv"), show_col_types=FALSE)

pca_vars <- c("multi_lambda1", "pol_lambda1", "sd_lambda1",
              "prop_dual_species", "prop_dual_interactions")

#--------------------------------------------------------------#
#-PCA of the five descriptors----------------------------------#
#--------------------------------------------------------------#

# Scale the descriptors to unit variance. They differ in units and range.
res <- PCA(macro %>% select(all_of(pca_vars)), scale.unit=TRUE, graph=FALSE)

loadings <- res$var$coord[, 1:2]
scores <- res$ind$coord[, 1:2]

# The sign of a PCA axis is arbitrary. Orient PC1 so the descriptors load
# positively (positive sum of loadings) and PC2 so the proportion of
# dual-role species loads negatively. Higher scores then mean stronger global
# network integration and stronger pollination topological dominance.
axis_sign <- c(sign(sum(loadings[, 1])), -sign(loadings["prop_dual_species", 2]))
loadings <- sweep(loadings, 2, axis_sign, "*")
scores <- sweep(scores, 2, axis_sign, "*")

#--------------------------------------------------------------#
#-Loadings and island scores, written to data/main_summary-----#
#--------------------------------------------------------------#

# Loadings of the five descriptors and percentage of the variance explained
# by each axis
variance <- res$eig[1:2, "percentage of variance"]

pca_loadings <- tibble(descriptor=rownames(loadings),
                       PC1=loadings[, 1], PC2=loadings[, 2]) %>%
  add_row(descriptor="Variance explained (%)", PC1=variance[1], PC2=variance[2])
print(pca_loadings)

# Score of each island on the two oriented axes
pca_scores <- tibble(network=macro$network,
                     global_integration=scores[, 1],
                     pollination_dominance=scores[, 2])

write_csv(pca_loadings, here("data", "main_summary", "pca_loadings.csv"))
write_csv(pca_scores, here("data", "main_summary", "pca_scores.csv"))
message("Saved: pca_loadings.csv and pca_scores.csv")
