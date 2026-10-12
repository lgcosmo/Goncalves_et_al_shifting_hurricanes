#--------------------------------------------------------------#
#-PCA of the topology descriptors------------------------------#
#--------------------------------------------------------------#

# Summarize the topology of the 11 island networks with a principal component
# analysis (PCA) of the five descriptors of calculate_macro_structure.R. Three
# measure the integration of the network, namely the leading eigenvalue of
# the multilayer network, of the pollination layer and of the seed-dispersal
# layer. The other two, the proportion of dual-role species and the
# proportion of interactions involving dual-role species, measure how
# strongly dual-role species couple the layers. The first two axes are kept
# and oriented here, once for the whole pipeline. All five descriptors load
# positively on the first axis, global network integration. On the second
# axis the dual-role descriptors load negatively, so that higher scores mean
# stronger pollination topological dominance. The scores are used in the Fig. 2
# inset, the regressions of linear_model_pca.R and the pSEM of Fig. 4
# through island_data.csv.
#
# Inputs:  data/main_summary/macro_structure.csv (calculate_macro_structure.R)
# Outputs: data/main_summary/pca_scores.csv, the score of each island on
#          global_integration and pollination_dominance (build_island_data.R)
#          data/main_summary/pca_loadings.csv, the loadings of the five
#          descriptors and the variance explained (Extended Data Table 2)

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
#-PCA of the five descriptors (Extended Data Table 2)----------#
#--------------------------------------------------------------#

# The descriptors are scaled to unit variance, as they differ in units and
# range
res <- PCA(macro %>% select(all_of(pca_vars)), scale.unit=TRUE, graph=FALSE)

loadings <- res$var$coord[, 1:2]
scores <- res$ind$coord[, 1:2]

# The sign of a PCA axis is arbitrary. The first axis is oriented so that the
# descriptors load positively, and the second so that the proportion of
# dual-role species loads negatively. Higher scores then mean stronger global
# network integration and stronger pollination topological dominance.
axis_sign <- c(sign(sum(loadings[, 1])), -sign(loadings["prop_dual_species", 2]))
loadings <- sweep(loadings, 2, axis_sign, "*")
scores <- sweep(scores, 2, axis_sign, "*")

#--------------------------------------------------------------#
#-Loadings and island scores, written to data/main_summary-----#
#--------------------------------------------------------------#

# Loadings of the five descriptors and the share of the variation explained
# by each axis, the rows of Extended Data Table 2
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
