#--------------------------------------------------------------#
#-Island biogeography and the PCA axes-------------------------#
#--------------------------------------------------------------#

# Test whether island biogeography explains the topology of the 11 island
# networks, summarized by the two PCA axes of calculate_pca.R. Each axis is
# regressed on the four biogeographic variables of the islands. They are
# isolation (distance to the continent), area (log-transformed), elevation
# and island type (continental or oceanic). The continuous predictors are
# standardized, so that their coefficients are comparable. With 11 islands,
# each model has six residual degrees of freedom. The R-squared and P value
# of each model are the ones reported in the Results, and the island-level
# submodels of the pSEM (structural_equation_model.R, Fig. 4) regress the
# same axes on the same variables.
#
# Inputs:  data/main_summary/island_data.csv, written by build_island_data.R
# Outputs: the summaries of the two models, printed to the console

library(here)
library(tidyverse)

#--------------------------------------------------------------#
#-Island biogeography and PCA scores of the 11 islands---------#
#--------------------------------------------------------------#

island_data <- read_csv(here("data", "main_summary", "island_data.csv"), show_col_types=FALSE)

#--------------------------------------------------------------#
#-Each PCA axis against island biogeography (R-squared and P)--#
#--------------------------------------------------------------#

# Global network integration against island biogeography, for the Results
lm1 <- lm(global_integration ~ scale(isolation) + scale(log(area)) + scale(elevation) + type,
          data=island_data)
summary(lm1)

# Pollination topological dominance against island biogeography, the
# relationship of the Fig. 2 inset
lm2 <- lm(pollination_dominance ~ scale(isolation) + scale(log(area)) + scale(elevation) + type,
          data=island_data)
summary(lm2)
