#--------------------------------------------------------------#
#-Island biogeography and the PCA axes-------------------------#
#--------------------------------------------------------------#

# Test whether island biogeography explains the topology of the 11 island
# networks. Each of the two PCA axes is regressed on isolation (distance to
# the continent), area (log-transformed), elevation and island type
# (continental or oceanic). The continuous predictors are standardized to
# make their coefficients comparable. With 11 islands, each model has six
# residual degrees of freedom.
#
# Inputs:  data/main_summary/island_data.csv
# Outputs: none. The summaries of the two models are printed to the console.

library(here)
library(tidyverse)

#--------------------------------------------------------------#
#-Island biogeography and PCA scores of the 11 islands---------#
#--------------------------------------------------------------#

island_data <- read_csv(here("data", "main_summary", "island_data.csv"), show_col_types=FALSE)

#--------------------------------------------------------------#
#-Each PCA axis against island biogeography--------------------#
#--------------------------------------------------------------#

# Global network integration (PC1) against island biogeography
lm1 <- lm(global_integration ~ scale(isolation) + scale(log(area)) + scale(elevation) + type,
          data=island_data)
summary(lm1)

# Pollination topological dominance (PC2) against island biogeography
lm2 <- lm(pollination_dominance ~ scale(isolation) + scale(log(area)) + scale(elevation) + type,
          data=island_data)
summary(lm2)
