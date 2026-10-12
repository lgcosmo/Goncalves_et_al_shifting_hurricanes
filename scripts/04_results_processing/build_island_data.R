#--------------------------------------------------------------#
#-Island biogeography, hurricane regime and network topology---#
#--------------------------------------------------------------#

# Join the island-level descriptors of the 11 islands into one table with one
# row per island. The table supplies the island biogeography of the pSEM
# (Fig. 4), the regressions of the PCA axes on island biogeography
# (linear_model_pca.R) and the island descriptors drawn by figure_2.R and
# figure_3.R. Island biogeography covers island type (continental or
# oceanic), area, maximum elevation and isolation (distance to the
# continent). The hurricane regime holds the empirical hurricane frequency
# (lambda, estimated from the historical hurricane record). Network topology
# covers species richness, connectance, the leading eigenvalue, nestedness
# (NODF) and modularity (Barber's Q) of each layer and of the multilayer
# network, the proportions of dual-role species and of their interactions,
# and the two oriented PCA axes of calculate_pca.R, global network
# integration and pollination topological dominance.
#
# Inputs:  src/Island_index.csv (island biogeography)
#          src/empirical_lambdas.csv (empirical_lambda_pmf.R)
#          data/main_summary/macro_structure.csv (calculate_macro_structure.R)
#          data/main_summary/pca_scores.csv (calculate_pca.R)
# Outputs: data/main_summary/island_data.csv
# Used by: linear_model_pca.R, structural_equation_model.R, figure_2.R and
#          figure_3.R

library(here)
library(tidyverse)

source(here("src", "island_labels.R"))

#--------------------------------------------------------------#
#-Island descriptors, keyed by the island code-----------------#
#--------------------------------------------------------------#

island_index <- read_csv(here("src", "Island_index.csv"), show_col_types=FALSE)

# Empirical hurricane frequency of each island, with the number of hurricanes
# and the years of record behind it
emp_lambdas <- read_csv(here("src", "empirical_lambdas.csv"), show_col_types=FALSE) %>%
  select(network, lambda, n_events, observation_period_years)

macro <- read_csv(here("data", "main_summary", "macro_structure.csv"), show_col_types=FALSE)

pca_scores <- read_csv(here("data", "main_summary", "pca_scores.csv"), show_col_types=FALSE)

#--------------------------------------------------------------#
#-One row per island, written to island_data.csv---------------#
#--------------------------------------------------------------#

# The frequency table defines the 11 islands, named after the manuscript;
# every other table joins it by the island code
island_data <- emp_lambdas %>%
  mutate(island=island_names[network], .after=network) %>%
  left_join(island_index %>% select(-island), by="network") %>%
  left_join(macro, by="network") %>%
  left_join(pca_scores, by="network")

write_csv(island_data, here("data", "main_summary", "island_data.csv"))
message("Saved: island_data.csv")
