#--------------------------------------------------------------#
#-Island biogeography, hurricane regime and network topology---#
#--------------------------------------------------------------#

# Join the descriptors of the 11 islands into one table with one row per
# island: island biogeography (type, area, maximum elevation and isolation),
# empirical hurricane frequency (lambda), network topology (richness,
# connectance, leading eigenvalue, NODF and Q of each layer and of the
# multilayer network, proportions of dual-role species and of their
# interactions) and the two PCA axes (global network integration and
# pollination topological dominance).
#
# Inputs:  src/Island_index.csv
#          src/empirical_lambdas.csv
#          data/main_summary/macro_structure.csv
#          data/main_summary/pca_scores.csv
# Outputs: data/main_summary/island_data.csv

library(here)
library(tidyverse)

source(here("src", "island_labels.R"))

#--------------------------------------------------------------#
#-Loading the island descriptors-------------------------------#
#--------------------------------------------------------------#

island_index <- read_csv(here("src", "Island_index.csv"), show_col_types=FALSE)

# Empirical lambda of each island, with the number of hurricanes and the
# years of record
emp_lambdas <- read_csv(here("src", "empirical_lambdas.csv"), show_col_types=FALSE) %>%
  select(network, lambda, n_events, observation_period_years)

macro <- read_csv(here("data", "main_summary", "macro_structure.csv"), show_col_types=FALSE)

pca_scores <- read_csv(here("data", "main_summary", "pca_scores.csv"), show_col_types=FALSE)

#--------------------------------------------------------------#
#-One row per island-------------------------------------------#
#--------------------------------------------------------------#

# Start from the table of lambdas, add the island names and join the other
# tables by the island code
island_data <- emp_lambdas %>%
  mutate(island=island_names[network], .after=network) %>%
  left_join(island_index %>% select(-island), by="network") %>%
  left_join(macro, by="network") %>%
  left_join(pca_scores, by="network")

write_csv(island_data, here("data", "main_summary", "island_data.csv"))
message("Saved: island_data.csv")
