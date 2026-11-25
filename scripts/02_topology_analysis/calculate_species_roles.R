#--------------------------------------------------------------#
#-Role metrics of species in the multilayer networks-----------#
#--------------------------------------------------------------#

# Compute three role metrics for every plant and bird of each island. The
# within-layer degree z-score compares the degree of a species in one layer
# with the degrees of its guild (plants or birds) on the same island, for the
# pollination layer (z_pollination) and the seed-dispersal layer
# (z_dispersal). The multilayer participation coefficient (MPC) measures how
# evenly a species splits its interactions between the two layers, from 0
# (one layer only) to 0.5 (an equal split).
#
# Inputs:  src/multi_pol/B_<island>.csv
#          src/multi_sd/S_<island>.csv
# Outputs: data/main_summary/species_roles.csv

library(here)
library(tidyverse)

#--------------------------------------------------------------#
#-Role metrics of the species of one island--------------------#
#--------------------------------------------------------------#

# Role metrics of every species of one island. Degrees are counted on the
# binarized layers. The z-scores are computed within guilds: plants with
# plants and birds with birds. Returns a tibble with one row per species.
calculate_network_roles <- function(network) {

  df_pol <- read_csv(here("src", "multi_pol", paste0("B_", network, ".csv")), show_col_types=FALSE)
  df_sd  <- read_csv(here("src", "multi_sd",  paste0("S_", network, ".csv")), show_col_types=FALSE)

  mat_pol <- (as.matrix(df_pol[, -1]) > 0) * 1
  mat_sd  <- (as.matrix(df_sd[, -1]) > 0) * 1

  # Degree of each plant (rows) and bird (columns) in each layer
  species <- bind_rows(
    tibble(species_id=df_pol[[1]], species_set="Plant",
           degree_pol=rowSums(mat_pol), degree_sd=rowSums(mat_sd)),
    tibble(species_id=colnames(df_pol)[-1], species_set="Animal",
           degree_pol=colSums(mat_pol), degree_sd=colSums(mat_sd))
  )

  # z-scores within guilds, and MPC from the share of each layer in the total
  # degree of the species
  species %>%
    group_by(species_set) %>%
    mutate(
      z_pollination = (degree_pol - mean(degree_pol)) / sd(degree_pol),
      z_dispersal   = (degree_sd - mean(degree_sd)) / sd(degree_sd)
    ) %>%
    ungroup() %>%
    mutate(
      k_total = degree_pol + degree_sd,
      mpc = 1 - ((degree_pol / k_total)^2 + (degree_sd / k_total)^2),
      network = network
    )
}

#--------------------------------------------------------------#
#-Role metrics of the species of all islands-------------------#
#--------------------------------------------------------------#

networks <- list.files(here("src", "multi_pol"), pattern="^B_.*\\.csv$") %>%
  str_remove("^B_") %>%
  str_remove("\\.csv$")

all_roles <- map_dfr(networks, calculate_network_roles)

dir.create(here("data", "main_summary"), showWarnings=FALSE, recursive=TRUE)
write_csv(all_roles, here("data", "main_summary", "species_roles.csv"))
message("Saved: species_roles.csv")
