#--------------------------------------------------------------#
#-Functional roles of species in the multilayer networks-------#
#--------------------------------------------------------------#

# Classify every plant and bird of each island's multilayer network into a
# functional role from how its interactions are distributed within and across
# the pollination and seed-dispersal layers. Three role metrics describe each
# species. The within-layer degree z-score compares its degree in a layer
# with the degrees of its guild (plants or birds) on the same island,
# separately for the pollination layer (z_pollination) and the seed-dispersal
# layer (z_dispersal). The multilayer participation coefficient (MPC)
# measures how evenly it splits its interactions between the two layers,
# from 0 (one layer only) to 0.5 (an equal split). A z-score above 0.75
# marks high within-layer connectivity (Olesen et al. 2007) and an MPC above
# 0.3 marks participation in both layers. The thresholds sort the species
# into Multilayer Hubs, Pollination Hubs, Dispersal Hubs, Multilayer
# Connectors and Peripheral species, the role counts of the Results. The
# pSEM of Fig. 4 uses the two z-scores and the MPC of the bird species as
# predictors of species fragility.
#
# Inputs:  src/multi_pol/B_<island>.csv and src/multi_sd/S_<island>.csv,
#          written by build_networks_multilayer.R
# Outputs: data/main_summary/species_roles.csv, one row per species and
#          island, read by calculate_auc_extinction.R on its way to the pSEM
#          (structural_equation_model.R)

library(here)
library(tidyverse)

# A z-score above 0.75 marks high within-layer connectivity and an MPC above
# 0.3 marks participation in both layers
z_threshold <- 0.75
p_threshold <- 0.3

#--------------------------------------------------------------#
#-Role metrics and roles of the species of one island----------#
#--------------------------------------------------------------#

# Role metrics and role of every species of one island. Degrees are counted
# on the binarized layers. The z-scores are computed within guilds, so plants
# are compared with plants and birds with birds. Multilayer Hubs are highly
# connected in both layers, or in one layer while participating in both.
# Pollination Hubs and Dispersal Hubs are highly connected in one layer and
# concentrate their interactions there. Multilayer Connectors are highly
# connected in neither layer but participate in both, and Peripheral species
# in neither layer with interactions mostly in one. The rules are tested in
# that order, so each species takes the first role that fits.
# Arguments: network is the island code shared by the two matrix files in
# src/multi_pol and src/multi_sd. Returns a tibble with one row per species,
# holding its guild, its degree in each layer and in total, the two z-scores,
# the MPC, its role and the island.
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

  # Within-layer degree z-scores by guild, the MPC from the share of each
  # layer in the species' total degree, and the role from the thresholds
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
      role = case_when(
        (z_pollination > z_threshold & z_dispersal > z_threshold) |
          ((z_pollination > z_threshold | z_dispersal > z_threshold) & mpc > p_threshold) ~ "Multilayer Hub",
        z_pollination > z_threshold ~ "Pollination Hub",
        z_dispersal > z_threshold   ~ "Dispersal Hub",
        mpc > p_threshold           ~ "Multilayer Connector",
        TRUE                        ~ "Peripheral"
      ),
      network = network
    )
}

#--------------------------------------------------------------#
#-Every island, written to species_roles.csv-------------------#
#--------------------------------------------------------------#

networks <- list.files(here("src", "multi_pol"), pattern="^B_.*\\.csv$") %>%
  str_remove("^B_") %>%
  str_remove("\\.csv$")

all_roles <- map_dfr(networks, calculate_network_roles)

dir.create(here("data", "main_summary"), showWarnings=FALSE, recursive=TRUE)
write_csv(all_roles, here("data", "main_summary", "species_roles.csv"))
message("Saved: species_roles.csv")
