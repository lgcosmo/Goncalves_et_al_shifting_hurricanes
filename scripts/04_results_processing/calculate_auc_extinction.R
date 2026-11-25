#--------------------------------------------------------------#
#-Species persistence from the AUC of the extinction curves----#
#--------------------------------------------------------------#

# Compute the persistence of each species on each island under the two
# hurricane scenarios, the two responses of the pSEM. The extinction curve of
# a species is the fraction of simulations in which it went extinct along
# the values of hurricane frequency (lambda) or of the intensity shift. The
# area under this curve (AUC) is divided by the range of the x axis, so it
# lies between 0 and 1. Persistence is 1 - AUC, the mean probability that the
# species persists across the simulated values. The 95% confidence interval
# comes from the binomial variance of each point of the curve.
#
# Inputs:  data/main_summary/main_summary_metrics_pmf.csv
#          data/main_summary/main_summary_metrics_lambda.csv
#          data/main_summary/species_roles.csv
# Outputs: data/main_summary/main_summary_auc.csv

library(here)
library(tidyverse)

source(here("src", "aux_functions.R"))

#--------------------------------------------------------------#
#-AUC of the extinction curves in the two scenarios------------#
#--------------------------------------------------------------#

df_pmf    <- read_csv(here("data", "main_summary", "main_summary_metrics_pmf.csv"),    show_col_types=FALSE)
df_lambda <- read_csv(here("data", "main_summary", "main_summary_metrics_lambda.csv"), show_col_types=FALSE)

# x is lambda in the frequency scenario and the shift in the intensity
# scenario
results_pmf    <- compute_auc(df_pmf,    x_col="intensity", scenario_label="PMF")
results_lambda <- compute_auc(df_lambda, x_col="lambda",    scenario_label="Lambda")

#--------------------------------------------------------------#
#-Persistence and the species-level predictors of the pSEM-----#
#--------------------------------------------------------------#

# Guild and role metrics of each species
species_roles <- read_csv(here("data", "main_summary", "species_roles.csv"), show_col_types=FALSE) %>%
  select(network, species_id, species_set, mpc, z_pollination, z_dispersal)

# 95% confidence interval of the AUC from its standard error. Persistence has
# the same standard error, so its interval is one minus the interval of the
# AUC, with the limits exchanged.
results_all <- bind_rows(results_pmf, results_lambda) %>%
  mutate(
    auc_lower=auc - 1.96 * auc_se,
    auc_upper=auc + 1.96 * auc_se,
    persistence=1 - auc,
    persistence_lower=1 - auc_upper,
    persistence_upper=1 - auc_lower
  ) %>%
  left_join(species_roles, by=c("network", "species_id"))

write_csv(results_all, here("data", "main_summary", "main_summary_auc.csv"))
message("Saved: main_summary_auc.csv")
