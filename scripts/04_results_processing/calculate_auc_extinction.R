#--------------------------------------------------------------#
#-Species fragility as the AUC of the extinction curves--------#
#--------------------------------------------------------------#

# Compute the fragility of every species in every island network to hurricane
# frequency and to hurricane intensity, the two responses of the pSEM (Fig. 4
# and the SEM tables). Species fragility is the normalized area under the
# extinction curve (AUC), the curve that relates the fraction of simulations
# in which the species went extinct to hurricane frequency (lambda) or to
# hurricane intensity (the shift s of probability mass toward categories 4
# and 5). The area is computed separately for the two scenarios and divided
# by the range of the x axis, so that higher AUC means greater fragility. A
# 95% confidence interval follows from the binomial variance of each point of
# the curve. Each species then receives its guild (plant or bird), its
# within-layer degree z-scores, its multilayer participation coefficient
# (MPC) and its functional role, the species-level predictors of the pSEM.
#
# Inputs:  data/main_summary/main_summary_metrics_pmf.csv    (survival_curves.R)
#          data/main_summary/main_summary_metrics_lambda.csv (survival_curves.R)
#          data/main_summary/species_roles.csv (calculate_species_roles.R)
# Outputs: data/main_summary/main_summary_auc.csv
# Helpers: compute_auc and calculate_auc_se in src/aux_functions.R
# Used by: structural_equation_model.R

library(here)
library(tidyverse)

source(here("src", "aux_functions.R"))

#--------------------------------------------------------------#
#-Species fragility in the frequency and intensity scenarios---#
#--------------------------------------------------------------#

df_pmf    <- read_csv(here("data", "main_summary", "main_summary_metrics_pmf.csv"),    show_col_types=FALSE)
df_lambda <- read_csv(here("data", "main_summary", "main_summary_metrics_lambda.csv"), show_col_types=FALSE)

# The x axis of each curve is lambda in the frequency sweep and the intensity
# shift in the intensity sweep. The scenario labels PMF and Lambda follow the
# names of the simulation folders
results_pmf    <- compute_auc(df_pmf,    x_col="intensity", scenario_label="PMF")
results_lambda <- compute_auc(df_lambda, x_col="lambda",    scenario_label="Lambda")

#--------------------------------------------------------------#
#-Species-level predictors of the pSEM in the AUC table--------#
#--------------------------------------------------------------#

# Guild, role metrics and role of every species, from calculate_species_roles.R
species_roles <- read_csv(here("data", "main_summary", "species_roles.csv"), show_col_types=FALSE) %>%
  select(network, species_id, species_set, mpc, z_pollination, z_dispersal, role)

# 95% confidence interval of each AUC from its analytical standard error. The
# network and type rows of the curves keep an empty species id and gain no
# species-level predictors.
results_all <- bind_rows(results_pmf, results_lambda) %>%
  mutate(
    auc_lower=auc - 1.96 * auc_se,
    auc_upper=auc + 1.96 * auc_se
  ) %>%
  left_join(species_roles, by=c("network", "species_id"))

write_csv(results_all, here("data", "main_summary", "main_summary_auc.csv"))
message("Saved: main_summary_auc.csv")
