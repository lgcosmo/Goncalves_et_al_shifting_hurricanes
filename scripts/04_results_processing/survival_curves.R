#--------------------------------------------------------------#
#-Extinction curves of the frequency and intensity sweeps------#
#--------------------------------------------------------------#

# Summarize the event-driven hurricane simulations of the empirical networks
# into the extinction curves of the two scenarios of Fig. 3. Under increased
# hurricane frequency, lambda rises from 0.01 to 0.80 hurricanes per year
# (data/main_sims_lambda/). Under increased hurricane intensity, a share s of
# the probability mass of categories 1 to 3 moves to categories 4 and 5
# (data/main_sims_pmf/). Each simulation file holds one island under one
# parameter value, with 1,000 replicates of 300 years. The extinctions are
# summarized at three levels. The species level holds the fraction of
# simulations in which each species went extinct. The type level averages
# that fraction over the species of each layer membership. The network level
# holds the fraction of species that went extinct, averaged over replicates.
# The frequency curves are then interpolated at the empirical frequency of
# each island, for the points that mark the islands in Fig. 3c.
#
# Inputs:  data/main_sims_pmf/*.csv    (main_simulations_empirical_pmf.jl)
#          data/main_sims_lambda/*.csv (main_simulations_empirical.jl)
#          src/empirical_lambdas.csv   (empirical_lambda_pmf.R)
# Outputs: data/main_summary/main_summary_metrics_pmf.csv
#          data/main_summary/main_summary_metrics_lambda.csv
#          data/main_summary/main_summary_empirical_interpolated.csv
# Helpers: process_sim_dir and interpolate_at_lambda in src/aux_functions.R
# Used by: calculate_auc_extinction.R and figure_3.R

library(here)
library(future.apply)

source(here("src", "aux_functions.R"))

# Workers read the simulation files in parallel, one file at a time
plan(multisession)

dir_pmf    <- here("data", "main_sims_pmf")
dir_lambda <- here("data", "main_sims_lambda")
dir_out    <- here("data", "main_summary")
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)

#--------------------------------------------------------------#
#-Extinction curves of both scenarios (Fig. 3c and 3d)---------#
#--------------------------------------------------------------#

# The labels PMF (intensity) and Lambda (frequency) name the two sweeps in
# the progress messages, as in the scenario column of main_summary_auc.csv
metrics_pmf <- process_sim_dir(dir_pmf, "PMF")
fwrite(metrics_pmf, file.path(dir_out, "main_summary_metrics_pmf.csv"))

metrics_lambda <- process_sim_dir(dir_lambda, "Lambda")
fwrite(metrics_lambda, file.path(dir_out, "main_summary_metrics_lambda.csv"))

#--------------------------------------------------------------#
#-Extinctions at the empirical frequency of each island (Fig. 3c)#
#--------------------------------------------------------------#

# The sweep raises lambda in steps of 0.01, whereas the empirical frequency
# of an island can fall between two simulated values
emp_lambdas <- fread(here("src", "empirical_lambdas.csv"))

interp <- interpolate_at_lambda(metrics_lambda, emp_lambdas)
fwrite(interp, file.path(dir_out, "main_summary_empirical_interpolated.csv"))

message("Saved: main_summary_metrics_pmf.csv, main_summary_metrics_lambda.csv ",
        "and main_summary_empirical_interpolated.csv")
