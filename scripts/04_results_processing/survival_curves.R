#--------------------------------------------------------------#
#-Extinction curves of the two hurricane scenarios-------------#
#--------------------------------------------------------------#

# Build the extinction curves of the empirical networks from the raw
# simulation files, at the species, type and network levels. In one scenario
# hurricane frequency (lambda) goes from 0.01 to 0.80 hurricanes per year. In
# the other, hurricane intensity increases as a probability s (the shift)
# moves from categories 1 to 3 to categories 4 and 5. Also interpolate the
# frequency curves at the empirical lambda of each island.
#
# Inputs:  data/main_sims_pmf/*.csv
#          data/main_sims_lambda/*.csv
#          src/empirical_lambdas.csv
# Outputs: data/main_summary/main_summary_metrics_pmf.csv
#          data/main_summary/main_summary_metrics_lambda.csv
#          data/main_summary/main_summary_empirical_interpolated.csv

library(here)
library(future.apply)

source(here("src", "aux_functions.R"))

# Read the simulation files in parallel
plan(multisession)

dir_pmf    <- here("data", "main_sims_pmf")
dir_lambda <- here("data", "main_sims_lambda")
dir_out    <- here("data", "main_summary")
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)

#--------------------------------------------------------------#
#-Extinction curves of both scenarios--------------------------#
#--------------------------------------------------------------#

# PMF is the intensity scenario and Lambda the frequency scenario
metrics_pmf <- process_sim_dir(dir_pmf, "PMF")
fwrite(metrics_pmf, file.path(dir_out, "main_summary_metrics_pmf.csv"))

metrics_lambda <- process_sim_dir(dir_lambda, "Lambda")
fwrite(metrics_lambda, file.path(dir_out, "main_summary_metrics_lambda.csv"))

#--------------------------------------------------------------#
#-Extinctions at the empirical frequency of each island--------#
#--------------------------------------------------------------#

# Lambda was simulated in steps of 0.01, and the empirical lambda of an
# island can fall between two simulated values
emp_lambdas <- fread(here("src", "empirical_lambdas.csv"))

interp <- interpolate_at_lambda(metrics_lambda, emp_lambdas)
fwrite(interp, file.path(dir_out, "main_summary_empirical_interpolated.csv"))

message("Saved: main_summary_metrics_pmf.csv, main_summary_metrics_lambda.csv ",
        "and main_summary_empirical_interpolated.csv")
