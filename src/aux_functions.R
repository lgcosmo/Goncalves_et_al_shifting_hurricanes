#--------------------------------------------------------------#
#-Helpers of the results-processing scripts--------------------#
#--------------------------------------------------------------#

# Functions shared by the scripts in scripts/04_results_processing: counting
# the extinctions in the simulation files, building the extinction curves
# and computing the area under them (AUC).
#
# Inputs:  simulation files from scripts/03_simulations, each with one
#          network under one hurricane regime and one row per species and
#          replicate
# Outputs: none. The functions return tables that the scripts save.

library(data.table)
library(tidyverse)
library(future.apply)
library(here)

#--------------------------------------------------------------#
#-Extinctions in one simulation file---------------------------#
#--------------------------------------------------------------#

# Count the extinctions of one simulation file. At the network level, the
# fraction of species extinct in each replicate, summarized as mean and sd
# over replicates. At the species level, the fraction of replicates in which
# each species went extinct.
# Returns a list with the data.tables network (one row) and species (one row
# per species).
summarize_extinctions <- function(file) {

  sims <- fread(file, select=c("network", "sim_id", "species_id", "type",
                               "is_extinct", "lambda", "intensity"))
  sims[, is_extinct := as.logical(is_extinct)]

  # The island and the hurricane regime are the same in every row of a file
  regime <- sims[1, .(network, lambda, intensity)]

  # Network level: fraction of species extinct per replicate, then mean and sd
  per_replicate <- sims[, .(fraction_extinct=mean(is_extinct)), by=sim_id]
  network <- data.table(
    regime,
    n_species=uniqueN(sims$species_id),
    n_replicates=nrow(per_replicate),
    mean_fraction_extinct=mean(per_replicate$fraction_extinct),
    sd_fraction_extinct=sd(per_replicate$fraction_extinct)
  )

  # Species level: fraction of replicates in which each species went extinct
  species <- sims[, .(fraction_extinct=mean(is_extinct), n_replicates=.N),
                  by=.(species_id, type)]
  species <- cbind(species, regime)

  list(network=network, species=species)
}

#--------------------------------------------------------------#
#-Extinction curves of the two hurricane scenarios-------------#
#--------------------------------------------------------------#

# Arrange the extinctions of one simulation file as the rows of an extinction
# curve at three levels: each species, the mean over the species of each type
# (Pollination, Seed-dispersal or Multilayer) and the whole network.
# Returns a data.table with network, lambda, intensity, level, name,
# extinction_fraction, type and species_id.
calculate_survival_curve <- function(file) {

  extinctions <- summarize_extinctions(file)
  regime <- extinctions$network[, .(network, lambda, intensity)]

  species <- extinctions$species[, .(species_id, type, extinction_fraction=fraction_extinct)]
  species[, `:=`(level="Species", name=species_id)]

  types <- species[, .(extinction_fraction=mean(extinction_fraction)), by=type]
  types[, `:=`(level="Type", name=type, species_id=NA_character_)]

  network <- data.table(extinction_fraction=extinctions$network$mean_fraction_extinct,
                        level="Network", name=regime$network, type="All",
                        species_id=NA_character_)

  # One table with the three levels and the hurricane regime of the file
  metrics <- rbindlist(list(network, types, species), use.names=TRUE)
  metrics <- cbind(regime, metrics)
  setcolorder(metrics, c("network", "lambda", "intensity", "level", "name",
                         "extinction_fraction", "type", "species_id"))

  return(metrics)
}

# Build the extinction curves of all simulation files in a folder, in
# parallel. label ("PMF" or "Lambda") identifies the scenario in the progress
# message. Returns one data.table with all curves.
process_sim_dir <- function(dir, label) {

  files <- list.files(dir, pattern="\\.csv$", full.names=TRUE)
  if (length(files) == 0) stop(paste("No simulation files found in:", dir))
  message(label, " scenario: ", length(files), " simulation files")

  curves <- future_lapply(files, calculate_survival_curve, future.seed=TRUE)

  rbindlist(curves)
}

#--------------------------------------------------------------#
#-Extinctions at the empirical hurricane frequency-------------#
#--------------------------------------------------------------#

# Interpolate the extinction fraction at the empirical hurricane frequency of
# each island. Lambda was simulated in steps of 0.01, so the value at the
# empirical lambda is interpolated linearly between the two closest values.
# emp_lams has the columns network and lambda. Returns a data.table with one
# row per species, level and type.
interpolate_at_lambda <- function(metrics, emp_lams) {

  results <- lapply(seq_len(nrow(emp_lams)), function(i) {

    net <- emp_lams$network[i]
    lam <- emp_lams$lambda[i]
    sub <- metrics[network == net]
    setorder(sub, lambda)

    sub[, .(
      lambda              = lam,
      extinction_fraction = approx(x=lambda, y=extinction_fraction, xout=lam)$y,
      network             = net
    ), by = .(species_id, level, type)]
  })

  rbindlist(results)
}

#--------------------------------------------------------------#
#-Area under the extinction curves (AUC)-----------------------#
#--------------------------------------------------------------#

# Normalized area under one extinction curve (trapezoidal rule) and its
# standard error. The area is divided by the range of x, so it lies between 0
# and 1. Species persistence is one minus this area. Each point of the curve
# is a binomial proportion from n_reps replicates, with variance
# y (1 - y) / n_reps, and these variances are propagated through the
# trapezoid coefficients.
# Returns a list with auc and auc_se, both NA when the curve has fewer than
# two points or no range in x.
calculate_auc_se <- function(x, y, n_reps=1000) {

  ord <- order(x)
  x   <- x[ord]
  y   <- y[ord]
  n   <- length(x)

  if (n < 2) return(list(auc=NA, auc_se=NA))

  dx  <- diff(x)
  auc <- sum((y[-n] + y[-1]) / 2 * dx)

  # Variance of each point of the curve and its trapezoid coefficient
  var_y           <- (y * (1 - y)) / n_reps
  c_coeff         <- numeric(n)
  c_coeff[1]      <- 0.5 * dx[1]
  c_coeff[n]      <- 0.5 * dx[n - 1]
  if (n > 2) c_coeff[2:(n-1)] <- 0.5 * (x[3:n] - x[1:(n-2)])

  auc_se  <- sqrt(sum(c_coeff^2 * var_y))
  x_range <- max(x) - min(x)

  if (x_range == 0) return(list(auc=NA, auc_se=NA))
  list(auc=auc / x_range, auc_se=auc_se / x_range)
}

# Apply calculate_auc_se to the extinction curves of one scenario, one curve
# per network, type and species. x_col is the column with the values of
# lambda or of the intensity shift. Returns a tibble with auc and auc_se.
compute_auc <- function(df, x_col, scenario_label) {
  df %>%
    group_by(network, type, level, species_id) %>%
    summarise(
      calc = list(calculate_auc_se(.data[[x_col]], extinction_fraction)),
      .groups = "drop"
    ) %>%
    mutate(
      auc      = map_dbl(calc, "auc"),
      auc_se   = map_dbl(calc, "auc_se"),
      scenario = scenario_label
    ) %>%
    select(-calc)
}
