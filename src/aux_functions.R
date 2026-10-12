#--------------------------------------------------------------#
#-Helpers of the results-processing scripts--------------------#
#--------------------------------------------------------------#

# Functions shared by the scripts of scripts/04_results_processing.
# summarize_extinctions reads one simulation file and is the single place
# where extinctions are counted. survival_curves.R builds the extinction
# curves with calculate_survival_curve and process_sim_dir and interpolates
# them with interpolate_at_lambda; process_null_simulations.R and
# null_model_comparison.R summarize the null and empirical networks with
# summarize_extinctions; calculate_auc_extinction.R computes species
# fragility with compute_auc and calculate_auc_se. Each script sources this
# file after loading its packages.
#
# Inputs:  simulation files written by the scripts of 03_simulations, each
#          holding one network under one hurricane regime with one row per
#          species and replicate
# Outputs: none. The functions return tables that the scripts write.

library(data.table)
library(tidyverse)
library(future.apply)
library(here)

#--------------------------------------------------------------#
#-Extinctions in one simulation file---------------------------#
#--------------------------------------------------------------#

# Count the extinctions of one simulation file at the network level and at
# the species level. At the network level, the fraction of species that went
# extinct in each replicate is averaged over replicates, and its standard
# deviation measures the stochastic uncertainty of the hurricane simulations.
# At the species level, the fraction of replicates in which each species went
# extinct. The file holds one network under one hurricane regime, and every
# species appears in every replicate.
# Arguments: file is the simulation file. Returns a list with two
# data.tables. network has one row with network, lambda, intensity,
# n_species, n_replicates, mean_fraction_extinct and sd_fraction_extinct.
# species has one row per species with species_id, type (the layers in which
# the species interacts: Pollination, Seed-dispersal or Multilayer),
# fraction_extinct, n_replicates, network, lambda and intensity.
summarize_extinctions <- function(file) {

  sims <- fread(file, select=c("network", "sim_id", "species_id", "type",
                               "is_extinct", "lambda", "intensity"))
  sims[, is_extinct := as.logical(is_extinct)]

  # The island and the hurricane regime are the same in every row of a file
  regime <- sims[1, .(network, lambda, intensity)]

  # Network level: fraction of species lost in each replicate, then its mean
  # and standard deviation over replicates
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
#-Extinction curves of the frequency and intensity sweeps------#
#--------------------------------------------------------------#

# Arrange the extinctions of one simulation file as the rows of an extinction
# curve, at three levels. The species level holds the fraction of simulations
# in which each species went extinct. The type level averages that fraction
# over the species of each layer membership. The network level holds the
# fraction of species that went extinct, averaged over replicates, which
# equals the mean over species because every species appears in every
# replicate. figure_3.R draws the network level, and
# calculate_auc_extinction.R computes species fragility at all three levels.
# Arguments: file is the simulation file. Returns a data.table with one row
# per network, per type and per species holding network, lambda, intensity,
# level ("Network", "Type" or "Species"), name (the island, the type or the
# species), extinction_fraction, type ("All" at the network level) and
# species_id (empty above the species level).
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

  # One table with the three levels, the network row first, and the island
  # and hurricane regime of the file on every row
  metrics <- rbindlist(list(network, types, species), use.names=TRUE)
  metrics <- cbind(regime, metrics)
  setcolorder(metrics, c("network", "lambda", "intensity", "level", "name",
                         "extinction_fraction", "type", "species_id"))

  return(metrics)
}

# Build the extinction curves of every simulation file of one sweep, in
# parallel, and bind them.
# Arguments: dir is the folder of simulation files; label names the sweep in
# the progress message ("PMF" or "Lambda"). Returns one data.table with the
# network, type and species rows of every simulation file.
process_sim_dir <- function(dir, label) {

  files <- list.files(dir, pattern="\\.csv$", full.names=TRUE)
  if (length(files) == 0) stop(paste("No simulation files found in:", dir))
  message(label, " sweep: ", length(files), " simulation files")

  curves <- future_lapply(files, calculate_survival_curve, future.seed=TRUE)

  rbindlist(curves)
}

#--------------------------------------------------------------#
#-Extinctions at the empirical hurricane frequency-------------#
#--------------------------------------------------------------#

# Interpolate the extinction fraction of each species at the empirical
# hurricane frequency of its island, for the points that place each island
# on its extinction curve in Fig. 3c. The frequency sweep simulates lambda on
# a grid, so the value at the empirical lambda is interpolated linearly
# between the two neighboring grid values, for each species, level and type.
# Arguments: metrics is the table of extinction curves of the frequency
# sweep; emp_lams holds the columns network and lambda with the empirical
# frequency of each island. Returns a data.table with species_id, level,
# type, network, the empirical lambda and the interpolated
# extinction_fraction.
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
#-Species fragility as the AUC of the extinction curves--------#
#--------------------------------------------------------------#

# Compute the normalized area under one extinction curve and its analytical
# standard error with the trapezoidal rule. The area is the species fragility
# that the pSEM (Fig. 4) explains. It is divided by the range of x, so that it
# lies between 0 and 1 and higher values mean greater fragility. Each point
# of the curve is a binomial proportion estimated from n_reps replicates. Its
# variance, y (1 - y) / n_reps, propagates through the trapezoid coefficients
# of the area. These coefficients are half the width of the first and last
# intervals at the ends of the curve and half the distance between the two
# neighbors at every interior point.
# Arguments: x is the hurricane frequency or intensity at each point; y is
# the extinction fraction at each x; n_reps is the number of replicates
# behind each y. Returns a list with auc and auc_se, both NA when the curve
# has fewer than two points or no range in x.
calculate_auc_se <- function(x, y, n_reps=1000) {

  ord <- order(x)
  x   <- x[ord]
  y   <- y[ord]
  n   <- length(x)

  if (n < 2) return(list(auc=NA, auc_se=NA))

  dx  <- diff(x)
  auc <- sum((y[-n] + y[-1]) / 2 * dx)

  # Variance of each point of the curve, and the trapezoid coefficient that
  # weights it in the area
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

# Apply calculate_auc_se to every extinction curve of one scenario, one curve
# per network, per type and per species in each island network.
# Arguments: df is the table of extinction curves of the scenario; x_col
# names the column of the sweep ("lambda" or "intensity"); scenario_label is
# written to the scenario column ("PMF" or "Lambda"). Returns a tibble with
# network, type, level, species_id, auc, auc_se and scenario.
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
