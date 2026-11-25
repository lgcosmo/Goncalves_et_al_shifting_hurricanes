#--------------------------------------------------------------#
#-Comparing the empirical networks with the null models--------#
#--------------------------------------------------------------#

# Compare the extinctions of the empirical network of each island with null
# models I to II-B under the empirical hurricane regime of the island, at
# the network level (fraction of species extinct) and at the species level
# (fraction of replicates in which a species went extinct). The z-score is
# the empirical value minus the mean of the null networks, divided by their
# sd. It gets large when the null networks differ little from each other.
# The difference (empirical minus null mean) has a 95% confidence interval
# that combines the noise of the 1,000 empirical replicates and of the 100
# null networks. Negative values mean fewer extinctions in the empirical
# network.
#
# Inputs:  data/main_sims_pmf/*_pmf0.0_tyears300.csv
#          data/main_summary/null_models_network_extinction.csv
#          data/main_summary/null_models_species_extinction.csv
# Outputs: data/main_summary/null_models_comparison.csv

library(data.table)
library(here)

source(here("src", "aux_functions.R"))

dir_empirical <- here("data", "main_sims_pmf")
dir_summary   <- here("data", "main_summary")

#--------------------------------------------------------------#
#-1. Extinctions in the empirical networks---------------------#
#--------------------------------------------------------------#

# An intensity shift of 0 is the empirical hurricane regime of each island
empirical_files <- list.files(dir_empirical, pattern="pmf0\\.0_tyears300\\.csv$",
                              full.names=TRUE)

# One file per island. Standard error at the network level from the sd over
# replicates, and at the species level from the binomial distribution.
empirical <- lapply(empirical_files, function(file) {

  extinctions <- summarize_extinctions(file)

  network <- extinctions$network[, .(island=network,
                                     empirical=mean_fraction_extinct,
                                     empirical_se=sd_fraction_extinct / sqrt(n_replicates))]

  species <- extinctions$species[, .(island=network, species_id,
                                     empirical_species_type=type,
                                     empirical=fraction_extinct,
                                     empirical_se=sqrt(fraction_extinct * (1 - fraction_extinct) / n_replicates))]

  list(network=network, species=species)
})

empirical_network <- rbindlist(lapply(empirical, function(x) x$network))
empirical_species <- rbindlist(lapply(empirical, function(x) x$species))

#--------------------------------------------------------------#
#-2. Extinctions in the null models----------------------------#
#--------------------------------------------------------------#

# Mean and sd of each null model on each island over its null networks. The
# mean of each null network carries its own simulation noise, so the standard
# error of the null mean (sd / sqrt(n)) includes both the differences among
# null networks and the noise of the simulations.
null_networks <- fread(file.path(dir_summary, "null_models_network_extinction.csv"))
null_network <- null_networks[, .(
  n_null_networks=.N,
  null_mean=mean(mean_fraction_extinct),
  null_sd=sd(mean_fraction_extinct)
), by=.(island, model)]

null_species_all <- fread(file.path(dir_summary, "null_models_species_extinction.csv"))
null_species <- null_species_all[, .(
  n_null_networks=.N,
  null_mean=mean(fraction_replicates_extinct),
  null_sd=sd(fraction_replicates_extinct)
), by=.(island, model, species_id)]

#--------------------------------------------------------------#
#-3. Compare the empirical networks with each null model-------#
#--------------------------------------------------------------#

# Compute the two measures for a table with the empirical and null values.
# Returns the table with the new columns.
compare <- function(comparison) {
  comparison[, null_se := null_sd / sqrt(n_null_networks)]

  # z-score: is the empirical network extreme among the null networks?
  comparison[, z_score := fifelse(null_sd > 0, (empirical - null_mean) / null_sd,
                                  NA_real_)]

  # Difference and its 95% confidence interval: do the empirical and null
  # extinctions differ by more than the noise of the simulations?
  comparison[, difference := empirical - null_mean]
  comparison[, se_difference := sqrt(empirical_se^2 + null_se^2)]
  comparison[, `:=`(ci_lower=difference - 1.96 * se_difference,
                    ci_upper=difference + 1.96 * se_difference)]
  comparison[, differs := ci_lower > 0 | ci_upper < 0]
  comparison
}

network_comparison <- compare(merge(null_network, empirical_network, by="island"))
network_comparison[, `:=`(level="network", species_id=NA_character_,
                          empirical_species_type=NA_character_)]

species_comparison <- compare(merge(null_species, empirical_species,
                                    by=c("island", "species_id")))
species_comparison[, level := "species"]

null_comparison <- rbindlist(list(network_comparison, species_comparison), use.names=TRUE)
setcolorder(null_comparison, c(
  "model", "island", "level", "species_id", "empirical_species_type",
  "empirical", "empirical_se", "null_mean", "null_sd", "null_se", "n_null_networks",
  "z_score", "difference", "se_difference", "ci_lower", "ci_upper", "differs"
))
setorder(null_comparison, model, island, level, species_id)

fwrite(null_comparison, file.path(dir_summary, "null_models_comparison.csv"))

#--------------------------------------------------------------#
#-4. Print the network-level results---------------------------#
#--------------------------------------------------------------#

print(null_comparison[level == "network",
                      .(model, island, empirical=round(empirical, 3),
                        null_mean=round(null_mean, 3), z_score=round(z_score, 2),
                        difference=round(difference, 3), ci_lower=round(ci_lower, 3),
                        ci_upper=round(ci_upper, 3), differs)])

message("Saved: null_models_comparison.csv")
