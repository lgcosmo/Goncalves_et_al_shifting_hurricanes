#--------------------------------------------------------------#
#-Null Model Comparison----------------------------------------#
#--------------------------------------------------------------#

# Compare the extinctions of each island's empirical network with those of
# null models I to II-B, all simulated under the island's empirical hurricane
# regime. The comparison works at the network level (the fraction of species
# that went extinct) and at the species level (the fraction of replicates in
# which a species went extinct). Two measures compare the empirical value with
# each null model. The z-score, drawn in Extended Data Fig. 1, is the
# empirical value minus the mean of the null networks, divided by their
# standard deviation. Small differences among null networks make this
# standard deviation small and the z-score large. The difference with
# simulation uncertainty is the empirical value minus the null-model mean,
# with a 95% confidence interval that combines the noise of the 1,000
# empirical replicates and of the 100 null networks, each estimated from
# 1,000 replicates. Negative values mean that the empirical network lost
# fewer species than the null model.
#
# Inputs:  data/main_sims_pmf/*_pmf0.0_tyears300.csv (empirical networks under
#          the empirical regime)
#          data/main_summary/null_models_network_extinction.csv and
#          data/main_summary/null_models_species_extinction.csv, both written by
#          process_null_simulations.R
# Outputs: data/main_summary/null_models_comparison.csv
# Helpers: summarize_extinctions in src/aux_functions.R

library(data.table)
library(here)

source(here("src", "aux_functions.R"))

dir_empirical <- here("data", "main_sims_pmf")
dir_summary   <- here("data", "main_summary")

#--------------------------------------------------------------#
#-1. Extinctions in the empirical networks---------------------#
#--------------------------------------------------------------#

# The intensity shift of 0 reproduces each island's empirical hurricane regime
empirical_files <- list.files(dir_empirical, pattern="pmf0\\.0_tyears300\\.csv$",
                              full.names=TRUE)

# Only 11 files, one per island, so they are read one after the other. The
# standard error of the network level follows from the spread across
# replicates, and that of the species level from the binomial distribution.
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

# Summarize each null model on each island across its null networks. The
# standard error of the null-model mean divides the standard deviation across
# null networks by the square root of their number. Each network mean already
# carries its own simulation noise, so this standard error includes both the
# differences among null networks and the noise of the simulations.
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

# Both measures, for one table holding the empirical and null values
compare <- function(comparison) {
  comparison[, null_se := null_sd / sqrt(n_null_networks)]

  # z-score: is the empirical network extreme among the null networks?
  comparison[, z_score := fifelse(null_sd > 0, (empirical - null_mean) / null_sd,
                                  NA_real_)]

  # Difference with simulation uncertainty: do the empirical and null
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
