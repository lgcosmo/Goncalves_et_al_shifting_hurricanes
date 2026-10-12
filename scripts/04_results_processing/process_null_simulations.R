#--------------------------------------------------------------#
#-Processing the simulations of the decoupling null models-----#
#--------------------------------------------------------------#

# Summarize the hurricane simulations of null models I to II-B into the
# extinctions that null_model_comparison.R compares with the empirical
# networks (Extended Data Fig. 1). Each simulation file holds one null
# network under its island's empirical hurricane regime, with 1,000
# replicates. Every null network yields two summaries. The network level
# holds the fraction of species that went extinct in each replicate, with its
# mean and standard deviation across replicates. The species level holds the
# fraction of simulations in which each species went extinct.
#
# Inputs:  data/sims/<island>_null_<model>_<id>_multilayer_*.csv
#          (main_simulations_null_decoupling.jl)
# Outputs: data/main_summary/null_models_network_extinction.csv
#          data/main_summary/null_models_species_extinction.csv
# Helpers: summarize_extinctions in src/aux_functions.R
# Used by: null_model_comparison.R

library(data.table)
library(here)
library(parallel)

source(here("src", "aux_functions.R"))

#--------------------------------------------------------------#
#-1. Settings--------------------------------------------------#
#--------------------------------------------------------------#

# Null models to process. The simulation files of null model X are named
# <island>_null_X_<id>_multilayer_*.csv
null_models <- c("I", "I-A", "II", "II-B")

# Number of parallel workers; each reads and summarizes one file at a time
n_workers <- max(1, detectCores() - 1)

dir_sims <- here("data", "sims")
dir_out  <- here("data", "main_summary")
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)

#--------------------------------------------------------------#
#-2. Simulation files------------------------------------------#
#--------------------------------------------------------------#

# The network name at the start of each file name identifies the island, the
# null model and the replicate, for example 05_dominica_null_I-A_0042
files <- data.table(file=list.files(dir_sims, pattern="_null_[A-Z0-9-]+_[0-9]{4}_multilayer_.*\\.csv$",
                                    full.names=TRUE))
files[, network_id := sub("_multilayer_.*$", "", basename(file))]
files[, island     := sub("_null_[A-Z0-9-]+_[0-9]{4}$", "", network_id)]
files[, model      := sub("^.*_null_([A-Z0-9-]+)_[0-9]{4}$", "\\1", network_id)]
files[, null_id    := as.integer(sub("^.*_([0-9]{4})$", "\\1", network_id))]
files <- files[model %in% null_models]

message("Simulation files found per null model:")
print(files[, .N, by=model])

#--------------------------------------------------------------#
#-3. Extinctions in each null network--------------------------#
#--------------------------------------------------------------#

# Summarize one simulation file at the network and species levels with
# summarize_extinctions, and tag both tables with the island, null model and
# null network of the file.
# Arguments: i is the row of the file in files, the table of simulation
# files. Returns a list with the network and species tables of that null
# network.
summarize_null_network <- function(i, files) {

  extinctions <- summarize_extinctions(files$file[i])
  ids <- files[i, .(island, model, null_id, network_id)]

  network <- cbind(ids, extinctions$network[, .(n_species, n_replicates,
                                                mean_fraction_extinct,
                                                sd_fraction_extinct)])

  # The species type gives the layers in which the species interacts in this
  # null network (Pollination, Seed-dispersal or Multilayer) and can differ
  # from its empirical type
  species <- cbind(ids, extinctions$species[, .(species_id, null_species_type=type,
                                                n_replicates,
                                                fraction_replicates_extinct=fraction_extinct)])

  list(network=network, species=species)
}

# Each worker uses a single data.table thread, so the workers do not compete
# for the same cores, and receives the summary function
cluster <- makeCluster(n_workers)
invisible(clusterEvalQ(cluster, {
  library(data.table)
  setDTthreads(1)
}))
clusterExport(cluster, "summarize_extinctions")
summaries <- parLapply(cluster, seq_len(nrow(files)), summarize_null_network, files=files)
stopCluster(cluster)

#--------------------------------------------------------------#
#-4. Save the summaries----------------------------------------#
#--------------------------------------------------------------#

# One row per null network, and one row per species and null network
network_extinction <- rbindlist(lapply(summaries, function(x) x$network))
species_extinction <- rbindlist(lapply(summaries, function(x) x$species))

setorder(network_extinction, model, island, null_id)
setcolorder(species_extinction, c("island", "model", "null_id", "network_id", "species_id",
                                  "null_species_type", "n_replicates",
                                  "fraction_replicates_extinct"))
setorder(species_extinction, model, island, null_id, species_id)

fwrite(network_extinction, file.path(dir_out, "null_models_network_extinction.csv"))
fwrite(species_extinction, file.path(dir_out, "null_models_species_extinction.csv"))

message("Saved: null_models_network_extinction.csv and null_models_species_extinction.csv")
