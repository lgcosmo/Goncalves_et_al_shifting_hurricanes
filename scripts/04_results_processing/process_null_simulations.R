#--------------------------------------------------------------#
#-Processing the simulations of the decoupling null models-----#
#--------------------------------------------------------------#

# Summarize the hurricane simulations of null models I to II-B. Each
# simulation file has one null network under the empirical hurricane regime
# of its island, with 1,000 replicates. For every null network, compute the
# fraction of species extinct per replicate (mean and sd over replicates)
# and the fraction of replicates in which each species went extinct.
#
# Inputs:  data/sims/<island>_null_<model>_<id>_multilayer_*.csv
# Outputs: data/main_summary/null_models_network_extinction.csv
#          data/main_summary/null_models_species_extinction.csv

library(data.table)
library(here)
library(parallel)

source(here("src", "aux_functions.R"))

#--------------------------------------------------------------#
#-1. Settings--------------------------------------------------#
#--------------------------------------------------------------#

# Null models to process
null_models <- c("I", "I-A", "II", "II-B")

# Number of parallel workers
n_workers <- max(1, detectCores() - 1)

dir_sims <- here("data", "sims")
dir_out  <- here("data", "main_summary")
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)

#--------------------------------------------------------------#
#-2. Simulation files------------------------------------------#
#--------------------------------------------------------------#

# The file name starts with the island, the null model and the id of the
# null network, for example 05_dominica_null_I-A_0042
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

# Summarize the simulation file in row i of files at the network and species
# levels. Returns a list with the two tables, both identified by island,
# null model and null network.
summarize_null_network <- function(i, files) {

  extinctions <- summarize_extinctions(files$file[i])
  ids <- files[i, .(island, model, null_id, network_id)]

  network <- cbind(ids, extinctions$network[, .(n_species, n_replicates,
                                                mean_fraction_extinct,
                                                sd_fraction_extinct)])

  # Layers in which the species interacts in this null network (Pollination,
  # Seed-dispersal or Multilayer). Can differ from its empirical type.
  species <- cbind(ids, extinctions$species[, .(species_id, null_species_type=type,
                                                n_replicates,
                                                fraction_replicates_extinct=fraction_extinct)])

  list(network=network, species=species)
}

# One data.table thread per worker, to avoid competition for the same cores.
# The workers also need summarize_extinctions.
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
