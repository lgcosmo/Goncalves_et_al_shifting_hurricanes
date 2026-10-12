#--------------------------------------------------------------#
#-Building the Multilayer Networks of the Lesser Antilles------#
#--------------------------------------------------------------#

# Build the multilayer bird-pollination and seed-dispersal network of each of
# the 11 islands of the Lesser Antilles from the field records. The
# pollination layer records the plant species whose pollen each bird carried,
# and the seed-dispersal layer the plant species whose seeds appeared in the
# fecal samples of each bird. Each layer becomes a binary incidence matrix
# with plants in rows and birds in columns. The two layers of an island share
# the same rows and columns, with an empty row or column for a species
# recorded in one layer only, so that every later script can tell dual-role
# birds apart from single-role birds. The matrices in src/ at present hold
# the number of records of each interaction, from an earlier run without
# binarization; every script that reads them keeps only presence or absence,
# except the proportion of interactions involving dual-role species of
# calculate_macro_structure.R, which a rebuilt binary set would change.
#
# Inputs:  src/pl_barcode/<island>.csv (one file per island, one row per
#          sampled bird and one column per plant, holding the pollen records)
#          src/dataseed_carib_noweight.csv (one row per island, bird and
#          plant, holding the seed records)
# Outputs: src/multi_pol/B_<island>.csv (pollination layer)
#          src/multi_sd/S_<island>.csv (seed-dispersal layer)
# Used by: scripts/02_topology_analysis, build_null_decoupling.jl,
#          scripts/03_simulations and figure_2.R

library(here)
library(data.table)
library(stringr)

#--------------------------------------------------------------#
#-Pollen and seed records of the 11 islands--------------------#
#--------------------------------------------------------------#

# One file of pollen records per island, named by the island code that labels
# every later file (01_trinidad to 11_saba). The seed records of all islands
# share one file.
pollen_files <- list.files(here("src", "pl_barcode"), pattern="\\.csv$", full.names=TRUE)
pollen_records <- lapply(pollen_files, fread)
names(pollen_records) <- sub("\\.csv$", "", basename(pollen_files))

seed_records <- fread(here("src", "dataseed_carib_noweight.csv"))

# Species names are squished and set to title case before matching, so that
# a bird spelled differently in the two records still receives one column
# and stays a dual-role bird. Two pollen files are encoded in Latin-1 and
# carry non-breaking spaces after some plant names, which become ordinary
# spaces before the squish.
tidy_name <- function(x) {
  x <- as.character(x)
  latin1 <- !validUTF8(x)
  x[latin1] <- iconv(x[latin1], from="latin1", to="UTF-8")
  str_to_title(str_squish(x))
}

#--------------------------------------------------------------#
#-Aligned binary layers, written to src/multi_pol and multi_sd-#
#--------------------------------------------------------------#

# Pollination layer of one island: the pollen records of each plant, summed
# over the sampled birds of each species. Every column other than the island,
# the sample code and the bird species is a plant. Plants whose names
# coincide after tidying are summed into one row.
# Arguments: records is the table of pollen records of the island. Returns a
# matrix of record counts with plants in rows and birds in columns.
pollination_layer <- function(records) {

  plants <- setdiff(names(records), c("island", "code", "species"))
  counts <- records[, lapply(.SD, sum, na.rm=TRUE), by=.(bird=tidy_name(species)),
                    .SDcols=plants]

  layer <- t(as.matrix(counts[, ..plants]))
  colnames(layer) <- counts$bird
  layer <- rowsum(layer, tidy_name(plants))

  return(layer[rownames(layer) != "", , drop=FALSE])
}

# Seed-dispersal layer of one island: the seeds of each plant found in the
# fecal samples of each bird species.
# Arguments: records is the table of seed records of the island. Returns a
# matrix of seed counts with plants in rows and birds in columns.
dispersal_layer <- function(records) {

  counts <- records[, .(seeds=sum(seeds)), by=.(plant=tidy_name(plant), bird=tidy_name(bird))]
  wide <- dcast(counts, plant ~ bird, value.var="seeds", fill=0)

  layer <- as.matrix(wide[, -1])
  rownames(layer) <- wide$plant

  return(layer)
}

# Place one layer on the shared rows and columns of its island, filling with
# zeros the plants and birds recorded in the other layer only.
# Arguments: layer is a matrix with plants in rows and birds in columns;
# plants and birds are the shared row and column names. Returns the aligned
# matrix.
align_layer <- function(layer, plants, birds) {

  aligned <- matrix(0, nrow=length(plants), ncol=length(birds), dimnames=list(plants, birds))
  aligned[rownames(layer), colnames(layer)] <- layer

  return(aligned)
}

# Write one layer with the plant names in the first column, the format that
# every later script reads.
# Arguments: layer is the aligned matrix; path is the CSV file. Returns
# nothing.
write_layer <- function(layer, path) {
  write.csv(cbind(plant=rownames(layer), as.data.frame(layer)), path, row.names=FALSE)
}

#--------------------------------------------------------------#
#-Aligned binary layers, written to src/multi_pol and multi_sd-#
#--------------------------------------------------------------#

dir.create(here("src", "multi_pol"), showWarnings=FALSE)
dir.create(here("src", "multi_sd"), showWarnings=FALSE)

for (net in names(pollen_records)) {

  B <- pollination_layer(pollen_records[[net]])
  S <- dispersal_layer(seed_records[island == net])

  # Every plant and bird recorded in either layer, in alphabetical order
  plants <- sort(unique(c(rownames(B), rownames(S))))
  birds  <- sort(unique(c(colnames(B), colnames(S))))

  # Presence or absence of each interaction on the shared rows and columns
  B <- (align_layer(B, plants, birds) > 0) * 1
  S <- (align_layer(S, plants, birds) > 0) * 1

  # A sampled bird without pollen or seeds, or a plant recorded with no bird,
  # takes no part in the network and leaves the matrices
  keep_plants <- rowSums(B) + rowSums(S) > 0
  keep_birds  <- colSums(B) + colSums(S) > 0
  B <- B[keep_plants, keep_birds, drop=FALSE]
  S <- S[keep_plants, keep_birds, drop=FALSE]

  write_layer(B, here("src", "multi_pol", paste0("B_", net, ".csv")))
  write_layer(S, here("src", "multi_sd", paste0("S_", net, ".csv")))
}

message("Saved: the two layers of ", length(pollen_records), " islands in src/multi_pol and src/multi_sd")
