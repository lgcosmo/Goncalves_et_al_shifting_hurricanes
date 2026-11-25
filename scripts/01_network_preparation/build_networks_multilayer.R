#--------------------------------------------------------------#
#-Building the multilayer networks of the Lesser Antilles------#
#--------------------------------------------------------------#

# Build the multilayer network of each of the 11 islands from the field
# records. The pollination layer has the plants whose pollen each bird
# carried. The seed-dispersal layer has the plants whose seeds were found in
# the fecal samples of each bird. Each layer is a binary incidence matrix
# with plants in rows and birds in columns. The two layers of an island share
# the same rows and columns.
# Note: the matrices stored in src/ have the number of records of each
# interaction instead of 0 and 1. Only the proportion of interactions
# involving dual-role species uses these counts.
#
# Inputs:  src/pl_barcode/<island>.csv
#          src/dataseed_carib_noweight.csv
# Outputs: src/multi_pol/B_<island>.csv
#          src/multi_sd/S_<island>.csv

library(here)
library(data.table)
library(stringr)

#--------------------------------------------------------------#
#-Pollen and seed records of the 11 islands--------------------#
#--------------------------------------------------------------#

# One file of pollen records per island, named by the island code
# (01_trinidad to 11_saba). The seed records of all islands are in one file.
pollen_files <- list.files(here("src", "pl_barcode"), pattern="\\.csv$", full.names=TRUE)
pollen_records <- lapply(pollen_files, fread)
names(pollen_records) <- sub("\\.csv$", "", basename(pollen_files))

seed_records <- fread(here("src", "dataseed_carib_noweight.csv"))

# Standardize species names (single spaces, title case). A bird spelled
# differently in the pollen and seed records then gets one column and stays
# a dual-role bird. Two pollen files are encoded in Latin-1 and have
# non-breaking spaces after some plant names. They are converted to UTF-8
# before the squish.
tidy_name <- function(x) {
  x <- as.character(x)
  latin1 <- !validUTF8(x)
  x[latin1] <- iconv(x[latin1], from="latin1", to="UTF-8")
  str_to_title(str_squish(x))
}

#--------------------------------------------------------------#
#-Functions to build and align the two layers------------------#
#--------------------------------------------------------------#

# Pollination layer of one island: pollen records of each plant, summed over
# the sampled birds of each species. Every column other than island, code
# and species is a plant. Plants with the same name after tidying are summed
# into one row. Returns a matrix of record counts with plants in rows and
# birds in columns.
pollination_layer <- function(records) {

  plants <- setdiff(names(records), c("island", "code", "species"))
  counts <- records[, lapply(.SD, sum, na.rm=TRUE), by=.(bird=tidy_name(species)),
                    .SDcols=plants]

  layer <- t(as.matrix(counts[, ..plants]))
  colnames(layer) <- counts$bird
  layer <- rowsum(layer, tidy_name(plants))

  return(layer[rownames(layer) != "", , drop=FALSE])
}

# Seed-dispersal layer of one island: seeds of each plant found in the fecal
# samples of each bird species. Returns a matrix of seed counts with plants
# in rows and birds in columns.
dispersal_layer <- function(records) {

  counts <- records[, .(seeds=sum(seeds)), by=.(plant=tidy_name(plant), bird=tidy_name(bird))]
  wide <- dcast(counts, plant ~ bird, value.var="seeds", fill=0)

  layer <- as.matrix(wide[, -1])
  rownames(layer) <- wide$plant

  return(layer)
}

# Place one layer on the shared rows and columns of its island. Plants and
# birds recorded only in the other layer are filled with zeros.
# Returns the aligned matrix.
align_layer <- function(layer, plants, birds) {

  aligned <- matrix(0, nrow=length(plants), ncol=length(birds), dimnames=list(plants, birds))
  aligned[rownames(layer), colnames(layer)] <- layer

  return(aligned)
}

# Write one layer to a CSV file, with the plant names in the first column.
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

  # Binarize each layer on the shared rows and columns
  B <- (align_layer(B, plants, birds) > 0) * 1
  S <- (align_layer(S, plants, birds) > 0) * 1

  # Drop birds sampled without pollen or seeds and plants recorded with no
  # bird
  keep_plants <- rowSums(B) + rowSums(S) > 0
  keep_birds  <- colSums(B) + colSums(S) > 0
  B <- B[keep_plants, keep_birds, drop=FALSE]
  S <- S[keep_plants, keep_birds, drop=FALSE]

  write_layer(B, here("src", "multi_pol", paste0("B_", net, ".csv")))
  write_layer(S, here("src", "multi_sd", paste0("S_", net, ".csv")))
}

message("Saved: the two layers of ", length(pollen_records), " islands in src/multi_pol and src/multi_sd")
