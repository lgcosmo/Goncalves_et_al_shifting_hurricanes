#--------------------------------------------------------------#
#-Island codes, names and order--------------------------------#
#--------------------------------------------------------------#

# Code of each island network (01_trinidad to 11_saba, the prefix of the
# network files) and the island name used in the manuscript, ordered from
# south to north.
#
# Inputs:  none
# Outputs: island_names (island code to island name) and island_order (the
#          island names from south to north)

island_names <- c(
  "01_trinidad"  = "Trinidad",
  "02_tobago"    = "Tobago",
  "03_grenada"   = "Grenada",
  "04_stvincent" = "St. Vincent",
  "05_dominica"  = "Dominica",
  "06_mtserrat"  = "Montserrat",
  "07_antigua"   = "Antigua",
  "08_nevis"     = "Nevis",
  "09_stkitts"   = "St. Kitts",
  "10_eustatios" = "St. Eustatius",
  "11_saba"      = "Saba"
)

island_order <- unname(island_names)
