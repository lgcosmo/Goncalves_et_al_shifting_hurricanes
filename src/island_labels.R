#--------------------------------------------------------------#
#-Island codes, names and order--------------------------------#
#--------------------------------------------------------------#

# The code of each island network (01_trinidad to 11_saba, the prefix of
# every network file) and the island name used in the manuscript, in the
# south-to-north order of the islands. Every script that labels islands
# sources this file, so that the names and their order match across the
# figures and tables.
#
# Inputs:  none
# Outputs: island_names, a named vector from island code to island name;
#          island_order, the island names from south to north
# Used by: build_island_data.R, figure_2.R, figure_3.R,
#          extended_data_figure_1.R and diagnostic/plot_null_networks.R

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
