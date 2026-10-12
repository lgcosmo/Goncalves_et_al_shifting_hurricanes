#--------------------------------------------------------------#
#-Helpers shared by the figure scripts-------------------------#
#--------------------------------------------------------------#

# Functions that every figure script of scripts/05_visualization sources
# after loading its packages.
#
# Inputs:  none
# Outputs: none. save_figure writes the files of one figure.
# Used by: figure_2.R, figure_3.R, figure_4.R and extended_data_figure_1.R

library(here)
library(ggplot2)

# Save one figure as PNG, PDF and SVG in figures/, the three formats the
# manuscript and its submission use.
# Arguments: figure is the ggplot object; name is the file name without
# extension (figure_2, extended_data_figure_1); width and height are in
# inches; dpi is the resolution of the PNG; bg is the background color.
# Returns nothing.
save_figure <- function(figure, name, width, height, dpi=300, bg="white") {

  dir.create(here("figures"), showWarnings=FALSE, recursive=TRUE)

  for (format in c("png", "pdf", "svg")) {
    ggsave(here("figures", paste0(name, ".", format)), figure,
           width=width, height=height, dpi=dpi, bg=bg)
  }

  message("Saved: figures/", name, ".png, .pdf and .svg")
}
