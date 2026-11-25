#--------------------------------------------------------------#
#-Helpers shared by the figure scripts-------------------------#
#--------------------------------------------------------------#

# Functions shared by the scripts in scripts/05_visualization.
#
# Inputs:  none
# Outputs: none. save_figure saves the files of one figure.

library(here)
library(ggplot2)

# Save one figure as PNG, PDF and SVG in figures/. name is the file name
# without extension. width and height are in inches.
save_figure <- function(figure, name, width, height, dpi=300, bg="white") {

  dir.create(here("figures"), showWarnings=FALSE, recursive=TRUE)

  for (format in c("png", "pdf", "svg")) {
    ggsave(here("figures", paste0(name, ".", format)), figure,
           width=width, height=height, dpi=dpi, bg=bg)
  }

  message("Saved: figures/", name, ".png, .pdf and .svg")
}
