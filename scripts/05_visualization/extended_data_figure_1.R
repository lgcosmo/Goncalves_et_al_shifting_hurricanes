#--------------------------------------------------------------#
#-Extended Data Figure 1: empirical versus null extinctions----#
#--------------------------------------------------------------#

# Extended Data Figure 1: extinctions of the empirical network of each island
# compared with its null networks (null models I, I-A, II and II-B, 100 null
# networks each), all under the empirical hurricane regime of the island.
# The network figure shows the z-score of the fraction of species extinct.
# A second figure, not in the manuscript, shows the z-scores of each species
# with the network-level z-score of each null model. Negative z-scores mean
# fewer extinctions in the empirical network than in the null networks.
#
# Inputs:  data/main_summary/null_models_comparison.csv
# Outputs: figures/extended_data_figure_1.png, .pdf and .svg
#          figures/extended_data_figure_1_species.png, .pdf and .svg

library(here)
library(data.table)
library(ggplot2)
library(ggpubr)

source(here("src", "island_labels.R"))
source(here("src", "figure_helpers.R"))

#--------------------------------------------------------------#
#-1. Load and prepare the data---------------------------------#
#--------------------------------------------------------------#

null_results <- fread(here("data", "main_summary", "null_models_comparison.csv"))

# Names of the null models in the legend
model_labels <- c(
  "I"    = "Null I (bird species focused, full decoupling)",
  "I-A"  = "Null I-A (bird species focused, half decoupling)",
  "II"   = "Null II (interaction focused, full decoupling)",
  "II-B" = "Null II-B (interaction focused, half decoupling)"
)

# One row per island, from Trinidad (south) to Saba (north). The four null
# models are offset within the row, from null I at the top to null II-B at
# the bottom.
model_offsets <- c("I"=0.27, "I-A"=0.09, "II"=-0.09, "II-B"=-0.27)

null_results[, model := factor(model, levels=names(model_labels))]
null_results[, y := match(island, names(island_names)) +
                    model_offsets[as.character(model)]]

network_results <- null_results[level == "network"]

# Label with the z-score to one decimal. Values within 0.05 of zero are
# written as 0.0 instead of -0.0.
network_results[, z_label := fifelse(abs(z_score) < 0.05, "0.0", sprintf("%.1f", z_score))]

species_results <- null_results[level == "species"]

#--------------------------------------------------------------#
#-2. Figure settings-------------------------------------------#
#--------------------------------------------------------------#

# Green for the bird species focused null models and orange for the
# interaction focused ones. Circles for full and triangles for half
# decoupling.
pal_model <- c("I"="#1f7a5c", "I-A"="#1f7a5c", "II"="#d8843a", "II-B"="#d8843a")
shape_model <- c("I"=21, "I-A"=24, "II"=21, "II-B"=24)
ink <- "#243238"

theme_null <- theme_pubr(base_size=13) +
  theme(
    panel.grid.major.y=element_line(color="#e4e7e6", linewidth=0.35),
    panel.grid.minor=element_blank(),
    axis.text=element_text(color=ink),
    axis.title=element_text(color=ink),
    legend.position="bottom",
    legend.title=element_blank(),
    legend.text=element_text(color=ink),
    plot.margin=margin(8, 12, 8, 8)
  )

# Elements shared by both figures: island rows, zero line, and the colors
# and shapes of the null models with one legend in two rows. Returns a list
# of ggplot components.
null_model_layout <- function() {
  list(
    geom_vline(xintercept=0, color=ink, linewidth=0.55),
    scale_y_continuous(breaks=seq_along(island_order), labels=island_order,
                       limits=c(0.55, length(island_order) + 0.45), expand=c(0, 0)),
    scale_color_manual(values=pal_model, guide="none"),
    scale_fill_manual(values=pal_model, labels=model_labels, breaks=names(model_labels)),
    scale_shape_manual(values=shape_model, labels=model_labels, breaks=names(model_labels)),
    guides(fill=guide_legend(ncol=2, byrow=TRUE), shape=guide_legend(ncol=2, byrow=TRUE)),
    labs(x="Empirical departure from null model (z-score)", y=NULL),
    theme_null
  )
}

#--------------------------------------------------------------#
#-3. Extended Data Figure 1: network-level z-scores------------#
#--------------------------------------------------------------#

# Each line goes from the null expectation (zero) to the empirical z-score,
# with the value written next to the point
p1 <- ggplot(network_results) +
  null_model_layout() +
  geom_segment(aes(x=0, xend=z_score, y=y, yend=y, color=model),
               linewidth=0.9, alpha=0.75) +
  geom_point(aes(x=z_score, y=y, fill=model, shape=model),
             color=ink, size=3.2, stroke=0.6) +
  geom_text(aes(x=z_score, y=y, label=z_label, hjust=ifelse(z_score < 0, 1.5, -0.5)),
            color=ink, size=2.8) +
  scale_x_continuous(breaks=seq(-12, 2, by=2)) +
  coord_cartesian(xlim=c(-12, 1.5), clip="off")

save_figure(p1, "extended_data_figure_1", width=10, height=10, dpi=400)

#--------------------------------------------------------------#
#-4. Species figure: species z-scores, not in the manuscript---#
#--------------------------------------------------------------#

# Boxes and small points are the z-scores of the species, in the color of
# their null model. The large markers are the network-level z-scores.
p2 <- ggplot() +
  null_model_layout() +
  geom_boxplot(data=species_results,
               aes(x=z_score, y=y, group=interaction(island, model), fill=model),
               orientation="y", width=0.14, outlier.shape=NA,
               alpha=0.55, color=ink, linewidth=0.3, show.legend=FALSE) +
  geom_point(data=species_results, aes(x=z_score, y=y),
             color=ink, size=0.5, alpha=0.25,
             position=position_jitter(width=0, height=0.025, seed=42)) +
  geom_point(data=network_results, aes(x=z_score, y=y, fill=model, shape=model),
             color=ink, size=2.8, stroke=0.6) +
  scale_x_continuous(breaks=seq(-12, 12, by=4)) +
  coord_cartesian(xlim=c(-14.5, 12.5))

save_figure(p2, "extended_data_figure_1_species", width=10, height=10, dpi=400)
