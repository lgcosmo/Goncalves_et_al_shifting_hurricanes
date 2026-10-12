#--------------------------------------------------------------#
#-Extended Data Figure 1: Empirical versus Null Extinctions----#
#--------------------------------------------------------------#

# Plot how the extinctions of each island's empirical network compare with
# those of its null networks, all simulated under the island's empirical
# hurricane regime. The null networks decouple the pollination and
# seed-dispersal layers while holding the remaining topological properties
# fixed, with 100 null networks per island and null model (Methods, null
# models I, I-A, II and II-B). The network figure is Extended Data Figure 1.
# It shows the z-score of the empirical fraction of species that went extinct
# among the null networks. The species figure shows the z-scores of the
# fraction of simulations in which each species went extinct, with the
# network-level z-score of each null model. That figure is not yet placed in
# the manuscript. Negative z-scores mean that the empirical network lost
# fewer species than its null counterparts.
#
# Inputs:  data/main_summary/null_models_comparison.csv
#          (null_model_comparison.R)
# Outputs: figures/extended_data_figure_1.png, .pdf and .svg
#          figures/extended_data_figure_1_species.png, .pdf and .svg

library(here)
library(data.table)
library(ggplot2)
library(ggpubr)

source(here("src", "island_labels.R"))
source(here("src", "figure_helpers.R"))

#--------------------------------------------------------------#
#-1. Load and Prepare Data-------------------------------------#
#--------------------------------------------------------------#

null_results <- fread(here("data", "main_summary", "null_models_comparison.csv"))

# Names of the null models in the figure legends, keyed by the names used in
# the simulation files
model_labels <- c(
  "I"    = "Null I (bird species focused, full decoupling)",
  "I-A"  = "Null I-A (bird species focused, half decoupling)",
  "II"   = "Null II (interaction focused, full decoupling)",
  "II-B" = "Null II-B (interaction focused, half decoupling)"
)

# Each island occupies one row, from Trinidad in the south to Saba in the
# north, and its four null models sit side by side within the row, from null
# I at the top to null II-B at the bottom
model_offsets <- c("I"=0.27, "I-A"=0.09, "II"=-0.09, "II-B"=-0.27)

null_results[, model := factor(model, levels=names(model_labels))]
null_results[, y := match(island, names(island_names)) +
                    model_offsets[as.character(model)]]

network_results <- null_results[level == "network"]

# Point labels show the z-score to one decimal. Values within 0.05 of zero
# read 0.0 rather than -0.0.
network_results[, z_label := fifelse(abs(z_score) < 0.05, "0.0", sprintf("%.1f", z_score))]

species_results <- null_results[level == "species"]

#--------------------------------------------------------------#
#-2. Figure Settings-------------------------------------------#
#--------------------------------------------------------------#

# Green marks the bird species focused null models and orange the interaction
# focused ones; circles mark full decoupling and triangles half decoupling
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

# Elements shared by both figures: the island rows, the zero line, and the
# colors and shapes of the null models with a single legend in two rows.
# Returns a list of ggplot components, added to a plot with +.
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
#-3. Extended Data Figure 1: Network-level z-scores------------#
#--------------------------------------------------------------#

# Each line begins at the null expectation and ends at the empirical z-score.
# Point labels make the magnitude and direction of the response explicit.
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
#-4. Species Figure: Species z-scores, not in the manuscript---#
#--------------------------------------------------------------#

# The boxes and points show the response of individual species. Each box
# takes the color of its null model, and the solid marker gives the
# whole-network response of that null model, shown in the network figure.
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
