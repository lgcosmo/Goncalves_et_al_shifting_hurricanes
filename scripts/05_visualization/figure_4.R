#--------------------------------------------------------------#
#-Figure 4: pSEM path diagram----------------------------------#
#--------------------------------------------------------------#

# Draw Figure 4, the path diagram of the piecewise structural equation model
# (pSEM) fitted by structural_equation_model.R. The figure shows the
# significant directed paths through which island biogeography, network
# topology and species functional roles shape the fragility of bird species
# to hurricane frequency and to hurricane intensity. Each variable occupies a
# card with its name and an icon, and the cards sit in four panels, one per
# level of organization and one for the two fragility responses. Blue arrows
# are positive standardized path coefficients and red arrows negative ones,
# so a blue arrow into a response means greater fragility. Arrow width grows
# with the absolute coefficient. Correlated errors are not drawn. Panels,
# cards, arrows and icons are placed by hand in plot coordinates, which
# coord_fixed draws with equal units on both axes.
#
# Inputs:  tables/sem_coefficients_clean.csv (structural_equation_model.R)
#          figures/icons/*.png
# Outputs: figures/figure_4.png, .pdf and .svg

library(here)
library(tidyverse)
library(ggforce)
library(grid)
library(png)

source(here("src", "figure_helpers.R"))

#--------------------------------------------------------------#
#-Load the significant pSEM paths------------------------------#
#--------------------------------------------------------------#

# Only the directed paths that piecewiseSEM marks as significant become
# arrows. The row that tests island type as a whole has no estimate and is
# dropped; its two island types keep their own estimates.
sem_coef <- read_csv(here("tables", "sem_coefficients_clean.csv"),
                     show_col_types=FALSE) %>%
  filter(path_type == "Directed path",
         !is.na(estimate),
         !is.na(significance),
         significance != "")

#--------------------------------------------------------------#
#-Palette, arrow widths and drawing helpers--------------------#
#--------------------------------------------------------------#

# One color per panel, in three shades for the panel border, the panel
# background and the cards. Blue and red mark the sign of the paths.
pal_sem <- list(
  ink="#000000",
  positive="#0072B2",
  negative="#B23A2F",
  island="#6FB27C",
  island_fill="#EAF7EA",
  island_card="#D9F3D8",
  network="#9ECBF4",
  network_fill="#E6F3FF",
  network_card="#D5E9FA",
  trait="#FFA55C",
  trait_fill="#FFF1E5",
  trait_card="#FFE0C4",
  hurricane="#F9D875",
  hurricane_fill="#FFF6D9",
  hurricane_card="#FFE9AC"
)

# Range of arrow widths onto which the absolute coefficients are rescaled
path_width_range <- c(0.45, 2.10)

# Expand a table of rectangles into the four corners that geom_shape needs to
# draw rounded boxes. Arguments: df holds one rectangle per row, with its id,
# xmin, xmax, ymin, ymax, fill, color, alpha and linewidth. Returns a tibble
# with four rows per rectangle, its corners in drawing order.
rect_shape <- function(df) {

  rect_list <- lapply(seq_len(nrow(df)), function(i) {

    row_i <- df[i, ]

    tibble(
      id=row_i$id,
      x=c(row_i$xmin, row_i$xmax, row_i$xmax, row_i$xmin),
      y=c(row_i$ymin, row_i$ymin, row_i$ymax, row_i$ymax),
      fill=row_i$fill,
      color=row_i$color,
      alpha=row_i$alpha,
      linewidth=row_i$linewidth
    )
  })

  return(bind_rows(rect_list))
}

# Place one icon on the figure. Arguments: p is the ggplot; filename names
# the PNG in figures/icons; xmin, xmax, ymin and ymax bound the icon in the
# coordinates of the plot. Returns the plot with the icon added as a raster
# annotation.
add_png_icon <- function(p, filename, xmin, xmax, ymin, ymax) {

  img <- png::readPNG(here("figures", "icons", filename))
  grob <- rasterGrob(img, interpolate=TRUE)

  p <- p +
    annotation_custom(grob, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

  return(p)
}

#--------------------------------------------------------------#
#-Panels and cards of the path diagram-------------------------#
#--------------------------------------------------------------#

# The four panels, one per level of organization and one for the two
# fragility responses
panel_boxes <- tribble(
  ~id, ~xmin, ~xmax, ~ymin, ~ymax, ~fill, ~color, ~alpha, ~linewidth,
  "island_panel",    0.45,  3.70, 1.25, 8.20, pal_sem$island_fill,    pal_sem$island,    0.78, 0.55,
  "fragility_panel", 4.40, 10.60, 7.85, 9.85, pal_sem$hurricane_fill, pal_sem$hurricane, 0.88, 0.55,
  "trait_panel",    11.25, 14.90, 1.25, 8.20, pal_sem$trait_fill,     pal_sem$trait,     0.78, 0.55,
  "network_panel",   4.40, 10.60, 0.15, 2.15, pal_sem$network_fill,   pal_sem$network,   0.78, 0.55
)

# One card per variable. The label sits on one side of the card and the icon
# on the other, so label_x and label_y are off-centre.
cards <- tribble(
  ~id, ~group, ~label, ~xmin, ~xmax, ~ymin, ~ymax, ~fill, ~color, ~label_x, ~label_y, ~size,
  "area",      "island",    "Island\narea",                             0.70, 3.45, 6.75, 7.80, pal_sem$island_card,    pal_sem$island,    1.22, 7.27, 4.3,
  "elevation", "island",    "Island\nelevation",                        0.70, 3.45, 5.35, 6.40, pal_sem$island_card,    pal_sem$island,    1.22, 5.88, 4.0,
  "isolation", "island",    "Island\nisolation",                        0.70, 3.45, 3.95, 5.00, pal_sem$island_card,    pal_sem$island,    1.24, 4.47, 4.0,
  "continental","island",   "Continental",                              0.92, 3.25, 2.65, 3.38, pal_sem$island_card,    pal_sem$island,    1.50, 3.02, 3.8,
  "oceanic",    "island",   "Oceanic",                                  0.92, 3.25, 1.85, 2.58, pal_sem$island_card,    pal_sem$island,    1.50, 2.22, 3.8,
  "frequency",  "outcome",  "Fragility to hurricane frequency",         4.60, 7.35, 8.15, 9.55, pal_sem$hurricane_card, pal_sem$hurricane, 5.98, 9.38, 4.2,
  "intensity",  "outcome",  "Fragility to hurricane intensity",         7.65,10.40, 8.15, 9.55, pal_sem$hurricane_card, pal_sem$hurricane, 9.02, 9.38, 4.2,
  "global",     "network",  "Global network\nintegration",              4.68, 7.25, 0.40, 1.90, pal_sem$network_card,   pal_sem$network,   5.96, 0.66, 4.0,
  "dominance",  "network",  "Pollination topological\ndominance",        7.75,10.32, 0.40, 1.90, pal_sem$network_card,   pal_sem$network,   9.04, 0.66, 3.8,
  "pol_degree", "trait",    "Pollination layer\ndegree (z-score)",       11.50,14.70, 6.75, 7.80, pal_sem$trait_card,     pal_sem$trait,     13.85, 7.28, 4.0,
  "sd_degree",  "trait",    "Dispersal layer\ndegree (z-score)",         11.50,14.70, 5.05, 6.10, pal_sem$trait_card,     pal_sem$trait,     13.85, 5.58, 4.0,
  "mpc",        "trait",    "Multilayer\nparticipation",                11.50,14.70, 3.35, 4.40, pal_sem$trait_card,     pal_sem$trait,     13.85, 3.88, 4.0,
  "body_mass",  "trait",    "Body mass",                                11.50,14.70, 1.65, 2.70, pal_sem$trait_card,     pal_sem$trait,     13.85, 2.17, 4.0
) %>%
  mutate(alpha=0.98, linewidth=0.50)

# A dashed box groups the continental and oceanic cards as island type
island_type_box <- tribble(
  ~id, ~xmin, ~xmax, ~ymin, ~ymax, ~fill, ~color, ~alpha, ~linewidth,
  "island_type", 0.70, 3.45, 1.65, 3.70, "#FFFFFF", pal_sem$island, 0.35, 0.38
)

# Names of the three levels along the sides of the figure, and the label of
# the island-type box
side_labels <- tribble(
  ~label, ~x, ~y, ~angle, ~color, ~size,
  "Island biogeography",       0.08, 4.65,  90, pal_sem$island, 6.2,
  "Species position and traits", 15.10, 4.65, -90, pal_sem$trait, 6.2,
  "Network topology",          7.50, -0.18,   0, "#6CC4D9", 6.0,
  "Island type",               2.05, 3.58,   0, pal_sem$ink, 3.4
)

#--------------------------------------------------------------#
#-Arrow layout of the significant paths------------------------#
#--------------------------------------------------------------#

# Point on the edge of one card where a path starts or ends: the middle of
# its left or right edge at the height of the label, or the middle of its top
# or bottom edge at the x of the label, so that a path stays attached to its
# card when the card moves.
# Arguments: id is the card id in cards; side is "left", "right", "top" or
# "bottom". Returns a named vector with x and y.
card_edge <- function(id, side) {
  card <- cards[cards$id == id, ]
  c(x=switch(side, left=card$xmin, right=card$xmax, card$label_x),
    y=switch(side, top=card$ymax, bottom=card$ymin, card$label_y))
}

# Each path leaves one side of its predictor card and enters one side of its
# response card. A path is drawn only when the coefficient table holds it as
# significant, so a path that becomes significant after a refit needs a new
# row here. The path from isolation enters the dominance card left of its
# label, clear of the path that leaves the card. The two paths from body mass
# to the roles above it are curved.
path_layout <- tribble(
  ~response, ~predictor, ~from, ~from_side, ~to, ~to_side, ~xend_offset, ~curvature,
  "Fragility to hurricane frequency", "Island area",
  "area", "right", "frequency", "bottom", 0.00, 0.00,
  "Fragility to hurricane frequency", "Elevation",
  "elevation", "right", "frequency", "bottom", 0.00, 0.00,
  "Fragility to hurricane frequency", "Isolation distance",
  "isolation", "right", "frequency", "bottom", 0.00, 0.00,
  "Fragility to hurricane frequency", "Global network integration",
  "global", "top", "frequency", "bottom", 0.00, 0.00,
  "Fragility to hurricane frequency", "Pollination topological dominance",
  "dominance", "top", "frequency", "bottom", 0.00, 0.00,
  "Fragility to hurricane frequency", "Multilayer participation coefficient",
  "mpc", "left", "frequency", "bottom", 0.00, 0.00,
  "Fragility to hurricane frequency", "Pollination layer degree (z-score)",
  "pol_degree", "left", "frequency", "bottom", 0.00, 0.00,
  "Fragility to hurricane intensity", "Island area",
  "area", "right", "intensity", "bottom", 0.00, 0.00,
  "Fragility to hurricane intensity", "Isolation distance",
  "isolation", "right", "intensity", "bottom", 0.00, 0.00,
  "Fragility to hurricane intensity", "Continental islands",
  "continental", "right", "intensity", "bottom", 0.00, 0.00,
  "Fragility to hurricane intensity", "Oceanic islands",
  "oceanic", "right", "intensity", "bottom", 0.00, 0.00,
  "Fragility to hurricane intensity", "Global network integration",
  "global", "top", "intensity", "bottom", 0.00, 0.00,
  "Fragility to hurricane intensity", "Pollination layer degree (z-score)",
  "pol_degree", "left", "intensity", "bottom", 0.00, 0.00,
  "Pollination topological dominance", "Isolation distance",
  "isolation", "right", "dominance", "top", -0.54, 0.00,
  "Pollination layer degree (z-score)", "Oceanic islands",
  "oceanic", "right", "pol_degree", "left", 0.00, 0.00,
  "Multilayer participation coefficient", "Body mass",
  "body_mass", "left", "mpc", "left", 0.00, -0.55,
  "Seed-dispersal layer degree (z-score)", "Body mass",
  "body_mass", "left", "sd_degree", "left", 0.00, -0.55
)

# Join the layout with the coefficients, keep the significant paths and place
# their ends on the card edges
arrow_paths <- path_layout %>%
  inner_join(sem_coef %>%
               select(response, predictor, estimate),
             by=c("response", "predictor")) %>%
  mutate(
    x=map2_dbl(from, from_side, ~card_edge(.x, .y)[["x"]]),
    y=map2_dbl(from, from_side, ~card_edge(.x, .y)[["y"]]),
    xend=map2_dbl(to, to_side, ~card_edge(.x, .y)[["x"]]) + xend_offset,
    yend=map2_dbl(to, to_side, ~card_edge(.x, .y)[["y"]]),
    # Blue for positive and red for negative coefficients. Width grows with
    # the absolute coefficient.
    color=if_else(estimate > 0, pal_sem$positive, pal_sem$negative),
    linewidth=scales::rescale(abs(estimate), to=path_width_range)
  )

arrow_paths_straight <- arrow_paths %>%
  filter(curvature == 0)

arrow_paths_curved <- arrow_paths %>%
  filter(curvature != 0)

#--------------------------------------------------------------#
#-Draw the path diagram----------------------------------------#
#--------------------------------------------------------------#

# Panels and the island-type box are drawn first, the arrows next and the
# cards last, so that every arrow starts and ends underneath a card
p_sem <- ggplot() +
  ggforce::geom_shape(data=rect_shape(panel_boxes),
                      aes(x=x, y=y, group=id, fill=fill, color=color,
                          alpha=alpha, linewidth=linewidth),
                      radius=unit(0.16, "cm"),
                      show.legend=FALSE) +
  ggforce::geom_shape(data=rect_shape(island_type_box),
                      aes(x=x, y=y, group=id, fill=fill, color=color,
                          alpha=alpha, linewidth=linewidth),
                      radius=unit(0.10, "cm"),
                      linetype="dashed",
                      show.legend=FALSE)

p_sem <- p_sem +
  geom_segment(data=arrow_paths_straight,
               aes(x=x, y=y, xend=xend, yend=yend,
                   color=color, linewidth=linewidth),
               lineend="round",
               arrow=arrow(length=unit(0.17, "cm"), type="closed"),
               show.legend=FALSE)

# geom_curve takes one curvature per layer, so the curved paths are added one
# at a time
for (i in seq_len(nrow(arrow_paths_curved))) {

  p_sem <- p_sem +
    geom_curve(data=arrow_paths_curved[i, ],
               aes(x=x, y=y, xend=xend, yend=yend),
               curvature=arrow_paths_curved$curvature[i],
               color=arrow_paths_curved$color[i],
               linewidth=arrow_paths_curved$linewidth[i],
               lineend="round",
               arrow=arrow(length=unit(0.17, "cm"), type="closed"))
}

p_sem <- p_sem +
  ggforce::geom_shape(data=rect_shape(cards),
                      aes(x=x, y=y, group=id, fill=fill, color=color,
                          alpha=alpha, linewidth=linewidth),
                      radius=unit(0.07, "cm"),
                      show.legend=FALSE) +
  geom_text(data=cards,
            aes(x=label_x, y=label_y, label=label, size=size),
            color=pal_sem$ink, lineheight=0.95,
            show.legend=FALSE) +
  geom_text(data=side_labels,
            aes(x=x, y=y, label=label, angle=angle, color=color, size=size),
            show.legend=FALSE) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_alpha_identity() +
  scale_linewidth_identity() +
  scale_size_identity() +
  coord_fixed(xlim=c(0, 15.65), ylim=c(-0.45, 10.05), clip="off") +
  theme_void() +
  theme(
    plot.background=element_rect(fill="white", color=NA),
    plot.margin=margin(10, 24, 10, 24)
  )

#--------------------------------------------------------------#
#-Place the card icons-----------------------------------------#
#--------------------------------------------------------------#

# Each icon sits beside the label of its card, within the box given in plot
# coordinates. Island icons go to the right of their labels, role icons to
# the left, fragility icons below and topology icons above.
icons <- tribble(
  ~file,                         ~xmin, ~xmax, ~ymin, ~ymax,
  "area.png",                     1.95,  3.22,  6.88,  7.65,
  "elevation.png",                2.05,  3.25,  5.58,  6.28,
  "distance.png",                 1.95,  3.35,  4.18,  4.85,
  "continental.png",              2.02,  3.10,  2.76,  3.24,
  "oceanic.png",                  2.02,  3.10,  1.96,  2.44,
  "hurr_freq.png",                5.45,  6.55,  8.32,  9.10,
  "hurr_int.png",                 8.40,  9.55,  8.28,  9.10,
  "global_integration.png",       5.15,  6.85,  0.95,  1.80,
  "pollination.png",              8.15,  9.55,  1.00,  1.86,
  "pollination_zscore.png",      11.70, 12.90,  6.86,  7.70,
  "dispersal_zscore.png",        11.70, 12.90,  5.16,  6.00,
  "multilayer_participation.png", 11.70, 12.95,  3.48,  4.28,
  "mass.png",                    11.70, 12.95,  1.78,  2.62
)

for (i in seq_len(nrow(icons))) {
  p_sem <- add_png_icon(p_sem, icons$file[i],
                        icons$xmin[i], icons$xmax[i], icons$ymin[i], icons$ymax[i])
}

#--------------------------------------------------------------#
#-Save Figure 4------------------------------------------------#
#--------------------------------------------------------------#

save_figure(p_sem, "figure_4", width=15, height=10)
print(p_sem)
