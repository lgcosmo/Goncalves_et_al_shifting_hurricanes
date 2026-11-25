#--------------------------------------------------------------#
#-Figure 4: pSEM path diagram----------------------------------#
#--------------------------------------------------------------#

# Draw Figure 4, the path diagram of the pSEM. It shows the significant
# directed paths from island biogeography, network topology and species
# roles to the persistence of bird species under hurricane frequency and
# under hurricane intensity. The diagram reads from bottom to top: islands
# at the bottom, network topology and species roles in the middle and the
# two persistence responses at the top. Blue arrows are positive path
# coefficients and red arrows negative ones, so a blue arrow into a response
# means greater persistence. Arrow width and arrowhead size increase with
# the absolute coefficient. Correlated errors are not drawn. Panels and
# cards are placed in plot coordinates, with equal units on both axes.
#
# Inputs:  tables/sem_coefficients_clean.csv
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

# Keep the significant directed paths. The row that tests island type as a
# whole has no estimate and is dropped.
sem_coef <- read_csv(here("tables", "sem_coefficients_clean.csv"),
                     show_col_types=FALSE) %>%
  filter(path_type == "Directed path",
         !is.na(estimate),
         !is.na(significance),
         significance != "")

#--------------------------------------------------------------#
#-Palette, arrow sizes and drawing helpers---------------------#
#--------------------------------------------------------------#

# One color per panel, in three shades: panel border, panel background and
# cards. Blue and red are the sign of the paths.
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

# Ranges for rescaling the absolute coefficients: arrow width, and length
# and half-width of the arrowheads (plot units)
path_width_range <- c(0.70, 2.60)
head_length_range <- c(0.17, 0.27)
head_halfwidth_range <- c(0.065, 0.125)

# Distance between neighboring arrows on the edge of a card, margin without
# arrows next to each corner, and gap between the arrow tip and the card
# (plot units)
slot_spacing <- 0.38
edge_margin <- 0.25
tip_gap <- 0.05

# Expand a table of rectangles (id, xmin, xmax, ymin, ymax and aesthetics)
# into their four corners, as geom_shape needs to draw rounded boxes.
# Returns a tibble with four rows per rectangle.
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

# Build the triangle of each arrowhead from the tip of the arrow (xtip,
# ytip), its unit vector (ux, uy), head_length and head_halfwidth. Returns a
# tibble with three rows per arrow.
head_shape <- function(df) {

  head_list <- lapply(seq_len(nrow(df)), function(i) {

    row_i <- df[i, ]

    # Middle of the base of the head, one head length behind the tip
    xbase <- row_i$xtip - row_i$head_length * row_i$ux
    ybase <- row_i$ytip - row_i$head_length * row_i$uy

    # The two corners, on either side of the arrow along its normal
    tibble(
      path=row_i$path,
      x=c(row_i$xtip,
          xbase - row_i$head_halfwidth * row_i$uy,
          xbase + row_i$head_halfwidth * row_i$uy),
      y=c(row_i$ytip,
          ybase + row_i$head_halfwidth * row_i$ux,
          ybase - row_i$head_halfwidth * row_i$ux),
      color=row_i$color
    )
  })

  return(bind_rows(head_list))
}

# Add one PNG icon from figures/icons to the plot p, inside the box given by
# xmin, xmax, ymin and ymax. The icon keeps its proportions.
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

# The four panels, from bottom to top: island biogeography, network topology
# (left) and species roles (right), and the two persistence responses. The
# space between the two middle panels is free for the paths from the islands
# to the responses.
panel_boxes <- tribble(
  ~id, ~xmin, ~xmax, ~ymin, ~ymax, ~fill, ~color, ~alpha, ~linewidth,
  "island_panel",     2.225, 12.375, 0.20,  2.65, pal_sem$island_fill,    pal_sem$island,    0.78, 0.55,
  "network_panel",    0.45,   3.55,  3.40,  7.70, pal_sem$network_fill,   pal_sem$network,   0.78, 0.55,
  "trait_panel",     11.05,  15.60,  3.40,  7.70, pal_sem$trait_fill,     pal_sem$trait,     0.78, 0.55,
  "persistence_panel", 4.20, 10.40,  9.25, 11.45, pal_sem$hurricane_fill, pal_sem$hurricane, 0.88, 0.55
)

# One card per variable, with its icon. The label is at the bottom of the
# card and the icon above it, except in the two responses (label at the
# top). Isolation is the first island card and oceanic the last one, each
# below the panel that its path enters.
cards <- tribble(
  ~id, ~group, ~label, ~icon, ~xmin, ~xmax, ~ymin, ~ymax, ~label_side, ~size,
  "isolation",   "island",    "Island\nisolation",                  "distance.png",                  2.475,  4.125, 0.80,  2.30, "bottom", 4.0,
  "area",        "island",    "Island\narea",                       "area.png",                      4.425,  6.075, 0.80,  2.30, "bottom", 4.0,
  "elevation",   "island",    "Island\nelevation",                  "elevation.png",                 6.375,  8.025, 0.80,  2.30, "bottom", 4.0,
  "continental", "island",    "Continental",                        "continental.png",               8.475, 10.125, 0.80,  2.30, "bottom", 4.0,
  "oceanic",     "island",    "Oceanic",                            "oceanic.png",                  10.325, 11.975, 0.80,  2.30, "bottom", 4.0,
  "dominance",   "network",   "Pollination topological\ndominance", "pollination.png",               0.70,   3.30,  3.65,  5.15, "bottom", 4.0,
  "global",      "network",   "Global network\nintegration",        "global_integration.png",        0.70,   3.30,  5.95,  7.45, "bottom", 4.0,
  "pol_degree",  "trait",     "Pollination layer\ndegree (z-score)", "pollination_zscore.png",       11.30,  13.20,  3.65,  5.15, "bottom", 4.0,
  "body_mass",   "trait",     "Body mass",                          "mass.png",                     13.45,  15.35,  3.65,  5.15, "bottom", 4.0,
  "mpc",         "trait",     "Multilayer\nparticipation",          "multilayer_participation.png", 11.30,  13.20,  5.95,  7.45, "bottom", 4.0,
  "sd_degree",   "trait",     "Dispersal layer\ndegree (z-score)",  "dispersal_zscore.png",         13.45,  15.35,  5.95,  7.45, "bottom", 4.0,
  "frequency",   "hurricane", "Persistence under\nhurricane frequency", "hurr_freq.png",            4.40,   7.15,  9.55, 11.15, "top",    4.2,
  "intensity",   "hurricane", "Persistence under\nhurricane intensity", "hurr_int.png",             7.45,  10.20,  9.55, 11.15, "top",    4.2
) %>%
  mutate(
    # Color of the panel, in the shade of the cards
    fill=unlist(pal_sem)[paste0(group, "_card")],
    color=unlist(pal_sem)[group],
    alpha=0.98,
    linewidth=0.50,
    # Label centered in the card
    label_x=(xmin + xmax) / 2,
    label_y=if_else(label_side == "bottom", ymin + 0.36, ymax - 0.38),
    # The icon uses the space left by the label
    icon_xmin=xmin + 0.15,
    icon_xmax=xmax - 0.15,
    icon_ymin=if_else(label_side == "bottom", ymin + 0.68, ymin + 0.14),
    icon_ymax=if_else(label_side == "bottom", ymax - 0.12, ymax - 0.74)
  )

# A dashed box groups the continental and oceanic cards as island type
island_type_box <- tribble(
  ~id, ~xmin, ~xmax, ~ymin, ~ymax, ~fill, ~color, ~alpha, ~linewidth,
  "island_type", 8.325, 12.125, 0.40, 2.45, "#FFFFFF", pal_sem$island, 0.35, 0.38
)

# Names of the three levels, on the outer side of their panels, and the label
# of the island-type box
side_labels <- tribble(
  ~label, ~x, ~y, ~angle, ~color, ~size,
  "Island biogeography",          7.30, -0.18,   0, pal_sem$island, 6.2,
  "Network topology",             0.08,  5.55,  90, "#6CC4D9",      6.2,
  "Species position and traits", 15.97,  5.55, -90, pal_sem$trait,  6.2,
  "Island type",                 10.225, 0.60,   0, pal_sem$ink,    3.4
)

#--------------------------------------------------------------#
#-Arrow layout of the significant paths------------------------#
#--------------------------------------------------------------#

# Side of the predictor card that each path leaves and side of the response
# card that it enters. Paths leave the top of an island card or the inner
# side of a middle card and enter the bottom of their response, so the
# diagram reads upward. Only the paths that are significant in the
# coefficient table are drawn. from_shift and to_shift move one end of a
# path along its edge (plot units): the path from body mass enters the
# multilayer participation card right of its middle, which keeps it steep.
path_layout <- tribble(
  ~response, ~predictor, ~from, ~from_side, ~to, ~to_side, ~from_shift, ~to_shift,
  "Persistence under hurricane frequency", "Island area",
  "area", "top", "frequency", "bottom", 0.00, 0.00,
  "Persistence under hurricane frequency", "Elevation",
  "elevation", "top", "frequency", "bottom", 0.00, 0.00,
  "Persistence under hurricane frequency", "Isolation distance",
  "isolation", "top", "frequency", "bottom", 0.00, 0.00,
  "Persistence under hurricane frequency", "Global network integration",
  "global", "right", "frequency", "bottom", 0.00, 0.00,
  "Persistence under hurricane frequency", "Pollination topological dominance",
  "dominance", "right", "frequency", "bottom", 0.00, 0.00,
  "Persistence under hurricane frequency", "Multilayer participation coefficient",
  "mpc", "left", "frequency", "bottom", 0.00, 0.00,
  "Persistence under hurricane frequency", "Pollination layer degree (z-score)",
  "pol_degree", "left", "frequency", "bottom", 0.00, 0.00,
  "Persistence under hurricane intensity", "Island area",
  "area", "top", "intensity", "bottom", 0.00, 0.00,
  "Persistence under hurricane intensity", "Isolation distance",
  "isolation", "top", "intensity", "bottom", 0.00, 0.00,
  "Persistence under hurricane intensity", "Continental islands",
  "continental", "top", "intensity", "bottom", 0.00, 0.00,
  "Persistence under hurricane intensity", "Oceanic islands",
  "oceanic", "top", "intensity", "bottom", 0.00, 0.00,
  "Persistence under hurricane intensity", "Global network integration",
  "global", "right", "intensity", "bottom", 0.00, 0.00,
  "Persistence under hurricane intensity", "Pollination layer degree (z-score)",
  "pol_degree", "left", "intensity", "bottom", 0.00, 0.00,
  "Pollination topological dominance", "Isolation distance",
  "isolation", "top", "dominance", "bottom", 0.00, 0.00,
  "Pollination layer degree (z-score)", "Oceanic islands",
  "oceanic", "top", "pol_degree", "bottom", 0.00, 0.00,
  "Multilayer participation coefficient", "Body mass",
  "body_mass", "top", "mpc", "bottom", 0.00, 0.55,
  "Seed-dispersal layer degree (z-score)", "Body mass",
  "body_mass", "top", "sd_degree", "bottom", 0.00, 0.00
)

# Join the layout with the coefficients. Only significant paths remain.
arrow_paths <- path_layout %>%
  inner_join(sem_coef %>%
               select(response, predictor, estimate),
             by=c("response", "predictor")) %>%
  mutate(path=row_number())

# Edges and middle of each card
card_geometry <- cards %>%
  transmute(card=id, xmin, xmax, ymin, ymax,
            cx=(xmin + xmax) / 2,
            cy=(ymin + ymax) / 2)

# The two ends of each path, one row per end, with its card and the card at
# the other end
path_ends <- bind_rows(
  arrow_paths %>%
    transmute(path, end="start", card=from, side=from_side, other=to,
              shift=from_shift),
  arrow_paths %>%
    transmute(path, end="end", card=to, side=to_side, other=from,
              shift=to_shift)
) %>%
  left_join(card_geometry, by="card") %>%
  left_join(card_geometry %>%
              transmute(other=card, other_cx=cx, other_cy=cy),
            by="other")

# Each end has its own point on the edge of its card. Ends that share an
# edge are spread around its middle, ordered by the direction of the card at
# the other end, so paths that share a card do not cross next to it. The
# direction is an angle that increases from left to right on a top or bottom
# edge and from bottom to top on a left or right edge.
path_ends <- path_ends %>%
  mutate(
    dx=other_cx - cx,
    dy=other_cy - cy,
    direction=case_when(
      side == "top"    ~ atan2(dx, dy),
      side == "bottom" ~ atan2(dx, -dy),
      side == "right"  ~ atan2(dy, dx),
      side == "left"   ~ atan2(dy, -dx)
    ),
    edge_length=if_else(side %in% c("top", "bottom"), xmax - xmin, ymax - ymin)
  ) %>%
  group_by(card, side) %>%
  mutate(
    n_ends=n(),
    # Smaller spacing when the edge cannot hold all the ends
    spacing=if_else(n_ends > 1,
                    pmin(slot_spacing,
                         (edge_length - 2 * edge_margin) / (n_ends - 1)),
                    0),
    offset=(rank(direction, ties.method="first") - (n_ends + 1) / 2) * spacing +
      shift
  ) %>%
  ungroup() %>%
  mutate(
    x=case_when(side == "left" ~ xmin, side == "right" ~ xmax, TRUE ~ cx + offset),
    y=case_when(side == "bottom" ~ ymin, side == "top" ~ ymax, TRUE ~ cy + offset)
  )

# Build each arrow from its two ends: color from the sign of the coefficient,
# width and head from its absolute value. The tip stops short of the card
# and the line ends inside the head.
arrow_paths <- arrow_paths %>%
  left_join(path_ends %>%
              filter(end == "start") %>%
              select(path, x, y),
            by="path") %>%
  left_join(path_ends %>%
              filter(end == "end") %>%
              select(path, xend=x, yend=y),
            by="path") %>%
  mutate(
    # Unit vector along the path
    path_length=sqrt((xend - x)^2 + (yend - y)^2),
    ux=(xend - x) / path_length,
    uy=(yend - y) / path_length,
    # Blue for positive and red for negative coefficients
    color=if_else(estimate > 0, pal_sem$positive, pal_sem$negative),
    linewidth=scales::rescale(abs(estimate), to=path_width_range),
    head_length=scales::rescale(abs(estimate), to=head_length_range),
    head_halfwidth=scales::rescale(abs(estimate), to=head_halfwidth_range),
    xtip=xend - tip_gap * ux,
    ytip=yend - tip_gap * uy,
    xline=xtip - 0.85 * head_length * ux,
    yline=ytip - 0.85 * head_length * uy
  ) %>%
  # Draw the widest arrows first, so the thin ones stay on top
  arrange(desc(linewidth))

arrow_heads <- head_shape(arrow_paths) %>%
  mutate(path=factor(path, levels=arrow_paths$path))

#--------------------------------------------------------------#
#-Draw the path diagram----------------------------------------#
#--------------------------------------------------------------#

# Panels and the island-type box first, then the arrows and the cards last,
# so every arrow starts underneath its card
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

# Line of each arrow, and its head as a filled triangle with a thin white
# outline
p_sem <- p_sem +
  geom_segment(data=arrow_paths,
               aes(x=x, y=y, xend=xline, yend=yline,
                   color=color, linewidth=linewidth),
               lineend="butt",
               show.legend=FALSE) +
  geom_polygon(data=arrow_heads,
               aes(x=x, y=y, group=path, fill=color),
               color="white", linewidth=0.25,
               show.legend=FALSE)

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
  coord_fixed(xlim=c(-0.15, 16.20), ylim=c(-0.50, 11.60), expand=FALSE,
              clip="off") +
  theme_void() +
  theme(
    plot.background=element_rect(fill="white", color=NA),
    plot.margin=margin(10, 24, 10, 24)
  )

#--------------------------------------------------------------#
#-Place the card icons-----------------------------------------#
#--------------------------------------------------------------#

# Each icon sits in its card, above the label or below it in the two
# responses
for (i in seq_len(nrow(cards))) {
  p_sem <- add_png_icon(p_sem, cards$icon[i],
                        cards$icon_xmin[i], cards$icon_xmax[i],
                        cards$icon_ymin[i], cards$icon_ymax[i])
}

#--------------------------------------------------------------#
#-Save Figure 4------------------------------------------------#
#--------------------------------------------------------------#

save_figure(p_sem, "figure_4", width=15, height=10.9)
print(p_sem)
