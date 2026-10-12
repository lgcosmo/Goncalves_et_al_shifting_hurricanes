#--------------------------------------------------------------#
#-Figure 2: Map of the islands and their multilayer networks---#
#--------------------------------------------------------------#

# Figure 2 maps the 11 islands of the Lesser Antilles and draws, beside each
# island, its multilayer bird-pollination and seed-dispersal network. Each
# network shows the pollination layer on the left, the birds in the centre and
# the seed-dispersal layer on the right, so that the dual-role birds are the
# species joining the two layers. An inset in the open sea plots pollination
# topological dominance, the second PCA axis of network topology, against
# distance to the continent, with island type given by point shape.
#
# Inputs:  src/multi_pol/B_*.csv and src/multi_sd/S_*.csv, the incidence
#          matrices of each island, built by build_networks_multilayer.R
#          data/main_summary/island_data.csv, written by build_island_data.R,
#          with the isolation, island type and PCA scores of each island
#          src/shape_islands/*.shp, the shapefiles of the islands
# Outputs: figures/figure_2.png, .pdf and .svg
# Helpers: island_names in src/island_labels.R and save_figure in
#          src/figure_helpers.R

library(here)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)
library(sf)
library(ggnetwork)
library(network)
library(grid)

source(here("src", "island_labels.R"))
source(here("src", "figure_helpers.R"))

#--------------------------------------------------------------#
#-Palette of Fig. 2--------------------------------------------#
#--------------------------------------------------------------#

# Green marks the pollination layer, rust the seed-dispersal layer and violet
# the dual-role birds; plants take the darker shade of their layer and
# single-role birds the lighter one. Connector colors the lines that join each
# island to its panel, and ink the text and outlines on the muted map.
pal <- list(
  plant_pol    = "#1f7a5c",
  plant_sd     = "#b9551f",
  edge_pol     = "#2f9d75",
  edge_sd      = "#d8843a",
  animal_pol   = "#7fc8a9",
  animal_sd    = "#f0b18a",
  animal_both  = "#7a3b8f",
  ocean        = "#eaf3f4",
  land_context = "#cfd6d0",
  island_fill  = "#f4ead9",
  island_line  = "#8f7668",
  connector    = "#6f8792",
  ink          = "#243238"
)

#--------------------------------------------------------------#
#-Island order and panel anchors of Fig. 2---------------------#
#--------------------------------------------------------------#

# Anchor (lon_target, lat_target) of each network panel, the name of the
# island in the shapefiles, and the latitude shift that separates the drawn
# island from a close neighbour. The panels form two columns, one on each
# side of the island arc, each ordered from north to south. The longitudes
# keep the line from each island to its panel short, so that it reads as a
# callout and not as an interaction. Trinidad moves a little south, away from
# Tobago, and St. Kitts and Nevis move apart, because the enlarged islands sit
# too close to these neighbours.
island_config <- tribble(
  ~network,        ~shapefile_name,    ~lon_target, ~lat_target, ~lat_shift,
  "11_saba",       "Saba",             -68.0,        21.2,         0.00,
  "10_eustatios",  "Saint Eustatius",  -68.0,        18.3,         0.00,
  "09_stkitts",    "Saint Kitts",      -68.0,        15.4,         0.18,
  "08_nevis",      "Nevis",            -68.0,        12.5,        -0.18,
  "06_mtserrat",   "Montserrat",       -68.0,         9.6,         0.00,
  "07_antigua",    "Antigua",          -56.8,        21.2,         0.00,
  "05_dominica",   "Dominica",         -56.8,        18.4,         0.00,
  "04_stvincent",  "Saint Vincent",    -56.8,        15.6,         0.00,
  "03_grenada",    "Grenada",          -56.8,        12.8,         0.00,
  "02_tobago",     "Tobago",           -56.8,        10.0,         0.00,
  "01_trinidad",   "Trinidad",         -56.8,         7.2,        -0.35
) %>%
  mutate(label=island_names[network])

#--------------------------------------------------------------#
#-Network panels of Fig. 2, one per island---------------------#
#--------------------------------------------------------------#

# Draw the multilayer network of one island as a panel of the map. Plants of
# the pollination layer form the left column, birds the central column and
# plants of the seed-dispersal layer the right column. Each interaction is a
# line from a plant to a bird, colored by its layer, so that the dual-role
# birds are the only species joined to both plant columns.
# Arguments: net_id is the island code of the matrix file names, for example
# "01_trinidad". Returns a ggplot grob to place on the map with
# annotation_custom.
build_network_panel <- function(net_id) {

  pol_mat <- read.csv(here("src", "multi_pol", paste0("B_", net_id, ".csv")),
                      check.names=FALSE, row.names=1)
  sd_mat  <- read.csv(here("src", "multi_sd", paste0("S_", net_id, ".csv")),
                      check.names=FALSE, row.names=1)

  # List the interactions of one incidence matrix, plants in rows and birds in
  # columns. Any positive entry draws a link.
  # Arguments: mat is the matrix; layer names the layer in the output. Returns
  # a data frame with one row per interaction.
  to_edges <- function(mat, layer) {
    mat %>%
      as.matrix() %>%
      as.table() %>%
      as.data.frame() %>%
      setNames(c("plant", "animal", "interaction")) %>%
      filter(interaction > 0) %>%
      mutate(layer=layer,
             plant=as.character(plant),
             animal=as.character(animal))
  }

  edges_pol <- to_edges(pol_mat, "Pollination")
  edges_sd  <- to_edges(sd_mat,  "Seed dispersal")

  # Dual-role birds sit at the top of the central column, above the single-role
  # pollinators and then the single-role dispersers, because they are the
  # species that join the two layers.
  animals_pol <- unique(edges_pol$animal)
  animals_sd  <- unique(edges_sd$animal)
  all_animals <- union(animals_pol, animals_sd)

  animal_role <- tibble(animal=all_animals) %>%
    mutate(
      role = case_when(
        animal %in% animals_pol & animal %in% animals_sd ~ "both",
        animal %in% animals_pol                          ~ "pollinator",
        TRUE                                             ~ "disperser"
      ),
      role_priority = case_when(
        role == "both"       ~ 1,
        role == "pollinator" ~ 2,
        TRUE                 ~ 3
      )
    )

  # A plant recorded in both layers appears once in each plant column. Within
  # a column, species are spaced evenly from the top to the bottom of the
  # panel, so that islands of different richness fill the same area.
  nodes <- bind_rows(
    tibble(id=paste0(unique(edges_pol$plant), "_Pol"),
           col=0, type="Plant_Pol"),
    tibble(id=paste0(unique(edges_sd$plant), "_SD"),
           col=2, type="Plant_SD"),
    animal_role %>%
      transmute(
        id=animal,
        col=1,
        type=case_when(
          role == "both"       ~ "Animal_Both",
          role == "pollinator" ~ "Animal_Pollinator",
          TRUE                 ~ "Animal_Disperser"
        )
      )
  ) %>%
    group_by(col) %>%
    mutate(
      sort_key = if_else(
        col == 1,
        as.numeric(animal_role$role_priority[match(id, animal_role$animal)]),
        as.numeric(factor(id))
      ),
      y = {
        n <- n()
        seq(1, 0, length.out=n)[rank(sort_key, ties.method="first")]
      }
    ) %>%
    ungroup()

  graph_edges <- bind_rows(
    edges_pol %>% transmute(from=paste0(plant, "_Pol"), to=animal, layer),
    edges_sd  %>% transmute(from=paste0(plant, "_SD"),  to=animal, layer)
  )

  net <- network(graph_edges, directed=FALSE)
  net_names <- network.vertex.names(net)
  nodes_ord <- nodes[match(net_names, nodes$id), ]

  set.vertex.attribute(net, "type", nodes_ord$type)
  layout_mat <- as.matrix(nodes_ord %>% select(col, y))

  gn <- ggnetwork(net, layout=layout_mat)

  edge_cols <- c("Pollination"=pal$edge_pol,
                 "Seed dispersal"=pal$edge_sd)
  node_cols <- c(
    "Plant_Pol"         = pal$plant_pol,
    "Plant_SD"          = pal$plant_sd,
    "Animal_Pollinator" = pal$animal_pol,
    "Animal_Disperser"  = pal$animal_sd,
    "Animal_Both"       = pal$animal_both
  )

  # A transparent background lets the panel sit directly on the sea of the map
  p <- ggplot(gn, aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_edges(aes(color=layer), alpha=0.5, linewidth=0.25) +
    geom_nodes(aes(fill=type), shape=21, color="white", size=1.8, stroke=0.3) +
    scale_color_manual(values=edge_cols, guide="none") +
    scale_fill_manual(values=node_cols, guide="none") +
    theme_void() +
    theme(
      plot.background=element_rect(fill=NA, color=NA),
      panel.background=element_rect(fill=NA, color=NA),
      plot.margin=margin(2, 2, 2, 2)
    )

  return(ggplotGrob(p))
}

#--------------------------------------------------------------#
#-Inset of Fig. 2: pollination dominance and isolation---------#
#--------------------------------------------------------------#

# Draw the inset, which plots pollination topological dominance against
# distance to the continent, one point per island, with a linear fit. The
# dominance is the second PCA axis of network topology, oriented so that
# larger values mean a multilayer structure dominated by the pollination
# layer. Point shape gives island type; point color repeats isolation, so the
# inset has no color legend. Returns a ggplot grob to place in the open sea
# between the northern network panels.
build_topology_inset <- function() {

  island_data <- read.csv(here("data", "main_summary", "island_data.csv"),
                          stringsAsFactors=FALSE)

  # The inset stays compact, with no title, no color legend and thin axis
  # lines, so that it reads as a secondary panel inside the map.
  p <- ggplot(data=island_data,
              aes(x=isolation, y=pollination_dominance)) +
    geom_point(aes(shape=type, color=isolation), size=2.0, alpha=0.9) +
    geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.45) +
    scale_color_viridis_c(option="viridis", guide="none") +
    scale_shape_manual(
      values=c("Continental"=16, "Oceanic"=17),
      name="Island type:"
    ) +
    labs(x="Distance to the continent",
         y="Pollination topological dominance") +
    theme_pubr(base_size=14) +
    theme(
      plot.background=element_rect(fill="white", color=NA),
      panel.background=element_rect(fill="white", color=NA),
      axis.title.x=element_text(color=pal$ink, size=11),
      axis.title.y=element_text(color=pal$ink, size=9.5,
                                margin=margin(r=3)),
      axis.text=element_text(color=pal$ink, size=10),
      axis.line=element_line(linewidth=0.25, color=pal$ink),
      axis.ticks=element_line(linewidth=0.25, color=pal$ink),
      legend.position="top",
      legend.direction="horizontal",
      legend.justification="center",
      legend.background=element_rect(fill="white", color=NA),
      legend.title=element_text(color=pal$ink, size=10),
      legend.text=element_text(color=pal$ink, size=10),
      legend.key.size=unit(0.40, "cm"),
      legend.margin=margin(0, 0, 0, 0),
      plot.margin=margin(5, 6, 5, 6)
    )

  return(ggplotGrob(p))
}

#--------------------------------------------------------------#
#-Map of Fig. 2: island shapefiles-----------------------------#
#--------------------------------------------------------------#

shape_dir <- here("src", "shape_islands")

# Spherical geometry (S2) is switched off because some of the Caribbean
# shapefiles hold duplicate vertices. They draw correctly, but they stop the
# sf operations that repair, scale and shift the islands below.
sf_use_s2(FALSE)

# The Caribbean shapefile names its islands in the column caribbean1 and sets
# the coordinate reference system of the map
carib_sf <- st_read(file.path(shape_dir, "select_caribbean_islands.shp"), quiet=TRUE)
target_crs <- st_crs(carib_sf)

# Read a single-island shapefile, bring it to the reference system of the map,
# repair its geometry and name the island.
# Arguments: file is the shapefile; name is the island name. Returns an sf
# object with the columns Name and geometry.
read_island <- function(file, name) {
  island <- st_read(file.path(shape_dir, file), quiet=TRUE) %>%
    st_transform(target_crs) %>%
    st_make_valid()
  island$Name <- name
  select(island, Name, geometry)
}

all_islands_sf <- rbind(
  carib_sf %>% st_make_valid() %>% select(Name=caribbean1, geometry) %>% filter(!is.na(Name)),
  read_island("Trinidad.shp", "Trinidad"),
  read_island("Tobago.shp", "Tobago")
) %>%
  mutate(is_study=Name %in% island_config$shapefile_name)

# The 11 study islands are enlarged by a quarter about their centroids, so
# that they stay legible at the scale of the whole map, and shifted by the
# latitude given in island_config. The other islands stay as drawn.
scale_study <- 1.25

geom_list <- lapply(seq_len(nrow(all_islands_sf)), function(i) {
  g <- st_geometry(all_islands_sf[i, ])
  if (!all_islands_sf$is_study[i]) return(g)

  ctr <- st_centroid(g)
  lat_shift <- island_config$lat_shift[match(all_islands_sf$Name[i], island_config$shapefile_name)]
  g <- (g - ctr) * scale_study + ctr + c(0, lat_shift)
  st_crs(g) <- target_crs
  g
})

st_geometry(all_islands_sf) <- do.call(c, geom_list)
st_crs(all_islands_sf) <- target_crs

# The centroid of each study island, as drawn after scaling and shifting,
# anchors the line to its network panel
study_sf <- all_islands_sf %>% filter(is_study)
centroids <- st_centroid(study_sf) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(shapefile_name=study_sf$Name) %>%
  rename(lon_cent=X, lat_cent=Y)

layout_final <- island_config %>%
  left_join(centroids, by="shapefile_name")

#--------------------------------------------------------------#
#-Draw the network panels and the inset------------------------#
#--------------------------------------------------------------#

# The eleven network panels and the inset are drawn once, before the map, so
# that the loop that places them only positions grobs
panel_list <- lapply(setNames(layout_final$network, layout_final$network), build_network_panel)
topology_grob <- build_topology_inset()

#--------------------------------------------------------------#
#-Base map and panel positions---------------------------------#
#--------------------------------------------------------------#

# The map spans the Lesser Antilles from Trinidad to Saba, with open sea on
# both sides for the network panels and at the top for the topology inset
xlim <- c(-70.3, -54.4)
ylim <- c(5.7, 23.3)

# Width, height and inward shift of each network panel, in degrees of the map
inset_w <- 2.60
inset_h <- 2.20
inset_shift <- 0.70

# Each panel is pulled toward its island, and the line from the island ends at
# the inner edge of the panel, so that the callout stays short
layout_final <- layout_final %>%
  mutate(
    inset_side = if_else(lon_target < lon_cent, "left", "right"),
    inset_center_x = if_else(inset_side == "left",
                             lon_target + inset_shift,
                             lon_target - inset_shift),
    connector_x = if_else(inset_side == "left",
                          inset_center_x + inset_w / 2,
                          inset_center_x - inset_w / 2),
    connector_y = lat_target
  )

p_map <- ggplot() +
  geom_sf(data=all_islands_sf %>% filter(!is_study),
          fill=pal$land_context, color="white", linewidth=0.18, alpha=0.55) +
  geom_curve(data=layout_final,
             aes(x=lon_cent, y=lat_cent, xend=connector_x, yend=connector_y),
             color=pal$connector, linewidth=0.30, alpha=0.60,
             curvature=0.08, lineend="round",
             arrow=arrow(length=unit(0.11, "cm"), type="closed")) +
  geom_sf(data=all_islands_sf %>% filter(is_study),
          fill=pal$island_fill, color=pal$island_line, linewidth=0.42) +
  coord_sf(xlim=xlim, ylim=ylim, expand=FALSE) +
  theme_void() +
  theme(
    panel.background=element_rect(fill=pal$ocean, color=NA),
    plot.background=element_rect(fill=pal$ocean, color=NA),
    plot.margin=margin(8, 8, 8, 8)
  )

#--------------------------------------------------------------#
#-Place the island names, the network panels and the inset-----#
#--------------------------------------------------------------#

# The island name sits just above its network panel, clear of the topology
# inset in the north
for (i in seq_len(nrow(layout_final))) {

  row <- layout_final[i, ]

  p_map <- p_map +
    annotate("text",
      x=row$inset_center_x, y=row$lat_target + inset_h / 2 + 0.20,
      label=row$label, color=pal$ink, size=3.2, fontface="bold",
      hjust=0.5, vjust=0
    ) +
    annotation_custom(
      grob=panel_list[[row$network]],
      xmin=row$inset_center_x - inset_w / 2,
      xmax=row$inset_center_x + inset_w / 2,
      ymin=row$lat_target - inset_h / 2,
      ymax=row$lat_target + inset_h / 2
    )
}

# The inset fills the open sea at the top of the map, between the two columns
# of network panels and clear of the Saba and Antigua panels. Its 5:4 bounds
# are close to square, so that the regression is not flattened.
p_map <- p_map +
  annotation_custom(
    grob=topology_grob,
    xmin=-65.00, xmax=-60.00,
    ymin=19.10, ymax=23.10
  )

#--------------------------------------------------------------#
#-Legend of the network panels---------------------------------#
#--------------------------------------------------------------#

# The legend sits in the sea south of the island arc. Green and rust segments
# stand for the pollination and seed-dispersal links, and a violet point for
# the dual-role birds, one row each from the top.
legend_x <- -62.1
legend_y <- 6.25
legend_gap <- 0.42

legend_items <- tribble(
  ~row, ~kind,     ~color,          ~label,
  2,    "segment", pal$edge_pol,    "Pollination links",
  1,    "segment", pal$edge_sd,     "Seed-dispersal links",
  0,    "point",   pal$animal_both, "Animals in both layers"
) %>%
  mutate(y=legend_y + row * legend_gap)

p_map <- p_map +
  geom_segment(data=filter(legend_items, kind == "segment"),
               aes(x=legend_x - 1.15, xend=legend_x - 0.65, y=y, yend=y, color=color),
               linewidth=1.0, alpha=0.85) +
  geom_point(data=filter(legend_items, kind == "point"),
             aes(x=legend_x - 0.90, y=y, fill=color),
             color="white", shape=21, size=3.0, stroke=0.35) +
  geom_text(data=legend_items,
            aes(x=legend_x - 0.48, y=y, label=label),
            color=pal$ink, size=2.8, hjust=0, vjust=0.5) +
  scale_color_identity() +
  scale_fill_identity()

#--------------------------------------------------------------#
#-Save Fig. 2--------------------------------------------------#
#--------------------------------------------------------------#

# The sea color fills the margins, so the saved figure matches the map
save_figure(p_map, "figure_2", width=10, height=12, bg=pal$ocean)
print(p_map)
