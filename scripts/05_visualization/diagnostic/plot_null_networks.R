#--------------------------------------------------------------#
#-Empirical and null multilayer networks-----------------------#
#--------------------------------------------------------------#

# Diagnostic figures, not in the manuscript. Draw the empirical network of
# each island in the center of a page, surrounded by four of its null
# networks, to check by eye what each null model rewires. Species keep the
# positions of the empirical panel in the null panels. Bird colors show the
# role of each bird in the network of its panel, so violet marks the birds
# that remain dual-role after the decoupling. Layout and colors follow the
# network insets of Figure 2.
#
# Inputs:  src/multi_pol/ and src/multi_sd/
#          src/null_<model>/
# Outputs: scripts/05_visualization/diagnostic/
#          null_<model>_networks_<island>.png, .pdf and .svg

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

source(here("src", "island_labels.R"))

#--------------------------------------------------------------#
#-Null models, palette and panel layout------------------------#
#--------------------------------------------------------------#

# Null models to draw, with files in src/null_X/ ending in _null_X_<id>.
# null_ids are the four null networks (out of 100) drawn for each island.
null_models <- c("I", "I-A", "II", "II-B")
null_ids <- c(1, 2, 3, 4)

# Short description of each null model for the figure subtitle
null_descriptions <- c(
  "I"="all connectors retain a single layer; the layers are separated",
  "I-A"="half of the connectors retain a single layer",
  "II"="interactions of dual-role birds are rewired until the layers are separated",
  "II-B"="interactions of dual-role birds are rewired until they retain half of them"
)
output_dir <- here("scripts", "05_visualization", "diagnostic")
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

# Colors of Figure 2: green for the pollination layer, rust for the
# seed-dispersal layer and violet for the dual-role birds
pal <- list(
  plant_pol   = "#1f7a5c",
  plant_sd    = "#b9551f",
  edge_pol    = "#2f9d75",
  edge_sd     = "#d8843a",
  animal_pol  = "#7fc8a9",
  animal_sd   = "#f0b18a",
  animal_both = "#7a3b8f",
  ocean       = "#eaf3f4",
  ink         = "#243238"
)

edge_cols <- c("Pollination"=pal$edge_pol,
               "Seed dispersal"=pal$edge_sd)
node_cols <- c("Plant_Pol"=pal$plant_pol,
               "Plant_SD"=pal$plant_sd,
               "Animal_Pollinator"=pal$animal_pol,
               "Animal_Disperser"=pal$animal_sd,
               "Animal_Both"=pal$animal_both)
node_labels <- c("Pollination plant", "Seed-dispersal plant", "Pollinator",
                 "Seed disperser", "Bird in both layers")

# Folder and file suffix of the five networks of each null model. The
# empirical network comes first because its panel sets the species positions
# of the other four.
network_configs <- lapply(setNames(null_models, null_models), function(null_model) {
  data.table(
    folder=c("", rep(paste0("null_", null_model), 4)),
    suffix=c("", paste0("_null_", null_model, "_", sprintf("%04d", null_ids))),
    label=c("Empirical", paste0("Null model ", null_model, " - ", sprintf("%04d", null_ids)))
  )
})

# A-D are the four null panels and E the empirical panel. The empty cells
# above and below E place the empirical network in the center.
panel_design <- "
A#C
AEC
BED
B#D
"

#--------------------------------------------------------------#
#-Build the five panels for each null model and island---------#
#--------------------------------------------------------------#

# One figure per null model and island
for (null_model in null_models) for (net_id in names(island_names)) {

  network_config <- network_configs[[null_model]]
  panel_list <- vector("list", nrow(network_config))

  for (i in seq_len(nrow(network_config))) {

    # Binary incidence matrices of the two layers (plants in rows, birds in
    # columns), with the same rows and columns
    net_dir <- file.path(here("src"), network_config$folder[i])
    net_file <- paste0(net_id, network_config$suffix[i], ".csv")
    pol_df <- fread(file.path(net_dir, "multi_pol", paste0("B_", net_file)))
    sd_df <- fread(file.path(net_dir, "multi_sd", paste0("S_", net_file)))
    stopifnot(identical(pol_df$plant, sd_df$plant),
              identical(names(pol_df), names(sd_df)))

    pol_mat <- as.matrix(pol_df[, -1, with=FALSE])
    sd_mat <- as.matrix(sd_df[, -1, with=FALSE])
    rownames(pol_mat) <- pol_df$plant
    rownames(sd_mat) <- sd_df$plant

    # Every positive entry is one interaction. Plants appear once per layer,
    # so a plant with both flower visits and fruit consumption is drawn in
    # both plant columns.
    pol_links <- which(pol_mat > 0, arr.ind=TRUE)
    sd_links <- which(sd_mat > 0, arr.ind=TRUE)
    edges <- rbindlist(list(
      data.table(from=paste0(rownames(pol_mat)[pol_links[, 1]], "_Pol"),
                 to=colnames(pol_mat)[pol_links[, 2]], layer="Pollination"),
      data.table(from=paste0(rownames(sd_mat)[sd_links[, 1]], "_SD"),
                 to=colnames(sd_mat)[sd_links[, 2]], layer="Seed dispersal")
    ))

    # Bird roles are computed for every network. Violet marks the dual-role
    # birds of the panel, not those of the empirical network.
    birds_pol <- colSums(pol_mat > 0) > 0
    birds_sd <- colSums(sd_mat > 0) > 0
    bird_types <- fifelse(birds_pol & birds_sd, "Animal_Both",
                         fifelse(birds_pol, "Animal_Pollinator", "Animal_Disperser"))
    # Proportion of interactions of dual-role birds over both layers, shown
    # in the panel subtitle. Every null model lowers it, to zero when the
    # layers are separated.
    dual_birds <- birds_pol & birds_sd
    n_dual_interactions <- sum(pol_mat[, dual_birds] > 0) + sum(sd_mat[, dual_birds] > 0)
    n_interactions <- sum(pol_mat > 0) + sum(sd_mat > 0)

    plants_pol <- rownames(pol_mat)[rowSums(pol_mat > 0) > 0]
    plants_sd <- rownames(sd_mat)[rowSums(sd_mat > 0) > 0]
    nodes <- rbindlist(list(
      data.table(id=paste0(plants_pol, "_Pol"), col=0, type="Plant_Pol"),
      data.table(id=colnames(pol_mat)[birds_pol | birds_sd], col=1,
                 type=bird_types[birds_pol | birds_sd]),
      data.table(id=paste0(plants_sd, "_SD"), col=2, type="Plant_SD")
    ))

    # Order of the species as in Figure 2: dual-role birds of the empirical
    # network first, then pollinators and seed dispersers. Plants are ordered
    # by name within each layer. The null panels reuse these positions.
    if (i == 1) {
      nodes[, sort_key := fifelse(col == 1,
        match(type, c("Animal_Both", "Animal_Pollinator", "Animal_Disperser")),
        as.integer(factor(id)))]
      nodes[, y := seq(1, 0, length.out=.N)[rank(sort_key, ties.method="first")],
            by=col]
      layout_nodes <- nodes[, .(id, x=col, y)]
    }

    stopifnot(setequal(nodes$id, layout_nodes$id))
    nodes <- merge(nodes[, .(id, type)], layout_nodes, by="id", sort=FALSE)
    edges <- merge(edges, layout_nodes, by.x="from", by.y="id", sort=FALSE)
    edges <- merge(edges, layout_nodes, by.x="to", by.y="id", sort=FALSE,
                   suffixes=c("", "end"))

    # Translucent links and white-outlined species markers, as in Figure 2.
    # The legend keeps all roles even when a null network lacks one, such as
    # the dual-role birds in null model I.
    panel_list[[i]] <- ggplot() +
      geom_segment(data=edges, aes(x=x, y=y, xend=xend, yend=yend, color=layer),
                   alpha=0.5, linewidth=0.25,
                   show.legend=c(colour=TRUE, fill=FALSE)) +
      geom_point(data=nodes, aes(x=x, y=y, fill=type),
                 shape=21, color="white", size=1.8, stroke=0.3,
                 show.legend=c(colour=FALSE, fill=TRUE)) +
      scale_color_manual(values=edge_cols, limits=names(edge_cols),
                         labels=c("Pollination links", "Seed-dispersal links")) +
      scale_fill_manual(values=node_cols, limits=names(node_cols),
                        labels=node_labels) +
      coord_cartesian(xlim=c(-0.08, 2.08), ylim=c(-0.03, 1.03), expand=FALSE) +
      labs(title=network_config$label[i],
           subtitle=sprintf("Interactions of dual-role birds: %.2f (%d of %d)",
                            n_dual_interactions / n_interactions,
                            n_dual_interactions, n_interactions),
           color=NULL, fill=NULL) +
      guides(color=guide_legend(order=1, nrow=1,
                               override.aes=list(alpha=1, linewidth=0.8)),
             fill=guide_legend(order=2, nrow=1, override.aes=list(size=3))) +
      theme_void() +
      theme(
        plot.title=element_text(size=15, face="bold", color=pal$ink, hjust=0.5),
        plot.subtitle=element_text(size=11, color=pal$ink, hjust=0.5,
                                   margin=margin(b=10)),
        plot.background=element_rect(fill=NA, color=NA),
        panel.background=element_rect(fill=NA, color=NA),
        plot.margin=margin(10, 14, 10, 14)
      )
  }

  #--------------------------------------------------------------#
  #-Arrange panels and save--------------------------------------#
  #--------------------------------------------------------------#

  p1 <- wrap_plots(A=panel_list[[2]], B=panel_list[[3]],
                   C=panel_list[[4]], D=panel_list[[5]], E=panel_list[[1]],
                   design=panel_design, guides="collect") +
    plot_annotation(
      title=island_names[[net_id]],
      subtitle=paste0("Empirical network and four null model ", null_model,
                      " counterparts (", null_descriptions[[null_model]], ")"),
      caption=paste("Species positions are fixed to the empirical network.",
                    "Bird colors indicate their role in each panel."),
      theme=theme(
        plot.title=element_text(size=23, face="bold", color=pal$ink, hjust=0.5),
        plot.subtitle=element_text(size=13, color=pal$ink, hjust=0.5),
        plot.caption=element_text(size=11, color=pal$ink, hjust=0.5),
        plot.background=element_rect(fill=pal$ocean, color=NA),
        plot.margin=margin(18, 18, 14, 18)
      )
    ) &
    theme(legend.position="bottom", legend.box="vertical",
          legend.text=element_text(size=11, color=pal$ink),
          legend.background=element_rect(fill=NA, color=NA))

  # PNG for quick inspection, PDF and SVG for zooming into the links
  for (format in c("png", "pdf", "svg")) {
    ggsave(file.path(output_dir, paste0("null_", null_model, "_networks_", net_id, ".", format)),
           p1, width=16, height=11, dpi=300, bg=pal$ocean)
  }
}
