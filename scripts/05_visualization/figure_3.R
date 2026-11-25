#--------------------------------------------------------------#
#-Figure 3: hurricane regimes and network extinction-----------#
#--------------------------------------------------------------#

# Draw Figure 3. Panels a and b use Trinidad as an example: the Poisson
# probabilities of 0 to 4 hurricanes in a year at its empirical frequency
# and at an increased frequency (a), and its PMF of hurricane categories 1
# to 5, empirical and after probability moves from categories 1 to 3 to
# categories 4 and 5 (b). Panels c and d show the fraction of species
# extinct on each island under increasing hurricane frequency (lambda from
# 0.01 to 0.80) and under increasing hurricane intensity (shift from 0 to
# 1). The point on each curve is the empirical regime of the island. Colors
# separate continental from oceanic islands.
#
# Inputs:  src/empirical_lambdas.csv and src/empirical_pmfs.csv
#          data/main_summary/island_data.csv
#          data/main_summary/main_summary_metrics_lambda.csv
#          data/main_summary/main_summary_metrics_pmf.csv
#          data/main_summary/main_summary_empirical_interpolated.csv
# Outputs: figures/figure_3.png, .pdf and .svg

library(here)
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)
library(tidyr)
library(viridisLite)

source(here("src", "island_labels.R"))
source(here("src", "figure_helpers.R"))

#--------------------------------------------------------------#
#-Hurricane regimes and extinction curves of the islands-------#
#--------------------------------------------------------------#

# Replace the island codes of one column by the island names, as a factor
# ordered from Trinidad (south) to Saba (north). Returns df.
relabel_islands <- function(df, column) {
  df[[column]] <- factor(island_names[as.character(df[[column]])], levels=island_order)
  return(df)
}

# Fraction of species extinct on each island at each value of lambda (panel
# c) and of the intensity shift (panel d), and the same fraction interpolated
# at the empirical lambda of each island (points of panel c)
results_lambda <- fread(here("data", "main_summary", "main_summary_metrics_lambda.csv")) %>%
  relabel_islands("network")
results_interp <- fread(here("data", "main_summary", "main_summary_empirical_interpolated.csv")) %>%
  relabel_islands("network")
results_pmf    <- fread(here("data", "main_summary", "main_summary_metrics_pmf.csv")) %>%
  relabel_islands("network")

# Empirical hurricane regimes. Lambda is the expected number of hurricanes
# per year and the PMF is the probability of each Saffir-Simpson category.
lambda_estimates <- fread(here("src", "empirical_lambdas.csv")) %>%
  relabel_islands("network") %>%
  rename(island=network)
intensity_pmf <- fread(here("src", "empirical_pmfs.csv")) %>%
  relabel_islands("network") %>%
  rename(island=network)

# Type of each island (continental or oceanic) and its distance to the
# continent, used for the island colors
island_types <- fread(here("data", "main_summary", "island_data.csv")) %>%
  relabel_islands("network") %>%
  transmute(island=network, type, isolation=as.numeric(isolation))

#--------------------------------------------------------------#
#-Colors and theme of the four panels--------------------------#
#--------------------------------------------------------------#

# Trinidad is the example island of panels a and b
example_island <- "Trinidad"

# Dark neutral tone of the axis text, tags and legend
pal_fig3 <- list(ink="#253238")

# Colors separate the two island types and follow the distance to the
# continent within each type. Continental islands are warm: red for Trinidad
# and orange for Tobago. Oceanic islands follow the viridis scale, from dark
# violet (closest) to yellow green (most isolated). The oceanic colors are
# spaced by the rank of isolation, because the northern islands are at
# similar distances, and the scale stops before its yellow end, which is
# hard to read on white.
continental_islands <- island_types %>%
  filter(type == "Continental") %>%
  arrange(isolation) %>%
  pull(island) %>%
  as.character()

oceanic_islands <- island_types %>%
  filter(type == "Oceanic") %>%
  arrange(isolation) %>%
  pull(island) %>%
  as.character()

continental_cols <- c("#B2182B", "#EE8532")
names(continental_cols) <- continental_islands

oceanic_cols <- viridis(length(oceanic_islands), begin=0.10, end=0.92)
names(oceanic_cols) <- oceanic_islands

island_cols <- c(continental_cols, oceanic_cols)[island_order]

# Sort a results table so the continental islands come last and are drawn
# over the oceanic ones. Returns df with the factor draw_order.
order_for_drawing <- function(df) {
  df$draw_order <- factor(as.character(df$network),
                          levels=c(oceanic_islands, continental_islands))
  df <- df[order(df$draw_order), ]
  return(df)
}

results_lambda <- order_for_drawing(results_lambda)
results_interp <- order_for_drawing(results_interp)
results_pmf    <- order_for_drawing(results_pmf)

# Colors of the two regimes in panels a and b: blue for the empirical regime
# and the red of Trinidad for the escalated one
pal_regime <- c(
  "Empirical" = "#2166AC",
  "Escalated" = unname(island_cols[example_island])
)

# Color scale of the islands in panels c and d. Each panel shows the legend
# of one island type (continental in c, oceanic in d), and patchwork
# collects the two legends below the figure.
scale_color_island <- function(islands, legend_title) {
  scale_color_manual(
    values=island_cols,
    breaks=islands,
    name=legend_title,
    guide=guide_legend(
      nrow=2,
      byrow=TRUE,
      override.aes=list(linewidth=1.6, alpha=1)
    )
  )
}

# Fill scale for the points of panels c and d, same palette and no legend
scale_fill_island <- function() {
  scale_fill_manual(
    values=island_cols,
    breaks=island_order,
    guide="none"
  )
}

# Theme of the four panels, based on theme_pubr
theme_fig3 <- function(base_size=14) {
  theme_pubr(base_size=base_size) +
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.title=element_text(color=pal_fig3$ink, size=base_size),
      axis.text=element_text(color=pal_fig3$ink, size=base_size),
      legend.title=element_text(color=pal_fig3$ink),
      legend.text=element_text(color=pal_fig3$ink),
      plot.tag=element_text(face="bold", color=pal_fig3$ink),
      plot.margin=margin(8, 8, 8, 16)
    )
}

#--------------------------------------------------------------#
#-Probability distributions of the two regimes (a, b)----------#
#--------------------------------------------------------------#

# Outline of a probability distribution drawn as touching bars of width 1:
# the staircase along the top of the bars, closed on the x axis at both
# ends. Returns a data frame with the corners x and y.
distribution_outline <- function(values, probability) {

  # Two upper corners per bar, half a unit to each side of its value. The
  # outline starts and ends at a probability of 0.
  outline <- data.frame(
    x=c(min(values) - 0.5, rep(values, each=2) + c(-0.5, 0.5), max(values) + 0.5),
    y=c(0, rep(probability, each=2), 0)
  )

  return(outline)
}

# Draw the empirical and the escalated distribution of the example island,
# one over the other, as transparent shapes with a solid outline. Three
# lines of text name the island and the two regimes in the colors of the
# distributions, so the panel has no legend. data has the columns value,
# probability and regime. label_x and label_hjust place the text.
# Returns a ggplot object.
plot_regime <- function(data, regime_labels, x_lab, y_max, label_x, label_hjust) {

  # Outline of the distribution of each regime
  outlines <- lapply(regime_levels, function(regime_name) {
    regime_data <- data[data$regime == regime_name, ]
    regime_data <- regime_data[order(regime_data$value), ]
    outline <- distribution_outline(regime_data$value, regime_data$probability)
    outline$regime <- factor(regime_name, levels=regime_levels)
    return(outline)
  })
  outlines <- do.call(rbind, outlines)

  # The three lines of text, from the island name downwards
  labels <- data.frame(
    label=c(example_island, regime_labels[["Empirical"]], regime_labels[["Escalated"]]),
    color=c(pal_fig3$ink, pal_regime[["Empirical"]], pal_regime[["Escalated"]]),
    face=c("bold", "plain", "plain"),
    y=y_max * c(0.96, 0.87, 0.78)
  )

  p <- ggplot() +
    geom_polygon(data=outlines, aes(x=x, y=y, fill=regime), alpha=0.40) +
    geom_path(data=outlines, aes(x=x, y=y, color=regime), linewidth=0.8) +
    geom_text(data=labels,
              aes(x=label_x, y=y, label=label, fontface=face),
              color=labels$color, hjust=label_hjust, size=4.4) +
    scale_fill_manual(values=pal_regime, guide="none") +
    scale_color_manual(values=pal_regime, guide="none") +
    scale_x_continuous(breaks=sort(unique(data$value))) +
    scale_y_continuous(limits=c(0, y_max), expand=expansion(mult=c(0, 0.02))) +
    labs(x=x_lab, y="Probability") +
    theme_fig3()

  return(p)
}

# Order of the regimes: the escalated distribution is drawn over the
# empirical one
regime_levels <- names(pal_regime)

#--------------------------------------------------------------#
#-Panel a: number of hurricanes in a year----------------------#
#--------------------------------------------------------------#

# Hurricanes follow a Poisson process, so lambda gives the probability of 0
# to 4 hurricanes in a year. The escalated regime is 0.60 hurricanes per
# year, close to the highest lambda of panel c.
k_values <- 0:4
lambda_empirical <- lambda_estimates$lambda[lambda_estimates$island == example_island]
lambda_escalated <- 0.60

poisson_example <- rbind(
  data.frame(regime="Empirical", value=k_values,
             probability=dpois(k_values, lambda_empirical)),
  data.frame(regime="Escalated", value=k_values,
             probability=dpois(k_values, lambda_escalated))
)
poisson_example$regime <- factor(poisson_example$regime, levels=regime_levels)

frequency_labels <- c(
  "Empirical" = sprintf("Empirical frequency, %.2f per year", lambda_empirical),
  "Escalated" = sprintf("Increased frequency, %.2f per year", lambda_escalated)
)

p3.1 <- plot_regime(poisson_example, frequency_labels,
                    x_lab="Number of hurricanes in a year",
                    y_max=1, label_x=4.5, label_hjust=1)

#--------------------------------------------------------------#
#-Panel b: hurricane categories--------------------------------#
#--------------------------------------------------------------#

# Move probability from weak to strong hurricanes, as in the intensity
# scenario of the simulations. A shift s is removed from categories 1 to 3
# and added to categories 4 and 5. Categories keep their relative
# probabilities within each group, and the shift cannot exceed the
# probability of categories 1 to 3. Returns a data frame with category and
# probability.
shift_intensity_pmf <- function(pmf_df, shift_amount) {

  p_low  <- sum(pmf_df$prob[pmf_df$category <= 3])
  p_high <- sum(pmf_df$prob[pmf_df$category >= 4])

  actual_shift <- min(shift_amount, p_low)

  scale_low <- ifelse(p_low > 0,
                      (p_low - actual_shift) / p_low,
                      0)
  scale_high <- ifelse(p_high > 0,
                       (p_high + actual_shift) / p_high,
                       0)

  shifted_pmf <- pmf_df %>%
    mutate(
      probability=case_when(
        category <= 3 ~ prob * scale_low,
        TRUE          ~ prob * scale_high
      )
    ) %>%
    select(category, probability)

  return(shifted_pmf)
}

category_values <- 1:5

# The escalated regime is a shift of 0.30 on the x axis of panel d, where 0
# is the empirical PMF
intensity_shift <- 0.30

# Empirical PMF of Trinidad, with probability 0 for the categories absent
# from its record
intensity_empirical <- data.frame(category=category_values) %>%
  left_join(
    intensity_pmf %>%
      filter(island == example_island) %>%
      mutate(category=as.integer(gsub("Cat_", "", Category))) %>%
      select(category, prob),
    by="category"
  ) %>%
  mutate(prob=replace_na(prob, 0))

intensity_escalated <- shift_intensity_pmf(intensity_empirical, intensity_shift)

intensity_example <- rbind(
  data.frame(regime="Empirical", value=intensity_empirical$category,
             probability=intensity_empirical$prob),
  data.frame(regime="Escalated", value=intensity_escalated$category,
             probability=intensity_escalated$probability)
)
intensity_example$regime <- factor(intensity_example$regime, levels=regime_levels)

intensity_labels <- c(
  "Empirical" = "Empirical intensity",
  "Escalated" = sprintf("Shift in intensity of %.2f", intensity_shift)
)

p3.2 <- plot_regime(intensity_example, intensity_labels,
                    x_lab="Hurricane category",
                    y_max=0.6, label_x=0.5, label_hjust=0)

#--------------------------------------------------------------#
#-Panel c: extinction under increasing frequency---------------#
#--------------------------------------------------------------#

# One curve per island: fraction of species extinct for lambda from 0.01 to
# 0.80 hurricanes per year. The point is the empirical lambda of the island
# (interpolated).
p3.3 <- ggplot() +
  geom_line(data=results_lambda %>% filter(level == "Network"),
            aes(x=lambda, y=extinction_fraction, group=draw_order, color=network),
            linewidth=0.9, alpha=0.95) +
  geom_point(data=results_interp %>% filter(level == "Network"),
             aes(x=lambda, y=extinction_fraction, fill=network),
             color="black", shape=21, size=3.0, stroke=0.4) +
  scale_color_island(continental_islands, "Continental") +
  scale_fill_island() +
  scale_x_continuous(limits=c(0, 0.8),
                     breaks=seq(0, 0.8, by=0.2),
                     expand=expansion(mult=c(0, 0.02))) +
  scale_y_continuous(limits=c(0, 1),
                     breaks=seq(0, 1, by=0.25),
                     expand=expansion(mult=c(0, 0.02))) +
  labs(x="Average number of hurricanes per year",
       y="Extinction") +
  coord_cartesian(clip="off") +
  theme_fig3()

#--------------------------------------------------------------#
#-Panel d: extinction under increasing intensity---------------#
#--------------------------------------------------------------#

# One curve per island: fraction of species extinct for shifts from 0 to 1.
# The point at a shift of 0 is the empirical PMF of the island. The points
# are on the y axis, so the panel is not clipped at its edges.
p3.4 <- ggplot(data=results_pmf %>% filter(level == "Network"),
               aes(x=intensity, y=extinction_fraction, group=draw_order, color=network)) +
  geom_line(linewidth=0.9, alpha=0.95) +
  geom_point(data=results_pmf %>% filter(level == "Network", intensity == 0),
             aes(x=intensity, y=extinction_fraction, fill=network),
             color="black", shape=21, size=3.0, stroke=0.4) +
  scale_color_island(oceanic_islands, "Oceanic") +
  scale_fill_island() +
  scale_x_continuous(limits=c(0, 1),
                     breaks=seq(0, 1, by=0.25),
                     expand=expansion(mult=c(0, 0.02))) +
  scale_y_continuous(limits=c(0, 1),
                     breaks=seq(0, 1, by=0.25),
                     expand=expansion(mult=c(0, 0.02))) +
  labs(x="Shift in hurricane intensity", y="Extinction") +
  coord_cartesian(clip="off") +
  theme_fig3()

#--------------------------------------------------------------#
#-Assemble and save Figure 3-----------------------------------#
#--------------------------------------------------------------#

# Regimes of Trinidad (a, b) above the extinction curves of the 11 islands
# (c, d), with the two legends side by side below the figure
figure_3 <- ((p3.1 + p3.2) / (p3.3 + p3.4)) +
  plot_layout(guides="collect", heights=c(1.0, 1.05)) +
  plot_annotation(tag_levels="a") &
  theme(
    legend.position="bottom",
    legend.direction="horizontal",
    legend.title=element_text(size=11, face="bold"),
    legend.text=element_text(size=10),
    legend.key.width=unit(0.70, "cm"),
    legend.key.height=unit(0.35, "cm"),
    legend.margin=margin(0, 0, 0, 0),
    legend.box.margin=margin(0, 0, 0, 0),
    legend.box="horizontal",
    legend.spacing.x=unit(1.0, "cm"),
    legend.key.spacing.x=unit(0.25, "cm"),
    legend.key.spacing.y=unit(0.04, "cm")
  )

save_figure(figure_3, "figure_3", width=11, height=10)
print(figure_3)
