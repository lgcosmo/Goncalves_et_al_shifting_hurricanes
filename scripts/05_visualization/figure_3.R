#--------------------------------------------------------------#
#-Figure 3: Hurricane regimes and network extinction-----------#
#--------------------------------------------------------------#

# Figure 3 shows the hurricane regimes of the 11 islands and the extinctions
# they cause in the simulations. Panel a draws, for each island, the Poisson
# probabilities of 0 to 4 hurricanes in a year at its empirical hurricane
# frequency, and panel b the empirical probability mass function (PMF) of
# hurricane categories 1 to 5. Both panels add, for Trinidad, a dashed line
# with an escalated regime, an increased hurricane frequency in a and
# probability mass shifted from categories 1 to 3 to categories 4 and 5 in b.
# Panels c and d show the fraction of species that went extinct on each island
# under increasing hurricane frequency (lambda from 0.01 to 0.80) and under
# increasing hurricane intensity (shift from 0 to 1). The point on each curve
# marks the island's empirical regime.
#
# Inputs:  src/empirical_lambdas.csv and src/empirical_pmfs.csv, the empirical
#          hurricane frequency and category PMF of each island, written by
#          empirical_lambda_pmf.R
#          data/main_summary/island_data.csv, written by build_island_data.R,
#          with the isolation of each island
#          data/main_summary/main_summary_metrics_lambda.csv and
#          main_summary_metrics_pmf.csv, the extinctions of the frequency and
#          intensity sweeps, and main_summary_empirical_interpolated.csv, the
#          extinctions at each island's empirical frequency, all written by
#          survival_curves.R
# Outputs: figures/figure_3.png, .pdf and .svg
# Helpers: island_names and island_order in src/island_labels.R and
#          save_figure in src/figure_helpers.R

library(here)
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)
library(tidyr)
library(grid)

source(here("src", "island_labels.R"))
source(here("src", "figure_helpers.R"))

#--------------------------------------------------------------#
#-Hurricane regimes (a, b) and extinction sweeps (c, d)--------#
#--------------------------------------------------------------#

# Replace the island codes of one column by the island names, ordered from
# Trinidad in the south to Saba in the north, so that colors and legend
# entries match across the four panels.
# Arguments: df is a data frame; column names the column with island codes.
# Returns df with that column as a factor whose levels follow island_order.
relabel_islands <- function(df, column) {
  df[[column]] <- factor(island_names[as.character(df[[column]])], levels=island_order)
  return(df)
}

# The two sweeps give the fraction of species that went extinct on each island
# at every hurricane frequency (panel c) and every intensity shift (panel d).
# The interpolated table gives the same fraction at each island's empirical
# frequency, which places the island's point on its curve in panel c.
results_lambda <- fread(here("data", "main_summary", "main_summary_metrics_lambda.csv")) %>%
  relabel_islands("network")
results_interp <- fread(here("data", "main_summary", "main_summary_empirical_interpolated.csv")) %>%
  relabel_islands("network")
results_pmf    <- fread(here("data", "main_summary", "main_summary_metrics_pmf.csv")) %>%
  relabel_islands("network")

# The empirical hurricane regimes of panels a and b, estimated from the
# historical hurricane record. Lambda is the expected number of hurricanes per
# year, and the PMF gives the probability of each Saffir-Simpson category from
# 1 to 5.
lambda_estimates <- fread(here("src", "empirical_lambdas.csv")) %>%
  relabel_islands("network") %>%
  rename(island=network)
intensity_pmf <- fread(here("src", "empirical_pmfs.csv")) %>%
  relabel_islands("network") %>%
  rename(island=network)

# Distance to the continent of each island, which orders the colors of the
# oceanic islands below
island_isolation <- fread(here("data", "main_summary", "island_data.csv")) %>%
  relabel_islands("network") %>%
  transmute(island=network, isolation=as.numeric(isolation))

#--------------------------------------------------------------#
#-Island colors and theme of the four panels-------------------#
#--------------------------------------------------------------#

# Trinidad is the island drawn under the escalated regimes in panels a and b.
# Its deep violet keeps the pair of solid and dashed lines visible against the
# empirical range and the curves of the other islands.
example_island <- "Trinidad"

# Neutral tones for the band spanning the empirical range of the islands, its
# label, the Trinidad annotations and the axis text
pal_fig3 <- list(
  empirical_fill = "#edf2ef",
  empirical_line = "#5d716d",
  example        = "#2f3136",
  ink            = "#253238"
)

# Trinidad and Tobago, the continental islands, take violet and coral. The nine
# oceanic islands follow a gradient from blue to yellow with increasing
# distance to the continent, spaced evenly along the rank of isolation rather
# than distance itself, because the northern islands lie at similar distances.
continental_cols <- c(
  "Trinidad" = "#4B1D6B",
  "Tobago"   = "#D95F59"
)

oceanic_isolation <- island_isolation %>%
  filter(!island %in% names(continental_cols)) %>%
  arrange(isolation)

oceanic_gradient <- grDevices::colorRampPalette(
  c("#2C7FB8", "#2BB3A3", "#8BC34A", "#D6C500")
)(256)

oceanic_cols <- oceanic_gradient[round(seq(1, 256, length.out=nrow(oceanic_isolation)))]
names(oceanic_cols) <- as.character(oceanic_isolation$island)

island_cols <- c(continental_cols, oceanic_cols)[island_order]

# Color scale shared by the four panels. Only panel c shows the island
# legend, which patchwork then collects below the assembled figure.
# Arguments: show_guide is TRUE to draw the legend. Returns a
# scale_color_manual with the island palette.
scale_color_island <- function(show_guide=FALSE) {

  island_guide <- if (show_guide) {
    guide_legend(
      nrow=2,
      byrow=TRUE,
      override.aes=list(linewidth=1.4, alpha=1)
    )
  } else {
    "none"
  }

  scale_color_manual(
    values=island_cols,
    breaks=island_order,
    name="Island",
    guide=island_guide
  )
}

# Fill scale for the points of panels c and d, with the same palette and no
# legend. Returns a scale_fill_manual.
scale_fill_island <- function() {
  scale_fill_manual(
    values=island_cols,
    breaks=island_order,
    guide="none"
  )
}

# Theme shared by the four panels, built on theme_pubr without grid lines.
# Arguments: base_size is the font size in points. Returns a ggplot theme.
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
#-Panel a: Poisson probabilities at the empirical frequency----#
#--------------------------------------------------------------#

# Hurricanes follow a Poisson process, so the hurricane frequency lambda of an
# island gives the probability of 0 to 4 hurricanes in a year. A band spans the
# range of the 11 islands, one line per island shows its empirical frequency,
# and the dashed line shows Trinidad at an increased frequency of 0.60
# hurricanes per year, near the top of the frequency sweep.
k_values <- 0:4
lambda_example <- 0.60

# Probabilities at the empirical frequency of each island, and their range
# across islands at each number of hurricanes
poisson_islands <- crossing(lambda_estimates, k=k_values) %>%
  mutate(probability=dpois(k, lambda))

poisson_summary <- poisson_islands %>%
  group_by(k) %>%
  summarise(prob_min=min(probability), prob_max=max(probability), .groups="drop")

# Trinidad under its empirical frequency, for the solid line and the
# annotation, and under the increased frequency, for the dashed line
poisson_example_baseline <- poisson_islands %>%
  filter(island == example_island)

example_lambda <- unique(poisson_example_baseline$lambda)

poisson_example <- tibble(island=example_island, k=k_values,
                          probability=dpois(k_values, lambda_example))

p3.1 <- ggplot() +
  geom_ribbon(data=poisson_summary,
              aes(x=k, ymin=prob_min, ymax=prob_max),
              fill=pal_fig3$empirical_fill, alpha=0.95) +
  geom_line(data=poisson_islands,
            aes(x=k, y=probability, group=island, color=island),
            alpha=0.72, linewidth=0.55) +
  geom_line(data=poisson_example_baseline,
            aes(x=k, y=probability, color=island),
            linewidth=1.15, alpha=1) +
  geom_line(data=poisson_example,
            aes(x=k, y=probability, color=island),
            linewidth=1.15, alpha=1, linetype="dashed") +
  annotate("curve",
           x=0.25, y=0.83, xend=1.05, yend=0.34,
           curvature=-0.28, color=pal_fig3$example, linewidth=0.45,
           arrow=arrow(length=unit(0.16, "cm"), type="closed")) +
  annotate("text",
           x=0.10, y=0.93, label="Empirical island range",
           color=pal_fig3$empirical_line, hjust=0, size=4.0) +
  annotate("text",
           x=1.32, y=0.42,
           label=paste0(example_island, " frequency\n",
                        "solid = ", sprintf("%.2f", example_lambda),
                        " hurricanes/year\n",
                        "dashed = ", sprintf("%.2f", lambda_example),
                        " hurricanes/year"),
           color=pal_fig3$example, hjust=0, size=4.0) +
  scale_x_continuous(breaks=k_values) +
  scale_color_island(show_guide=FALSE) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(mult=c(0, 0.04))) +
  labs(x="Number of hurricanes in a year", y="Probability") +
  theme_fig3()

#--------------------------------------------------------------#
#-Panel b: Empirical PMF of hurricane categories---------------#
#--------------------------------------------------------------#

# Shift probability mass from weak to strong hurricanes for the dashed line of
# Trinidad, as shift_pmf in main_simulations_empirical_pmf.jl does for the
# intensity scenario. A shift s removes mass s from categories 1 to 3 and adds
# it to categories 4 and 5. Within each group, the categories keep their
# relative probabilities, and the shift cannot exceed the mass of categories
# 1 to 3.
# Arguments: pmf_df holds the empirical PMF of one island, with the columns
# category and prob; shift_amount is the shift s. Returns a data frame with
# the columns category and probability.
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

# The shift drawn for Trinidad is 0.30, on the scale of the x axis of panel d,
# where 0 is the empirical PMF and larger values move more probability mass
# from categories 1 to 3 to categories 4 and 5.
intensity_shift <- 0.30

# Fill the five categories of every island, with probability 0 for any
# category absent from its hurricane record
intensity_complete <- crossing(
  island=factor(island_order, levels=island_order),
  category=category_values
) %>%
  left_join(
    intensity_pmf %>%
      mutate(category=as.integer(gsub("Cat_", "", Category))) %>%
      select(island, category, prob),
    by=c("island", "category")
  ) %>%
  mutate(prob=replace_na(prob, 0))

# Range of the 11 islands in each category, and Trinidad before and after
# the shift
intensity_summary <- intensity_complete %>%
  group_by(category) %>%
  summarise(prob_min=min(prob), prob_max=max(prob), .groups="drop")

intensity_example_baseline <- intensity_complete %>%
  filter(island == example_island)

intensity_example <- shift_intensity_pmf(intensity_example_baseline, intensity_shift) %>%
  mutate(island=example_island)

p3.2 <- ggplot() +
  geom_linerange(data=intensity_summary,
                 aes(x=category, ymin=prob_min, ymax=prob_max),
                 color=pal_fig3$empirical_fill, linewidth=7.5,
                 alpha=0.95, lineend="round") +
  geom_line(data=intensity_complete,
            aes(x=category, y=prob, group=island, color=island),
            alpha=0.55, linewidth=0.55) +
  geom_point(data=intensity_complete,
             aes(x=category, y=prob, color=island),
             alpha=0.70, size=1.9,
             position=position_jitter(width=0.07, height=0, seed=10)) +
  geom_line(data=intensity_example_baseline,
            aes(x=category, y=prob, color=island),
            linewidth=1.15, alpha=1) +
  geom_line(data=intensity_example,
            aes(x=category, y=probability, color=island),
            linewidth=1.15, alpha=1, linetype="dashed") +
  annotate("segment",
           x=2.15, xend=4.35, y=0.51, yend=0.51,
           color=pal_fig3$example, linewidth=0.45,
           arrow=arrow(length=unit(0.16, "cm"), type="closed")) +
  annotate("text",
           x=3.25, y=0.60,
           label=paste0(example_island, " intensity\n",
                        "solid = empirical; dashed = shift ",
                        sprintf("%.2f", intensity_shift)),
           color=pal_fig3$example, hjust=0.5, size=4.0) +
  annotate("text",
           x=1.00, y=0.45, label="Empirical island range",
           color=pal_fig3$empirical_line, hjust=0, size=4.0) +
  scale_x_continuous(breaks=category_values) +
  scale_color_island(show_guide=FALSE) +
  scale_y_continuous(limits=c(0, 0.68), expand=expansion(mult=c(0, 0.04))) +
  labs(x="Hurricane category", y="Probability") +
  theme_fig3()

#--------------------------------------------------------------#
#-Panel c: Extinction under increasing frequency---------------#
#--------------------------------------------------------------#

# One curve per island gives the fraction of species that went extinct across
# the frequency sweep, from 0.01 to 0.80 hurricanes per year. The point on
# each curve marks the island's empirical frequency, interpolated from the
# sweep.
p3.3 <- ggplot() +
  geom_line(data=results_lambda %>% filter(level == "Network"),
            aes(x=lambda, y=extinction_fraction, group=network, color=network),
            linewidth=0.75, alpha=0.92) +
  geom_point(data=results_interp %>% filter(level == "Network"),
             aes(x=lambda, y=extinction_fraction, fill=network),
             color="white", shape=21, size=3.0, stroke=0.35) +
  scale_color_island(show_guide=TRUE) +
  scale_fill_island() +
  scale_x_continuous(limits=c(0, 0.8),
                     breaks=seq(0, 0.8, by=0.2),
                     expand=expansion(mult=c(0, 0.02))) +
  scale_y_continuous(limits=c(0, 1),
                     breaks=seq(0, 1, by=0.25),
                     expand=expansion(mult=c(0, 0.02))) +
  labs(x="Average number of hurricanes per year",
       y="Extinction") +
  theme_fig3()

#--------------------------------------------------------------#
#-Panel d: Extinction under increasing intensity---------------#
#--------------------------------------------------------------#

# One curve per island gives the fraction of species that went extinct across
# the intensity sweep, from a shift of 0 to a shift of 1. The point at a shift
# of 0 marks the island's empirical PMF, which is its empirical hurricane
# regime.
p3.4 <- ggplot(data=results_pmf %>% filter(level == "Network"),
               aes(x=intensity, y=extinction_fraction, group=network, color=network)) +
  geom_line(linewidth=0.75, alpha=0.92) +
  geom_point(data=results_pmf %>% filter(level == "Network", intensity == 0),
             aes(x=intensity, y=extinction_fraction, fill=network),
             color="white", shape=21, size=3.0, stroke=0.35) +
  scale_color_island(show_guide=FALSE) +
  scale_fill_island() +
  scale_x_continuous(limits=c(0, 1),
                     breaks=seq(0, 1, by=0.25),
                     expand=expansion(mult=c(0, 0.02))) +
  scale_y_continuous(limits=c(0, 1),
                     breaks=seq(0, 1, by=0.25),
                     expand=expansion(mult=c(0, 0.02))) +
  labs(x="Shift in hurricane intensity", y="Extinction") +
  theme_fig3()

#--------------------------------------------------------------#
#-Assemble and save Fig. 3-------------------------------------#
#--------------------------------------------------------------#

# The hurricane regimes (a, b) sit above the extinction curves (c, d), and the
# island legend of panel c is collected below the figure
figure_3 <- ((p3.1 + p3.2) / (p3.3 + p3.4)) +
  plot_layout(guides="collect", heights=c(1.0, 1.05)) +
  plot_annotation(tag_levels="a") &
  theme(
    legend.position="bottom",
    legend.direction="horizontal",
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    legend.key.width=unit(0.70, "cm"),
    legend.key.height=unit(0.35, "cm"),
    legend.margin=margin(0, 0, 0, 0),
    legend.box.margin=margin(0, 0, 0, 0),
    legend.spacing.x=unit(0.16, "cm"),
    legend.spacing.y=unit(0.04, "cm")
  )

save_figure(figure_3, "figure_3", width=11, height=10)
print(figure_3)
