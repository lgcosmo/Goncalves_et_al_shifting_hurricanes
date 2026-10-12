#--------------------------------------------------------------#
#-Empirical hurricane regime of each island--------------------#
#--------------------------------------------------------------#

# Estimate the hurricane regime of each of the 11 islands from the historical
# hurricane record. The hurricane frequency lambda is the number of recorded
# hurricanes per year of the island's record, the rate of the Poisson process
# of the simulations. The hurricane intensity is the probability mass
# function (PMF) of Saffir-Simpson categories 1 to 5 among the recorded
# hurricanes. Both appear in Fig. 3a and 3b (figure_3.R) and define the
# empirical regime of every simulation script.
#
# Inputs:  src/hurricane_crono.csv, the historical hurricane record, with one
#          row per hurricane and island (network, name, year, and one
#          indicator column per category)
# Outputs: src/empirical_lambdas.csv, one row per island with its frequency
#          src/empirical_pmfs.csv, one row per island and category
# Used by: the three scripts of 03_simulations, survival_curves.R,
#          build_island_data.R and figure_3.R

library(here)
library(tidyverse)

#--------------------------------------------------------------#
#-Historical hurricane record----------------------------------#
#--------------------------------------------------------------#

hurricane_data <- read_csv(here("src", "hurricane_crono.csv"), show_col_types=FALSE)

#--------------------------------------------------------------#
#-Hurricane frequency (Fig. 3a)--------------------------------#
#--------------------------------------------------------------#

# Hurricanes per year of each island's record, the expected number of
# hurricanes per year of the Poisson process. The record of an island runs
# from its first to its last recorded hurricane, so the observation period
# differs among islands (150 to 173 years).
lambda_estimates <- hurricane_data %>%
  group_by(network) %>%
  summarise(
    first_year=min(year),
    last_year=max(year),
    n_events=n(),
    .groups="drop"
  ) %>%
  mutate(
    observation_period_years=last_year - first_year + 1,
    lambda=round(n_events / observation_period_years, 3)
  )

#--------------------------------------------------------------#
#-Hurricane intensity (Fig. 3b)--------------------------------#
#--------------------------------------------------------------#

# Probability of each Saffir-Simpson category among the recorded hurricanes
# of each island. Each hurricane carries one category indicator.
intensity_pmf <- hurricane_data %>%
  pivot_longer(starts_with("Cat_"), names_to="Category", values_to="is_present") %>%
  filter(is_present == 1) %>%
  count(network, Category, name="n_category") %>%
  group_by(network) %>%
  mutate(prob=n_category / sum(n_category)) %>%
  ungroup()

#--------------------------------------------------------------#
#-Hurricane regimes, written to src/---------------------------#
#--------------------------------------------------------------#

write_csv(lambda_estimates, here("src", "empirical_lambdas.csv"))
write_csv(intensity_pmf, here("src", "empirical_pmfs.csv"))
message("Saved: empirical_lambdas.csv and empirical_pmfs.csv")
