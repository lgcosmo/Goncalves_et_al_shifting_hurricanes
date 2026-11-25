#--------------------------------------------------------------#
#-Empirical hurricane regime of each island--------------------#
#--------------------------------------------------------------#

# Estimate the hurricane regime of each of the 11 islands from the
# historical hurricane record. Hurricane frequency (lambda) is the number of
# recorded hurricanes per year of record, the rate of the Poisson process of
# the simulations. Hurricane intensity is the probability mass function
# (PMF) of Saffir-Simpson categories 1 to 5 among the recorded hurricanes.
#
# Inputs:  src/hurricane_crono.csv
# Outputs: src/empirical_lambdas.csv
#          src/empirical_pmfs.csv

library(here)
library(tidyverse)

#--------------------------------------------------------------#
#-Historical hurricane record----------------------------------#
#--------------------------------------------------------------#

hurricane_data <- read_csv(here("src", "hurricane_crono.csv"), show_col_types=FALSE)

#--------------------------------------------------------------#
#-Hurricane frequency (lambda)---------------------------------#
#--------------------------------------------------------------#

# Hurricanes per year of record of each island. The record of an island goes
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
#-Hurricane intensity (PMF of categories)----------------------#
#--------------------------------------------------------------#

# Probability of each Saffir-Simpson category among the recorded hurricanes
# of each island. Each hurricane has one category indicator.
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
