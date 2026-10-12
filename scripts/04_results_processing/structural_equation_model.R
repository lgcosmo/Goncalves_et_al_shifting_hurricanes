#--------------------------------------------------------------#
#-Piecewise Structural Equation Model of Species Fragility-----#
#--------------------------------------------------------------#

# Fit the piecewise structural equation model (pSEM) of the fragility of bird
# species to hurricanes, reported in the SEM tables and drawn in Figure 4. The
# two responses are species fragility to hurricane frequency and to hurricane
# intensity, the normalized area under the extinction curve (AUC) of each bird
# species in the two scenarios. Higher AUC means greater fragility. Predictors
# enter at three levels of organization. Island biogeography holds island type
# (continental or oceanic), area, elevation and isolation. Network topology
# enters through the two oriented PCA axes of calculate_pca.R, global network
# integration and pollination topological dominance. The species level adds
# the multilayer participation coefficient
# (MPC), the two within-layer degree z-scores and body mass. Directed paths
# run from all three levels to each response, from island biogeography to
# each topology axis, and from body mass to each functional role. Every
# variable is standardized before fitting, so the estimates are standardized
# path coefficients. Fisher's C tests the overall fit.
#
# Inputs:  data/main_summary/main_summary_auc.csv (calculate_auc_extinction.R)
#          data/main_summary/island_data.csv (build_island_data.R)
#          src/bodymass.csv
# Outputs: tables/sem_coefficients_clean.csv, sem_covariances_clean.csv,
#          sem_r_squared_clean.csv, sem_model_fit.csv,
#          sem_directed_separation_clean.csv and sem_variable_legend.csv,
#          collected in tables/sem_results_clean.xlsx. The pSEM summary is
#          printed to the console.
# Used by: figure_4.R draws the significant directed paths

library(here)
library(tidyverse)
library(piecewiseSEM)
library(lme4)

#--------------------------------------------------------------#
#-Load Fragility, Island and Body-Mass Data--------------------#
#--------------------------------------------------------------#

# Species names are lower-cased on both sides of the body-mass join, so the
# join ignores their capitalization
auc_data    <- read_csv(here("data", "main_summary", "main_summary_auc.csv"), show_col_types=FALSE)
island_data <- read_csv(here("data", "main_summary", "island_data.csv"),      show_col_types=FALSE)
bodymass    <- read_csv(here("src", "bodymass.csv"),                          show_col_types=FALSE) %>%
  mutate(species_clean=str_to_lower(species))

#--------------------------------------------------------------#
#-Build the pSEM Dataset, One Row per Bird and Island----------#
#--------------------------------------------------------------#

# Only bird species enter the pSEM, one row per species and island, with its
# fragility to hurricane intensity (auc_pmf) and to hurricane frequency
# (auc_lambda) side by side and its three role metrics
auc_wide <- auc_data %>%
  filter(level == "Species", species_set == "Animal") %>%
  select(network, species_id, scenario, auc, mpc, z_pollination, z_dispersal) %>%
  pivot_wider(names_from=scenario, values_from=auc, names_prefix="auc_") %>%
  rename_with(tolower, starts_with("auc_")) %>%
  mutate(species_clean=str_to_lower(species_id))

# Each species row receives the biogeography and the two PCA axes of its
# island, and the body mass of its species
sem_data <- auc_wide %>%
  left_join(
    island_data %>% select(network, island, type, area, elevation, isolation,
                           global_integration, pollination_dominance),
    by="network"
  ) %>%
  left_join(bodymass, by="species_clean")

#--------------------------------------------------------------#
#-Standardize the pSEM Variables-------------------------------#
#--------------------------------------------------------------#

# Every variable is z-scored, so the estimates are standardized path
# coefficients. Island area and body mass are log-transformed first because
# both span orders of magnitude.
sem_df <- sem_data %>%
  mutate(
    area_sc        = as.numeric(scale(log(area))),
    elev_sc        = as.numeric(scale(elevation)),
    isol_sc        = as.numeric(scale(isolation)),
    PC1_sc         = as.numeric(scale(global_integration)),
    PC2_sc         = as.numeric(scale(pollination_dominance)),
    mpc            = as.numeric(scale(mpc)),
    z_pollination  = as.numeric(scale(z_pollination)),
    z_dispersal    = as.numeric(scale(z_dispersal)),
    bm_avg         = as.numeric(scale(log(bm_avg))),
    auc_pmf        = as.numeric(scale(auc_pmf)),
    auc_lambda     = as.numeric(scale(auc_lambda)),
    island         = as.factor(island),
    type           = as.factor(type),
    species_id     = as.factor(species_id)
  ) %>%
  as.data.frame()

# The two topology axes are island-level variables, so their models are
# fitted to one row per island (N = 11)
island_level_df <- sem_df %>%
  distinct(island, type, area_sc, elev_sc, isol_sc, PC1_sc, PC2_sc) %>%
  as.data.frame()

#--------------------------------------------------------------#
#-Fit the Piecewise SEM----------------------------------------#
#--------------------------------------------------------------#

# The mixed models use the bobyqa optimizer with a large evaluation budget so
# that every one of them converges
ctrl <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=200000))

model_list <- psem(

  # Species fragility to hurricane intensity (auc_pmf) and to hurricane
  # frequency (auc_lambda). A random intercept per species accounts for the
  # same bird species occurring on several islands.
  lmer(auc_pmf    ~ PC1_sc + PC2_sc +
                    mpc + z_pollination + z_dispersal + bm_avg +
                    type + area_sc + elev_sc + isol_sc + (1 | species_id),
       data=sem_df, control=ctrl),

  lmer(auc_lambda ~ PC1_sc + PC2_sc +
                    mpc + z_pollination + z_dispersal + bm_avg +
                    type + area_sc + elev_sc + isol_sc + (1 | species_id),
       data=sem_df, control=ctrl),

  # The two topology axes respond to island biogeography, one row per island
  lm(PC1_sc ~ type + area_sc + elev_sc + isol_sc, data=island_level_df),
  lm(PC2_sc ~ type + area_sc + elev_sc + isol_sc, data=island_level_df),

  # The three functional roles respond to body mass and to island
  # biogeography, with random intercepts for island and for species
  lmer(mpc           ~ bm_avg + type + area_sc + elev_sc + isol_sc + (1 | island) + (1 | species_id), data=sem_df, control=ctrl),
  lmer(z_pollination ~ bm_avg + type + area_sc + elev_sc + isol_sc + (1 | island) + (1 | species_id), data=sem_df, control=ctrl),
  lmer(z_dispersal   ~ bm_avg + type + area_sc + elev_sc + isol_sc + (1 | island) + (1 | species_id), data=sem_df, control=ctrl),

  # Correlated errors among the three role metrics, between the two fragility
  # responses and between the two topology axes. These pairs co-vary without
  # a directed path between them.
  z_pollination %~~% z_dispersal,
  z_pollination %~~% mpc,
  z_dispersal   %~~% mpc,
  auc_pmf       %~~% auc_lambda,
  PC1_sc        %~~% PC2_sc,

  data=sem_df
)

#--------------------------------------------------------------#
#-Print the pSEM Summary to the Console------------------------#
#--------------------------------------------------------------#

cat("\n==============================================\n")
cat("          RUNNING GLOBAL SEM MODEL            \n")
cat("==============================================\n")

# If the summary of the whole pSEM fails, each component model is summarized
# on its own
tryCatch(
  print(summary(model_list, .progressBar=FALSE)),
  error = function(e) {
    cat("\n[WARN] Global summary failed. Printing individual model summaries.\n\n")
    for (mod in model_list) {
      if (is.null(mod) || is.atomic(mod)) next
      tryCatch(print(summary(mod)), error = function(e) {})
      cat("---\n")
    }
  }
)

#--------------------------------------------------------------#
#-Export the SEM Tables----------------------------------------#
#--------------------------------------------------------------#

# piecewiseSEM reports the variables by their code names. One table gives,
# for every variable, its code, its name in the manuscript and its
# definition. Island type is categorical, so piecewiseSEM reports one
# coefficient per island type, and those two rows come first so that
# "type = Oceanic" is translated before "type". The table is also the legend
# written beside the SEM tables.
sem_variables <- tribble(
  ~code,                ~name,                                   ~definition,
  "type = Oceanic",     "Oceanic islands",                       "Coefficient for oceanic islands in categorical island-type effects.",
  "type = Continental", "Continental islands",                   "Coefficient for continental islands in categorical island-type effects.",
  "auc_pmf",            "Fragility to hurricane intensity",      "Species fragility in response to shifts in hurricane intensity.",
  "auc_lambda",         "Fragility to hurricane frequency",      "Species fragility in response to increases in hurricane frequency.",
  "PC1_sc",             "Global network integration",            "First PCA axis of island-level multilayer structure; larger values indicate globally more integrated networks.",
  "PC2_sc",             "Pollination topological dominance",     "Second PCA axis of island-level multilayer structure, oriented so that larger values indicate stronger pollination topological dominance.",
  "mpc",                "Multilayer participation coefficient",  "Species-level participation across pollination and seed-dispersal layers.",
  "z_pollination",      "Pollination layer degree (z-score)",    "Standardized degree in the pollination layer.",
  "z_dispersal",        "Seed-dispersal layer degree (z-score)", "Standardized degree in the seed-dispersal layer.",
  "bm_avg",             "Body mass",                             "Average species body mass.",
  "area_sc",            "Island area",                           "Scaled log island area.",
  "elev_sc",            "Elevation",                             "Scaled maximum island elevation.",
  "isol_sc",            "Isolation distance",                    "Scaled distance from the continent.",
  "type",               "Island type",                           "Island type, continental or oceanic."
)

sem_names <- setNames(sem_variables$name, sem_variables$code)

# Translate the variable names reported by piecewiseSEM into their manuscript
# names. Correlated errors carry a leading "~~", which is removed before the
# lookup. Arguments: x is a vector of names. Returns the manuscript names, and
# the original name where no translation exists.
clean_sem_name <- function(x) {
  code <- str_replace(as.character(x), "^~~", "")
  coalesce(unname(sem_names[code]), code)
}

# Translate the variable names inside an independence claim, which
# piecewiseSEM writes as a formula. Arguments: x is a vector of claims.
# Returns the claims with every code name replaced by its manuscript name.
clean_sem_claim <- function(x) {
  str_replace_all(as.character(x), sem_names)
}

# Convert a column of the piecewiseSEM tables to numbers. The test of a
# categorical predictor has no estimate, which piecewiseSEM stores as "-", so
# the conversion turns it into NA without a warning. Arguments: x is the
# column. Returns a numeric vector.
as_sem_numeric <- function(x) {
  return(suppressWarnings(as.numeric(as.character(x))))
}

# piecewiseSEM leaves the column of significance stars unnamed
raw_coefs <- coefs(model_list)
names(raw_coefs)[names(raw_coefs) == ""] <- "Significance"

# Directed paths and correlated errors are told apart by the "~~" prefix of
# the response
sem_coefficients_all <- raw_coefs %>%
  as_tibble() %>%
  mutate(
    path_type=if_else(str_detect(Response, "^~~"),
                      "Correlated error", "Directed path"),
    response=clean_sem_name(Response),
    predictor=clean_sem_name(Predictor),
    estimate=as_sem_numeric(Estimate),
    std_error=as_sem_numeric(Std.Error),
    df=as_sem_numeric(DF),
    statistic=as_sem_numeric(Crit.Value),
    p_value=as_sem_numeric(P.Value),
    std_estimate=as_sem_numeric(Std.Estimate),
    significance=na_if(str_trim(as.character(Significance)), "")
  ) %>%
  transmute(
    path_type,
    response,
    predictor,
    estimate,
    std_error,
    df,
    statistic,
    p_value,
    std_estimate,
    significance
  )

sem_coefficients <- sem_coefficients_all %>%
  filter(path_type == "Directed path")

sem_covariances <- sem_coefficients_all %>%
  filter(path_type == "Correlated error") %>%
  transmute(
    path_type,
    variable_1=response,
    variable_2=predictor,
    estimate,
    df,
    statistic,
    p_value,
    std_estimate,
    significance
  )

# For the mixed models, marginal R-squared covers the fixed effects and
# conditional R-squared adds the random intercepts
sem_r_squared <- rsquared(model_list) %>%
  as_tibble() %>%
  transmute(
    response=clean_sem_name(Response),
    method,
    marginal_r2=Marginal,
    conditional_r2=Conditional
  )

# Overall fit. Fisher's C combines the tests of directed separation below. A
# high P-value means that the paths left out of the model are consistent with
# the data.
sem_fit_aic <- AIC(model_list)
sem_fit_fisher <- fisherC(model_list)

sem_model_fit <- tibble(
  statistic=c("AIC", "Number of parameters", "AIC sample size",
              "Fisher's C", "Fisher's C df", "Fisher's C p-value"),
  value=c(sem_fit_aic$AIC, sem_fit_aic$K, sem_fit_aic$n,
          sem_fit_fisher$Fisher.C, sem_fit_fisher$df,
          sem_fit_fisher$P.Value)
)

# Tests of directed separation, one per pair of variables without a directed
# path between them. The unnamed column again holds the significance stars.
raw_directed_separation <- dSep(model_list, .progressBar=FALSE)
names(raw_directed_separation)[names(raw_directed_separation) == ""] <- "Significance"

sem_directed_separation <- raw_directed_separation %>%
  as_tibble() %>%
  transmute(
    independence_claim=clean_sem_claim(Independ.Claim),
    test_type=Test.Type,
    df=DF,
    statistic=Crit.Value,
    p_value=P.Value,
    significance=na_if(str_trim(as.character(Significance)), "")
  )

# The legend of the code names, written beside the tables
sem_variable_legend <- sem_variables %>%
  rename(descriptive_name=name)

# The six tables, written one by one as CSV and together as one workbook
sem_tables <- list(
  "Coefficients"=sem_coefficients,
  "Covariances"=sem_covariances,
  "R-squared"=sem_r_squared,
  "Model fit"=sem_model_fit,
  "Directed separation"=sem_directed_separation,
  "Variable legend"=sem_variable_legend
)

table_files <- c(
  "Coefficients"="sem_coefficients_clean.csv",
  "Covariances"="sem_covariances_clean.csv",
  "R-squared"="sem_r_squared_clean.csv",
  "Model fit"="sem_model_fit.csv",
  "Directed separation"="sem_directed_separation_clean.csv",
  "Variable legend"="sem_variable_legend.csv"
)

dir.create(here("tables"), showWarnings=FALSE, recursive=TRUE)

for (name in names(sem_tables)) {
  write_csv(sem_tables[[name]], here("tables", table_files[[name]]), na="")
}

writexl::write_xlsx(sem_tables, here("tables", "sem_results_clean.xlsx"))

message("Saved: the six SEM tables in tables/ and tables/sem_results_clean.xlsx")
