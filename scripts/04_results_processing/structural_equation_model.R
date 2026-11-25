#--------------------------------------------------------------#
#-Piecewise structural equation model of species persistence---#
#--------------------------------------------------------------#

# Fit the piecewise structural equation model (pSEM) of the persistence of
# bird species under hurricanes. The two responses are persistence under
# hurricane frequency and under hurricane intensity (1 - AUC of the
# extinction curve), so a positive path into a response favors the species.
# Predictors enter at three levels: island biogeography (type, area,
# elevation and isolation), network topology (the two PCA axes, global
# network integration and pollination topological dominance) and species
# (MPC, the two within-layer degree z-scores and body mass). Directed paths
# run from the three levels to each response, from island biogeography to
# each PCA axis, and from body mass and island biogeography to each role
# metric. All variables are standardized, so the estimates are standardized
# path coefficients.
#
# Inputs:  data/main_summary/main_summary_auc.csv
#          data/main_summary/island_data.csv
#          src/bodymass.csv
# Outputs: the six tables/sem_*.csv files (coefficients, covariances,
#          R-squared, model fit, directed separation and variable legend)
#          and tables/sem_results_clean.xlsx with all of them

library(here)
library(tidyverse)
library(piecewiseSEM)
library(lme4)

#--------------------------------------------------------------#
#-Load persistence, island and body-mass data------------------#
#--------------------------------------------------------------#

# Species names in lower case on both sides of the body-mass join
auc_data    <- read_csv(here("data", "main_summary", "main_summary_auc.csv"), show_col_types=FALSE)
island_data <- read_csv(here("data", "main_summary", "island_data.csv"),      show_col_types=FALSE)
bodymass    <- read_csv(here("src", "bodymass.csv"),                          show_col_types=FALSE) %>%
  mutate(species_clean=str_to_lower(species))

#--------------------------------------------------------------#
#-Build the pSEM dataset, one row per bird and island----------#
#--------------------------------------------------------------#

# Only birds enter the pSEM. One row per species and island, with persistence
# under intensity (persistence_pmf) and under frequency (persistence_lambda)
# side by side and the three role metrics.
persistence_wide <- auc_data %>%
  filter(level == "Species", species_set == "Animal") %>%
  select(network, species_id, scenario, persistence, mpc, z_pollination, z_dispersal) %>%
  pivot_wider(names_from=scenario, values_from=persistence, names_prefix="persistence_") %>%
  rename_with(tolower, starts_with("persistence_")) %>%
  mutate(species_clean=str_to_lower(species_id))

# Add the biogeography and PCA axes of the island and the body mass of the
# species
sem_data <- persistence_wide %>%
  left_join(
    island_data %>% select(network, island, type, area, elevation, isolation,
                           global_integration, pollination_dominance),
    by="network"
  ) %>%
  left_join(bodymass, by="species_clean")

#--------------------------------------------------------------#
#-Standardize the pSEM variables-------------------------------#
#--------------------------------------------------------------#

# z-score every variable. Island area and body mass are log-transformed
# first because they span orders of magnitude.
sem_df <- sem_data %>%
  mutate(
    area_sc            = as.numeric(scale(log(area))),
    elev_sc            = as.numeric(scale(elevation)),
    isol_sc            = as.numeric(scale(isolation)),
    PC1_sc             = as.numeric(scale(global_integration)),
    PC2_sc             = as.numeric(scale(pollination_dominance)),
    mpc                = as.numeric(scale(mpc)),
    z_pollination      = as.numeric(scale(z_pollination)),
    z_dispersal        = as.numeric(scale(z_dispersal)),
    bm_avg             = as.numeric(scale(log(bm_avg))),
    persistence_pmf    = as.numeric(scale(persistence_pmf)),
    persistence_lambda = as.numeric(scale(persistence_lambda)),
    island             = as.factor(island),
    type               = as.factor(type),
    species_id         = as.factor(species_id)
  ) %>%
  as.data.frame()

# The two PCA axes are island-level variables, so their models use one row
# per island (N = 11)
island_level_df <- sem_df %>%
  distinct(island, type, area_sc, elev_sc, isol_sc, PC1_sc, PC2_sc) %>%
  as.data.frame()

#--------------------------------------------------------------#
#-Fit the piecewise SEM----------------------------------------#
#--------------------------------------------------------------#

# bobyqa optimizer with a large number of evaluations, so all mixed models
# converge
ctrl <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=200000))

model_list <- psem(

  # Persistence under hurricane intensity and under hurricane frequency. The
  # random intercept accounts for bird species that occur on several islands.
  lmer(persistence_pmf    ~ PC1_sc + PC2_sc +
                            mpc + z_pollination + z_dispersal + bm_avg +
                            type + area_sc + elev_sc + isol_sc + (1 | species_id),
       data=sem_df, control=ctrl),

  lmer(persistence_lambda ~ PC1_sc + PC2_sc +
                            mpc + z_pollination + z_dispersal + bm_avg +
                            type + area_sc + elev_sc + isol_sc + (1 | species_id),
       data=sem_df, control=ctrl),

  # The two PCA axes as a function of island biogeography
  lm(PC1_sc ~ type + area_sc + elev_sc + isol_sc, data=island_level_df),
  lm(PC2_sc ~ type + area_sc + elev_sc + isol_sc, data=island_level_df),

  # The three role metrics as a function of body mass and island
  # biogeography, with random intercepts for island and species
  lmer(mpc           ~ bm_avg + type + area_sc + elev_sc + isol_sc + (1 | island) + (1 | species_id), data=sem_df, control=ctrl),
  lmer(z_pollination ~ bm_avg + type + area_sc + elev_sc + isol_sc + (1 | island) + (1 | species_id), data=sem_df, control=ctrl),
  lmer(z_dispersal   ~ bm_avg + type + area_sc + elev_sc + isol_sc + (1 | island) + (1 | species_id), data=sem_df, control=ctrl),

  # Correlated errors among the three role metrics, between the two
  # persistence responses and between the two PCA axes
  z_pollination   %~~% z_dispersal,
  z_pollination   %~~% mpc,
  z_dispersal     %~~% mpc,
  persistence_pmf %~~% persistence_lambda,
  PC1_sc          %~~% PC2_sc,

  data=sem_df
)

#--------------------------------------------------------------#
#-Print the pSEM summary to the console------------------------#
#--------------------------------------------------------------#

cat("\n==============================================\n")
cat("          RUNNING GLOBAL SEM MODEL            \n")
cat("==============================================\n")

# If the summary of the whole pSEM fails, summarize each model on its own
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
#-Export the SEM tables----------------------------------------#
#--------------------------------------------------------------#

# Code, manuscript name and definition of each variable. piecewiseSEM
# reports one coefficient per island type, and these two rows come first so
# "type = Oceanic" is translated before "type".
sem_variables <- tribble(
  ~code,                ~name,                                   ~definition,
  "type = Oceanic",     "Oceanic islands",                       "Coefficient for oceanic islands in categorical island-type effects.",
  "type = Continental", "Continental islands",                   "Coefficient for continental islands in categorical island-type effects.",
  "persistence_pmf",    "Persistence under hurricane intensity", "Species persistence across shifts in hurricane intensity, one minus the normalized area under the extinction curve.",
  "persistence_lambda", "Persistence under hurricane frequency", "Species persistence across increases in hurricane frequency, one minus the normalized area under the extinction curve.",
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

# Translate the piecewiseSEM variable names into their manuscript names.
# Correlated errors have a leading "~~", removed before the lookup. Names
# without a translation are returned as they are.
clean_sem_name <- function(x) {
  code <- str_replace(as.character(x), "^~~", "")
  coalesce(unname(sem_names[code]), code)
}

# Translate the variable names inside an independence claim (a formula)
clean_sem_claim <- function(x) {
  str_replace_all(as.character(x), sem_names)
}

# Convert a column of the piecewiseSEM tables to numbers. The test of a
# categorical predictor has no estimate ("-"), which becomes NA.
as_sem_numeric <- function(x) {
  return(suppressWarnings(as.numeric(as.character(x))))
}

# The column of significance stars has no name in piecewiseSEM
raw_coefs <- coefs(model_list)
names(raw_coefs)[names(raw_coefs) == ""] <- "Significance"

# Correlated errors have the "~~" prefix in the response
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

# Marginal R-squared: fixed effects. Conditional R-squared: fixed effects
# plus random intercepts.
sem_r_squared <- rsquared(model_list) %>%
  as_tibble() %>%
  transmute(
    response=clean_sem_name(Response),
    method,
    marginal_r2=Marginal,
    conditional_r2=Conditional
  )

# Overall fit. Fisher's C combines the tests of directed separation. A high
# P-value means the paths left out of the model are consistent with the data.
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
# path between them
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

# Legend of the variable codes
sem_variable_legend <- sem_variables %>%
  rename(descriptive_name=name)

# Save the six tables as CSV and together in one workbook
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
