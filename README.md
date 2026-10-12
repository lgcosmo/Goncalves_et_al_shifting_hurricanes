# Hurricane regimes and multilayer mutualistic networks in the Lesser Antilles

This repository contains the data, simulations and figures for the manuscript:

**Shifting hurricane regimes threaten island mutualisms**

by Fernando Gonçalves, Leandro G. Cosmo, Max Vollstädt, Fabio Tarazona, Seokmin Kim, Tianying Zhang, Hanna Welzel, Emmeli Claré, Andreas Norrild, Zhao Li, Brandon Whitley, Laura Jones, Natasha de Vere, Mark Hulme, Christopher Searcy, Mauro Galetti, Benno Simmons, Christopher Kaiser-Bunbury, Bo Dalsgaard, and Jordi Bascompte.

Pollination and seed dispersal by birds sustain the plant communities of islands. Hurricanes disturb these mutualisms, and shifting climate regimes are expected to bring more frequent and more intense storms to the Caribbean. This project asks how shifting hurricane regimes affect the joint persistence of pollination and seed dispersal on islands. In this study, we combined field-sampled multilayer bird-pollination and seed-dispersal networks from 11 islands of the Lesser Antilles with simulations of their population dynamics under empirical and escalated hurricane regimes.

## Biological system

We study one multilayer network per island, built from the field records of birds and plants, as follows:

| Community component | Representation in the data |
|---|---|
| Islands | 11 islands, coded `01_trinidad` to `11_saba` in every file |
| Pollination layer | Binary plant-bird matrix from the pollen carried by each sampled bird (`src/multi_pol/`) |
| Seed-dispersal layer | Binary plant-bird matrix from the seeds in the fecal samples of each bird (`src/multi_sd/`) |
| Dual-role birds | Birds that pollinate and disperse seeds, the connectors that couple the two layers |
| Hurricane regime | Hurricane frequency and category probability distribution of each island, estimated from the historical record |

The two layers of an island share the same rows and columns, with an empty row or column (row/column of only zeros) for a species recorded in one layer only. This alignment lets every script tell dual-role birds apart from single-role birds. The stored matrices hold the number of records of each interaction, and every analysis reads them as presence or absence, except the proportion of interactions involving dual-role species, which weights each interaction by its records.

## Hurricane regimes

The empirical regime of each island combines two quantities estimated from the historical hurricane record (`src/hurricane_crono.csv`). The hurricane frequency, $\lambda_i$, is the number of recorded hurricanes per year of island $i$'s record (150 to 173 years, depending on the island), the rate of the Poisson process of the simulations. The hurricane intensity is the probability mass function (PMF) of Saffir-Simpson categories among the recorded hurricanes, with $p_{c,i}$ the probability that a hurricane striking island $i$ belongs to category $c$. Each category carries its own damage fraction, $d_c$, the share of abundance that one hurricane of that category removes from every species:

| Saffir-Simpson category | Damage fraction $d_c$ |
|---|---:|
| Category 1 | 0.1 |
| Category 2 | 0.2 |
| Category 3 | 0.3 |
| Category 4 | 0.6 |
| Category 5 | 0.8 |

The manuscript examines two escalation scenarios. Under increased hurricane frequency, every island keeps its empirical PMF and $\lambda_i$ varies from 0.01 to 0.80 hurricanes per year in steps of 0.01. Under increased hurricane intensity, every island keeps its empirical $\lambda_i$ and a shift $s$ moves probability mass from weak (categories 1 to 3) to strong (categories 4 and 5) hurricanes, from 0 to 1 in steps of 0.05. The shifted mass is $\Delta s = \min(s, p_{\mathrm{low}})$, categories 1 to 3 are rescaled by $(p_{\mathrm{low}} - \Delta s)/p_{\mathrm{low}}$ and categories 4 and 5 by $(p_{\mathrm{high}} + \Delta s)/p_{\mathrm{high}}$, where $p_{\mathrm{low}}$ and $p_{\mathrm{high}}$ sum the empirical probabilities of the two groups.

## Population dynamics

Species abundances follow the continuous-time model of the Methods. Briefly, the abundance of species $i$ changes with its intrinsic growth rate, intraspecific competition, intra-guild competition, and the mutualistic benefits of its pollination and seed-dispersal partners. The benefits saturate through a Type II functional response, so abundant partners bring diminishing returns:

$$
\frac{dN_i}{dt} = N_i \left( r_i - N_i + \sum_{j \neq i} c_{ij} N_j + \frac{\sum_{j} \gamma_{ij}^{(P)} N_j + \sum_{k} \gamma_{ik}^{(D)} N_k}{1 + h \left( \sum_{j} \gamma_{ij}^{(P)} N_j + \sum_{k} \gamma_{ik}^{(D)} N_k \right)} \right).
$$

Here, $r_i$ is the intrinsic growth rate of species $i$, $c_{ij}$ is the strength of intra-guild competition, and $h$ is the handling time of the functional response. In the simulations, $c_{ij}$ equals the base strength $c$ divided by the number of other guild members present in the layer, summed over the two layers. The terms $\gamma_{ij}^{(P)}$ and $\gamma_{ik}^{(D)}$ are the specific mutualistic benefits species $i$ receives from partner $j$ via pollination and from partner $k$ via seed dispersal:

$$
\gamma_{ij}^{(P)} = \frac{\gamma_0^{(P)} a_{ij}^{(P)}}{\left(k_i^{(P)}\right)^{\delta}}, \qquad
\gamma_{ik}^{(D)} = \frac{\gamma_0^{(D)} a_{ik}^{(D)}}{\left(k_i^{(D)}\right)^{\delta}},
$$

where $a_{ij}^{(P)}$ and $a_{ik}^{(D)}$ are the binary adjacency matrices of the two layers, $k_i^{(P)}$ and $k_i^{(D)}$ are the numbers of partners of species $i$ in each layer, and $\gamma_0^{(P)}$ and $\gamma_0^{(D)}$ are the base mutualistic strengths. The exponent $\delta$ sets the generalist trade-off, so that the higher the value, the more specialists receive a greater per-interaction benefit than generalists.

Each simulation run consists of two steps. First, the population dynamics is simulated until the abundances of all species reach a steady state, $N_i^{\ast}$, without any hurricanes striking. Second, hurricanes strike at random times for $t_{\mathrm{years}} = 300$ years. The waiting time between hurricanes is exponential with mean $1/\lambda_i$, as in a Poisson process, and each hurricane removes a fraction $d_c$ of the abundance of each species, such that $N_i(t^{+}) = N_i(t)(1 - d_c)$. The fraction removed depends on the category of the hurricane, which occurs with probability given by the PMF of hurricanes in each island. A species whose abundance falls below its extinction threshold, defined as $\epsilon_i = 0.01 N_i^{\ast}$, goes extinct. Every scenario runs 1,000 replicates with the following parameters:

| Parameter | Value | Meaning |
|---|---:|---|
| $\gamma_0^{(P)}$ | 0.5 | Base mutualistic strength of pollination (`m` in the code) |
| $\gamma_0^{(D)}$ | 0.5 | Base mutualistic strength of seed dispersal (`d` in the code) |
| $c$ | -0.1 | Intra-guild competition strength (`c` in the code) |
| $h$ | 0.5 | Handling time of the Type II functional response (`ht` in the code) |
| $\delta$ | 0.5 | Exponent of the generalist trade-off |
| $r_i$ | $U[0.0, 0.5]$ | Intrinsic growth rate, drawn at random for each replicate |
| $\epsilon_i$ | $0.01 N_i^{\ast}$ | Extinction threshold of species $i$ |

## Network structure and species roles

We summarized the structure of each island network with five descriptors. Three of those measure the integration of the networks, namely the leading eigenvalue of the multilayer network, of the pollination layer and of the seed-dispersal layer. The other two, the proportion of dual-role species and the proportion of interactions involving dual-role species, measure how strongly dual-role species couple the layers. Then, a principal component analysis reduces the five descriptors to two axes. All five descriptors load positively on the first axis, which we defined as the global network integration axis. On the second axis the dual-role descriptors load negatively, whereas the leading eigenvalue of the pollination layer loads positively. As such, higher scores in this PC can be interpreted as networks that are species-rich and connected around the pollination layer (the pollination topological dominance axis). Linear models then regress each PC axis on island biogeography descriptors, specifically isolation, area, elevation and island type (continental or oceanic).

At the species level, three role metrics describe how each plant and bird distributes its interactions. The within-layer degree z-scores compare the degree of species $i$ in each layer with the degrees of its guild on the same island:

$$
z_i^{(P)} = \frac{k_i^{(P)} - \bar{k}^{(P)}}{\sigma_k^{(P)}}, \qquad
z_i^{(D)} = \frac{k_i^{(D)} - \bar{k}^{(D)}}{\sigma_k^{(D)}},
$$

where $\bar{k}$ and $\sigma_k$ are the mean and standard deviation of the degrees of the guild in that layer. The multilayer participation coefficient measures how evenly species $i$ splits its interactions between the two layers:

$$
\mathrm{MPC}_i = 1 - \left[ \left( \frac{k_i^{(P)}}{k_i^{(P)} + k_i^{(D)}} \right)^2 + \left( \frac{k_i^{(D)}}{k_i^{(P)} + k_i^{(D)}} \right)^2 \right],
$$

from 0, when a species interacts in one layer only, to 0.5, when it splits its interactions equally between the two layers. Thresholds on the z-scores ($z > 0.75$) and the MPC ($\mathrm{MPC} > 0.3$) sort the species into Multilayer Hubs, Pollination Hubs, Dispersal Hubs, Multilayer Connectors and Peripheral species.

## Species fragility and the structural equation model

Species fragility is the normalized area under the extinction curve (AUC), the curve that relates the fraction of simulations in which a species went extinct, $y_t$, to the perturbation level, $x_t$, hurricane frequency in one scenario and the intensity shift in the other:

$$
\mathrm{AUC} = \frac{1}{x_n - x_1} \sum_{t=1}^{n-1} \frac{y_t + y_{t+1}}{2} \, \Delta x_t,
$$

where $\Delta x_t = x_{t+1} - x_t$ and the normalization by the range of the perturbation axis makes the AUC comparable across scenarios. Higher AUC means greater fragility. A piecewise structural equation model (pSEM) then connects the fragility of bird species to predictors at three levels of organization. The first one consists of island biogeography, namely island type, area, elevation and isolation. The second level is at the network-scale, within which network structure is encoded through the two PCA axes (global integration and pollination topological dominance). Finally, the species level, containing the MPC, the two within-layer degree z-scores and the body mass of birds. Directed paths run from all three levels to each fragility response, from island biogeography to each PCA axis, and from body mass to each role metric. Every variable is standardized before fitting, so the estimates are standardized path coefficients, and Fisher's C tests the overall fit.

## Null models that deliberately decouple the two layers

We used four null models to test whether the coupling of the two layers by dual-role birds shapes the extinctions that hurricanes cause. Every null network progressively decouples the pollination and seed-dispersal layers of an empirical network while holding its remaining structure fixed. Specifically, species richness, the number of interactions in each layer, and the identity and degree of the plants in each layer are held fixed. A rewired interaction keeps its plant and its layer, and the receiving bird interacts exclusively in that layer, does not yet interact with that plant, and is never a hummingbird receiving a fruit.

| Null model | Decoupling level | Share decoupled | Outcome |
|---|---|---|---|
| I | Bird species | All connectors | The layers separate |
| I-A | Bird species | Half of the connectors | Dual-role birds remain |
| II | Interactions | All dual-role interactions | The layers separate |
| II-B | Interactions | Half of the dual-role interactions | Dual-role birds remain |

Every island receives 100 null networks per null model, each built from its own seed and recorded in a manifest (`src/null_<model>/network_manifest.csv`). The null networks are then simulated under the empirical hurricane regime of their island, and the extinctions of the empirical networks are compared with those of their null counterparts through z-scores, $z = (x - \mu)/\sigma$, where $x$ is the empirical fraction of extinctions and $\mu$ and $\sigma$ are its mean and standard deviation across the 100 null networks.

## Manuscript figures

Figure 1 of the manuscript is a conceptual figure produced outside this repository, so the figure scripts follow the manuscript numbering from Figure 2 on.

| Figure | Biological question | Active script | Main output |
|---|---|---|---|
| Figure 2 | How are the pollination and seed-dispersal layers of the 11 islands structured, and how does pollination topological dominance change with isolation? | `scripts/05_visualization/figure_2.R` | `figures/figure_2.pdf` |
| Figure 3 | How do hurricane regimes differ among islands, and how many species go extinct as hurricane frequency and intensity escalate? | `scripts/05_visualization/figure_3.R` | `figures/figure_3.pdf` |
| Figure 4 | Which island, network and species properties shape the fragility of bird species to hurricanes? | `scripts/05_visualization/figure_4.R` | `figures/figure_4.pdf` |
| Extended Data Figure 1 | Does the coupling of the two layers by dual-role birds change how many species hurricanes drive extinct? | `scripts/05_visualization/extended_data_figure_1.R` | `figures/extended_data_figure_1.pdf` |

## Repository structure

```text
Goncalves_et_al_shifting_hurricanes/
|-- src/
|   |-- pl_barcode/                           # Pollination records of each island (field data)
|   |-- dataseed_carib.csv                    # Fecal sample records of seeds with counts (field data)
|   |-- dataseed_carib_noweight.csv           # Fecal sample records of seeds, binarized
|   |-- multi_pol/                            # Pollination layer of each island (B_<island>.csv)
|   |-- multi_sd/                             # Seed-dispersal layer of each island (S_<island>.csv)
|   |-- null_<model>/                         # Null networks and summary data of models I, I-A, II, II-B
|   |-- hurricane_crono.csv                   # Historical hurricane record of the 11 islands
|   |-- empirical_lambdas.csv                 # Empirical hurricane frequency of each island
|   |-- empirical_pmfs.csv                    # Empirical category PMF of each island
|   |-- Island_index.csv                      # Island biogeography
|   |-- bodymass.csv                          # Body mass of the bird species
|   |-- hummingbirds.csv                      # Hummingbird list of the null models
|   |-- shape_islands/                        # Shapefiles of the islands (Figure 2)
|   |-- island_labels.R                       # Island names and order of the figures
|   |-- figure_helpers.R                      # Figure saver
|   |-- aux_functions.R                       # Simulation summaries and AUC helpers
|   |-- main_functions.jl                     # Population dynamics, hurricanes, simulation wrapper
|   |-- aux_functions.jl                      # Network loader, regimes and interaction matrix
|   |-- null_io_functions.jl                  # Reading and writing of the null matrices
|   |-- null_decoupling_functions.jl          # Rewiring of the four null models
|   `-- null_decoupling_checks.jl             # Checks of every null network
|-- scripts/
|   |-- 01_network_preparation/
|   |   |-- build_networks_multilayer.R       # Binary layers of each island from the field records
|   |   `-- build_null_decoupling.jl          # Null networks of the four null models
|   |-- 02_topology_analysis/
|   |   |-- calculate_macro_structure.R       # Structure descriptors of each island network
|   |   |-- calculate_pca.R                   # PCA axes of network structure
|   |   |-- calculate_species_roles.R         # Role metrics and functional roles of every species
|   |   `-- linear_model_pca.R                # PCA axes against island biogeography
|   |-- 03_simulations/
|   |   |-- empirical_lambda_pmf.R            # Hurricane regime of each island from the record
|   |   |-- main_simulations_empirical.jl     # Frequency sweep of the empirical networks
|   |   |-- main_simulations_empirical_pmf.jl # Intensity sweep of the empirical networks
|   |   `-- main_simulations_null_decoupling.jl # Null networks under the empirical regimes
|   |-- 04_results_processing/
|   |   |-- survival_curves.R                 # Extinction curves of the two sweeps
|   |   |-- build_island_data.R               # Island-level table of the pSEM and figures
|   |   |-- calculate_auc_extinction.R        # Species fragility as the AUC of the curves
|   |   |-- process_null_simulations.R        # Extinctions of the null networks
|   |   |-- null_model_comparison.R           # Empirical versus null extinctions
|   |   `-- structural_equation_model.R       # pSEM of species fragility
|   `-- 05_visualization/
|       |-- figure_2.R                        # Map, networks and PCA inset
|       |-- figure_3.R                        # Hurricane regimes and extinction curves
|       |-- figure_4.R                        # pSEM path diagram
|       |-- extended_data_figure_1.R          # Empirical versus null z-scores
|       `-- diagnostic/plot_null_networks.R   # Diagnostic plots of the null networks
|-- data/                                     # Summaries of the simulations (main_summary/)
|-- tables/                                   # pSEM tables (sem_results_clean.xlsx and its CSVs)
`-- figures/                                  # Manuscript figures
```

## Running the analyses

The workflow combines R and Julia. The R scripts uses the library `here` to define paths, so run these scripts from the project root folder (opening the .Rproj file). THe required R packages are as follows:

```r
install.packages(c(
  "here",
  "tidyverse",
  "data.table",
  "stringr",
  "future.apply",
  "FactoMineR",
  "piecewiseSEM",
  "lme4",
  "ggpubr",
  "patchwork",
  "sf",
  "ggnetwork",
  "network",
  "ggforce",
  "png"
))
```

The Julia scripts activate the project with DrWatson, and `Project.toml` and `Manifest.toml` pin their packages (Julia 1.12.5). Instantiate the environment once:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Within this project, we already provide the network matrices, the hurricane regimes, the summarized outputs in `data/main_summary/` and the pSEM tables in `tables/`, which together reproduce every figure without rerunning the simulations. The full pipeline rebuilds all of them from the field records, in the following order.

First, build the networks, the hurricane regimes and the null networks:

```r
source("scripts/01_network_preparation/build_networks_multilayer.R")
source("scripts/03_simulations/empirical_lambda_pmf.R")
```

```bash
julia scripts/01_network_preparation/build_null_decoupling.jl
```

Next, compute the structure of the networks and the roles of their species:

```r
source("scripts/02_topology_analysis/calculate_macro_structure.R")
source("scripts/02_topology_analysis/calculate_pca.R")
source("scripts/02_topology_analysis/calculate_species_roles.R")
```

Then run the simulations. This can take some time, so we already provide summary files of simulation results to reproduce the figures. Nevertheless, simulations can be run through the following scripts:

```bash
julia scripts/03_simulations/main_simulations_empirical.jl
julia scripts/03_simulations/main_simulations_empirical_pmf.jl
julia scripts/03_simulations/main_simulations_null_decoupling.jl
```

The three scripts start their own worker processes (`addprocs` at the top of each script; adjust the number to the computer number of cores) and run 880, 231 and 4,400 scenarios of 1,000 replicates. The raw simulation files are saved within `data/main_sims_lambda/`, `data/main_sims_pmf/` and `data/sims/`. As we mentioned above, they are too large for the repository, so only their summaries in `data/main_summary/` are provided here.

Next, summarize the simulations and fit the models. `linear_model_pca.R` uses the island dataset, so it runs after `build_island_data.R`:

```r
source("scripts/04_results_processing/survival_curves.R")
source("scripts/04_results_processing/build_island_data.R")
source("scripts/02_topology_analysis/linear_model_pca.R")
source("scripts/04_results_processing/calculate_auc_extinction.R")
source("scripts/04_results_processing/process_null_simulations.R")
source("scripts/04_results_processing/null_model_comparison.R")
source("scripts/04_results_processing/structural_equation_model.R")
```

Finally, generate the manuscript figures:

```r
source("scripts/05_visualization/figure_2.R")
source("scripts/05_visualization/figure_3.R")
source("scripts/05_visualization/figure_4.R")
source("scripts/05_visualization/extended_data_figure_1.R")
```

The summaries are saved in `data/main_summary/`, the pSEM tables in `tables/`, and every figure in `figures/` as `.png`, `.pdf` and `.svg`.
