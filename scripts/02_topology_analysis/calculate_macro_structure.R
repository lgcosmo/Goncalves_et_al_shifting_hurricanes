#--------------------------------------------------------------#
#-Topology of the multilayer networks--------------------------#
#--------------------------------------------------------------#

# Compute the topology descriptors of each island's multilayer
# bird-pollination and seed-dispersal network. Five of them enter the PCA of
# calculate_pca.R (Extended Data Table 2) and, through its two axes, the Fig.
# 2 inset and the pSEM of Fig. 4. Three are the leading eigenvalue of the
# pollination layer, of the seed-dispersal layer and of the multilayer
# network, and two are the proportion of dual-role species and the
# proportion of interactions involving dual-role species. The leading
# eigenvalue predicts the resilience of a network to pulse perturbations such
# as hurricanes, and it correlates with nestedness, modularity, connectance
# and species richness. To document this correlation (Extended Data), the
# script also computes species richness, connectance, nestedness (NODF) and
# modularity (Barber's Q) for each layer and for the multilayer network, with
# z-scores of NODF and Q against null networks that keep the degree
# distributions of the empirical network. The multilayer network is the union
# of the two layers. A plant and a bird interact in it when they interact in
# either layer. The proportion of interactions involving dual-role species
# weights each interaction by its number of records, as stored in the
# matrices; every other descriptor counts interactions as presence or
# absence.
#
# Inputs:  src/multi_pol/B_<island>.csv and src/multi_sd/S_<island>.csv,
#          written by build_networks_multilayer.R
# Outputs: data/main_summary/macro_structure.csv, one row per island, read by
#          calculate_pca.R and build_island_data.R

library(here)
library(tidyverse)

#--------------------------------------------------------------#
#-Computing Nestedness (NODF)----------------------------------#
#--------------------------------------------------------------#

# Nestedness (NODF) of a binary incidence matrix, following Almeida-Neto et
# al. (2008). For each pair of plants, and for each pair of birds, the species
# with fewer partners contributes the share of its partners that it shares
# with the species with more partners. Pairs with the same degree contribute
# zero, so that NODF rewards decreasing fill and not only overlap. The metric
# is the mean contribution over all pairs, scaled from 0 (no specialist
# shares partners with a generalist) to 100 (the partners of every specialist
# form a subset of those of every generalist).
# Arguments: matrix is the incidence matrix with plants in rows and birds in
# columns, binarized before use; sort_for_nodf orders rows and columns by
# decreasing degree before the pairs are compared. Returns NODF between 0 and
# 100, or 0 for a matrix without interactions.
compute_nodf <- function(matrix, sort_for_nodf = TRUE) {
  B <- (matrix > 0) * 1
  if (sum(B) == 0) {
    return(0)
  }

  nrow_val <- nrow(B)
  ncol_val <- ncol(B)

  if (sort_for_nodf) {
    row_deg <- rowSums(B)
    col_deg <- colSums(B)
    B <- B[order(row_deg, decreasing = TRUE), order(col_deg, decreasing = TRUE), drop = FALSE]
  }

  row_deg <- rowSums(B)
  col_deg <- colSums(B)

  row_sum <- 0
  if (nrow_val > 1) {
    for (i in 1:(nrow_val - 1)) {
      for (j in (i + 1):nrow_val) {
        ki <- row_deg[i]
        kj <- row_deg[j]
        min_deg <- min(ki, kj)
        if (ki != kj && min_deg > 0) {
          overlap <- sum(B[i, ] & B[j, ])
          row_sum <- row_sum + overlap / min_deg
        }
      }
    }
  }

  col_sum <- 0
  if (ncol_val > 1) {
    for (i in 1:(ncol_val - 1)) {
      for (j in (i + 1):ncol_val) {
        ki <- col_deg[i]
        kj <- col_deg[j]
        min_deg <- min(ki, kj)
        if (ki != kj && min_deg > 0) {
          overlap <- sum(B[, i] & B[, j])
          col_sum <- col_sum + overlap / min_deg
        }
      }
    }
  }

  num_row_pairs <- if (nrow_val > 1) nrow_val * (nrow_val - 1) / 2 else 0
  num_col_pairs <- if (ncol_val > 1) ncol_val * (ncol_val - 1) / 2 else 0
  denom <- num_row_pairs + num_col_pairs

  if (denom == 0) {
    return(0)
  }

  return(100 * (row_sum + col_sum) / denom)
}

#--------------------------------------------------------------#
#-Bascompte Null Model 2---------------------------------------#
#--------------------------------------------------------------#

# One null network under null model 2 of Bascompte et al. (2003). Each plant
# and bird interact with a probability equal to the mean of the plant's fill
# (its degree over the number of birds) and the bird's fill (its degree over
# the number of plants). The null networks therefore keep, on average, the
# degree of each plant and of each bird, and with them the connectance, so
# that the z-scores of NODF and Q measure structure beyond what the degree
# distributions impose.
# Arguments: matrix is the empirical incidence matrix, binarized before use.
# Returns a binary incidence matrix of the same dimensions.
bascompte_null_model_2 <- function(matrix) {
  B <- (matrix > 0) * 1
  rows <- nrow(B)
  cols <- ncol(B)

  if (sum(B) == 0) {
    return(matrix(0, nrow = rows, ncol = cols))
  }

  p_row <- rowSums(B) / cols
  p_col <- colSums(B) / rows

  P <- 0.5 * (matrix(p_row, nrow = rows, ncol = cols) +
    matrix(p_col, nrow = rows, ncol = cols, byrow = TRUE))

  U <- matrix(runif(rows * cols), nrow = rows, ncol = cols)
  null_matrix <- (U < P) * 1
  return(null_matrix)
}

#--------------------------------------------------------------#
#-Barber Modularity for Bipartite Networks---------------------#
#--------------------------------------------------------------#

# Barber's modularity Q of a bipartite network under a given division of its
# plants and birds into modules. Standard modularity compares the interactions
# inside modules with those expected among all pairs of species, but a plant
# never interacts with another plant, nor a bird with another bird. Barber's Q
# therefore takes the expectation only between plants and birds. The expected
# number of interactions inside a module is the product of the summed degree
# of its plants and the summed degree of its birds, divided by the number of
# interactions in the network.
# Arguments: B is the binary incidence matrix; row_labels and col_labels give
# the module of each plant and of each bird as positive integers. Returns Q,
# or 0 for a network without interactions.
barber_Q <- function(B, row_labels, col_labels) {
  L <- sum(B)
  if (L == 0) {
    return(0)
  }

  m <- nrow(B)
  n <- ncol(B)

  max_label <- max(c(row_labels, col_labels))

  # Summed degree of the plants (aR) and of the birds (aC) in each module, and
  # the number of interactions inside each module (e)
  aR <- numeric(max_label)
  aC <- numeric(max_label)
  e <- numeric(max_label)

  row_deg <- rowSums(B)
  col_deg <- colSums(B)

  aR_agg <- tapply(row_deg, row_labels, sum)
  aR[as.numeric(names(aR_agg))] <- aR_agg

  aC_agg <- tapply(col_deg, col_labels, sum)
  aC[as.numeric(names(aC_agg))] <- aC_agg

  # An interaction lies inside a module when its plant and its bird share the
  # module
  indices <- which(B > 0, arr.ind = TRUE)
  rows_idx <- indices[, 1]
  cols_idx <- indices[, 2]

  r_mod <- row_labels[rows_idx]
  c_mod <- col_labels[cols_idx]

  same <- (r_mod == c_mod)

  e_agg <- tapply(rep(1, sum(same)), r_mod[same], sum)
  if (length(e_agg) > 0) {
    e[as.numeric(names(e_agg))] <- e_agg
  }

  # Interactions inside each module minus their expectation, summed over the
  # modules and divided by the number of interactions
  Q <- sum(e - (aR * aC) / L) / L
  return(Q)
}

# Renumber the modules as 1, 2, 3 and so on, dropping the labels that no plant
# or bird holds any more, so that the module vectors stay short after moves.
# Arguments: row_labels and col_labels are the current module labels of the
# plants and of the birds. Returns the renumbered labels and the number of
# modules.
compact_relabel <- function(row_labels, col_labels) {
  all_labels <- sort(unique(c(row_labels, col_labels)))
  mapping <- setNames(seq_along(all_labels), all_labels)

  return(list(
    row_labels = mapping[as.character(row_labels)],
    col_labels = mapping[as.character(col_labels)],
    count = length(all_labels)
  ))
}

# Greedy search for the division that maximizes Barber's Q. The search starts
# from a given division, or from every species in its own module. Each plant
# then moves to the module of one of its partners whenever the move increases
# Q, and the birds do the same. Both passes repeat until no species moves or
# max_iter passes have run. The gain of a move follows the bipartite form of
# the modularity increment, which needs only the interactions of the moving
# species with each module and the summed degrees of the two modules
# involved, so that Q is never recomputed from scratch during the search.
# Arguments: matrix is the incidence matrix, binarized before use; init_row
# and init_col give the starting division of the plants and of the birds
# (NULL puts every species in its own module); max_iter caps the number of
# passes. Returns Q and the module labels of the plants and of the birds.
local_move_barber <- function(matrix, init_row = NULL, init_col = NULL, max_iter = 50) {
  B <- (matrix > 0) * 1
  m <- nrow(B)
  n <- ncol(B)
  L <- sum(B) # number of interactions

  # A network without interactions has no modules and Q = 0
  if (L == 0) {
    return(list(Q = 0, r = 1:m, c = (m + 1):(m + n)))
  }

  # Without a starting division, every plant and every bird begins in its own
  # module, with the bird labels following the plant labels
  if (is.null(init_row)) {
    row_labels <- 1:m
  } else {
    row_labels <- init_row
  }

  if (is.null(init_col)) {
    col_labels <- (m + 1):(m + n)
  } else {
    col_labels <- init_col
  }

  row_deg <- rowSums(B)
  col_deg <- colSums(B)

  iter <- 0
  improved <- TRUE

  # Renumber the modules so that the label vectors index them directly
  cres <- compact_relabel(row_labels, col_labels)
  row_labels <- cres$row_labels
  col_labels <- cres$col_labels

  # Each pass moves the plants and then the birds. The search stops when a
  # pass moves no species.
  while (improved && iter < max_iter) {
    iter <- iter + 1
    improved <- FALSE

    max_lab <- max(c(row_labels, col_labels))
    aR <- numeric(max_lab)
    aC <- numeric(max_lab)

    # Summed degree of the plants (aR) and of the birds (aC) in each module
    aR_agg <- tapply(row_deg, row_labels, sum)
    aR[as.numeric(names(aR_agg))] <- aR_agg

    aC_agg <- tapply(col_deg, col_labels, sum)
    aC[as.numeric(names(aC_agg))] <- aC_agg

    # First pass. Each plant considers the modules of its partners.
    for (i in 1:m) {
      ki <- row_deg[i] # degree of plant i
      s <- row_labels[i] # current module of plant i

      neigh_cols <- which(B[i, ] > 0)
      neigh_mods <- unique(c(s, col_labels[neigh_cols]))

      best_t <- s
      best_delta <- 0

      k_i_to_mod <- numeric(max_lab)

      # Interactions of plant i with the birds of each module
      if (length(neigh_cols) > 0) {
        tab <- table(col_labels[neigh_cols])
        k_i_to_mod[as.numeric(names(tab))] <- as.numeric(tab)
      }

      sumA_s <- k_i_to_mod[s]

      # Gain in Q from moving plant i to module t. The gain is the change in
      # interactions inside modules minus the change in their expectation,
      # which for a plant depends on the summed bird degrees of the two
      # modules.
      for (t in neigh_mods) {
        if (t == s) next

        sumA_t <- k_i_to_mod[t]

        delta_num <- (sumA_t - sumA_s) - ki * (aC[t] - aC[s]) / L
        delta <- delta_num / L

        # A small tolerance keeps ties from moving species back and forth
        if (delta > best_delta + 1e-10) {
          best_delta <- delta
          best_t <- t
        }
      }

      if (best_t != s) {
        row_labels[i] <- best_t

        # The summed plant degrees of the two modules follow the move
        aR[s] <- aR[s] - ki
        aR[best_t] <- aR[best_t] + ki
        improved <- TRUE
      }
    }

    # Drop the modules left empty and renumber
    cres <- compact_relabel(row_labels, col_labels)
    row_labels <- cres$row_labels
    col_labels <- cres$col_labels

    # Summed degrees per module under the new labels
    max_lab <- max(c(row_labels, col_labels))
    aR <- numeric(max_lab)
    aC <- numeric(max_lab)
    aR_agg <- tapply(row_deg, row_labels, sum)
    aR[as.numeric(names(aR_agg))] <- aR_agg
    aC_agg <- tapply(col_deg, col_labels, sum)
    aC[as.numeric(names(aC_agg))] <- aC_agg

    # Second pass. Each bird considers the modules of its partners, with the
    # same rules as the plants.
    for (j in 1:n) {
      kj <- col_deg[j]
      s <- col_labels[j]

      neigh_rows <- which(B[, j] > 0)
      neigh_mods <- unique(c(s, row_labels[neigh_rows]))

      best_t <- s
      best_delta <- 0

      k_j_to_mod <- numeric(max_lab)
      if (length(neigh_rows) > 0) {
        tab <- table(row_labels[neigh_rows])
        k_j_to_mod[as.numeric(names(tab))] <- as.numeric(tab)
      }

      sumA_s <- k_j_to_mod[s]

      for (t in neigh_mods) {
        if (t == s) next

        sumA_t <- k_j_to_mod[t]
        # For a bird, the expectation depends on the summed plant degrees of
        # the two modules
        delta_num <- (sumA_t - sumA_s) - kj * (aR[t] - aR[s]) / L
        delta <- delta_num / L

        if (delta > best_delta + 1e-10) {
          best_delta <- delta
          best_t <- t
        }
      }

      if (best_t != s) {
        col_labels[j] <- best_t
        aC[s] <- aC[s] - kj
        aC[best_t] <- aC[best_t] + kj
        improved <- TRUE
      }
    }

    cres <- compact_relabel(row_labels, col_labels)
    row_labels <- cres$row_labels
    col_labels <- cres$col_labels
  }

  # Q of the final division
  Q <- barber_Q(B, row_labels, col_labels)
  return(list(Q = Q, row_labels = row_labels, col_labels = col_labels))
}

# Simulated annealing to leave the local optimum found by the greedy search.
# At each step one plant or bird, drawn at random, moves to a random module,
# either an existing one or a new one. A move that lowers Q is accepted with
# probability exp(delta / T), where the temperature T shrinks by the cooling
# factor at every step, so that the search explores early and settles late.
# The best division met along the way is kept.
# Arguments: matrix is the incidence matrix, binarized before use; init_row
# and init_col give the starting division of the plants and of the birds;
# steps is the number of proposed moves; T0 is the starting temperature;
# cooling multiplies the temperature at each step. Returns the highest Q
# found and its module labels.
sa_barber <- function(matrix, init_row, init_col, steps = 2000, T0 = 0.01, cooling = 0.999) {
  B <- (matrix > 0) * 1
  row_labels <- init_row
  col_labels <- init_col

  current_Q <- barber_Q(B, row_labels, col_labels)
  best_Q <- current_Q
  best_row <- row_labels
  best_col <- col_labels

  T_curr <- T0
  all_labels <- unique(c(row_labels, col_labels))

  for (i in 1:steps) {
    is_row <- runif(1) < 0.5

    if (is_row) {
      idx <- sample(length(row_labels), 1)
      old_lab <- row_labels[idx]
      possible_labels <- c(all_labels, max(all_labels) + 1)
      new_lab <- sample(possible_labels, 1)

      row_labels[idx] <- new_lab
      new_Q <- barber_Q(B, row_labels, col_labels)
    } else {
      idx <- sample(length(col_labels), 1)
      old_lab <- col_labels[idx]
      possible_labels <- c(all_labels, max(all_labels) + 1)
      new_lab <- sample(possible_labels, 1)

      col_labels[idx] <- new_lab
      new_Q <- barber_Q(B, row_labels, col_labels)
    }

    delta <- new_Q - current_Q

    # Metropolis rule. A gain is always accepted, and a loss with a
    # probability that falls as the temperature cools.
    if (delta >= 0 || runif(1) < exp(delta / T_curr)) {
      current_Q <- new_Q
      if (current_Q > best_Q) {
        best_Q <- current_Q
        best_row <- row_labels
        best_col <- col_labels
      }
      all_labels <- unique(c(row_labels, col_labels))
    } else {
      if (is_row) {
        row_labels[idx] <- old_lab
      } else {
        col_labels[idx] <- old_lab
      }
    }

    T_curr <- T_curr * cooling
  }

  return(list(Q = best_Q, row_labels = best_row, col_labels = best_col))
}

# Modularity of a network, as the highest Q found by a greedy search, a
# simulated annealing started from its division, and a second greedy search
# started from the annealed division. The annealing frees the search from the
# local optimum of the first pass, and the last pass refines the division
# that the annealing found.
# Arguments: matrix is the incidence matrix; sa_steps is the number of
# annealing steps. Returns Q.
hybrid_barber_modularity <- function(matrix, sa_steps = 2000) {
  res1 <- local_move_barber(matrix, max_iter = 50)
  res2 <- sa_barber(matrix, res1$row_labels, res1$col_labels, steps = sa_steps)
  res3 <- local_move_barber(matrix, init_row = res2$row_labels, init_col = res2$col_labels, max_iter = 50)

  best_Q <- max(res1$Q, res2$Q, res3$Q)
  return(best_Q)
}

#--------------------------------------------------------------#
#-Descriptors of one network-----------------------------------#
#--------------------------------------------------------------#

# Descriptors of one network, either a layer or the multilayer network.
# Species richness counts the plants and birds of the network, connectance
# divides its interactions by the possible plant and bird pairs, and the
# leading eigenvalue is the largest singular value of the incidence matrix,
# which equals the leading eigenvalue of the adjacency matrix of the
# bipartite network. Nestedness (NODF) and modularity (Barber's Q) are also
# standardized as z-scores against null_reps null networks of Bascompte's
# null model 2, which keep the degree distributions of the empirical
# network, so that nestedness and modularity can be compared across islands
# that differ in size and connectance. The annealing of each null network
# runs fewer steps than that of the empirical network to keep the running
# time low.
# Arguments: matrix is the incidence matrix of the species with at least one
# interaction, binarized before use; null_reps is the number of null
# networks. Returns a list with the richness, connectance, lambda1, the raw
# and standardized NODF and the raw and standardized Q, or NULL for a network
# without interactions.
analyze_network <- function(matrix, null_reps = 100) {
  B <- (matrix > 0) * 1
  if (sum(B) == 0) {
    return(NULL)
  }

  nodf_obs <- compute_nodf(B)
  q_obs <- hybrid_barber_modularity(B, sa_steps = 1000)

  # Leading eigenvalue, as the largest singular value of the incidence matrix
  sv <- svd(B)$d
  lambda1 <- sv[1]

  richness <- nrow(B) + ncol(B)
  connectance <- sum(B) / (nrow(B) * ncol(B))

  # NODF and Q of each null network, with a shorter annealing
  nodf_nulls <- numeric(null_reps)
  q_nulls <- numeric(null_reps)

  for (i in 1:null_reps) {
    B_null <- bascompte_null_model_2(B)
    nodf_nulls[i] <- compute_nodf(B_null)
    q_nulls[i] <- hybrid_barber_modularity(B_null, sa_steps = 200)
  }

  # z-scores of the empirical values among the null networks, set to zero
  # when the null networks do not vary
  nodf_mean <- mean(nodf_nulls)
  nodf_sd <- sd(nodf_nulls)
  nodf_z <- if (nodf_sd > 0) (nodf_obs - nodf_mean) / nodf_sd else 0

  q_mean <- mean(q_nulls)
  q_sd <- sd(q_nulls)
  q_z <- if (q_sd > 0) (q_obs - q_mean) / q_sd else 0

  return(list(
    richness = richness,
    connectance = connectance,
    lambda1 = lambda1,
    nodf_raw = nodf_obs,
    nodf_z = nodf_z,
    q_raw = q_obs,
    q_z = q_z
  ))
}

#--------------------------------------------------------------#
#-Dual-role species and their interactions---------------------#
#--------------------------------------------------------------#

# Coupling of the two layers of an island by its dual-role species. A
# dual-role species interacts in both layers. Such a bird both pollinates
# flowers and disperses seeds, and such a plant is both pollinated and
# dispersed by birds. The two descriptors are the proportion of dual-role
# species among the species with at least one interaction, and the proportion
# of all interactions, over both layers, that involve a dual-role plant or
# bird. The second sums the records of each interaction, as stored in the
# matrices, so that an interaction recorded many times weighs more. Both
# enter the PCA of calculate_pca.R.
# Arguments: mat_pol and mat_sd are the incidence matrices of the two layers
# as stored, with the number of records of each interaction, the same plants
# in rows and the same birds in columns. Returns a list with
# prop_dual_species and prop_dual_interactions.
dual_role_descriptors <- function(mat_pol, mat_sd) {

  # Species present in each layer, and in both
  plants_pol <- rowSums(mat_pol) > 0
  plants_sd  <- rowSums(mat_sd) > 0
  birds_pol  <- colSums(mat_pol) > 0
  birds_sd   <- colSums(mat_sd) > 0

  dual_plants <- plants_pol & plants_sd
  dual_birds  <- birds_pol & birds_sd

  # An interaction involves a dual-role species when its plant or its bird is
  # one
  involves_dual <- outer(dual_plants, dual_birds, FUN="|")

  list(
    prop_dual_species = (sum(dual_plants) + sum(dual_birds)) /
      (sum(plants_pol | plants_sd) + sum(birds_pol | birds_sd)),
    prop_dual_interactions = (sum(mat_pol[involves_dual]) + sum(mat_sd[involves_dual])) /
      (sum(mat_pol) + sum(mat_sd))
  )
}

#--------------------------------------------------------------#
#-Running every island-----------------------------------------#
#--------------------------------------------------------------#

# Read one incidence matrix as stored, with the number of records of each
# interaction, the plants in rows and the birds in columns. analyze_network
# binarizes it; dual_role_descriptors uses the records.
# Arguments: path is the CSV file with the plant names in its first column.
# Returns the matrix with the plant names as row names.
read_layer <- function(path) {
  df <- read_csv(path, show_col_types=FALSE)
  mat <- as.matrix(df[, -1])
  rownames(mat) <- df[[1]]
  return(mat)
}

# Keep the plants and birds with at least one interaction in a matrix, so
# that a bird that only disperses seeds does not count as an isolated species
# of the pollination layer.
# Arguments: mat is the binary matrix. Returns the matrix restricted to its
# connected plants and birds.
active_matrix <- function(mat) {
  mat[rowSums(mat) > 0, colSums(mat) > 0, drop=FALSE]
}

# Compute every descriptor for every island and write them to one table, one
# row per island. The multilayer network is the union of the two layers,
# which share their plants and birds. Takes no arguments, reads src/multi_pol
# and src/multi_sd and writes data/main_summary/macro_structure.csv.
run_all <- function() {

  networks <- list.files(here("src", "multi_pol"), pattern="^B_.*\\.csv$") %>%
    str_remove("^B_") %>%
    str_remove("\\.csv$")

  results <- list()

  for (net in networks) {
    cat("Analyzing network:", net, "\n")

    mat_pol <- read_layer(here("src", "multi_pol", paste0("B_", net, ".csv")))
    mat_sd  <- read_layer(here("src", "multi_sd",  paste0("S_", net, ".csv")))

    # Descriptors of each layer and of the multilayer network, in which a
    # plant and a bird interact when they interact in either layer
    res_pol   <- analyze_network(active_matrix(mat_pol))
    res_sd    <- analyze_network(active_matrix(mat_sd))
    res_multi <- analyze_network(active_matrix((mat_pol | mat_sd) * 1))

    # One row per island, with the descriptors of each network prefixed by
    # its name (pol, sd, multi) and the two dual-role descriptors
    results[[net]] <- bind_cols(
      tibble(network=net),
      as_tibble(res_pol)   %>% rename_with(~paste0("pol_", .x)),
      as_tibble(res_sd)    %>% rename_with(~paste0("sd_", .x)),
      as_tibble(res_multi) %>% rename_with(~paste0("multi_", .x)),
      as_tibble(dual_role_descriptors(mat_pol, mat_sd))
    )
  }

  final_df <- bind_rows(results)
  dir.create(here("data", "main_summary"), showWarnings=FALSE, recursive=TRUE)
  write_csv(final_df, here("data", "main_summary", "macro_structure.csv"))
  message("Saved: macro_structure.csv")
}

run_all()
