# ==============================================================================
# Sensitivity 1: Conservation budget sweep.
# S1 configuration (30% CN floor, no cap) with the 2030 carbon floor imposed.
# Budgets swept at 1.0x, 1.3x, 1.5x, 2.0x, 2.5x current spend (WTP).
# Low budgets may be infeasible; tryCatch records a feasible = FALSE row.
# Data loaded once; one results table saved as xlsx.
# Author: Tijmen Hennekes
# ==============================================================================
budget_multipliers <- c(1.0, 1.3, 1.5, 2.0, 2.5)   # budgets to sweep
# ----------------- 1. LOAD DATA -----------------------------------------
cat("Initializing core optimization environment...\n")
library(terra); library(sf); library(gurobi); library(prioritizr)
library(tibble); library(readxl); library(writexl)
set.seed(24)
out_dir <- "C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/04-outputs"
master_mask <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/forest_binary.tif")
zone_names <- c("Close_to_Nature", "Integrated", "Managed_Forest")
cost_stack <- c(
  terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/cost_cn.tif"),
  terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/cost_ifm.tif"),
  terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/cost_mf.tif"))
names(cost_stack) <- zone_names
cost_stack <- cost_stack * terra::cellSize(cost_stack[[1]], unit = "ha") / 1e6
ifm_cost_factor <- 0.25
cost_stack[["Integrated"]] <- ifm_cost_factor * cost_stack[["Close_to_Nature"]]
scherpenhuijzen_cn_mask <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/locked_in_cn.tif")
carbon_rast <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/carbon_inverse_adjusted.tif")
names(carbon_rast) <- zone_names
cat("Processing ecological distribution layers...\n")
individual_files_cn  <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_cn",  pattern = "\\.tif$", full.names = TRUE)
individual_files_ifm <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_ifm", pattern = "\\.tif$", full.names = TRUE)
individual_files_mf  <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_mf",  pattern = "\\.tif$", full.names = TRUE)
num_features <- length(individual_files_cn)
cat("Species layers found:", num_features, "\n")
stopifnot(length(individual_files_ifm) == num_features,
          length(individual_files_mf)  == num_features)
species_metadata       <- readxl::read_excel("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/IUCN_table_features_filtered.xlsx")
file_codes             <- gsub("_mean\\.tif$", "", basename(individual_files_cn))
species_metadata       <- species_metadata[match(file_codes, species_metadata$Code), ]
explicit_feature_names <- species_metadata$Code
process_zone_stack <- function(file_list, template_mask, names_vector) {
  raw_stack     <- terra::rast(file_list)
  aligned_stack <- terra::resample(raw_stack, template_mask, method = "bilinear")
  aligned_stack <- terra::mask(aligned_stack, template_mask)
  names(aligned_stack) <- names_vector
  return(aligned_stack)
}
species_cn_stack  <- process_zone_stack(individual_files_cn,  master_mask, explicit_feature_names)
species_ifm_stack <- process_zone_stack(individual_files_ifm, master_mask, explicit_feature_names)
species_mf_stack  <- process_zone_stack(individual_files_mf,  master_mask, explicit_feature_names)
cat("Species stacks ready:", num_features, "species x 3 zones\n")
clip_thresh <- 1e-6
species_cn_stack[species_cn_stack   < clip_thresh] <- 0
species_ifm_stack[species_ifm_stack < clip_thresh] <- 0
species_mf_stack[species_mf_stack   < clip_thresh] <- 0
carbon_cn_stack_raw  <- terra::mask(terra::resample(carbon_rast[["Close_to_Nature"]], master_mask, method = "bilinear"), master_mask)
carbon_ifm_stack_raw <- terra::mask(terra::resample(carbon_rast[["Integrated"]],      master_mask, method = "bilinear"), master_mask)
carbon_mf_stack_raw  <- terra::mask(terra::resample(carbon_rast[["Managed_Forest"]],  master_mask, method = "bilinear"), master_mask)
names(carbon_cn_stack_raw) <- names(carbon_ifm_stack_raw) <- names(carbon_mf_stack_raw) <- "Carbon_Sequestration"
carbon_constraint_stack <- c(carbon_cn_stack_raw, carbon_ifm_stack_raw, carbon_mf_stack_raw) / 1e6
names(carbon_constraint_stack) <- zone_names
fm_rast    <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/01_raw_data/fm_swe.tif")
fm_aligned <- terra::mask(terra::resample(fm_rast, master_mask, method = "near"), master_mask)
current_carbon <- terra::ifel(
  fm_aligned == 1, carbon_cn_stack_raw,
  terra::ifel(fm_aligned == 2, carbon_ifm_stack_raw,
              terra::ifel(fm_aligned == 3, carbon_mf_stack_raw, NA)))
current_mgmt_carbon_Mt <- terra::global(current_carbon, "sum", na.rm = TRUE)$sum / 1e6
carbon_target_Mt       <- current_mgmt_carbon_Mt * (47.321 / 43.366)   # 2030 target (+9%)
current_cost_map <- terra::ifel(
  fm_aligned == 1, cost_stack[["Close_to_Nature"]],
  terra::ifel(fm_aligned == 2, cost_stack[["Integrated"]],
              terra::ifel(fm_aligned == 3, cost_stack[["Managed_Forest"]], NA)))
current_mgmt_cost <- terra::global(current_cost_map, "sum", na.rm = TRUE)$sum
feature_biodiversity <- prioritizr::zones(
  species_cn_stack, species_ifm_stack, species_mf_stack,
  zone_names = zone_names, feature_names = explicit_feature_names)
targets_matrix <- matrix(0, nrow = num_features, ncol = 3)
feature_weights <- rep(1, num_features)
for (i in 1:num_features) {
  status <- species_metadata$status[i]
  if (status %in% c("CR", "EN", "VU")) {
    targets_matrix[i, 1] <- 1.00; targets_matrix[i, 2] <- 0.50; feature_weights[i] <- 10
  } else if (status %in% c("NT", "DD")) {
    targets_matrix[i, 1] <- 0.90; targets_matrix[i, 2] <- 0.45; feature_weights[i] <- 5
  } else {
    targets_matrix[i, 1] <- 0.80; targets_matrix[i, 2] <- 0.40; feature_weights[i] <- 1
  }
}
feature_weights_matrix      <- matrix(0, nrow = num_features, ncol = 3)
feature_weights_matrix[, 1] <- feature_weights
feature_weights_matrix[, 2] <- feature_weights * 0.5
scherpenhuijzen_aligned <- terra::resample(scherpenhuijzen_cn_mask, master_mask, method = "near")
n_cells  <- terra::ncell(master_mask)
n_forest <- sum(terra::values(master_mask, mat = FALSE) == 1, na.rm = TRUE)
locked_in_matrix      <- matrix(FALSE, nrow = n_cells, ncol = 3)
locked_in_matrix[, 1] <- terra::values(scherpenhuijzen_aligned, mat = FALSE) == 1
locked_in_matrix[is.na(locked_in_matrix)] <- FALSE
constraint_matrix_cn      <- matrix(0, nrow = n_cells, ncol = 3)
constraint_matrix_cn[, 1] <- 1 / n_forest
cat(sprintf("Data loaded. Carbon target %.2f Mt (+9%%). WTP baseline %.1f. Sweeping %d budgets.\n",
            carbon_target_Mt, current_mgmt_cost, length(budget_multipliers)))
# ============================= 2. SWEEP LOOP ==================================
sweep <- data.frame()
for (bm in budget_multipliers) {
  budget_tag    <- sub("\\.", "p", paste0(bm, "x"))
  solver_budget <- current_mgmt_cost * bm
  cat(sprintf("\n================  BUDGET %s  (%.1f)  ================\n", budget_tag, solver_budget))
  policy_problem <- prioritizr::problem(cost_stack, feature_biodiversity) %>%
    prioritizr::add_min_shortfall_objective(budget = solver_budget) %>%
    prioritizr::add_relative_targets(targets_matrix) %>%
    prioritizr::add_feature_weights(feature_weights_matrix) %>%
    prioritizr::add_locked_in_constraints(locked_in_matrix) %>%
    prioritizr::add_mandatory_allocation_constraints() %>%
    prioritizr::add_linear_constraints(threshold = 0.30, sense = ">=", data = constraint_matrix_cn) %>%
    prioritizr::add_linear_constraints(threshold = carbon_target_Mt, sense = ">=", data = carbon_constraint_stack)
  policy_solution <- tryCatch(
    policy_problem %>%
      prioritizr::add_gurobi_solver(gap = 0.05, time_limit = 10800, verbose = TRUE,
                                    numeric_focus = TRUE, threads = 4, first_feasible = FALSE) %>%
      solve(force = TRUE),
    error = function(e) { cat("  INFEASIBLE at budget", budget_tag, "-", conditionMessage(e), "\n"); NULL })
  if (is.null(policy_solution)) {
    sweep <- rbind(sweep, data.frame(
      budget_multiplier = bm, budget = round(solver_budget, 1), feasible = FALSE,
      CN_pct = NA, IFM_pct = NA, MF_pct = NA,
      mean_conservation_held = NA, pct_of_target = NA, carbon_Mt = NA,
      cost = NA, utilisation_pct = NA))
    next
  }
  solution_map <- prioritizr::category_layer(policy_solution)
  cn_cells  <- sum(terra::values(solution_map, mat = FALSE) == 1, na.rm = TRUE)
  ifm_cells <- sum(terra::values(solution_map, mat = FALSE) == 2, na.rm = TRUE)
  mf_cells  <- sum(terra::values(solution_map, mat = FALSE) == 3, na.rm = TRUE)
  n_pu_sol  <- cn_cells + ifm_cells + mf_cells
  target_cov <- as.data.frame(prioritizr::eval_target_coverage_summary(policy_problem, policy_solution))
  target_cov$zone <- as.character(target_cov$zone)
  cn_rep  <- target_cov[target_cov$zone == "Close_to_Nature", c("feature","absolute_target","absolute_held")]
  ifm_rep <- target_cov[target_cov$zone == "Integrated",      c("feature","absolute_target","absolute_held")]
  names(cn_rep)  <- c("species","abs_target_cn","abs_held_cn")
  names(ifm_rep) <- c("species","abs_target_ifm","abs_held_ifm")
  rep_species <- merge(cn_rep, ifm_rep, by = "species")
  rep_species$conservation_held <- rep_species$abs_held_cn + rep_species$abs_held_ifm
  rep_species$total_target      <- rep_species$abs_target_cn + rep_species$abs_target_ifm
  rep_species$pct_of_target     <- ifelse(rep_species$total_target > 0,
                                          100 * rep_species$conservation_held / rep_species$total_target, NA)
  c_cn  <- terra::global(terra::mask(carbon_cn_stack_raw,  terra::ifel(solution_map == 1, 1, NA)), "sum", na.rm = TRUE)$sum
  c_ifm <- terra::global(terra::mask(carbon_ifm_stack_raw, terra::ifel(solution_map == 2, 1, NA)), "sum", na.rm = TRUE)$sum
  c_mf  <- terra::global(terra::mask(carbon_mf_stack_raw,  terra::ifel(solution_map == 3, 1, NA)), "sum", na.rm = TRUE)$sum
  c_tot <- (c_cn + c_ifm + c_mf) / 1e6
  solution_cost <- terra::global(
    cost_stack[["Close_to_Nature"]] * (solution_map == 1) +
      cost_stack[["Integrated"]]    * (solution_map == 2) +
      cost_stack[["Managed_Forest"]] * (solution_map == 3),
    "sum", na.rm = TRUE)$sum
  sweep <- rbind(sweep, data.frame(
    budget_multiplier = bm, budget = round(solver_budget, 1), feasible = TRUE,
    CN_pct  = round(100 * cn_cells  / n_pu_sol, 1),
    IFM_pct = round(100 * ifm_cells / n_pu_sol, 1),
    MF_pct  = round(100 * mf_cells  / n_pu_sol, 1),
    mean_conservation_held = round(mean(rep_species$conservation_held), 3),
    pct_of_target          = round(mean(rep_species$pct_of_target, na.rm = TRUE), 1),
    carbon_Mt = round(c_tot, 2),
    cost = round(solution_cost, 1),
    utilisation_pct = round(100 * solution_cost / solver_budget, 1)))
  cat(sprintf("  -> CN %.1f / IFM %.1f / MF %.1f | secured %.1f%% | carbon %.2f Mt | util %.1f%%\n",
              100*cn_cells/n_pu_sol, 100*ifm_cells/n_pu_sol, 100*mf_cells/n_pu_sol,
              mean(rep_species$pct_of_target, na.rm = TRUE), c_tot, 100*solution_cost/solver_budget))
}
# ============================= 3. WRITE SWEEP =================================
cat("\n================  BUDGET SWEEP RESULTS (+9% floor, no cap)  ================\n")
print(sweep, row.names = FALSE)
tryCatch(
  writexl::write_xlsx(sweep, file.path(out_dir, "sensitivity1_budget_sweep.xlsx")),
  error = function(e) cat("WARNING: xlsx not written -", conditionMessage(e), "- close it in Excel.\n"))
cat("\nResults table saved to:", file.path(out_dir, "sensitivity1_budget_sweep.xlsx"), "\n")