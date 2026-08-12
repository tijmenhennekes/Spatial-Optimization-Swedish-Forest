# ==============================================================================
# Scenario 1b: Minimum-feasible policy solution
#   = S1, but at 1.3x the WTP budget, CN limited to 30% (floor only),
#     and the CORRECTED 2030 LULUCF carbon floor.
# Purpose: the smallest budget at which BOTH the 30% protection target and the
#          corrected 2030 carbon target can be met simultaneously. Supported by
#          the budget sweep (target reached ~1.25-1.3x current spend).
# Author: Tijmen Hennekes
# ==============================================================================
cat("Initializing scenario 1b (minimum-feasible) environment...\n")
library(terra); library(sf); library(gurobi); library(prioritizr)
library(tibble); library(readxl); library(writexl)
set.seed(24)
out_dir <- "C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/04-outputs"
master_mask <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/forest_binary.tif")
zone_names <- c("Close_to_Nature", "Integrated", "Managed_Forest")
cost_stack <- c(
  terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/cost_cn.tif"),
  terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/cost_ifm.tif"),
  terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/cost_mf.tif")
)
names(cost_stack) <- zone_names
cost_stack <- cost_stack * terra::cellSize(cost_stack[[1]], unit = "ha") / 1e6
ifm_cost_factor <- 0.25
cost_stack[["Integrated"]] <- ifm_cost_factor * cost_stack[["Close_to_Nature"]]
scherpenhuijzen_cn_mask <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/locked_in_cn.tif")
carbon_rast <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/carbon_inverse_adjusted.tif")
names(carbon_rast) <- zone_names
# --------------------------- FEATURE INPUTS -----------------------------------
cat("Processing ecological distribution layers...\n")
individual_files_cn  <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_cn",  pattern = "\\.tif$", full.names = TRUE)
individual_files_ifm <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_ifm", pattern = "\\.tif$", full.names = TRUE)
individual_files_mf  <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_mf",  pattern = "\\.tif$", full.names = TRUE)
num_features <- length(individual_files_cn)
species_metadata       <- readxl::read_excel("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/IUCN_table_features_filtered.xlsx")
file_codes             <- gsub("_mean\\.tif$", "", basename(individual_files_cn))
species_metadata       <- species_metadata[match(file_codes, species_metadata$Code), ]
explicit_feature_names <- species_metadata$Code
process_zone_stack <- function(file_list, template_mask, names_vector) {
  s <- terra::mask(terra::resample(terra::rast(file_list), template_mask, "bilinear"), template_mask)
  names(s) <- names_vector; s
}
species_cn_stack  <- process_zone_stack(individual_files_cn,  master_mask, explicit_feature_names)
species_ifm_stack <- process_zone_stack(individual_files_ifm, master_mask, explicit_feature_names)
species_mf_stack  <- process_zone_stack(individual_files_mf,  master_mask, explicit_feature_names)
clip_thresh <- 1e-6
species_cn_stack[species_cn_stack   < clip_thresh] <- 0
species_ifm_stack[species_ifm_stack < clip_thresh] <- 0
species_mf_stack[species_mf_stack   < clip_thresh] <- 0
carbon_cn_stack_raw  <- terra::mask(terra::resample(carbon_rast[["Close_to_Nature"]], master_mask, "bilinear"), master_mask)
carbon_ifm_stack_raw <- terra::mask(terra::resample(carbon_rast[["Integrated"]],      master_mask, "bilinear"), master_mask)
carbon_mf_stack_raw  <- terra::mask(terra::resample(carbon_rast[["Managed_Forest"]],  master_mask, "bilinear"), master_mask)
names(carbon_cn_stack_raw) <- names(carbon_ifm_stack_raw) <- names(carbon_mf_stack_raw) <- "Carbon_Sequestration"
carbon_constraint_stack <- c(carbon_cn_stack_raw, carbon_ifm_stack_raw, carbon_mf_stack_raw) / 1e6
names(carbon_constraint_stack) <- zone_names
# ---- Current management overlay: carbon reference and WTP baseline ----
fm_rast    <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/01_raw_data/fm_swe.tif")
fm_aligned <- terra::mask(terra::resample(fm_rast, master_mask, method = "near"), master_mask)
current_carbon <- terra::ifel(fm_aligned == 1, carbon_cn_stack_raw,
  terra::ifel(fm_aligned == 2, carbon_ifm_stack_raw,
              terra::ifel(fm_aligned == 3, carbon_mf_stack_raw, NA)))
current_mgmt_carbon_Mt <- terra::global(current_carbon, "sum", na.rm = TRUE)$sum / 1e6
# CORRECTED 2030 LULUCF ratio (Annex IIa, internally consistent pair). +3.955 Mt
# increment over the 2016-2018 baseline => ~+9%. MUST match scenario1_policy.R.
carbon_target_Mt <- current_mgmt_carbon_Mt * (47.321 / 43.366)
cat(sprintf("Current mgmt carbon: %.2f Mt  ->  corrected carbon target: %.2f Mt\n",
            current_mgmt_carbon_Mt, carbon_target_Mt))
current_cost_map <- terra::ifel(fm_aligned == 1, cost_stack[["Close_to_Nature"]],
  terra::ifel(fm_aligned == 2, cost_stack[["Integrated"]],
              terra::ifel(fm_aligned == 3, cost_stack[["Managed_Forest"]], NA)))
current_mgmt_cost <- terra::global(current_cost_map, "sum", na.rm = TRUE)$sum
# ------------------------ ZONES / TARGETS / WEIGHTS ---------------------------
feature_biodiversity <- prioritizr::zones(species_cn_stack, species_ifm_stack, species_mf_stack,
  zone_names = zone_names, feature_names = explicit_feature_names)
targets_matrix <- matrix(0, nrow = num_features, ncol = 3)
for (i in 1:num_features) {
  status <- species_metadata$status[i]
  if (status %in% c("CR","EN","VU")) { targets_matrix[i,1] <- 1.00; targets_matrix[i,2] <- 0.50 }
  else if (status %in% c("NT","DD")) { targets_matrix[i,1] <- 0.90; targets_matrix[i,2] <- 0.45 }
  else                                { targets_matrix[i,1] <- 0.80; targets_matrix[i,2] <- 0.40 }
}
feature_weights <- rep(1, num_features)
for (i in 1:num_features) {
  status <- species_metadata$status[i]
  if      (status %in% c("CR","EN","VU")) feature_weights[i] <- 10
  else if (status %in% c("NT","DD"))      feature_weights[i] <- 5
}
feature_weights_matrix      <- matrix(0, nrow = num_features, ncol = 3)
feature_weights_matrix[, 1] <- feature_weights
feature_weights_matrix[, 2] <- feature_weights * 0.5
# ------------------------------ CONSTRAINTS -----------------------------------
scherpenhuijzen_aligned <- terra::resample(scherpenhuijzen_cn_mask, master_mask, method = "near")
n_cells  <- terra::ncell(master_mask)
n_forest <- sum(terra::values(master_mask, mat = FALSE) == 1, na.rm = TRUE)
locked_in_matrix      <- matrix(FALSE, nrow = n_cells, ncol = 3)
locked_in_matrix[, 1] <- terra::values(scherpenhuijzen_aligned, mat = FALSE) == 1
locked_in_matrix[is.na(locked_in_matrix)] <- FALSE
constraint_matrix_cn      <- matrix(0, nrow = n_cells, ncol = 3)
constraint_matrix_cn[, 1] <- 1 / n_forest
total_cn_cost  <- terra::global(cost_stack[["Close_to_Nature"]], "sum", na.rm = TRUE)$sum
# ---- CHANGE 1: budget = 1.3x current management spend (minimum feasible) ----
solver_budget   <- 1.3 * current_mgmt_cost
budget_fraction <- solver_budget / total_cn_cost
cat(sprintf("Budget = 1.3x WTP: %.1f (%.1f%% of full-CN ceiling)\n", solver_budget, 100*budget_fraction))
cat(sprintf("Corrected carbon target (Mt): %.2f\n", carbon_target_Mt))
# ------------------------- PROBLEM (MIN-SHORTFALL) ----------------------------
# CHANGE 2: CN cap removed. 30% is a FLOOR only (EU minimum); CN may exceed it.
policy_problem <- prioritizr::problem(cost_stack, feature_biodiversity) %>%
  prioritizr::add_min_shortfall_objective(budget = solver_budget) %>%
  prioritizr::add_relative_targets(targets_matrix) %>%
  prioritizr::add_feature_weights(feature_weights_matrix) %>%
  prioritizr::add_locked_in_constraints(locked_in_matrix) %>%
  prioritizr::add_mandatory_allocation_constraints() %>%
  prioritizr::add_linear_constraints(threshold = 0.30, sense = ">=", data = constraint_matrix_cn) %>%
  prioritizr::add_linear_constraints(threshold = carbon_target_Mt, sense = ">=", data = carbon_constraint_stack)
cat("Running pre-solve feasibility check...\n")
prioritizr::presolve_check(policy_problem)
cat("\n========  LAUNCHING GUROBI - SCENARIO 1b (minimum-feasible)  ========\n")
policy_solution <- policy_problem %>%
  prioritizr::add_gurobi_solver(gap = 0.05, time_limit = 10800, verbose = TRUE,
                                numeric_focus = TRUE, threads = 4, first_feasible = FALSE) %>%
  solve(force = TRUE)
# ==============================================================================
# RESULTS
# ==============================================================================
solution_map <- prioritizr::category_layer(policy_solution)
cn_cells  <- sum(terra::values(solution_map, mat = FALSE) == 1, na.rm = TRUE)
ifm_cells <- sum(terra::values(solution_map, mat = FALSE) == 2, na.rm = TRUE)
mf_cells  <- sum(terra::values(solution_map, mat = FALSE) == 3, na.rm = TRUE)
n_pu_sol  <- cn_cells + ifm_cells + mf_cells
cat("\n--- ZONE ALLOCATION ---\n")
cat(sprintf("CN:  %d (%.1f%%)\nIFM: %d (%.1f%%)\nMF:  %d (%.1f%%)\n",
            cn_cells, 100*cn_cells/n_pu_sol, ifm_cells, 100*ifm_cells/n_pu_sol, mf_cells, 100*mf_cells/n_pu_sol))

# --- BIODIVERSITY: SECURED (CN + IFM) ---
target_cov <- as.data.frame(prioritizr::eval_target_coverage_summary(policy_problem, policy_solution))
target_cov$zone <- as.character(target_cov$zone)
sec_cov <- target_cov[target_cov$zone %in% c("Close_to_Nature", "Integrated"), ]
cat("\n--- BIODIVERSITY: TARGET COVERAGE (SECURED = CN + IFM) ---\n")
cat(sprintf("Secured zone-targets met (CN+IFM): %d / %d (%.1f%%)\n",
            sum(sec_cov$met), nrow(sec_cov), 100 * sum(sec_cov$met) / nrow(sec_cov)))
cat(sprintf("Mean relative shortfall (CN+IFM):  %.1f%%\n", 100 * mean(sec_cov$relative_shortfall, na.rm = TRUE)))
cn_rep  <- target_cov[target_cov$zone == "Close_to_Nature", c("feature","absolute_target","absolute_held")]
ifm_rep <- target_cov[target_cov$zone == "Integrated",      c("feature","absolute_target","absolute_held")]
mf_rep  <- target_cov[target_cov$zone == "Managed_Forest",  c("feature","absolute_held")]
names(cn_rep)  <- c("species","abs_target_cn","abs_held_cn")
names(ifm_rep) <- c("species","abs_target_ifm","abs_held_ifm")
names(mf_rep)  <- c("species","abs_held_mf")
rep_species <- merge(merge(cn_rep, ifm_rep, by = "species"), mf_rep, by = "species")
rep_species$status <- species_metadata$status[match(rep_species$species, species_metadata$Code)]
rep_species$total_held        <- rep_species$abs_held_cn + rep_species$abs_held_ifm + rep_species$abs_held_mf
rep_species$conservation_held <- rep_species$abs_held_cn + rep_species$abs_held_ifm
rep_species$total_target      <- rep_species$abs_target_cn + rep_species$abs_target_ifm
rep_species$pct_of_target <- ifelse(rep_species$total_target > 0,
                                    100 * rep_species$conservation_held / rep_species$total_target, NA)
rep_species$represented   <- rep_species$conservation_held >= rep_species$total_target
cat("\n--- BIODIVERSITY REPRESENTATION (SECURED; 264 species) ---\n")
cat(sprintf("Mean secured coverage (%% of target):        %.1f%%\n", mean(rep_species$pct_of_target, na.rm = TRUE)))
cat(sprintf("Species meeting target (secured >= target): %d / %d\n", sum(rep_species$represented), nrow(rep_species)))
cat(sprintf("Species with habitat present (any zone):    %d / %d\n", sum(rep_species$total_held > 0), nrow(rep_species)))
cat(sprintf("Species with secured habitat (CN+IFM > 0):  %d / %d\n", sum(rep_species$conservation_held > 0), nrow(rep_species)))
cat(sprintf("Effectively unprotected (<1%% secured):      %d\n", sum(rep_species$conservation_held < 0.01)))
if (sum(rep_species$conservation_held < 0.01) > 0) print(table(rep_species$status[rep_species$conservation_held < 0.01]))
rep_species$tier <- ifelse(rep_species$status %in% c("CR","EN","VU"), "1 CR/EN/VU",
                    ifelse(rep_species$status %in% c("NT","DD"),      "2 NT/DD", "3 LC/other"))
rep_species$tier <- factor(rep_species$tier, levels = c("1 CR/EN/VU","2 NT/DD","3 LC/other"))
tier_tbl <- do.call(rbind, lapply(levels(rep_species$tier), function(t) {
  s <- rep_species[rep_species$tier == t, ]
  data.frame(tier = t, species = nrow(s), avg_pct_secured = round(mean(s$pct_of_target, na.rm = TRUE), 1))
}))
cat("\n--- SECURED COVERAGE BY IUCN TIER (S1b) ---\n"); print(tier_tbl, row.names = FALSE)

# --- CARBON ---
cat("\n--- CARBON (relative) ---\n")
c_cn  <- terra::global(terra::mask(carbon_cn_stack_raw,  terra::ifel(solution_map == 1, 1, NA)), "sum", na.rm = TRUE)$sum
c_ifm <- terra::global(terra::mask(carbon_ifm_stack_raw, terra::ifel(solution_map == 2, 1, NA)), "sum", na.rm = TRUE)$sum
c_mf  <- terra::global(terra::mask(carbon_mf_stack_raw,  terra::ifel(solution_map == 3, 1, NA)), "sum", na.rm = TRUE)$sum
c_tot <- (c_cn + c_ifm + c_mf) / 1e6
cat(sprintf("S1b carbon: %.2f Mt  |  target %.2f Mt  |  met: %s  |  overshoot: %+.1f%%\n",
            c_tot, carbon_target_Mt, c_tot >= carbon_target_Mt - 1e-6, 100*(c_tot/carbon_target_Mt - 1)))

# --- COST ---
solution_cost <- terra::global(
  cost_stack[["Close_to_Nature"]] * (solution_map == 1) +
    cost_stack[["Integrated"]]      * (solution_map == 2) +
    cost_stack[["Managed_Forest"]]  * (solution_map == 3), "sum", na.rm = TRUE)$sum
cat(sprintf("\nS1b cost: %.1f  (%.1f%% of full-CN ceiling)  |  utilisation: %.1f%%\n",
            solution_cost, 100 * solution_cost / total_cn_cost, 100 * solution_cost / solver_budget))

# --- EXPORT ---
terra::writeRaster(solution_map, file.path(out_dir, "scenario1b_feasible_solution_zones.tif"), overwrite = TRUE)
tryCatch(writexl::write_xlsx(rep_species, file.path(out_dir, "scenario1b_feasible_representation.xlsx")),
         error = function(e) cat("WARNING: xlsx not written -", conditionMessage(e), "\n"))
cat("\nScenario 1b results saved to:", out_dir, "\n")
