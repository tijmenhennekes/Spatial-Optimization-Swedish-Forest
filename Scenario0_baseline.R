# ==============================================================================
# Scenario 0: Current Forest Baseline (Business-as-Usual reference)
# Author: Tijmen Hennekes
# Evaluates the CURRENT allocation (fm_swe.tif) - NOT an optimisation.
# Secured framing: targets are attainable only in CN + IFM (MF = production).
# ==============================================================================
cat("Initializing baseline evaluation environment...\n")
library(terra); library(sf); library(prioritizr); library(tibble); library(readxl); library(writexl)
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
# EULOCC1K cost is EUR/ha/yr; x cell area -> EUR/cell in EUR million. Uniform rescale.
cost_stack <- cost_stack * terra::cellSize(cost_stack[[1]], unit = "ha") / 1e6
ifm_cost_factor <- 0.25
cost_stack[["Integrated"]] <- ifm_cost_factor * cost_stack[["Close_to_Nature"]]
carbon_rast <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/carbon_inverse_adjusted.tif")
names(carbon_rast) <- zone_names
fm_rast    <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/01_raw_data/fm_swe.tif")
fm_aligned <- terra::mask(terra::resample(fm_rast, master_mask, method = "near"), master_mask)
cat("Processing species layers...\n")
individual_files_cn  <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_cn",  pattern = "\\.tif$", full.names = TRUE)
individual_files_ifm <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_ifm", pattern = "\\.tif$", full.names = TRUE)
individual_files_mf  <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_mf",  pattern = "\\.tif$", full.names = TRUE)
num_features <- length(individual_files_cn)
species_metadata       <- readxl::read_excel("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/IUCN_table_features_filtered.xlsx")
file_codes             <- gsub("_mean\\.tif$", "", basename(individual_files_cn))
species_metadata       <- species_metadata[match(file_codes, species_metadata$Code), ]
explicit_feature_names <- species_metadata$Code
process_zone_stack <- function(file_list, template_mask, names_vector) {
  s <- terra::rast(file_list)
  s <- terra::resample(s, template_mask, method = "bilinear")
  s <- terra::mask(s, template_mask)
  names(s) <- names_vector
  s
}
species_cn_stack  <- process_zone_stack(individual_files_cn,  master_mask, explicit_feature_names)
species_ifm_stack <- process_zone_stack(individual_files_ifm, master_mask, explicit_feature_names)
species_mf_stack  <- process_zone_stack(individual_files_mf,  master_mask, explicit_feature_names)
clip_thresh <- 1e-6
species_cn_stack[species_cn_stack   < clip_thresh] <- 0
species_ifm_stack[species_ifm_stack < clip_thresh] <- 0
species_mf_stack[species_mf_stack   < clip_thresh] <- 0
cat("Species stacks ready:", num_features, "species x 3 zones\n")
carbon_cn_stack_raw  <- terra::mask(terra::resample(carbon_rast[["Close_to_Nature"]], master_mask, method = "bilinear"), master_mask)
carbon_ifm_stack_raw <- terra::mask(terra::resample(carbon_rast[["Integrated"]],      master_mask, method = "bilinear"), master_mask)
carbon_mf_stack_raw  <- terra::mask(terra::resample(carbon_rast[["Managed_Forest"]],  master_mask, method = "bilinear"), master_mask)
feature_biodiversity <- prioritizr::zones(
  species_cn_stack, species_ifm_stack, species_mf_stack,
  zone_names = zone_names, feature_names = explicit_feature_names)
targets_matrix <- matrix(0, nrow = num_features, ncol = 3)
for (i in 1:num_features) {
  status <- species_metadata$status[i]
  if (status %in% c("CR", "EN", "VU")) { targets_matrix[i, 1] <- 1.00; targets_matrix[i, 2] <- 0.50 }
  else if (status %in% c("NT", "DD"))  { targets_matrix[i, 1] <- 0.90; targets_matrix[i, 2] <- 0.45 }
  else                                  { targets_matrix[i, 1] <- 0.80; targets_matrix[i, 2] <- 0.40 }
}
feature_weights <- rep(1, num_features)
for (i in 1:num_features) {
  status <- species_metadata$status[i]
  if      (status %in% c("CR", "EN", "VU")) feature_weights[i] <- 10
  else if (status %in% c("NT", "DD"))       feature_weights[i] <- 5
}
feature_weights_matrix      <- matrix(0, nrow = num_features, ncol = 3)
feature_weights_matrix[, 1] <- feature_weights
feature_weights_matrix[, 2] <- feature_weights * 0.5
total_cn_cost  <- terra::global(cost_stack[["Close_to_Nature"]], "sum", na.rm = TRUE)$sum
total_ifm_cost <- terra::global(cost_stack[["Integrated"]],      "sum", na.rm = TRUE)$sum
cat("Constructing current-management allocation...\n")
fm0    <- terra::ifel(is.na(fm_aligned), 0, fm_aligned)
pu_cn  <- !is.na(cost_stack[["Close_to_Nature"]])
pu_ifm <- !is.na(cost_stack[["Integrated"]])
pu_mf  <- !is.na(cost_stack[["Managed_Forest"]])
sol_cn  <- terra::ifel(pu_cn,  terra::ifel(fm0 == 1, 1, 0), NA)
sol_ifm <- terra::ifel(pu_ifm, terra::ifel(fm0 == 2, 1, 0), NA)
sol_mf  <- terra::ifel(pu_mf,  terra::ifel(fm0 == 3, 1, 0), NA)
current_solution <- c(sol_cn, sol_ifm, sol_mf)
names(current_solution) <- zone_names
scoring_problem <- prioritizr::problem(cost_stack, feature_biodiversity) %>%
  prioritizr::add_min_shortfall_objective(budget = total_cn_cost + total_ifm_cost) %>%
  prioritizr::add_relative_targets(targets_matrix) %>%
  prioritizr::add_feature_weights(feature_weights_matrix)

# --- MANAGEMENT REGIMES ---
cat("\n--- CURRENT MANAGEMENT REGIMES ---\n")
cn_n  <- sum(terra::values(fm_aligned, mat = FALSE) == 1, na.rm = TRUE)
ifm_n <- sum(terra::values(fm_aligned, mat = FALSE) == 2, na.rm = TRUE)
mf_n  <- sum(terra::values(fm_aligned, mat = FALSE) == 3, na.rm = TRUE)
tot_n <- cn_n + ifm_n + mf_n
cat(sprintf("CN:  %d (%.1f%%)\n", cn_n,  100 * cn_n  / tot_n))
cat(sprintf("IFM: %d (%.1f%%)\n", ifm_n, 100 * ifm_n / tot_n))
cat(sprintf("MF:  %d (%.1f%%)\n", mf_n,  100 * mf_n  / tot_n))

# --- BIODIVERSITY: TARGET COVERAGE (SECURED = CN + IFM) ---
cat("\n--- BIODIVERSITY: TARGET COVERAGE ---\n")
target_cov <- as.data.frame(prioritizr::eval_target_coverage_summary(scoring_problem, current_solution))
target_cov$zone <- as.character(target_cov$zone)
sec_cov <- target_cov[target_cov$zone %in% c("Close_to_Nature", "Integrated"), ]
cat(sprintf("Secured zone-targets met (CN+IFM): %d / %d (%.1f%%)\n",
            sum(sec_cov$met), nrow(sec_cov), 100 * sum(sec_cov$met) / nrow(sec_cov)))
cat(sprintf("Mean relative shortfall (CN+IFM):  %.1f%%\n", 100 * mean(sec_cov$relative_shortfall, na.rm = TRUE)))

cn_rep  <- target_cov[target_cov$zone == "Close_to_Nature", c("feature","relative_target","absolute_target","absolute_held","met")]
ifm_rep <- target_cov[target_cov$zone == "Integrated",      c("feature","relative_target","absolute_target","absolute_held","met")]
mf_rep  <- target_cov[target_cov$zone == "Managed_Forest",  c("feature","absolute_held")]
names(cn_rep)  <- c("species","target_cn","abs_target_cn","abs_held_cn","met_cn")
names(ifm_rep) <- c("species","target_ifm","abs_target_ifm","abs_held_ifm","met_ifm")
names(mf_rep)  <- c("species","abs_held_mf")
rep_species <- merge(merge(cn_rep, ifm_rep, by = "species"), mf_rep, by = "species")
grp_map <- c(AMP = "Amphibians", AVE = "Birds", MAM = "Mammals", REP = "Reptiles")
rep_species$group  <- grp_map[substr(rep_species$species, 1, 3)]
rep_species$group[is.na(rep_species$group)] <- substr(rep_species$species[is.na(rep_species$group)], 1, 3)
rep_species$status <- species_metadata$status[match(rep_species$species, species_metadata$Code)]
rep_species$total_held        <- rep_species$abs_held_cn + rep_species$abs_held_ifm + rep_species$abs_held_mf
rep_species$conservation_held <- rep_species$abs_held_cn + rep_species$abs_held_ifm
rep_species$total_target      <- rep_species$abs_target_cn + rep_species$abs_target_ifm
rep_species$pct_of_target <- ifelse(rep_species$total_target > 0,
                                    100 * rep_species$conservation_held / rep_species$total_target, NA)
rep_species$represented   <- rep_species$conservation_held >= rep_species$total_target

cat("\n--- BIODIVERSITY REPRESENTATION (SECURED = CN + IFM; 264 species) ---\n")
cat(sprintf("Mean secured representation vs target:      %.1f%% of target\n",
            mean(rep_species$pct_of_target, na.rm = TRUE)))
cat(sprintf("Species meeting target (secured >= target): %d / %d\n",
            sum(rep_species$represented), nrow(rep_species)))
cat(sprintf("Species with habitat present (any zone):    %d / %d\n",
            sum(rep_species$total_held > 0), nrow(rep_species)))
cat(sprintf("Species with secured habitat (CN+IFM > 0):  %d / %d\n",
            sum(rep_species$conservation_held > 0), nrow(rep_species)))
gap <- rep_species[order(rep_species$conservation_held), ]
cat(sprintf("Gap species with < 1%% secured habitat: %d\n", sum(rep_species$conservation_held < 0.01)))
if (sum(rep_species$conservation_held < 0.01) > 0) print(table(gap$status[gap$conservation_held < 0.01]))

# --- CARBON (relative, uncalibrated) ---
cat("\n--- CARBON SEQUESTRATION (relative, uncalibrated) ---\n")
c_cn  <- terra::global(terra::mask(carbon_cn_stack_raw,  terra::ifel(fm_aligned == 1, 1, NA)), "sum", na.rm = TRUE)$sum
c_ifm <- terra::global(terra::mask(carbon_ifm_stack_raw, terra::ifel(fm_aligned == 2, 1, NA)), "sum", na.rm = TRUE)$sum
c_mf  <- terra::global(terra::mask(carbon_mf_stack_raw,  terra::ifel(fm_aligned == 3, 1, NA)), "sum", na.rm = TRUE)$sum
c_tot <- (c_cn + c_ifm + c_mf) / 1e6
carbon_ref_2030_Mt <- c_tot * (47.3 / 46.2)
cat(sprintf("CN: %.2f  IFM: %.2f  MF: %.2f  Total: %.2f Mt\n", c_cn/1e6, c_ifm/1e6, c_mf/1e6, c_tot))
cat(sprintf("Relative 2030 reference: %.2f Mt\n", carbon_ref_2030_Mt))
cat(sprintf("S0 vs 2030 reference: %.1f%%\n", 100 * c_tot / carbon_ref_2030_Mt))

# --- COST ---
cat("\n--- MANAGEMENT COST (current allocation) ---\n")
current_cost <- terra::global(
  cost_stack[["Close_to_Nature"]] * (fm_aligned == 1) +
    cost_stack[["Integrated"]]      * (fm_aligned == 2) +
    cost_stack[["Managed_Forest"]]  * (fm_aligned == 3),
  "sum", na.rm = TRUE)$sum
cat(sprintf("Current management cost: %.1f  (%.1f%% of full-CN ceiling)\n",
            current_cost, 100 * current_cost / total_cn_cost))

# --- EXPORT ---
writexl::write_xlsx(rep_species, file.path(out_dir, "scenario0_representation_summary.xlsx"))
writexl::write_xlsx(target_cov,  file.path(out_dir, "scenario0_target_coverage_summary.xlsx"))
terra::writeRaster(fm_aligned, file.path(out_dir, "scenario0_current_zones.tif"), overwrite = TRUE)
cat("\nBaseline results saved to:", out_dir, "\n")

rep_species$tier <- ifelse(rep_species$status %in% c("CR","EN","VU"), "1 CR/EN/VU",
                           ifelse(rep_species$status %in% c("NT","DD"),      "2 NT/DD", "3 LC/other"))
rep_species$tier <- factor(rep_species$tier, levels = c("1 CR/EN/VU","2 NT/DD","3 LC/other"))

tier_tbl <- do.call(rbind, lapply(levels(rep_species$tier), function(t) {
  s <- rep_species[rep_species$tier == t, ]
  data.frame(tier = t, species = nrow(s),
             avg_pct_secured    = round(mean(s$pct_of_target,   na.rm = TRUE), 1),
             median_pct_secured = round(median(s$pct_of_target, na.rm = TRUE), 1))
}))
tier_tbl <- rbind(tier_tbl, data.frame(tier = "All", species = nrow(rep_species),
                                       avg_pct_secured    = round(mean(rep_species$pct_of_target,   na.rm = TRUE), 1),
                                       median_pct_secured = round(median(rep_species$pct_of_target, na.rm = TRUE), 1)))
print(tier_tbl, row.names = FALSE)