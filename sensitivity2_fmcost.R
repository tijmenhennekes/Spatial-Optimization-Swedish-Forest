# ==============================================================================
# Sensitivity 2: IFM (integrated management) cost factor.
# IFM opportunity cost as a fraction of CN cost; swept at 0.10, 0.25, 0.50.
# S1 configuration (30% CN floor, no cap); carbon reported, not constrained.
# WTP budget computed once and held fixed. One results table saved as xlsx.
# Author: Tijmen Hennekes
# ==============================================================================
ifm_factors         <- c(0.10, 0.25, 0.50)
baseline_ifm_factor <- 0.25
# ----------------- 1. LOAD DATA ONCE -----------------------------------------
cat("Initializing...\n")
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
cost_cn <- cost_stack[["Close_to_Nature"]]
scherpenhuijzen_cn_mask <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/locked_in_cn.tif")
carbon_rast <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/carbon_inverse_adjusted.tif")
names(carbon_rast) <- zone_names
individual_files_cn  <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_cn",  pattern = "\\.tif$", full.names = TRUE)
individual_files_ifm <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_ifm", pattern = "\\.tif$", full.names = TRUE)
individual_files_mf  <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_mf",  pattern = "\\.tif$", full.names = TRUE)
num_features <- length(individual_files_cn)
species_metadata <- readxl::read_excel("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/IUCN_table_features_filtered.xlsx")
file_codes <- gsub("_mean\\.tif$", "", basename(individual_files_cn))
species_metadata <- species_metadata[match(file_codes, species_metadata$Code), ]
explicit_feature_names <- species_metadata$Code
pz <- function(fl, tm, nv) { s <- terra::mask(terra::resample(terra::rast(fl), tm, "bilinear"), tm); names(s) <- nv; s }
species_cn_stack  <- pz(individual_files_cn,  master_mask, explicit_feature_names)
species_ifm_stack <- pz(individual_files_ifm, master_mask, explicit_feature_names)
species_mf_stack  <- pz(individual_files_mf,  master_mask, explicit_feature_names)
clip_thresh <- 1e-6
species_cn_stack[species_cn_stack < clip_thresh] <- 0
species_ifm_stack[species_ifm_stack < clip_thresh] <- 0
species_mf_stack[species_mf_stack < clip_thresh] <- 0
carbon_cn_stack_raw  <- terra::mask(terra::resample(carbon_rast[["Close_to_Nature"]], master_mask, "bilinear"), master_mask)
carbon_ifm_stack_raw <- terra::mask(terra::resample(carbon_rast[["Integrated"]],      master_mask, "bilinear"), master_mask)
carbon_mf_stack_raw  <- terra::mask(terra::resample(carbon_rast[["Managed_Forest"]],  master_mask, "bilinear"), master_mask)
fm_rast <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/01_raw_data/fm_swe.tif")
fm_aligned <- terra::mask(terra::resample(fm_rast, master_mask, "near"), master_mask)
current_carbon <- terra::ifel(fm_aligned == 1, carbon_cn_stack_raw,
                              terra::ifel(fm_aligned == 2, carbon_ifm_stack_raw, terra::ifel(fm_aligned == 3, carbon_mf_stack_raw, NA)))
carbon_reference_Mt <- (terra::global(current_carbon, "sum", na.rm = TRUE)$sum / 1e6) * (47.321 / 43.366)  # +9%, reported only
cost_ifm_baseline <- baseline_ifm_factor * cost_cn
current_cost_map  <- terra::ifel(fm_aligned == 1, cost_cn,
                                 terra::ifel(fm_aligned == 2, cost_ifm_baseline, terra::ifel(fm_aligned == 3, cost_stack[["Managed_Forest"]], NA)))
solver_budget <- terra::global(current_cost_map, "sum", na.rm = TRUE)$sum   # FIXED across sweep
feature_biodiversity <- prioritizr::zones(species_cn_stack, species_ifm_stack, species_mf_stack,
                                          zone_names = zone_names, feature_names = explicit_feature_names)
targets_matrix <- matrix(0, nrow = num_features, ncol = 3); feature_weights <- rep(1, num_features)
for (i in 1:num_features) { s <- species_metadata$status[i]
if (s %in% c("CR","EN","VU")) { targets_matrix[i,1]<-1.00; targets_matrix[i,2]<-0.50; feature_weights[i]<-10 }
else if (s %in% c("NT","DD")) { targets_matrix[i,1]<-0.90; targets_matrix[i,2]<-0.45; feature_weights[i]<-5 }
else                          { targets_matrix[i,1]<-0.80; targets_matrix[i,2]<-0.40; feature_weights[i]<-1 } }
feature_weights_matrix <- matrix(0, nrow = num_features, ncol = 3)
feature_weights_matrix[,1] <- feature_weights; feature_weights_matrix[,2] <- feature_weights * 0.5
scherpenhuijzen_aligned <- terra::resample(scherpenhuijzen_cn_mask, master_mask, "near")
n_cells <- terra::ncell(master_mask); n_forest <- sum(terra::values(master_mask, mat=FALSE)==1, na.rm=TRUE)
locked_in_matrix <- matrix(FALSE, nrow=n_cells, ncol=3)
locked_in_matrix[,1] <- terra::values(scherpenhuijzen_aligned, mat=FALSE)==1
locked_in_matrix[is.na(locked_in_matrix)] <- FALSE
constraint_matrix_cn <- matrix(0, nrow=n_cells, ncol=3); constraint_matrix_cn[,1] <- 1/n_forest
cat(sprintf("Data loaded. WTP budget (fixed): %.1f | 2030 carbon reference %.2f Mt (reported) | IFM factors: %s\n",
            solver_budget, carbon_reference_Mt, paste(ifm_factors, collapse=", ")))
# ============================= 2. SWEEP LOOP ==================================
sweep <- data.frame()
for (f in ifm_factors) {
  cost_stack[["Integrated"]] <- f * cost_cn
  cat(sprintf("\n===============  IFM cost factor = %.2f  ===============\n", f))
  prob <- prioritizr::problem(cost_stack, feature_biodiversity) %>%
    prioritizr::add_min_shortfall_objective(budget = solver_budget) %>%
    prioritizr::add_relative_targets(targets_matrix) %>%
    prioritizr::add_feature_weights(feature_weights_matrix) %>%
    prioritizr::add_locked_in_constraints(locked_in_matrix) %>%
    prioritizr::add_mandatory_allocation_constraints() %>%
    prioritizr::add_linear_constraints(threshold=0.30, sense=">=", data=constraint_matrix_cn)
  sol <- prob %>% prioritizr::add_gurobi_solver(gap=0.05, time_limit=10800, verbose=TRUE,
                                                numeric_focus=TRUE, threads=4, first_feasible=FALSE) %>% solve(force=TRUE)
  smap <- prioritizr::category_layer(sol)
  cn_c <- sum(terra::values(smap,mat=FALSE)==1,na.rm=TRUE); ifm_c <- sum(terra::values(smap,mat=FALSE)==2,na.rm=TRUE)
  mf_c <- sum(terra::values(smap,mat=FALSE)==3,na.rm=TRUE); npu <- cn_c+ifm_c+mf_c
  tc <- as.data.frame(prioritizr::eval_target_coverage_summary(prob, sol)); tc$zone <- as.character(tc$zone)
  cn  <- tc[tc$zone=="Close_to_Nature", c("feature","absolute_target","absolute_held")]
  ifm <- tc[tc$zone=="Integrated",      c("feature","absolute_target","absolute_held")]
  names(cn)<-c("sp","atc","ahc"); names(ifm)<-c("sp","ati","ahi")
  rs <- merge(cn,ifm,by="sp")
  rs$conservation_held <- rs$ahc+rs$ahi; rs$total_target <- rs$atc+rs$ati
  rs$pct_of_target <- ifelse(rs$total_target>0, 100*rs$conservation_held/rs$total_target, NA)
  cc <- terra::global(terra::mask(carbon_cn_stack_raw, terra::ifel(smap==1,1,NA)),"sum",na.rm=TRUE)$sum
  ci <- terra::global(terra::mask(carbon_ifm_stack_raw,terra::ifel(smap==2,1,NA)),"sum",na.rm=TRUE)$sum
  cm <- terra::global(terra::mask(carbon_mf_stack_raw, terra::ifel(smap==3,1,NA)),"sum",na.rm=TRUE)$sum
  c_tot <- (cc+ci+cm)/1e6
  scost <- terra::global(cost_cn*(smap==1) + cost_stack[["Integrated"]]*(smap==2) +
                           cost_stack[["Managed_Forest"]]*(smap==3), "sum", na.rm=TRUE)$sum
  sweep <- rbind(sweep, data.frame(
    ifm_cost_factor = f, CN_pct=round(100*cn_c/npu,1), IFM_pct=round(100*ifm_c/npu,1), MF_pct=round(100*mf_c/npu,1),
    mean_conservation_held=round(mean(rs$conservation_held),3),
    pct_of_target=round(mean(rs$pct_of_target,na.rm=TRUE),1), carbon_Mt=round(c_tot,2),
    cost=round(scost,1), utilisation_pct=round(100*scost/solver_budget,1)))
  cat(sprintf("  -> CN %.1f / IFM %.1f / MF %.1f | secured %.1f%% | carbon %.2f Mt | util %.1f%%\n",
              100*cn_c/npu,100*ifm_c/npu,100*mf_c/npu, mean(rs$pct_of_target,na.rm=TRUE), c_tot, 100*scost/solver_budget))
}
# ============================= 3. WRITE SWEEP =================================
cat("\n===============  IFM-COST SENSITIVITY RESULTS (30% floor, no cap, carbon reported)  ===============\n")
print(sweep, row.names=FALSE)
tryCatch(writexl::write_xlsx(sweep, file.path(out_dir, "sensitivity2_fmcost_sweep.xlsx")),
         error=function(e) cat("WARNING: xlsx not written -", conditionMessage(e), "\n"))
cat("\nResults table saved to:", file.path(out_dir, "sensitivity2_fmcost_sweep.xlsx"), "\n")