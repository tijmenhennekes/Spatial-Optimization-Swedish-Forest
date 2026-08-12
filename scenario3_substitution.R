# ==============================================================================
# Scenario 3: Substitution and the displacement factor.
# Carbon redefined to include HWP substitution: seq = sink + harvest x DF.
# Author: Tijmen Hennekes
# ==============================================================================
cat("Initializing S3 (substitution DF test)...\n")
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
cost_stack[["Integrated"]] <- 0.25 * cost_stack[["Close_to_Nature"]]
scherpenhuijzen_cn_mask <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/locked_in_cn.tif")

# ---- Sink and harvest layers (both 3-band CN/IFM/MF) ----
align <- function(r) terra::mask(terra::resample(r, master_mask, "bilinear"), master_mask)
sink_rast <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/carbon_inverse_adjusted.tif")
names(sink_rast) <- zone_names
sink_cn <- align(sink_rast[["Close_to_Nature"]]); sink_ifm <- align(sink_rast[["Integrated"]]); sink_mf <- align(sink_rast[["Managed_Forest"]])
harv_rast <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/01_raw_data/harv_swe.tif")
names(harv_rast) <- zone_names
harv_rast <- terra::mask(terra::project(harv_rast, master_mask, method = "bilinear"), master_mask)
harv_cn <- harv_rast[["Close_to_Nature"]]; harv_ifm <- harv_rast[["Integrated"]]; harv_mf <- harv_rast[["Managed_Forest"]]
cat(sprintf("Zone mean harvest -> CN %.3g  IFM %.3g  MF %.3g  (CN should be ~0)\n",
            terra::global(harv_cn,"mean",na.rm=TRUE)[[1]], terra::global(harv_ifm,"mean",na.rm=TRUE)[[1]], terra::global(harv_mf,"mean",na.rm=TRUE)[[1]]))

# ---- Species + targets/weights ----
f_cn <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_cn",  pattern="\\.tif$", full.names=TRUE)
f_ifm<- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_ifm", pattern="\\.tif$", full.names=TRUE)
f_mf <- list.files("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/biodiversity_norm/species_norm_mf",  pattern="\\.tif$", full.names=TRUE)
nf <- length(f_cn)
meta <- readxl::read_excel("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/02_processed_data/IUCN_table_features_filtered.xlsx")
meta <- meta[match(gsub("_mean\\.tif$","",basename(f_cn)), meta$Code), ]; nm <- meta$Code
pzs <- function(fl){ s <- terra::mask(terra::resample(terra::rast(fl), master_mask, "bilinear"), master_mask); names(s)<-nm; s }
s_cn<-pzs(f_cn); s_ifm<-pzs(f_ifm); s_mf<-pzs(f_mf)
cl<-1e-6; s_cn[s_cn<cl]<-0; s_ifm[s_ifm<cl]<-0; s_mf[s_mf<cl]<-0
feat <- prioritizr::zones(s_cn, s_ifm, s_mf, zone_names=zone_names, feature_names=nm)
tm <- matrix(0,nf,3)
for(i in 1:nf){ st<-meta$status[i]
if(st%in%c("CR","EN","VU")){tm[i,1]<-1.00;tm[i,2]<-0.50} else if(st%in%c("NT","DD")){tm[i,1]<-0.90;tm[i,2]<-0.45} else {tm[i,1]<-0.80;tm[i,2]<-0.40} }
fw<-rep(1,nf); for(i in 1:nf){st<-meta$status[i]; if(st%in%c("CR","EN","VU"))fw[i]<-10 else if(st%in%c("NT","DD"))fw[i]<-5}
fwm<-matrix(0,nf,3); fwm[,1]<-fw; fwm[,2]<-fw*0.5

# ---- Budget + constraints ----
fm_rast <- terra::rast("C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/01_raw_data/fm_swe.tif")
fm_aligned <- terra::mask(terra::resample(fm_rast, master_mask, method="near"), master_mask)
current_cost_map <- terra::ifel(fm_aligned==1, cost_stack[["Close_to_Nature"]], terra::ifel(fm_aligned==2, cost_stack[["Integrated"]], terra::ifel(fm_aligned==3, cost_stack[["Managed_Forest"]], NA)))
solver_budget <- terra::global(current_cost_map, "sum", na.rm=TRUE)$sum   # WTP
carbon_target_Mt <- 91.03
scherpenhuijzen_aligned <- terra::resample(scherpenhuijzen_cn_mask, master_mask, method="near")
n_cells <- terra::ncell(master_mask); n_forest <- sum(terra::values(master_mask, mat=FALSE)==1, na.rm=TRUE)
locked_in_matrix <- matrix(FALSE, n_cells, 3); locked_in_matrix[,1] <- terra::values(scherpenhuijzen_aligned, mat=FALSE)==1
locked_in_matrix[is.na(locked_in_matrix)] <- FALSE
cn_cons <- matrix(0, n_cells, 3); cn_cons[,1] <- 1/n_forest

score <- function(sol_cat){  # secured-biodiversity summary for a category solution
  sol <- c(terra::ifel(!is.na(cost_stack[[1]]), terra::ifel(sol_cat==1,1,0), NA),
           terra::ifel(!is.na(cost_stack[[2]]), terra::ifel(sol_cat==2,1,0), NA),
           terra::ifel(!is.na(cost_stack[[3]]), terra::ifel(sol_cat==3,1,0), NA)); names(sol)<-zone_names
           tc <- as.data.frame(prioritizr::eval_target_coverage_summary(scoring, sol)); tc$zone<-as.character(tc$zone)
           cn<-tc[tc$zone=="Close_to_Nature",c("feature","absolute_target","absolute_held")]; ifm<-tc[tc$zone=="Integrated",c("feature","absolute_target","absolute_held")]
           names(cn)<-c("sp","atc","ahc"); names(ifm)<-c("sp","ati","ahi"); rs<-merge(cn,ifm,by="sp")
           rs$cons<-rs$ahc+rs$ahi; rs$tt<-rs$atc+rs$ati; rs$pct<-ifelse(rs$tt>0,100*rs$cons/rs$tt,NA)
           list(cov=mean(rs$pct,na.rm=TRUE), unprot=sum(rs$cons<0.01))
}
# scoring problem (biodiversity only) for eval_target_coverage_summary
scoring <- prioritizr::problem(cost_stack, feat) %>% prioritizr::add_min_shortfall_objective(budget=solver_budget) %>%
  prioritizr::add_relative_targets(tm) %>% prioritizr::add_feature_weights(fwm)

# ============================= DF SWEEP =======================================
DFs <- c(0.47, 0.72)   # DF = 0 is the sink-only reference (= S1, infeasible at 91.0)
sweep <- data.frame()
for (DF in DFs) {
  cat(sprintf("\n=================  DF = %.2f  =================\n", DF))
  seq_cn <- sink_cn + harv_cn*DF; seq_ifm <- sink_ifm + harv_ifm*DF; seq_mf <- sink_mf + harv_mf*DF
  carbon_stack <- c(seq_cn, seq_ifm, seq_mf) / 1e6; names(carbon_stack) <- zone_names
  prob <- prioritizr::problem(cost_stack, feat) %>%
    prioritizr::add_min_shortfall_objective(budget = solver_budget) %>%
    prioritizr::add_relative_targets(tm) %>% prioritizr::add_feature_weights(fwm) %>%
    prioritizr::add_locked_in_constraints(locked_in_matrix) %>%
    prioritizr::add_mandatory_allocation_constraints() %>%
    prioritizr::add_linear_constraints(threshold = 0.30, sense = ">=", data = cn_cons) %>%
    prioritizr::add_linear_constraints(threshold = carbon_target_Mt, sense = ">=", data = carbon_stack)
  sol <- tryCatch(
    prob %>% prioritizr::add_gurobi_solver(gap=0.05, time_limit=10800, verbose=TRUE, numeric_focus=TRUE, threads=4) %>% solve(force=TRUE),
    error = function(e) { cat("  INFEASIBLE at DF =", DF, "\n"); NULL })
  if (is.null(sol)) { sweep <- rbind(sweep, data.frame(DF=DF, feasible=FALSE, CN=NA,IFM=NA,MF=NA,coverage=NA,unprotected=NA,sink_Mt=NA,total_Mt=NA)); next }
  smap <- prioritizr::category_layer(sol)
  a<-sum(terra::values(smap)==1,na.rm=TRUE); b<-sum(terra::values(smap)==2,na.rm=TRUE); d<-sum(terra::values(smap)==3,na.rm=TRUE); n<-a+b+d
  bio <- score(smap)
  sink_tot <- (terra::global(terra::mask(sink_cn,terra::ifel(smap==1,1,NA)),"sum",na.rm=TRUE)$sum +
                 terra::global(terra::mask(sink_ifm,terra::ifel(smap==2,1,NA)),"sum",na.rm=TRUE)$sum +
                 terra::global(terra::mask(sink_mf,terra::ifel(smap==3,1,NA)),"sum",na.rm=TRUE)$sum)/1e6
  total <- (terra::global(terra::mask(seq_cn,terra::ifel(smap==1,1,NA)),"sum",na.rm=TRUE)$sum +
              terra::global(terra::mask(seq_ifm,terra::ifel(smap==2,1,NA)),"sum",na.rm=TRUE)$sum +
              terra::global(terra::mask(seq_mf,terra::ifel(smap==3,1,NA)),"sum",na.rm=TRUE)$sum)/1e6
  sweep <- rbind(sweep, data.frame(DF=DF, feasible=TRUE,
                                   CN=round(100*a/n,1), IFM=round(100*b/n,1), MF=round(100*d/n,1),
                                   coverage=round(bio$cov,1), unprotected=bio$unprot,
                                   sink_Mt=round(sink_tot,2), total_Mt=round(total,2)))
  terra::writeRaster(smap, file.path(out_dir, sprintf("scenario3_DF%02.0f_zones.tif", DF*100)), overwrite=TRUE)
}
cat("\n=================  S3 DF SWEEP RESULTS  =================\n")
cat("(DF = 0, sink only, is scenario S1: infeasible at the 91.0 floor)\n")
print(sweep, row.names = FALSE)
writexl::write_xlsx(sweep, file.path(out_dir, "scenario3_DF_sweep.xlsx"))
cat("\nSaved: scenario3_DF_sweep.xlsx\n")
