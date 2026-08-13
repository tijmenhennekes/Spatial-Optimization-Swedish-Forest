# ==============================================================================
# Management-intensity transition -> GeoTIFF (style in QGIS)
# Cell value = signed intensity of directional change between S0 and S1.
# Intensity rank: CN = 1 (protected) < IFM = 2 < MF = 3 (intensive).
# delta = S1 - S0:
#   -2  MF -> CN   (strong de-intensification / protection gained)
#   -1  MF -> IFM  or  IFM -> CN   (one-step de-intensification)
#    0  unchanged
#   +1  CN -> IFM  or  IFM -> MF   (one-step intensification)
#   +2  CN -> MF   (strong intensification)
# Author: Tijmen Hennekes
# ==============================================================================
library(terra)
base    <- "C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data"
out_dir <- file.path(base, "04-outputs")

master_mask <- terra::rast(file.path(base, "02_processed_data/forest_binary.tif"))
fm  <- terra::rast(file.path(base, "01_raw_data/fm_swe.tif"))            # S0: 1=CN 2=IFM 3=MF
s1  <- terra::rast(file.path(out_dir, "scenario1_solution_zones.tif"))   # S1: 1=CN 2=IFM 3=MF

s0r <- terra::mask(terra::resample(fm, master_mask, "near"), master_mask)
s1r <- terra::mask(terra::resample(s1, master_mask, "near"), master_mask)

# signed intensity of change (-2 .. +2); negative = de-intensified, positive = intensified
delta <- terra::mask(s1r - s0r, master_mask)

terra::writeRaster(delta, file.path(out_dir, "transition_intensity.tif"),
                   overwrite = TRUE, datatype = "INT2S")

cat("Saved: transition_intensity.tif in", out_dir, "\n")
print(terra::freq(delta))   # cell counts per intensity class (for the caption)