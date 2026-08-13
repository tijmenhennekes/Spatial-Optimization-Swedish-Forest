# ==============================================================================
# Management-intensity transition map (S0 -> S2), GeoTIFF for styling in QGIS.
# Cell value = signed intensity of directional change between the baseline (S0)
# and the optimal policy solution (S2). Intensity rank: CN = 1 < IFM = 2 < MF = 3.
# delta = S2 - S0:
#   -2  MF -> CN   (strong de-intensification / protection gained)
#   -1  MF -> IFM  or  IFM -> CN   (one-step de-intensification)
#    0  unchanged
#   +1  CN -> IFM  or  IFM -> MF   (one-step intensification)
#   +2  CN -> MF   (strong intensification)
# Author: Tijmen Hennekes
# ==============================================================================
library(terra)
out_dir <- "C:/Users/Tijme/OneDrive/Msc ERM/Thesis work/00 Prioritizr data/04-outputs"

s0 <- terra::rast(file.path(out_dir, "scenario0_current_zones.tif"))             # S0: 1=CN 2=IFM 3=MF
s2 <- terra::rast(file.path(out_dir, "scenario1b_feasible_solution_zones.tif"))  # S2: 1=CN 2=IFM 3=MF

s2 <- terra::resample(s2, s0, method = "near")   # ensure identical grid
delta <- s2 - s0                                 # ranks ordered by intensity, so the difference is valid

terra::writeRaster(delta, file.path(out_dir, "transition_intensity.tif"),
                   overwrite = TRUE, datatype = "INT2S")

cat("Saved: transition_intensity.tif in", out_dir, "\n")
print(terra::freq(delta))   # cell counts per intensity class (for the caption)
