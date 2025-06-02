# =============================================================================
# Single Tree Detection from 1m nDSM Data across Baden-Württemberg
# Goal: Estimate spatial tree distribution as a proxy for 10x10 m forest inventory maps
# =============================================================================

################################################################################
# Script Information
# Author: Eisnecker, Philipp
# Date: 2025/02/18
# Purpose: Large-scale single tree detection from photogrammetric nDSM data
#          across the federal state of Baden-Württemberg. Tree detection is 
#          performed using a dynamic local maxima filter (LMF) based on canopy 
#          height. The goal is to derive spatially averaged stand density as a 
#          proxy for 10x10 m forest inventory grid maps.
#
#          Detected trees are exported as shapefiles per year and region ("Los").
#          Window size for LMF is adaptively defined as a function of height.
#          The script processes data year-wise and handles multiple tiles.
#
################################################################################

# -----------------------------------------------------------------------------
# 1. Initialization
# -----------------------------------------------------------------------------

# Record the total processing start time
start_time_total <- Sys.time()
cat(paste0("Starting full processing at ", format(start_time_total), "\n\n"))

# Load required packages
library(terra)       # Raster data handling
library(lidR)        # Tree detection and CHM analysis
library(sf)          # Spatial vector data (e.g. Shapefiles)
library(supercells)  # (Optional) for segmentation; currently unused
library(dplyr)       # Data manipulation
library(fs)          # File path utilities

# -----------------------------------------------------------------------------
# 2. Define dynamic window size function for tree detection
# -----------------------------------------------------------------------------

# This function defines the moving window size based on raster cell height.
# Higher trees = larger crowns = larger search windows.
f <- function(x) {
  y <- 2.6 * (-(exp(-0.08 * (x - 2)) - 1)) + 3
  y[x < 2] <- 3     # Minimum window size for low vegetation
  y[x > 20] <- 5    # Cap window size for very tall trees
  return(y)
}

# Plot the window size function for visualization
heights <- seq(-5, 30, 0.5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0, 5))

# -----------------------------------------------------------------------------
# 3. Helper function to create output directories if missing
# -----------------------------------------------------------------------------

create_path <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    cat("Created path: ", path, "\n")
  } else {
    cat("Path already exists: ", path, "\n")
  }
}

# -----------------------------------------------------------------------------
# 4. Main processing function for one mission (one region/year)
# -----------------------------------------------------------------------------

process_mission <- function(mission_path, year) {
  los_folders <- list.files(
    mission_path,
    pattern = "^Los",
    full.names = TRUE,
    include.dirs = TRUE
  )
  
  # Loop over all "Los" folders (regional tiles)
  lapply(los_folders, function(folder) {
    ttops_list <- list()
    
    ndsm_path <- paste0(folder, "/SURE/nDSM/Clips")
    ndsm_files <- list.files(
      ndsm_path,
      pattern = "^nDSM_.*_1km\\.tif$",
      full.names = TRUE
    )
    
    # Tree detection for each 1x1 km nDSM tile
    ttops_located_list <- lapply(ndsm_files, function(ndsm_file) {
      chm <- rast(ndsm_file)           # Load nDSM as Canopy Height Model
      lmf_algorithm <- lmf(f)          # Local maxima filter with adaptive window
      ttops_located <- lidR::locate_trees(chm, lmf_algorithm)
      st_zm(ttops_located, drop = TRUE, what = "ZM")  # Drop Z/M from geometry
    })
    
    # Combine all detected tree tops for the region
    merged_sf <- bind_rows(ttops_located_list)
    
    # Extract mission number and prepare output path
    mission_no <- sub(".*Los(\\d{3})$", "\\1", folder)
    save_path <- paste0("S:/Projekte/VorratAktuell/Daten/Shapes/ttops/", year, "/Los", mission_no)
    create_path(save_path)
    
    # Write results to shapefile
    st_write(merged_sf, paste0(save_path, "/ttops_", year, "_Los", mission_no, ".shp"))
    return(NULL)
  })
  
  return(NULL)
}

# -----------------------------------------------------------------------------
# 5. Run the detection for all specified years
# -----------------------------------------------------------------------------

lapply(years, function(year) {
  start_time_year <- Sys.time()
  
  mission_path <- paste0("S:/Surface_UTM/", year, "/")
  process_mission(mission_path, year)
  
  # Logging processing time per year
  end_time_year <- Sys.time()
  duration_year <- end_time_year - start_time_year
  cat(paste0("Finished processing year: ", year, " at ", format(end_time_year), "\n"))
  cat(paste0("Processing time for year ", year, ": ", round(duration_year, 2), " ", units(duration_year), "\n\n"))
  
  return(NULL)
})

# -----------------------------------------------------------------------------
# 6. Finalize processing
# -----------------------------------------------------------------------------

end_time_total <- Sys.time()
duration_total <- end_time_total - start_time_total
cat(paste0("Total processing completed at ", format(end_time_total), "\n"))
cat(paste0("Total duration: ", round(duration_total, 2), " ", units(duration_total), "\n"))

# =============================================================================
# End of script
# =============================================================================
