#==============================================================================
# 1. Load libraries
#============================================ ==================================
library(raster)
library(caret)
library(dplyr)
library(stringr)

#==============================================================================
# 2. Define settings
#==============================================================================
esms <- c("MIROC-ES2L", "MRI-ESM2-0")

esms <- c("MIROC-ES2L")
# esms <- c("MRI-ESM2-0")

# scenarios <- c("past1000", "ssp126", "ssp585")
scenarios <- c("historical", "past1000", "ssp126", "ssp585")
scenarios <- c("ssp585")


# algorithms <- c("RF", "GBM", "SVM")

algorithms <- c("RF")


# Define base paths
stack_base <- "/Volumes/MyPassport/2024_SouthernRightWhaleLogbooks/02envStacks"
model_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/02modelling"
output_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/03predictions"


# Define base paths
stack_base <- "/Volumes/MyPassport/2024_SouthernRightWhaleLogbooks/03envStacksStandardized"
model_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/02modelling/01standardized/MIROC-ES2L/01FitModels-noDepth/"
output_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/04standardizedPredictions-noDepth"

# Define base paths
stack_base <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/03envStacksStandardized"
model_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/02modelling/01standardized/MIROC-ES2L/01FitModels-noDepth/"
output_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/04standardizedPredictions-noDepth"




# Define base paths
stack_base <- "/Volumes/MyPassport/2024_SouthernRightWhaleLogbooks/02envStacks"
model_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/02modelling/00raw/"
output_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/05RawPredictions"


# Define base paths
stack_base <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/02envStacks_Standardized_0s/"
model_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/02modelling/01standardized/"
output_base <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/03predict/01standardized"



predictors <- c(
  "SeaFloorSlope", "DistanceToShore", "SeaSurfaceTemperature", "TemperatureDepth",
  "TemperatureGradientSurface", "TemperatureGradientDepth", "TempDifferenceSurfaceDepth",
  "SeaSurfaceSalinity", "SalinityGradientSurface", "SalinityGradientDepth",
  "SeaIceConcentration", "SeaSurfaceHeight", "UCurrentSurface", "VCurrentSurface",
  "EddyKineticEnergy", "UCurrentDepth", "VCurrentDepth", 
  "WCurrentDepth", "UWindStress", "VWindStress", "WindStressMagnitude",
  "WindStressDirection"
)

predictors <- c("Bathymetry",
                "SeaFloorSlope", "DistanceToShore", "SeaSurfaceTemperature", 
                "TemperatureGradientSurface", "TempDifferenceSurfaceDepth",
                "SeaSurfaceSalinity", "SalinityGradientSurface", 
                "SeaSurfaceHeight", "UCurrentSurface", "VCurrentSurface",
                "EddyKineticEnergy", 
                "UWindStress", "VWindStress", "WindStressMagnitude",
                "WindStressDirection"
)


#==============================================================================
# 3. Loop over ESMs, scenarios, and algorithms
#==============================================================================
for (esm in esms) {
  for (scenario in scenarios) {
    # List all .grd stack files in the appropriate folder
    stack_dir <- file.path(stack_base, esm, scenario)
    stack_files <- list.files(stack_dir, pattern = "\\.grd$", full.names = TRUE)
    # Skip problematic file (20220516)
    stack_files <- stack_files[!grepl("20220516", stack_files)]
    
    if (length(stack_files) == 0) {
      message("No .grd files found in: ", stack_dir)
      next
    }
    
    for (algo in algorithms) {
      # Load trained model
      # model_path <- file.path(model_base, esm, "01FitModels", algo, paste0(algo, "_model.rds"))
      model_path <- file.path(model_base, esm, "01FitModels-noDepth", algo, paste0(algo, "_model.rds"))
      
      if (!file.exists(model_path)) {
        message("Model not found: ", model_path)
        next
      }
      model <- readRDS(model_path)
      
      # predictors <- model$finalModel$xNames
      
      # Create output directory
      pred_outdir <- file.path(output_base, esm, algo, scenario)
      if (!dir.exists(pred_outdir)) dir.create(pred_outdir, recursive = TRUE)
      
      # Loop through available monthly stacks
      for (stack_path in stack_files) {
        fname <- basename(stack_path)
        
        # Match by filename pattern
        # if (grepl(paste0("monthly_stack_stand_", esm, "_", scenario), fname)) {
        if (grepl(paste0("monthly_stack_stand_", esm, "_", scenario), fname)) {
          
          date_str <- str_extract(fname, "\\d{8}")
          
          message("Predicting for ", esm, " | ", scenario, " | ", algo, " | ", date_str)
          
          # Load raster stack and filter variables
          s <- stack(stack_path)
          available_vars <- names(s)
          common_vars <- predictors[predictors %in% available_vars]
          
          if (length(common_vars) == 0) {
            message("No matching variables found in ", fname)
            next
          }
          
          s <- subset(s, common_vars)
          
          # Predict probability of presence ("pres")
          p <- raster::predict(s, model, type = "prob", index = 1, progress = "text")
          
          # Save prediction raster
          # output_file <- file.path(pred_outdir, paste0(date_str, "_", esm, "_", algo, "_standardized_prediction",  ".tif"))
          output_file <- file.path(pred_outdir, paste0(date_str, "_", esm, "_", algo, "_prediction",  ".tif"))
          writeRaster(p, output_file, overwrite = TRUE)
        }
      }
    }
  }
}
