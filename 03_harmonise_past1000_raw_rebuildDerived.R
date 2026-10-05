# ==============================================================================
# SOUTHERN RIGHT WHALE
# CREATE HARMONISED past1000 RAW STACKS
#
# INPUT - READ ONLY:
#   01envStacks_Raw_0s
#
# MONTHLY DELTAS:
#   02_ESM_harmonisation
#
# OUTPUT - NEW FILES ONLY:
#   03envStacks_Raw_Harmonised_v1
#
# STRATEGY
# --------
# 1. Read original past1000 RAW stack.
# 2. Apply monthly spatial delta ONLY to base dynamic variables.
# 3. NEVER modify static variables.
# 4. Recalculate derived variables from the corrected base variables,
#    following the original stacksSRW.R formulas.
# 5. Preserve original layer names, layer count and layer order.
# 6. Write a completely new .grd/.gri stack.
#
# ORIGINAL RAW DIRECTORY IS NEVER WRITTEN TO.
# ==============================================================================


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(terra)
library(dplyr)
library(stringr)


# ==============================================================================
# 1. PATHS
# ==============================================================================

RAW_ROOT <- paste0(
  "/Volumes/23327-15324-35628/",
  "2024_SouthernRightWhaleLogbooks/",
  "StackRestoration/",
  "01envStacks_Raw_0s"
)


DELTA_ROOT <- paste0(
  "/Volumes/23327-15324-35628/",
  "2024_SouthernRightWhaleLogbooks/",
  "StackRestoration/",
  "02_ESM_harmonisation"
)


OUTPUT_ROOT <- paste0(
  "/Volumes/23327-15324-35628/",
  "2024_SouthernRightWhaleLogbooks/",
  "StackRestoration/",
  "03envStacks_Raw_Harmonised_v1"
)


# ==============================================================================
# 2. ABSOLUTE SAFETY FUNCTIONS
# ==============================================================================

clean_path <- function(path) {
  
  path <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )
  
  sub(
    "/+$",
    "",
    path
  )
}


RAW_ROOT_CLEAN <- clean_path(
  RAW_ROOT
)


OUTPUT_ROOT_CLEAN <- clean_path(
  OUTPUT_ROOT
)


is_inside_raw <- function(path) {
  
  path_clean <- clean_path(
    path
  )
  
  identical(
    path_clean,
    RAW_ROOT_CLEAN
  ) ||
    startsWith(
      path_clean,
      paste0(
        RAW_ROOT_CLEAN,
        "/"
      )
    )
}


assert_safe_output <- function(path) {
  
  if (is_inside_raw(path)) {
    
    stop(
      paste0(
        "\n\n",
        "============================================================\n",
        "CRITICAL SAFETY STOP\n",
        "============================================================\n",
        "Attempted write inside ORIGINAL RAW directory:\n\n",
        path,
        "\n\n",
        "01envStacks_Raw_0s IS READ-ONLY.\n",
        "Nothing will be written there.\n",
        "============================================================\n"
      )
    )
    
  }
  
  invisible(
    TRUE
  )
}


if (
  identical(
    RAW_ROOT_CLEAN,
    OUTPUT_ROOT_CLEAN
  ) ||
  startsWith(
    OUTPUT_ROOT_CLEAN,
    paste0(
      RAW_ROOT_CLEAN,
      "/"
    )
  )
) {
  
  stop(
    "OUTPUT_ROOT is inside RAW_ROOT. Refusing to continue."
  )
}


assert_safe_output(
  OUTPUT_ROOT
)


# ==============================================================================
# 3. SETTINGS
# ==============================================================================

ESMS <- c(
  "MIROC-ES2L",
  "MRI-ESM2-0"
)


# ------------------------------------------------------------------------------
# STATIC VARIABLES
#
# NEVER corrected.
# NEVER recalculated here.
# Copied exactly from each RAW stack.
# ------------------------------------------------------------------------------

STATIC_VARIABLES <- c(
  "Bathymetry",
  "SeaFloorSlope",
  "DistanceToShore"
)


# ------------------------------------------------------------------------------
# BASE DYNAMIC VARIABLES
#
# These receive the monthly spatial delta when that delta exists.
#
# Some variables do not exist in both ESMs. That is fine.
# ------------------------------------------------------------------------------

BASE_DYNAMIC_VARIABLES <- c(
  "SeaSurfaceTemperature",
  "TemperatureSurface",
  "TemperatureDepth",
  "MixedLayerDepth",
  "SeaSurfaceSalinity",
  "SalinitySurface",
  "SalinityDepth",
  "SeaIceConcentration",
  "SeaSurfaceHeight",
  "UCurrentSurface",
  "VCurrentSurface",
  "UCurrentDepth",
  "VCurrentDepth",
  "WCurrentSurface",
  "WCurrentDepth",
  "UWindStress",
  "VWindStress"
)


# ------------------------------------------------------------------------------
# DERIVED VARIABLES
#
# DO NOT receive delta directly.
# Recalculated from corrected base variables.
# ------------------------------------------------------------------------------

DERIVED_VARIABLES <- c(
  "TemperatureGradientSurface",
  "TemperatureGradientDepth",
  "TempDifferenceSurfaceDepth",
  "SalinityGradientSurface",
  "SalinityGradientDepth",
  "DistanceToIceEdge",
  "EddyKineticEnergy",
  "WindStressMagnitude",
  "WindStressDirection"
)


# ------------------------------------------------------------------------------
# Output behaviour
#
# FALSE = safest.
# If rerun, an already-created harmonised stack is skipped.
# ------------------------------------------------------------------------------

OVERWRITE_OUTPUTS <- FALSE


# ==============================================================================
# 4. TERRA SETTINGS
# ==============================================================================

terra_temp <- file.path(
  tempdir(),
  "SRW_harmonised_stacks_v1"
)


dir.create(
  terra_temp,
  recursive = TRUE,
  showWarnings = FALSE
)


terra::terraOptions(
  tempdir = terra_temp,
  progress = 1
)


# ==============================================================================
# 5. EXTRACT DATE FROM STACK FILENAME
# ==============================================================================

extract_date <- function(file) {
  
  date_string <- stringr::str_extract(
    basename(file),
    "[0-9]{8}(?=\\.grd$)"
  )
  
  
  as.Date(
    date_string,
    format = "%Y%m%d"
  )
}


# ==============================================================================
# 6. INDEX RAW STACKS
# ==============================================================================

build_file_index <- function(directory) {
  
  files <- list.files(
    directory,
    pattern = "\\.grd$",
    full.names = TRUE
  )
  
  
  if (length(files) == 0) {
    
    stop(
      paste0(
        "No .grd files found in:\n",
        directory
      )
    )
    
  }
  
  
  dates <- extract_date(
    files
  )
  
  
  data.frame(
    file = files,
    date = dates,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      !is.na(date)
    ) %>%
    dplyr::mutate(
      
      Year = as.integer(
        format(
          date,
          "%Y"
        )
      ),
      
      Month = as.integer(
        format(
          date,
          "%m"
        )
      )
      
    ) %>%
    dplyr::arrange(
      date
    )
}


# ==============================================================================
# 7. GET MONTHLY DELTA FILE
# ==============================================================================

get_delta_file <- function(
    esm,
    month
) {
  
  month_label <- sprintf(
    "%02d",
    month
  )
  
  
  delta_file <- file.path(
    DELTA_ROOT,
    esm,
    "delta_layers",
    paste0(
      esm,
      "_delta_historical_minus_past1000_month_",
      month_label,
      ".grd"
    )
  )
  
  
  if (!file.exists(delta_file)) {
    
    stop(
      paste0(
        "\nMissing monthly delta:\n",
        delta_file
      )
    )
    
  }
  
  
  delta_file
}


# ==============================================================================
# 8. HELPER - RECALCULATE DISTANCE TO ICE EDGE
#
# Original rule:
#
#   SeaIceConcentration >= 15
#
# Then calculate distance to ice.
#
# For masking, the original stack-building script used sample_layer derived
# from surface salinity. We therefore preferentially use SeaSurfaceSalinity
# as the valid-domain mask when available.
# ==============================================================================

recalculate_distance_to_ice <- function(
    sea_ice,
    domain_mask
) {
  
  sea_ice_binary <- terra::ifel(
    sea_ice >= 15,
    1,
    NA
  )
  
  
  n_ice_cells <- terra::global(
    !is.na(
      sea_ice_binary
    ),
    fun = "sum",
    na.rm = TRUE
  )[1, 1]
  
  
  if (
    is.na(n_ice_cells) ||
    n_ice_cells == 0
  ) {
    
    distance_to_edge <- sea_ice
    
    terra::values(
      distance_to_edge
    ) <- NA
    
  } else {
    
    distance_to_edge <- terra::distance(
      sea_ice_binary,
      unit = "km"
    )
    
  }
  
  
  distance_to_edge <- terra::mask(
    distance_to_edge,
    domain_mask
  )
  
  
  names(
    distance_to_edge
  ) <- "DistanceToIceEdge"
  
  
  distance_to_edge
}


# ==============================================================================
# 9. CREATE OUTPUT ROOT
# ==============================================================================

assert_safe_output(
  OUTPUT_ROOT
)


dir.create(
  OUTPUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 10. LOG CONTAINERS
# ==============================================================================

file_log_list <- list()

variable_log_list <- list()

file_log_counter <- 1

variable_log_counter <- 1


# ==============================================================================
# 11. PROCESS EACH ESM INDEPENDENTLY
# ==============================================================================

for (esm in ESMS) {
  
  cat("\n\n")
  cat("################################################################\n")
  cat("################################################################\n")
  cat("\n")
  cat("PROCESSING HARMONISED past1000: ", esm, "\n", sep = "")
  cat("\n")
  cat("################################################################\n")
  cat("################################################################\n")
  
  
  # ============================================================================
  # 11.1 RAW INPUT
  # ============================================================================
  
  raw_past_dir <- file.path(
    RAW_ROOT,
    esm,
    "past1000"
  )
  
  
  if (!dir.exists(raw_past_dir)) {
    
    stop(
      paste0(
        "RAW past1000 directory not found:\n",
        raw_past_dir
      )
    )
    
  }
  
  
  # ============================================================================
  # 11.2 NEW OUTPUT DIRECTORY
  # ============================================================================
  
  output_past_dir <- file.path(
    OUTPUT_ROOT,
    esm,
    "past1000"
  )
  
  
  output_log_dir <- file.path(
    OUTPUT_ROOT,
    esm,
    "logs"
  )
  
  
  assert_safe_output(
    output_past_dir
  )
  
  
  assert_safe_output(
    output_log_dir
  )
  
  
  dir.create(
    output_past_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    output_log_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ============================================================================
  # 11.3 INDEX ALL past1000 MONTHLY STACKS
  # ============================================================================
  
  raw_index <- build_file_index(
    raw_past_dir
  )
  
  
  cat(
    "\nNumber of RAW past1000 stacks: ",
    nrow(raw_index),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Date range: ",
    min(raw_index$date),
    " -> ",
    max(raw_index$date),
    "\n",
    sep = ""
  )
  
  
  # ============================================================================
  # 11.4 LOOP OVER ALL past1000 MONTHS
  # ============================================================================
  
  for (i in seq_len(
    nrow(raw_index)
  )) {
    
    raw_file <- raw_index$file[i]
    
    this_date <- raw_index$date[i]
    
    this_year <- raw_index$Year[i]
    
    this_month <- raw_index$Month[i]
    
    
    cat("\n")
    cat("============================================================\n")
    
    cat(
      esm,
      " | ",
      as.character(this_date),
      " | ",
      i,
      "/",
      nrow(raw_index),
      "\n",
      sep = ""
    )
    
    cat("============================================================\n")
    
    
    # --------------------------------------------------------------------------
    # Monthly delta
    # --------------------------------------------------------------------------
    
    delta_file <- get_delta_file(
      esm = esm,
      month = this_month
    )
    
    
    # --------------------------------------------------------------------------
    # Output has the SAME filename but lives in NEW directory
    # --------------------------------------------------------------------------
    
    output_file <- file.path(
      output_past_dir,
      basename(raw_file)
    )
    
    
    assert_safe_output(
      output_file
    )
    
    
    # --------------------------------------------------------------------------
    # Skip existing derived files
    # --------------------------------------------------------------------------
    
    if (
      file.exists(output_file) &&
      !OVERWRITE_OUTPUTS
    ) {
      
      cat(
        "Output already exists -> skipped.\n"
      )
      
      next
    }
    
    
    # ==========================================================================
    # 11.5 READ ORIGINAL RAW + DELTA
    # ==========================================================================
    
    raw_stack <- terra::rast(
      raw_file
    )
    
    
    delta_stack <- terra::rast(
      delta_file
    )
    
    
    # ==========================================================================
    # 11.6 GEOMETRY CHECK
    # ==========================================================================
    
    geometry_ok <- terra::compareGeom(
      raw_stack,
      delta_stack,
      stopOnError = FALSE
    )
    
    
    if (!geometry_ok) {
      
      stop(
        paste0(
          "\nGeometry mismatch:\n",
          raw_file,
          "\nvs\n",
          delta_file
        )
      )
      
    }
    
    
    raw_variables <- names(
      raw_stack
    )
    
    
    delta_variables <- names(
      delta_stack
    )
    
    
    # ==========================================================================
    # 11.7 CREATE WORKING COPY IN MEMORY
    #
    # This does NOT alter the disk RAW.
    # ==========================================================================
    
    working_stack <- raw_stack
    
    
    corrected_base_variables <- character()
    
    base_variables_without_delta <- character()
    
    
    # ==========================================================================
    # 11.8 CORRECT BASE DYNAMIC VARIABLES
    # ==========================================================================
    
    cat("\nBASE DYNAMIC VARIABLES\n")
    
    
    for (variable in BASE_DYNAMIC_VARIABLES) {
      
      # ------------------------------------------------------------------------
      # Variable does not exist in this RAW stack
      # ------------------------------------------------------------------------
      
      if (!(variable %in% raw_variables)) {
        
        next
      }
      
      
      # ------------------------------------------------------------------------
      # Variable exists, delta exists -> CORRECT
      # ------------------------------------------------------------------------
      
      if (variable %in% delta_variables) {
        
        cat(
          "  ",
          variable,
          " -> applying monthly delta\n",
          sep = ""
        )
        
        
        corrected_layer <- raw_stack[[variable]] +
          delta_stack[[variable]]
        
        
        names(
          corrected_layer
        ) <- variable
        
        
        working_stack[[variable]] <- corrected_layer
        
        
        corrected_base_variables <- c(
          corrected_base_variables,
          variable
        )
        
        
        variable_log_list[[variable_log_counter]] <- data.frame(
          ESM = esm,
          Date = as.character(this_date),
          Variable = variable,
          Type = "base_dynamic",
          Action = "monthly_delta_applied",
          stringsAsFactors = FALSE
        )
        
        
        variable_log_counter <- variable_log_counter + 1
        
      } else {
        
        # ----------------------------------------------------------------------
        # RAW variable exists but no delta exists.
        #
        # Keep original RAW layer.
        # ----------------------------------------------------------------------
        
        cat(
          "  ",
          variable,
          " -> NO delta available, keeping RAW unchanged\n",
          sep = ""
        )
        
        
        base_variables_without_delta <- c(
          base_variables_without_delta,
          variable
        )
        
        
        variable_log_list[[variable_log_counter]] <- data.frame(
          ESM = esm,
          Date = as.character(this_date),
          Variable = variable,
          Type = "base_dynamic",
          Action = "no_delta_raw_preserved",
          stringsAsFactors = FALSE
        )
        
        
        variable_log_counter <- variable_log_counter + 1
      }
    }
    
    
    # ==========================================================================
    # 11.9 STATIC VARIABLES
    # ==========================================================================
    
    cat("\nSTATIC VARIABLES\n")
    
    
    for (variable in STATIC_VARIABLES) {
      
      if (variable %in% raw_variables) {
        
        cat(
          "  ",
          variable,
          " -> unchanged\n",
          sep = ""
        )
        
        
        working_stack[[variable]] <- raw_stack[[variable]]
        
        
        variable_log_list[[variable_log_counter]] <- data.frame(
          ESM = esm,
          Date = as.character(this_date),
          Variable = variable,
          Type = "static",
          Action = "raw_preserved",
          stringsAsFactors = FALSE
        )
        
        
        variable_log_counter <- variable_log_counter + 1
      }
    }
    
    
    # ==========================================================================
    # 11.10 RECALCULATE DERIVED VARIABLES
    # ==========================================================================
    
    cat("\nDERIVED VARIABLES\n")
    
    
    recalculated_derived <- character()
    
    preserved_derived <- character()
    
    
    # --------------------------------------------------------------------------
    # EDDY KINETIC ENERGY
    #
    # EKE = (UCurrentSurface^2 + VCurrentSurface^2) / 2
    # --------------------------------------------------------------------------
    
    if ("EddyKineticEnergy" %in% raw_variables) {
      
      if (
        all(
          c(
            "UCurrentSurface",
            "VCurrentSurface"
          ) %in% raw_variables
        )
      ) {
        
        x <- (
          working_stack[["UCurrentSurface"]]^2 +
            working_stack[["VCurrentSurface"]]^2
        ) / 2
        
        
        names(
          x
        ) <- "EddyKineticEnergy"
        
        
        working_stack[["EddyKineticEnergy"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "EddyKineticEnergy"
        )
        
        
        cat(
          "  EddyKineticEnergy -> recalculated\n"
        )
        
      } else {
        
        working_stack[["EddyKineticEnergy"]] <- raw_stack[["EddyKineticEnergy"]]
        
        preserved_derived <- c(
          preserved_derived,
          "EddyKineticEnergy"
        )
        
        cat(
          "  EddyKineticEnergy -> dependencies missing, RAW preserved\n"
        )
      }
    }
    
    
    # --------------------------------------------------------------------------
    # TEMPERATURE GRADIENT SURFACE
    # --------------------------------------------------------------------------
    
    if ("TemperatureGradientSurface" %in% raw_variables) {
      
      if ("TemperatureSurface" %in% raw_variables) {
        
        x <- terra::terrain(
          working_stack[["TemperatureSurface"]],
          v = "slope",
          unit = "degrees"
        )
        
        
        names(
          x
        ) <- "TemperatureGradientSurface"
        
        
        working_stack[["TemperatureGradientSurface"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "TemperatureGradientSurface"
        )
        
        
        cat(
          "  TemperatureGradientSurface -> recalculated\n"
        )
        
      } else {
        
        working_stack[["TemperatureGradientSurface"]] <- raw_stack[["TemperatureGradientSurface"]]
        
        preserved_derived <- c(
          preserved_derived,
          "TemperatureGradientSurface"
        )
      }
    }
    
    
    # --------------------------------------------------------------------------
    # TEMPERATURE GRADIENT DEPTH
    # --------------------------------------------------------------------------
    
    if ("TemperatureGradientDepth" %in% raw_variables) {
      
      if ("TemperatureDepth" %in% raw_variables) {
        
        x <- terra::terrain(
          working_stack[["TemperatureDepth"]],
          v = "slope",
          unit = "degrees"
        )
        
        
        names(
          x
        ) <- "TemperatureGradientDepth"
        
        
        working_stack[["TemperatureGradientDepth"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "TemperatureGradientDepth"
        )
        
        
        cat(
          "  TemperatureGradientDepth -> recalculated\n"
        )
        
      } else {
        
        working_stack[["TemperatureGradientDepth"]] <- raw_stack[["TemperatureGradientDepth"]]
        
        preserved_derived <- c(
          preserved_derived,
          "TemperatureGradientDepth"
        )
      }
    }
    
    
    # --------------------------------------------------------------------------
    # SALINITY GRADIENT SURFACE
    # --------------------------------------------------------------------------
    
    if ("SalinityGradientSurface" %in% raw_variables) {
      
      if ("SalinitySurface" %in% raw_variables) {
        
        x <- terra::terrain(
          working_stack[["SalinitySurface"]],
          v = "slope",
          unit = "degrees"
        )
        
        
        names(
          x
        ) <- "SalinityGradientSurface"
        
        
        working_stack[["SalinityGradientSurface"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "SalinityGradientSurface"
        )
        
        
        cat(
          "  SalinityGradientSurface -> recalculated\n"
        )
        
      } else {
        
        working_stack[["SalinityGradientSurface"]] <- raw_stack[["SalinityGradientSurface"]]
        
        preserved_derived <- c(
          preserved_derived,
          "SalinityGradientSurface"
        )
      }
    }
    
    
    # --------------------------------------------------------------------------
    # SALINITY GRADIENT DEPTH
    # --------------------------------------------------------------------------
    
    if ("SalinityGradientDepth" %in% raw_variables) {
      
      if ("SalinityDepth" %in% raw_variables) {
        
        x <- terra::terrain(
          working_stack[["SalinityDepth"]],
          v = "slope",
          unit = "degrees"
        )
        
        
        names(
          x
        ) <- "SalinityGradientDepth"
        
        
        working_stack[["SalinityGradientDepth"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "SalinityGradientDepth"
        )
        
        
        cat(
          "  SalinityGradientDepth -> recalculated\n"
        )
        
      } else {
        
        working_stack[["SalinityGradientDepth"]] <- raw_stack[["SalinityGradientDepth"]]
        
        preserved_derived <- c(
          preserved_derived,
          "SalinityGradientDepth"
        )
      }
    }
    
    
    # --------------------------------------------------------------------------
    # TEMPERATURE DIFFERENCE SURFACE - DEPTH
    # --------------------------------------------------------------------------
    
    if ("TempDifferenceSurfaceDepth" %in% raw_variables) {
      
      if (
        all(
          c(
            "TemperatureSurface",
            "TemperatureDepth"
          ) %in% raw_variables
        )
      ) {
        
        x <- working_stack[["TemperatureSurface"]] -
          working_stack[["TemperatureDepth"]]
        
        
        names(
          x
        ) <- "TempDifferenceSurfaceDepth"
        
        
        working_stack[["TempDifferenceSurfaceDepth"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "TempDifferenceSurfaceDepth"
        )
        
        
        cat(
          "  TempDifferenceSurfaceDepth -> recalculated\n"
        )
        
      } else {
        
        working_stack[["TempDifferenceSurfaceDepth"]] <- raw_stack[["TempDifferenceSurfaceDepth"]]
        
        preserved_derived <- c(
          preserved_derived,
          "TempDifferenceSurfaceDepth"
        )
      }
    }
    
    
    # --------------------------------------------------------------------------
    # WIND STRESS MAGNITUDE
    # --------------------------------------------------------------------------
    
    if ("WindStressMagnitude" %in% raw_variables) {
      
      if (
        all(
          c(
            "UWindStress",
            "VWindStress"
          ) %in% raw_variables
        )
      ) {
        
        x <- sqrt(
          working_stack[["UWindStress"]]^2 +
            working_stack[["VWindStress"]]^2
        )
        
        
        names(
          x
        ) <- "WindStressMagnitude"
        
        
        working_stack[["WindStressMagnitude"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "WindStressMagnitude"
        )
        
        
        cat(
          "  WindStressMagnitude -> recalculated\n"
        )
        
      } else {
        
        working_stack[["WindStressMagnitude"]] <- raw_stack[["WindStressMagnitude"]]
        
        preserved_derived <- c(
          preserved_derived,
          "WindStressMagnitude"
        )
      }
    }
    
    
    # --------------------------------------------------------------------------
    # WIND STRESS DIRECTION
    #
    # atan2(V, U) -> degrees -> 0:360
    # --------------------------------------------------------------------------
    
    if ("WindStressDirection" %in% raw_variables) {
      
      if (
        all(
          c(
            "UWindStress",
            "VWindStress"
          ) %in% raw_variables
        )
      ) {
        
        x <- atan2(
          working_stack[["VWindStress"]],
          working_stack[["UWindStress"]]
        ) * (
          180 / pi
        )
        
        
        x <- terra::ifel(
          x < 0,
          x + 360,
          x
        )
        
        
        names(
          x
        ) <- "WindStressDirection"
        
        
        working_stack[["WindStressDirection"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "WindStressDirection"
        )
        
        
        cat(
          "  WindStressDirection -> recalculated\n"
        )
        
      } else {
        
        working_stack[["WindStressDirection"]] <- raw_stack[["WindStressDirection"]]
        
        preserved_derived <- c(
          preserved_derived,
          "WindStressDirection"
        )
      }
    }
    
    
    # --------------------------------------------------------------------------
    # DISTANCE TO ICE EDGE
    #
    # Ice = SeaIceConcentration >= 15
    #
    # Prefer SeaSurfaceSalinity as valid-domain mask because the original
    # workflow constructed sample_layer from SOS.
    # --------------------------------------------------------------------------
    
    if ("DistanceToIceEdge" %in% raw_variables) {
      
      if ("SeaIceConcentration" %in% raw_variables) {
        
        if ("SeaSurfaceSalinity" %in% raw_variables) {
          
          domain_mask <- working_stack[["SeaSurfaceSalinity"]]
          
        } else {
          
          domain_mask <- working_stack[["SeaIceConcentration"]]
          
        }
        
        
        x <- recalculate_distance_to_ice(
          sea_ice = working_stack[["SeaIceConcentration"]],
          domain_mask = domain_mask
        )
        
        
        working_stack[["DistanceToIceEdge"]] <- x
        
        
        recalculated_derived <- c(
          recalculated_derived,
          "DistanceToIceEdge"
        )
        
        
        cat(
          "  DistanceToIceEdge -> recalculated\n"
        )
        
      } else {
        
        working_stack[["DistanceToIceEdge"]] <- raw_stack[["DistanceToIceEdge"]]
        
        preserved_derived <- c(
          preserved_derived,
          "DistanceToIceEdge"
        )
      }
    }
    
    
    # ==========================================================================
    # 11.11 VARIABLE-LEVEL DERIVED LOG
    # ==========================================================================
    
    if (length(recalculated_derived) > 0) {
      
      for (variable in recalculated_derived) {
        
        variable_log_list[[variable_log_counter]] <- data.frame(
          ESM = esm,
          Date = as.character(this_date),
          Variable = variable,
          Type = "derived",
          Action = "recalculated_from_harmonised_base",
          stringsAsFactors = FALSE
        )
        
        
        variable_log_counter <- variable_log_counter + 1
      }
    }
    
    
    if (length(preserved_derived) > 0) {
      
      for (variable in preserved_derived) {
        
        variable_log_list[[variable_log_counter]] <- data.frame(
          ESM = esm,
          Date = as.character(this_date),
          Variable = variable,
          Type = "derived",
          Action = "dependencies_missing_raw_preserved",
          stringsAsFactors = FALSE
        )
        
        
        variable_log_counter <- variable_log_counter + 1
      }
    }
    
    
    # ==========================================================================
    # 11.12 FORCE ORIGINAL RAW ORDER
    #
    # Very important:
    #   same variables
    #   same layer count
    #   same order
    # ==========================================================================
    
    working_stack <- working_stack[[raw_variables]]
    
    
    names(
      working_stack
    ) <- raw_variables
    
    
    if (
      terra::nlyr(
        working_stack
      ) !=
      terra::nlyr(
        raw_stack
      )
    ) {
      
      stop(
        paste0(
          "Layer-count mismatch after harmonisation:\n",
          raw_file
        )
      )
    }
    
    
    if (!identical(
      names(working_stack),
      names(raw_stack)
    )) {
      
      stop(
        paste0(
          "Layer-name/order mismatch after harmonisation:\n",
          raw_file
        )
      )
    }
    
    
    # ==========================================================================
    # 11.13 FINAL WRITE SAFETY
    # ==========================================================================
    
    assert_safe_output(
      output_file
    )
    
    
    if (is_inside_raw(output_file)) {
      
      stop(
        "CRITICAL SAFETY STOP: output points to RAW directory."
      )
    }
    
    
    # ==========================================================================
    # 11.14 WRITE NEW HARMONISED STACK
    # ==========================================================================
    
    terra::writeRaster(
      working_stack,
      filename = output_file,
      overwrite = OVERWRITE_OUTPUTS,
      filetype = "RRASTER"
    )
    
    
    cat(
      "\nNEW harmonised file written:\n  ",
      output_file,
      "\n",
      sep = ""
    )
    
    
    # ==========================================================================
    # 11.15 FILE-LEVEL LOG
    # ==========================================================================
    
    file_log_list[[file_log_counter]] <- data.frame(
      
      ESM = esm,
      
      Date = as.character(
        this_date
      ),
      
      Year = this_year,
      
      Month = this_month,
      
      RawFile = basename(
        raw_file
      ),
      
      OutputFile = basename(
        output_file
      ),
      
      NRawLayers = terra::nlyr(
        raw_stack
      ),
      
      NBaseVariablesCorrected = length(
        corrected_base_variables
      ),
      
      CorrectedBaseVariables = paste(
        corrected_base_variables,
        collapse = "; "
      ),
      
      NBaseVariablesWithoutDelta = length(
        base_variables_without_delta
      ),
      
      BaseVariablesWithoutDelta = paste(
        base_variables_without_delta,
        collapse = "; "
      ),
      
      NDerivedRecalculated = length(
        recalculated_derived
      ),
      
      DerivedRecalculated = paste(
        recalculated_derived,
        collapse = "; "
      ),
      
      NDerivedPreserved = length(
        preserved_derived
      ),
      
      DerivedPreserved = paste(
        preserved_derived,
        collapse = "; "
      ),
      
      stringsAsFactors = FALSE
      
    )
    
    
    file_log_counter <- file_log_counter + 1
    
    
    # ==========================================================================
    # 11.16 MEMORY CLEANUP
    # ==========================================================================
    
    rm(
      raw_stack,
      delta_stack,
      working_stack
    )
    
    gc()
  }
  
  
  cat("\n")
  cat("============================================================\n")
  cat("FINISHED ESM: ", esm, "\n", sep = "")
  cat("============================================================\n")
}


# ==============================================================================
# 12. SAVE LOGS
# ==============================================================================

file_log <- dplyr::bind_rows(
  file_log_list
)


variable_log <- dplyr::bind_rows(
  variable_log_list
)


global_file_log <- file.path(
  OUTPUT_ROOT,
  "harmonisation_file_log.csv"
)


global_variable_log <- file.path(
  OUTPUT_ROOT,
  "harmonisation_variable_log.csv"
)


assert_safe_output(
  global_file_log
)


assert_safe_output(
  global_variable_log
)


write.csv(
  file_log,
  global_file_log,
  row.names = FALSE
)


write.csv(
  variable_log,
  global_variable_log,
  row.names = FALSE
)


# ==============================================================================
# 13. SAVE SEPARATE LOGS PER ESM
# ==============================================================================

for (esm in ESMS) {
  
  esm_log_dir <- file.path(
    OUTPUT_ROOT,
    esm,
    "logs"
  )
  
  
  assert_safe_output(
    esm_log_dir
  )
  
  
  dir.create(
    esm_log_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  esm_file_log <- file_log %>%
    dplyr::filter(
      ESM == esm
    )
  
  
  esm_variable_log <- variable_log %>%
    dplyr::filter(
      ESM == esm
    )
  
  
  write.csv(
    esm_file_log,
    file.path(
      esm_log_dir,
      paste0(
        esm,
        "_harmonisation_file_log.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    esm_variable_log,
    file.path(
      esm_log_dir,
      paste0(
        esm,
        "_harmonisation_variable_log.csv"
      )
    ),
    row.names = FALSE
  )
}


# ==============================================================================
# 14. FINAL REPORT
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("#                                                              #\n")
cat("#              HARMONISED past1000 COMPLETE                    #\n")
cat("#                                                              #\n")
cat("################################################################\n")


cat("\nORIGINAL RAW SOURCE - READ ONLY:\n")
cat(RAW_ROOT, "\n")


cat("\nNEW HARMONISED OUTPUT:\n")
cat(OUTPUT_ROOT, "\n")


cat("\nSTATIC VARIABLES LEFT UNCHANGED:\n")
print(
  STATIC_VARIABLES
)


cat("\nBASE VARIABLES CORRECTED WHEN DELTA AVAILABLE:\n")
print(
  BASE_DYNAMIC_VARIABLES
)


cat("\nDERIVED VARIABLES RECALCULATED:\n")
print(
  DERIVED_VARIABLES
)


cat("\nFILE LOG:\n")
cat(global_file_log, "\n")


cat("\nVARIABLE LOG:\n")
cat(global_variable_log, "\n")