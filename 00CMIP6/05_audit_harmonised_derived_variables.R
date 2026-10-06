# ==============================================================================
# 05_audit_harmonised_derived_variables.R
#
# SRW CMIP6 / PMIP harmonisation audit
#
# GOALS
# -----
# 1) Verify, cell-by-cell, that EVERY derived variable stored in the harmonised
#    past1000 stacks was recalculated from the harmonised base variables using
#    exactly the same formulas as the original stack-building workflow.
#
# 2) Quantify the past1000 -> historical discontinuity for ALL derived variables,
#    before and after harmonisation.
#
# 3) Show the behaviour of the base variables used by each derived variable, so
#    we can distinguish:
#       A. wrong recalculation
#       B. correct recalculation but altered derived climatology
#
# IMPORTANT
# ---------
# * READ-ONLY for both RAW and harmonised stacks.
# * Writes ONLY to a new audit folder.
# * Does NOT use WCurrentSurface / WCurrentDepth / wo.
# * Uses the original SRW derived-variable formulas.
# ==============================================================================


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


# ==============================================================================
# 1. PATHS
# ==============================================================================

RAW_ROOT <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/01envStacks_Raw_0s"

HARM_ROOT <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/03envStacks_Raw_Harmonised_v1"

OUT_ROOT <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/05_audit_harmonised_derived_variables"

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

raw_norm <- normalizePath(RAW_ROOT, winslash = "/", mustWork = TRUE)
harm_norm <- normalizePath(HARM_ROOT, winslash = "/", mustWork = TRUE)
out_norm <- normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE)

if (identical(out_norm, raw_norm) || startsWith(out_norm, paste0(raw_norm, "/"))) {
  stop("SAFETY STOP: OUT_ROOT is inside RAW_ROOT.")
}

if (identical(out_norm, harm_norm) || startsWith(out_norm, paste0(harm_norm, "/"))) {
  stop("SAFETY STOP: OUT_ROOT is inside HARM_ROOT.")
}


# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

ESMS <- c("MIROC-ES2L", "MRI-ESM2-0")

PAST_YEARS <- 1840:1849
HIST_YEARS <- 1850:1859

# Numerical tolerance for the cell-by-cell identity audit.
# GRD/FLT serialization can introduce very small floating-point differences.
ABS_TOL <- 1e-6

DERIVED_VARIABLES <- c(
  "EddyKineticEnergy",
  "TemperatureGradientSurface",
  "TemperatureGradientDepth",
  "SalinityGradientSurface",
  "SalinityGradientDepth",
  "TempDifferenceSurfaceDepth",
  "WindStressMagnitude",
  "WindStressDirection",
  "DistanceToIceEdge"
)

CIRCULAR_VARIABLES <- c("WindStressDirection")

BASE_VARIABLES_USED <- c(
  "UCurrentSurface",
  "VCurrentSurface",
  "TemperatureSurface",
  "TemperatureDepth",
  "SalinitySurface",
  "SalinityDepth",
  "UWindStress",
  "VWindStress",
  "SeaIceConcentration"
)

VARIABLE_LABELS <- c(
  EddyKineticEnergy = "Eddy kinetic energy",
  TemperatureGradientSurface = "Surface temperature gradient",
  TemperatureGradientDepth = "Depth temperature gradient",
  SalinityGradientSurface = "Surface salinity gradient",
  SalinityGradientDepth = "Depth salinity gradient",
  TempDifferenceSurfaceDepth = "Surface-depth temperature difference",
  WindStressMagnitude = "Wind-stress magnitude",
  WindStressDirection = "Wind-stress direction",
  DistanceToIceEdge = "Distance to ice edge",
  UCurrentSurface = "Zonal surface current",
  VCurrentSurface = "Meridional surface current",
  TemperatureSurface = "Surface temperature",
  TemperatureDepth = "Depth temperature",
  SalinitySurface = "Surface salinity",
  SalinityDepth = "Depth salinity",
  UWindStress = "Zonal wind stress",
  VWindStress = "Meridional wind stress",
  SeaIceConcentration = "Sea-ice concentration"
)


# ==============================================================================
# 3. DATE / FILE HELPERS
# ==============================================================================

parse_stack_date <- function(path) {

  x <- basename(path)

  hit <- regmatches(
    x,
    regexpr("[0-9]{8}(?=\\.grd$)", x, perl = TRUE)
  )

  if (length(hit) == 0 || identical(hit, "")) {

    hits <- regmatches(
      x,
      gregexpr("[0-9]{8}", x, perl = TRUE)
    )[[1]]

    if (length(hits) == 0) {
      return(as.Date(NA))
    }

    hit <- tail(hits, 1)
  }

  as.Date(hit, format = "%Y%m%d")
}


list_stack_files <- function(root, esm, experiment, years) {

  folder <- file.path(root, esm, experiment)

  if (!dir.exists(folder)) {
    stop("Directory does not exist: ", folder)
  }

  files <- list.files(
    folder,
    pattern = "\\.grd$",
    full.names = TRUE,
    recursive = TRUE
  )

  dates <- as.Date(vapply(files, parse_stack_date, as.Date(NA)))

  keep <- !is.na(dates) &
    as.integer(format(dates, "%Y")) %in% years

  files[keep]
}


build_file_inventory <- function(root, esm, experiment, years, source_name) {

  files <- list_stack_files(
    root = root,
    esm = esm,
    experiment = experiment,
    years = years
  )

  tibble(
    ESM = esm,
    Experiment = experiment,
    Source = source_name,
    File = files,
    Date = as.Date(vapply(files, parse_stack_date, as.Date(NA))),
    Year = as.integer(format(Date, "%Y")),
    Month = as.integer(format(Date, "%m"))
  ) %>%
    arrange(Date)
}


# ==============================================================================
# 4. ORIGINAL DERIVED-VARIABLE FORMULAS
# ==============================================================================

# These reproduce the formulas used in the original SRW stack-building script.
#
# EddyKineticEnergy:
#   (UCurrentSurface^2 + VCurrentSurface^2) / 2
#
# Temperature gradients:
#   terrain(..., v = "slope", unit = "degrees")
#
# Salinity gradients:
#   terrain(..., v = "slope", unit = "degrees")
#
# TempDifferenceSurfaceDepth:
#   TemperatureSurface - TemperatureDepth
#
# WindStressMagnitude:
#   sqrt(UWindStress^2 + VWindStress^2)
#
# WindStressDirection:
#   atan2(VWindStress, UWindStress) * 180/pi, converted to 0..360
#
# DistanceToIceEdge:
#   threshold SeaIceConcentration >= 15, then terra::distance(..., unit = "km")
#   and mask to the ocean grid.


prepare_sea_ice_edge_audit <- function(sea_ice_layer, ocean_reference) {

  sea_ice_binary <- ifel(
    sea_ice_layer >= 15,
    1,
    NA
  )

  all_na <- terra::global(
    !is.na(sea_ice_binary),
    fun = "sum",
    na.rm = TRUE
  )[1, 1] == 0

  if (isTRUE(all_na)) {

    distance_to_edge <- rast(sea_ice_layer)
    values(distance_to_edge) <- NA

  } else {

    distance_to_edge <- terra::distance(
      sea_ice_binary,
      unit = "km"
    )
  }

  names(distance_to_edge) <- "DistanceToIceEdge"

  # Original workflow masked DistanceToIceEdge with sample_layer.
  # Here we reproduce that ocean mask from a base ocean layer already present
  # in the same final stack.
  ocean_mask <- ocean_reference / ocean_reference

  distance_to_edge <- mask(
    distance_to_edge,
    ocean_mask
  )

  distance_to_edge
}


recalculate_all_derived <- function(st) {

  needed <- c(
    "UCurrentSurface",
    "VCurrentSurface",
    "TemperatureSurface",
    "TemperatureDepth",
    "SalinitySurface",
    "SalinityDepth",
    "UWindStress",
    "VWindStress",
    "SeaIceConcentration"
  )

  missing_needed <- setdiff(
    needed,
    names(st)
  )

  if (length(missing_needed) > 0) {
    stop(
      "Missing base variables: ",
      paste(missing_needed, collapse = ", ")
    )
  }

  EddyKineticEnergy <- (st[["UCurrentSurface"]]^2 + st[["VCurrentSurface"]]^2) / 2
  names(EddyKineticEnergy) <- "EddyKineticEnergy"

  TemperatureGradientSurface <- terrain(
    st[["TemperatureSurface"]],
    v = "slope",
    unit = "degrees"
  )
  names(TemperatureGradientSurface) <- "TemperatureGradientSurface"

  TemperatureGradientDepth <- terrain(
    st[["TemperatureDepth"]],
    v = "slope",
    unit = "degrees"
  )
  names(TemperatureGradientDepth) <- "TemperatureGradientDepth"

  SalinityGradientSurface <- terrain(
    st[["SalinitySurface"]],
    v = "slope",
    unit = "degrees"
  )
  names(SalinityGradientSurface) <- "SalinityGradientSurface"

  SalinityGradientDepth <- terrain(
    st[["SalinityDepth"]],
    v = "slope",
    unit = "degrees"
  )
  names(SalinityGradientDepth) <- "SalinityGradientDepth"

  TempDifferenceSurfaceDepth <- st[["TemperatureSurface"]] - st[["TemperatureDepth"]]
  names(TempDifferenceSurfaceDepth) <- "TempDifferenceSurfaceDepth"

  WindStressMagnitude <- sqrt(st[["UWindStress"]]^2 + st[["VWindStress"]]^2)
  names(WindStressMagnitude) <- "WindStressMagnitude"

  WindStressDirection <- atan2(st[["VWindStress"]], st[["UWindStress"]]) * (180 / pi)
  WindStressDirection <- ifel(
    WindStressDirection < 0,
    WindStressDirection + 360,
    WindStressDirection
  )
  names(WindStressDirection) <- "WindStressDirection"

  DistanceToIceEdge <- prepare_sea_ice_edge_audit(
    sea_ice_layer = st[["SeaIceConcentration"]],
    ocean_reference = st[["SeaSurfaceSalinity"]]
  )

  rast(
    list(
      EddyKineticEnergy,
      TemperatureGradientSurface,
      TemperatureGradientDepth,
      SalinityGradientSurface,
      SalinityGradientDepth,
      TempDifferenceSurfaceDepth,
      WindStressMagnitude,
      WindStressDirection,
      DistanceToIceEdge
    )
  )
}


# ==============================================================================
# 5. RASTER DIFFERENCE HELPERS
# ==============================================================================

angular_difference_raster <- function(a, b) {

  # Signed shortest angular difference a - b in degrees.
  ((a - b + 180) %% 360) - 180
}


difference_stats <- function(stored, recalculated, variable) {

  if (variable %in% CIRCULAR_VARIABLES) {

    diff_r <- angular_difference_raster(
      stored,
      recalculated
    )

  } else {

    diff_r <- stored - recalculated
  }

  abs_diff <- abs(diff_r)

  n_valid <- terra::global(
    !is.na(diff_r),
    fun = "sum",
    na.rm = TRUE
  )[1, 1]

  mean_abs <- terra::global(
    abs_diff,
    fun = "mean",
    na.rm = TRUE
  )[1, 1]

  max_abs <- terra::global(
    abs_diff,
    fun = "max",
    na.rm = TRUE
  )[1, 1]

  rmse <- sqrt(
    terra::global(
      diff_r^2,
      fun = "mean",
      na.rm = TRUE
    )[1, 1]
  )

  mean_signed <- terra::global(
    diff_r,
    fun = "mean",
    na.rm = TRUE
  )[1, 1]

  tibble(
    NValidCells = n_valid,
    MeanSignedResidual = mean_signed,
    MAE = mean_abs,
    RMSE = rmse,
    MaxAbsResidual = max_abs
  )
}


# ==============================================================================
# 6. CELL-BY-CELL AUDIT OF HARMONISED DERIVATIVES
# ==============================================================================

message("\n============================================================")
message("1. CELL-BY-CELL DERIVED-VARIABLE AUDIT")
message("============================================================")

identity_results <- list()
identity_i <- 1

for (esm in ESMS) {

  files <- list_stack_files(
    root = HARM_ROOT,
    esm = esm,
    experiment = "past1000",
    years = PAST_YEARS
  )

  message("\n", esm, ": ", length(files), " harmonised past1000 files")

  for (i in seq_along(files)) {

    f <- files[i]
    this_date <- parse_stack_date(f)

    if (i == 1 || i %% 20 == 0 || i == length(files)) {
      message(
        "  ",
        i,
        "/",
        length(files),
        " | ",
        format(this_date),
        " | ",
        basename(f)
      )
    }

    st <- rast(f)

    missing_derived <- setdiff(
      DERIVED_VARIABLES,
      names(st)
    )

    if (length(missing_derived) > 0) {

      for (v in missing_derived) {

        identity_results[[identity_i]] <- tibble(
          ESM = esm,
          Date = this_date,
          Variable = v,
          NValidCells = NA_real_,
          MeanSignedResidual = NA_real_,
          MAE = NA_real_,
          RMSE = NA_real_,
          MaxAbsResidual = NA_real_,
          Pass = FALSE,
          Note = "Stored derived layer missing"
        )

        identity_i <- identity_i + 1
      }
    }

    present_derived <- intersect(
      DERIVED_VARIABLES,
      names(st)
    )

    if (length(present_derived) == 0) {
      next
    }

    recalc <- recalculate_all_derived(st)

    for (v in present_derived) {

      stats_v <- difference_stats(
        stored = st[[v]],
        recalculated = recalc[[v]],
        variable = v
      )

      max_res <- stats_v$MaxAbsResidual[1]

      pass_v <- is.finite(max_res) &&
        max_res <= ABS_TOL

      identity_results[[identity_i]] <- bind_cols(
        tibble(
          ESM = esm,
          Date = this_date,
          Variable = v
        ),
        stats_v,
        tibble(
          Pass = pass_v,
          Note = ifelse(
            pass_v,
            "stored == recalculated within tolerance",
            "residual exceeds tolerance"
          )
        )
      )

      identity_i <- identity_i + 1
    }
  }
}

identity_audit <- bind_rows(identity_results) %>%
  arrange(ESM, Variable, Date)

write_csv(
  identity_audit,
  file.path(
    OUT_ROOT,
    "01_derived_cellwise_identity_audit.csv"
  )
)


# ==============================================================================
# 7. SUMMARISE IDENTITY AUDIT
# ==============================================================================

identity_summary <- identity_audit %>%
  group_by(ESM, Variable) %>%
  summarise(
    NFiles = n(),
    NPassed = sum(Pass, na.rm = TRUE),
    NFailed = sum(!Pass, na.rm = TRUE),
    WorstMAE = max(MAE, na.rm = TRUE),
    WorstRMSE = max(RMSE, na.rm = TRUE),
    WorstMaxAbsResidual = max(MaxAbsResidual, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    AllPassed = NFailed == 0
  ) %>%
  arrange(ESM, Variable)

write_csv(
  identity_summary,
  file.path(
    OUT_ROOT,
    "02_derived_identity_summary.csv"
  )
)

message("\nIDENTITY AUDIT SUMMARY:")
print(identity_summary, n = Inf)


# ==============================================================================
# 8. SPATIAL-MEAN EXTRACTION FOR BASE + DERIVED VARIABLES
# ==============================================================================

AUDIT_VARIABLES <- c(
  BASE_VARIABLES_USED,
  DERIVED_VARIABLES
)


circular_mean_deg <- function(x) {

  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  r <- x * pi / 180

  a <- atan2(
    mean(sin(r)),
    mean(cos(r))
  ) * 180 / pi

  if (a < 0) {
    a <- a + 360
  }

  a
}


angular_difference_deg <- function(a, b) {
  ((a - b + 180) %% 360) - 180
}


extract_stack_spatial_means <- function(path, esm, source_name) {

  st <- rast(path)
  this_date <- parse_stack_date(path)

  present <- intersect(
    AUDIT_VARIABLES,
    names(st)
  )

  if (length(present) == 0) {
    return(tibble())
  }

  linear_present <- setdiff(
    present,
    CIRCULAR_VARIABLES
  )

  out <- list()

  if (length(linear_present) > 0) {

    g <- terra::global(
      st[[linear_present]],
      fun = "mean",
      na.rm = TRUE
    )

    out[[length(out) + 1]] <- tibble(
      ESM = esm,
      Source = source_name,
      Date = this_date,
      Year = as.integer(format(this_date, "%Y")),
      Month = as.integer(format(this_date, "%m")),
      Variable = linear_present,
      Value = as.numeric(g[, 1])
    )
  }

  if ("WindStressDirection" %in% present) {

    wd <- st[["WindStressDirection"]]

    mean_sin <- terra::global(
      sin(wd * pi / 180),
      fun = "mean",
      na.rm = TRUE
    )[1, 1]

    mean_cos <- terra::global(
      cos(wd * pi / 180),
      fun = "mean",
      na.rm = TRUE
    )[1, 1]

    angle <- atan2(
      mean_sin,
      mean_cos
    ) * 180 / pi

    if (is.finite(angle) && angle < 0) {
      angle <- angle + 360
    }

    out[[length(out) + 1]] <- tibble(
      ESM = esm,
      Source = source_name,
      Date = this_date,
      Year = as.integer(format(this_date, "%Y")),
      Month = as.integer(format(this_date, "%m")),
      Variable = "WindStressDirection",
      Value = angle
    )
  }

  bind_rows(out)
}


extract_period <- function(root, esm, experiment, years, source_name) {

  files <- list_stack_files(
    root = root,
    esm = esm,
    experiment = experiment,
    years = years
  )

  message(
    "\n",
    esm,
    " | ",
    source_name,
    ": ",
    length(files),
    " files"
  )

  out <- vector(
    "list",
    length(files)
  )

  for (i in seq_along(files)) {

    if (i == 1 || i %% 20 == 0 || i == length(files)) {
      message(
        "  ",
        i,
        "/",
        length(files)
      )
    }

    out[[i]] <- extract_stack_spatial_means(
      path = files[i],
      esm = esm,
      source_name = source_name
    )
  }

  bind_rows(out)
}


message("\n============================================================")
message("2. EXTRACTING RAW / HARMONISED / HISTORICAL SERIES")
message("============================================================")

series_list <- list()
series_i <- 1

for (esm in ESMS) {

  series_list[[series_i]] <- extract_period(
    RAW_ROOT,
    esm,
    "past1000",
    PAST_YEARS,
    "raw_past1000"
  )
  series_i <- series_i + 1

  series_list[[series_i]] <- extract_period(
    HARM_ROOT,
    esm,
    "past1000",
    PAST_YEARS,
    "harmonised_past1000"
  )
  series_i <- series_i + 1

  series_list[[series_i]] <- extract_period(
    RAW_ROOT,
    esm,
    "historical",
    HIST_YEARS,
    "raw_historical"
  )
  series_i <- series_i + 1
}

spatial_means <- bind_rows(series_list) %>%
  arrange(ESM, Variable, Date, Source)

write_csv(
  spatial_means,
  file.path(
    OUT_ROOT,
    "03_base_and_derived_monthly_spatial_means.csv"
  )
)


# ==============================================================================
# 9. MONTHLY CLIMATOLOGIES
# ==============================================================================

monthly_climatology <- spatial_means %>%
  group_by(
    ESM,
    Source,
    Variable,
    Month
  ) %>%
  group_modify(
    ~ tibble(
      ClimMean = if (.y$Variable[[1]] %in% CIRCULAR_VARIABLES) {
        circular_mean_deg(.x$Value)
      } else {
        mean(.x$Value, na.rm = TRUE)
      },
      NYears = sum(is.finite(.x$Value))
    )
  ) %>%
  ungroup()

write_csv(
  monthly_climatology,
  file.path(
    OUT_ROOT,
    "04_base_and_derived_monthly_climatologies.csv"
  )
)


# ==============================================================================
# 10. DISCONTINUITY BEFORE / AFTER HARMONISATION
# ==============================================================================

clim_wide <- monthly_climatology %>%
  select(
    ESM,
    Variable,
    Month,
    Source,
    ClimMean
  ) %>%
  pivot_wider(
    names_from = Source,
    values_from = ClimMean
  )


jump_linear <- clim_wide %>%
  filter(!Variable %in% CIRCULAR_VARIABLES) %>%
  mutate(
    RawSignedJump = raw_historical - raw_past1000,
    HarmonisedSignedJump = raw_historical - harmonised_past1000,
    RawAbsJump = abs(RawSignedJump),
    HarmonisedAbsJump = abs(HarmonisedSignedJump)
  )


jump_circular <- clim_wide %>%
  filter(Variable %in% CIRCULAR_VARIABLES) %>%
  mutate(
    RawSignedJump = angular_difference_deg(
      raw_historical,
      raw_past1000
    ),
    HarmonisedSignedJump = angular_difference_deg(
      raw_historical,
      harmonised_past1000
    ),
    RawAbsJump = abs(RawSignedJump),
    HarmonisedAbsJump = abs(HarmonisedSignedJump)
  )


monthly_jump <- bind_rows(
  jump_linear,
  jump_circular
) %>%
  arrange(
    ESM,
    Variable,
    Month
  )

write_csv(
  monthly_jump,
  file.path(
    OUT_ROOT,
    "05_monthly_discontinuity_before_after.csv"
  )
)


jump_summary <- monthly_jump %>%
  group_by(
    ESM,
    Variable
  ) %>%
  summarise(
    NMonths = sum(
      is.finite(RawAbsJump) &
        is.finite(HarmonisedAbsJump)
    ),
    RawMeanAbsJump = mean(
      RawAbsJump,
      na.rm = TRUE
    ),
    HarmonisedMeanAbsJump = mean(
      HarmonisedAbsJump,
      na.rm = TRUE
    ),
    RawMaxAbsJump = max(
      RawAbsJump,
      na.rm = TRUE
    ),
    HarmonisedMaxAbsJump = max(
      HarmonisedAbsJump,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    AbsoluteImprovement = RawMeanAbsJump - HarmonisedMeanAbsJump,
    ReductionPct = ifelse(
      is.finite(RawMeanAbsJump) &
        RawMeanAbsJump > .Machine$double.eps,
      100 * (
        1 -
          HarmonisedMeanAbsJump /
          RawMeanAbsJump
      ),
      NA_real_
    ),
    Status = case_when(
      !is.finite(RawMeanAbsJump) |
        !is.finite(HarmonisedMeanAbsJump) ~ "insufficient_data",
      HarmonisedMeanAbsJump < RawMeanAbsJump ~ "improved",
      HarmonisedMeanAbsJump > RawMeanAbsJump ~ "worsened",
      TRUE ~ "no_change"
    ),
    VariableType = ifelse(
      Variable %in% DERIVED_VARIABLES,
      "derived",
      "base"
    ),
    Label = unname(VARIABLE_LABELS[Variable])
  ) %>%
  arrange(
    ESM,
    VariableType,
    Variable
  )

write_csv(
  jump_summary,
  file.path(
    OUT_ROOT,
    "06_discontinuity_summary_all_base_and_derived.csv"
  )
)


derived_jump_summary <- jump_summary %>%
  filter(
    VariableType == "derived"
  )

write_csv(
  derived_jump_summary,
  file.path(
    OUT_ROOT,
    "07_discontinuity_summary_derived_only.csv"
  )
)


# ==============================================================================
# 11. MAP EACH DERIVED VARIABLE TO ITS BASE VARIABLES
# ==============================================================================

dependency_table <- tibble(
  Derived = c(
    "EddyKineticEnergy",
    "TemperatureGradientSurface",
    "TemperatureGradientDepth",
    "SalinityGradientSurface",
    "SalinityGradientDepth",
    "TempDifferenceSurfaceDepth",
    "WindStressMagnitude",
    "WindStressDirection",
    "DistanceToIceEdge"
  ),
  Base1 = c(
    "UCurrentSurface",
    "TemperatureSurface",
    "TemperatureDepth",
    "SalinitySurface",
    "SalinityDepth",
    "TemperatureSurface",
    "UWindStress",
    "UWindStress",
    "SeaIceConcentration"
  ),
  Base2 = c(
    "VCurrentSurface",
    NA,
    NA,
    NA,
    NA,
    "TemperatureDepth",
    "VWindStress",
    "VWindStress",
    NA
  )
)

dependency_diagnostics <- derived_jump_summary %>%
  select(
    ESM,
    Derived = Variable,
    DerivedRawJump = RawMeanAbsJump,
    DerivedHarmonisedJump = HarmonisedMeanAbsJump,
    DerivedReductionPct = ReductionPct,
    DerivedStatus = Status
  ) %>%
  left_join(
    dependency_table,
    by = "Derived"
  )

base_jump_lookup <- jump_summary %>%
  filter(
    VariableType == "base"
  ) %>%
  select(
    ESM,
    Variable,
    RawMeanAbsJump,
    HarmonisedMeanAbsJump,
    ReductionPct,
    Status
  )

dependency_diagnostics <- dependency_diagnostics %>%
  left_join(
    base_jump_lookup %>%
      rename(
        Base1 = Variable,
        Base1RawJump = RawMeanAbsJump,
        Base1HarmonisedJump = HarmonisedMeanAbsJump,
        Base1ReductionPct = ReductionPct,
        Base1Status = Status
      ),
    by = c("ESM", "Base1")
  ) %>%
  left_join(
    base_jump_lookup %>%
      rename(
        Base2 = Variable,
        Base2RawJump = RawMeanAbsJump,
        Base2HarmonisedJump = HarmonisedMeanAbsJump,
        Base2ReductionPct = ReductionPct,
        Base2Status = Status
      ),
    by = c("ESM", "Base2")
  )

write_csv(
  dependency_diagnostics,
  file.path(
    OUT_ROOT,
    "08_derived_vs_base_dependency_diagnostics.csv"
  )
)


# ==============================================================================
# 12. FIGURE — EXACT RECALCULATION AUDIT
# ==============================================================================

p_identity <- identity_summary %>%
  mutate(
    Label = unname(VARIABLE_LABELS[Variable])
  ) %>%
  ggplot(
    aes(
      x = WorstMaxAbsResidual,
      y = reorder(Label, WorstMaxAbsResidual)
    )
  ) +
  geom_vline(
    xintercept = ABS_TOL,
    linetype = 2,
    linewidth = 0.45
  ) +
  geom_point(
    size = 2.6
  ) +
  facet_wrap(
    ~ ESM,
    scales = "free_y"
  ) +
  scale_x_log10() +
  labs(
    title = "Cell-by-cell audit of recalculated derived variables",
    subtitle = paste0(
      "Maximum absolute residual between stored and independently recalculated layers; ",
      "vertical line = tolerance (",
      format(ABS_TOL, scientific = TRUE),
      ")"
    ),
    x = "Worst maximum absolute residual across 1840–1849 (log scale)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    OUT_ROOT,
    "Fig01_derived_formula_identity_audit.png"
  ),
  p_identity,
  width = 11,
  height = 7,
  dpi = 300
)


# ==============================================================================
# 13. FIGURE — ALL DERIVED VARIABLES, BEFORE VS AFTER
# ==============================================================================

derived_plot_df <- derived_jump_summary %>%
  select(
    ESM,
    Variable,
    Label,
    RawMeanAbsJump,
    HarmonisedMeanAbsJump
  ) %>%
  pivot_longer(
    cols = c(
      RawMeanAbsJump,
      HarmonisedMeanAbsJump
    ),
    names_to = "State",
    values_to = "MeanAbsJump"
  ) %>%
  mutate(
    State = recode(
      State,
      RawMeanAbsJump = "Before harmonisation",
      HarmonisedMeanAbsJump = "After harmonisation"
    )
  )

p_derived <- ggplot(
  derived_plot_df,
  aes(
    x = MeanAbsJump,
    y = Label,
    shape = State
  )
) +
  geom_point(
    size = 2.6,
    position = position_dodge(width = 0.55)
  ) +
  facet_wrap(
    ~ ESM,
    scales = "free_x"
  ) +
  labs(
    title = "Derived-variable discontinuity before and after harmonisation",
    subtitle = "Mean absolute monthly climatological difference: 1840–1849 vs 1850–1859",
    x = "Mean absolute discontinuity",
    y = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    OUT_ROOT,
    "Fig02_all_derived_discontinuity_before_after.png"
  ),
  p_derived,
  width = 11,
  height = 8,
  dpi = 300
)


# ==============================================================================
# 14. FIGURE — PERCENT CHANGE FOR ALL DERIVED VARIABLES
# ==============================================================================

p_reduction <- derived_jump_summary %>%
  filter(
    is.finite(ReductionPct)
  ) %>%
  ggplot(
    aes(
      x = ReductionPct,
      y = reorder(Label, ReductionPct)
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45
  ) +
  geom_point(
    size = 2.7
  ) +
  facet_wrap(
    ~ ESM,
    scales = "free_y"
  ) +
  labs(
    title = "Effect of harmonisation on every derived variable",
    subtitle = "Positive = reduced discontinuity; negative = increased discontinuity",
    x = "Change in mean absolute discontinuity (%)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    OUT_ROOT,
    "Fig03_all_derived_discontinuity_change_percent.png"
  ),
  p_reduction,
  width = 10,
  height = 8,
  dpi = 300
)


# ==============================================================================
# 15. FIGURE — DERIVED VARIABLE VS BASE VARIABLES
# ==============================================================================

dependency_long <- dependency_diagnostics %>%
  select(
    ESM,
    Derived,
    DerivedReductionPct,
    Base1,
    Base1ReductionPct,
    Base2,
    Base2ReductionPct
  ) %>%
  pivot_longer(
    cols = c(
      DerivedReductionPct,
      Base1ReductionPct,
      Base2ReductionPct
    ),
    names_to = "Role",
    values_to = "ReductionPct"
  ) %>%
  mutate(
    Variable = case_when(
      Role == "DerivedReductionPct" ~ Derived,
      Role == "Base1ReductionPct" ~ Base1,
      Role == "Base2ReductionPct" ~ Base2,
      TRUE ~ NA_character_
    ),
    Role = recode(
      Role,
      DerivedReductionPct = "Derived",
      Base1ReductionPct = "Base variable 1",
      Base2ReductionPct = "Base variable 2"
    ),
    DerivedLabel = unname(VARIABLE_LABELS[Derived])
  ) %>%
  filter(
    !is.na(Variable),
    is.finite(ReductionPct)
  )

p_dependency <- ggplot(
  dependency_long,
  aes(
    x = ReductionPct,
    y = Role,
    shape = Role
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  geom_point(
    size = 2.3
  ) +
  facet_grid(
    ESM ~ DerivedLabel,
    scales = "free_x"
  ) +
  labs(
    title = "Derived-variable response versus its harmonised base variables",
    subtitle = "Positive values indicate a smaller 1840s–1850s climatological discontinuity",
    x = "Reduction in mean absolute discontinuity (%)",
    y = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 7)
  )

ggsave(
  file.path(
    OUT_ROOT,
    "Fig04_derived_vs_base_variables.png"
  ),
  p_dependency,
  width = 18,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 16. CONSOLE REPORT
# ==============================================================================

message("\n============================================================")
message("FINAL AUDIT REPORT")
message("============================================================")

message("\nA) FORMULA IDENTITY")
message("Stored harmonised derived layer vs independently recalculated layer:")

print(
  identity_summary %>%
    select(
      ESM,
      Variable,
      NFiles,
      NPassed,
      NFailed,
      WorstMaxAbsResidual,
      AllPassed
    ),
  n = Inf
)

message("\nB) DERIVED-VARIABLE DISCONTINUITY")
print(
  derived_jump_summary %>%
    select(
      ESM,
      Variable,
      RawMeanAbsJump,
      HarmonisedMeanAbsJump,
      ReductionPct,
      Status
    ),
  n = Inf
)

message("\nC) DERIVED VARIABLES AND THEIR BASE COMPONENTS")
print(
  dependency_diagnostics,
  n = Inf
)

message("\n------------------------------------------------------------")

n_failed <- sum(
  !identity_audit$Pass,
  na.rm = TRUE
)

if (n_failed == 0) {

  message(
    "FORMULA CHECK: PASS — every stored derived layer matched the ",
    "independent recalculation within tolerance."
  )

} else {

  message(
    "FORMULA CHECK: ATTENTION — ",
    n_failed,
    " file-variable combinations exceeded the tolerance."
  )
}

n_improved <- sum(
  derived_jump_summary$Status == "improved",
  na.rm = TRUE
)

n_worsened <- sum(
  derived_jump_summary$Status == "worsened",
  na.rm = TRUE
)

message(
  "BOUNDARY CHECK: ",
  n_improved,
  " derived ESM-variable combinations improved; ",
  n_worsened,
  " worsened."
)

message("\nOutputs:")
message(OUT_ROOT)

message("\nMost important CSVs:")
message(
  file.path(
    OUT_ROOT,
    "02_derived_identity_summary.csv"
  )
)
message(
  file.path(
    OUT_ROOT,
    "07_discontinuity_summary_derived_only.csv"
  )
)
message(
  file.path(
    OUT_ROOT,
    "08_derived_vs_base_dependency_diagnostics.csv"
  )
)

message("\n============================================================")
message("AUDIT COMPLETE")
message("============================================================")
