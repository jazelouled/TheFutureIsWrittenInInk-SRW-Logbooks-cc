# ==============================================================================
# 04_validate_1850_harmonisation.R
#
# Southern Right Whale (SRW) CMIP6 environmental stacks
# Validate the past1000 -> historical transition around 1850
#
# Purpose
# -------
# Quantify and visualise whether the harmonisation applied to past1000 reduces
# the discontinuity with historical conditions at the 1850 experiment boundary.
#
# Comparison:
#   RAW past1000        : 1840-1849
#   HARMONISED past1000: 1840-1849
#   RAW historical      : 1850-1859
#
# The script:
#   1. Reads stacks without modifying them.
#   2. Calculates spatial means for the final dynamic SDM predictors.
#   3. Builds monthly climatologies for the three comparison datasets.
#   4. Quantifies the 1850 discontinuity before and after harmonisation.
#   5. Produces annual transition plots and a summary figure of jump reduction.
#
# IMPORTANT
# ---------
# * RAW stacks are READ-ONLY.
# * Harmonised stacks are READ-ONLY.
# * All outputs are written to a NEW validation directory.
# * Static predictors are intentionally excluded from the transition analysis.
# * WindStressDirection is treated as a circular variable.
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

OUT_ROOT <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/04_harmonisation_validation_1850"

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# Absolute safety check: validation outputs must NEVER be written inside RAW_ROOT.
raw_norm <- normalizePath(RAW_ROOT, winslash = "/", mustWork = TRUE)
out_norm <- normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE)

if (startsWith(out_norm, paste0(raw_norm, "/")) || identical(out_norm, raw_norm)) {
  stop("SAFETY STOP: OUT_ROOT is inside RAW_ROOT.")
}


# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

ESMS <- c(
  "MIROC-ES2L",
  "MRI-ESM2-0"
)

PAST_YEARS <- 1840:1849
HIST_YEARS <- 1850:1859

# Final dynamic predictors used in the SDM workflow.
# Static predictors (Bathymetry, SeaFloorSlope, DistanceToShore) are excluded
# because they should not show an experiment-boundary discontinuity.
LINEAR_VARIABLES <- c(
  "SeaSurfaceTemperature",
  "TemperatureGradientSurface",
  "TempDifferenceSurfaceDepth",
  "SeaSurfaceSalinity",
  "SalinityGradientSurface",
  "SeaSurfaceHeight",
  "UCurrentSurface",
  "VCurrentSurface",
  "EddyKineticEnergy",
  "UWindStress",
  "VWindStress",
  "WindStressMagnitude"
)

CIRCULAR_VARIABLES <- c(
  "WindStressDirection"
)

VARIABLES <- c(LINEAR_VARIABLES, CIRCULAR_VARIABLES)

VARIABLE_LABELS <- c(
  SeaSurfaceTemperature = "Sea-surface temperature",
  TemperatureGradientSurface = "Surface temperature gradient",
  TempDifferenceSurfaceDepth = "Surface-depth temperature difference",
  SeaSurfaceSalinity = "Sea-surface salinity",
  SalinityGradientSurface = "Surface salinity gradient",
  SeaSurfaceHeight = "Sea-surface height",
  UCurrentSurface = "Zonal surface current",
  VCurrentSurface = "Meridional surface current",
  EddyKineticEnergy = "Eddy kinetic energy",
  UWindStress = "Zonal wind stress",
  VWindStress = "Meridional wind stress",
  WindStressMagnitude = "Wind-stress magnitude",
  WindStressDirection = "Wind-stress direction"
)


# ==============================================================================
# 3. HELPERS
# ==============================================================================

parse_stack_date <- function(path) {
  x <- basename(path)

  # Prefer an 8-digit YYYYMMDD token immediately before .grd.
  hit <- regmatches(x, regexpr("[0-9]{8}(?=\\.grd$)", x, perl = TRUE))

  if (length(hit) == 0 || identical(hit, "")) {
    # Fallback: take the final 8-digit token anywhere in the filename.
    all_hits <- regmatches(x, gregexpr("[0-9]{8}", x, perl = TRUE))[[1]]

    if (length(all_hits) == 0) {
      return(as.Date(NA))
    }

    hit <- tail(all_hits, 1)
  }

  as.Date(hit, format = "%Y%m%d")
}


list_stack_files <- function(root, esm, experiment, years) {

  folder <- file.path(root, esm, experiment)

  if (!dir.exists(folder)) {
    warning("Directory does not exist: ", folder)
    return(character(0))
  }

  files <- list.files(
    folder,
    pattern = "\\.grd$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(files) == 0) {
    warning("No .grd files found in: ", folder)
    return(character(0))
  }

  dates <- as.Date(vapply(files, parse_stack_date, as.Date(NA)))

  keep <- !is.na(dates) & as.integer(format(dates, "%Y")) %in% years

  files[keep]
}


circular_mean_deg <- function(x) {

  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  radians <- x * pi / 180

  angle <- atan2(
    mean(sin(radians)),
    mean(cos(radians))
  ) * 180 / pi

  if (angle < 0) {
    angle <- angle + 360
  }

  angle
}


angular_difference_deg <- function(a, b) {

  # Signed shortest difference a - b in [-180, 180].
  ((a - b + 180) %% 360) - 180
}


spatial_mean_one_stack <- function(path, esm, source_name) {

  stack_date <- parse_stack_date(path)

  if (is.na(stack_date)) {
    warning("Could not parse date from: ", path)
    return(NULL)
  }

  r <- rast(path)
  available <- names(r)

  present_linear <- intersect(LINEAR_VARIABLES, available)
  present_circular <- intersect(CIRCULAR_VARIABLES, available)
  missing_vars <- setdiff(VARIABLES, available)

  out <- list()

  if (length(present_linear) > 0) {

    vals <- terra::global(
      r[[present_linear]],
      fun = "mean",
      na.rm = TRUE
    )

    out[[length(out) + 1]] <- tibble(
      ESM = esm,
      Source = source_name,
      Date = stack_date,
      Year = as.integer(format(stack_date, "%Y")),
      Month = as.integer(format(stack_date, "%m")),
      Variable = present_linear,
      Value = as.numeric(vals[, 1])
    )
  }

  if (length(present_circular) > 0) {

    for (v in present_circular) {

      angle_radians <- r[[v]] * pi / 180

      mean_sin <- terra::global(
        sin(angle_radians),
        fun = "mean",
        na.rm = TRUE
      )[1, 1]

      mean_cos <- terra::global(
        cos(angle_radians),
        fun = "mean",
        na.rm = TRUE
      )[1, 1]

      angle <- atan2(mean_sin, mean_cos) * 180 / pi

      if (is.finite(angle) && angle < 0) {
        angle <- angle + 360
      }

      out[[length(out) + 1]] <- tibble(
        ESM = esm,
        Source = source_name,
        Date = stack_date,
        Year = as.integer(format(stack_date, "%Y")),
        Month = as.integer(format(stack_date, "%m")),
        Variable = v,
        Value = angle
      )
    }
  }

  if (length(missing_vars) > 0) {
    message(
      "Missing in ",
      esm,
      " | ",
      source_name,
      " | ",
      format(stack_date),
      ": ",
      paste(missing_vars, collapse = ", ")
    )
  }

  if (length(out) == 0) {
    return(NULL)
  }

  bind_rows(out)
}


extract_source <- function(root, esm, experiment, years, source_name) {

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

  if (length(files) == 0) {
    return(tibble())
  }

  result <- vector("list", length(files))

  for (i in seq_along(files)) {

    if (i == 1 || i %% 20 == 0 || i == length(files)) {
      message(
        "  ",
        i,
        "/",
        length(files),
        "  ",
        basename(files[i])
      )
    }

    result[[i]] <- spatial_mean_one_stack(
      path = files[i],
      esm = esm,
      source_name = source_name
    )
  }

  bind_rows(result)
}


summarise_temporal_mean <- function(df) {

  df %>%
    group_by(ESM, Source, Variable) %>%
    group_modify(
      ~ tibble(
        MeanValue = if (.y$Variable[[1]] %in% CIRCULAR_VARIABLES) {
          circular_mean_deg(.x$Value)
        } else {
          mean(.x$Value, na.rm = TRUE)
        }
      )
    ) %>%
    ungroup()
}


# ==============================================================================
# 4. EXTRACT SPATIAL MEANS
# ==============================================================================

message("\n============================================================")
message("EXTRACTING SPATIAL MEANS")
message("============================================================")

all_monthly <- list()
counter <- 1

for (esm in ESMS) {

  # RAW past1000: 1840-1849
  all_monthly[[counter]] <- extract_source(
    root = RAW_ROOT,
    esm = esm,
    experiment = "past1000",
    years = PAST_YEARS,
    source_name = "raw_past1000"
  )
  counter <- counter + 1

  # Harmonised past1000: 1840-1849
  all_monthly[[counter]] <- extract_source(
    root = HARM_ROOT,
    esm = esm,
    experiment = "past1000",
    years = PAST_YEARS,
    source_name = "harmonised_past1000"
  )
  counter <- counter + 1

  # RAW historical: 1850-1859
  all_monthly[[counter]] <- extract_source(
    root = RAW_ROOT,
    esm = esm,
    experiment = "historical",
    years = HIST_YEARS,
    source_name = "raw_historical"
  )
  counter <- counter + 1
}

monthly_spatial_means <- bind_rows(all_monthly) %>%
  arrange(ESM, Variable, Date, Source)

write_csv(
  monthly_spatial_means,
  file.path(OUT_ROOT, "01_monthly_spatial_means.csv")
)


# ==============================================================================
# 5. DATA COMPLETENESS AUDIT
# ==============================================================================

completeness <- monthly_spatial_means %>%
  count(ESM, Source, Variable, name = "NMonths") %>%
  complete(
    ESM = ESMS,
    Source = c(
      "raw_past1000",
      "harmonised_past1000",
      "raw_historical"
    ),
    Variable = VARIABLES,
    fill = list(NMonths = 0)
  ) %>%
  arrange(ESM, Variable, Source)

write_csv(
  completeness,
  file.path(OUT_ROOT, "02_data_completeness.csv")
)

message("\n============================================================")
message("DATA COMPLETENESS")
message("============================================================")

print(
  completeness %>%
    filter(NMonths != 120),
  n = Inf
)


# ==============================================================================
# 6. MONTHLY CLIMATOLOGIES
# ==============================================================================

monthly_climatology <- monthly_spatial_means %>%
  group_by(ESM, Source, Variable, Month) %>%
  group_modify(
    ~ tibble(
      ClimMean = if (.y$Variable[[1]] %in% CIRCULAR_VARIABLES) {
        circular_mean_deg(.x$Value)
      } else {
        mean(.x$Value, na.rm = TRUE)
      },
      N = sum(is.finite(.x$Value))
    )
  ) %>%
  ungroup()

write_csv(
  monthly_climatology,
  file.path(OUT_ROOT, "03_monthly_climatologies.csv")
)


# ==============================================================================
# 7. QUANTIFY THE DISCONTINUITY
# ==============================================================================

clim_wide <- monthly_climatology %>%
  select(ESM, Variable, Month, Source, ClimMean) %>%
  pivot_wider(
    names_from = Source,
    values_from = ClimMean
  )

linear_jump <- clim_wide %>%
  filter(!Variable %in% CIRCULAR_VARIABLES) %>%
  mutate(
    JumpRaw = raw_historical - raw_past1000,
    JumpHarmonised = raw_historical - harmonised_past1000,
    AbsJumpRaw = abs(JumpRaw),
    AbsJumpHarmonised = abs(JumpHarmonised)
  )

circular_jump <- clim_wide %>%
  filter(Variable %in% CIRCULAR_VARIABLES) %>%
  mutate(
    JumpRaw = angular_difference_deg(
      raw_historical,
      raw_past1000
    ),
    JumpHarmonised = angular_difference_deg(
      raw_historical,
      harmonised_past1000
    ),
    AbsJumpRaw = abs(JumpRaw),
    AbsJumpHarmonised = abs(JumpHarmonised)
  )

monthly_jump <- bind_rows(
  linear_jump,
  circular_jump
) %>%
  arrange(ESM, Variable, Month)

write_csv(
  monthly_jump,
  file.path(OUT_ROOT, "04_monthly_transition_differences.csv")
)


transition_metrics <- monthly_jump %>%
  group_by(ESM, Variable) %>%
  summarise(
    NMonths = sum(
      is.finite(AbsJumpRaw) &
        is.finite(AbsJumpHarmonised)
    ),
    RawDiscontinuity = mean(
      AbsJumpRaw,
      na.rm = TRUE
    ),
    HarmonisedDiscontinuity = mean(
      AbsJumpHarmonised,
      na.rm = TRUE
    ),
    RawMaxMonthlyDiscontinuity = max(
      AbsJumpRaw,
      na.rm = TRUE
    ),
    HarmonisedMaxMonthlyDiscontinuity = max(
      AbsJumpHarmonised,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Reduction = RawDiscontinuity - HarmonisedDiscontinuity,
    ReductionPct = ifelse(
      is.finite(RawDiscontinuity) &
        RawDiscontinuity > .Machine$double.eps,
      100 * (
        1 -
          HarmonisedDiscontinuity /
          RawDiscontinuity
      ),
      NA_real_
    ),
    Status = case_when(
      !is.finite(RawDiscontinuity) |
        !is.finite(HarmonisedDiscontinuity) ~ "insufficient_data",
      HarmonisedDiscontinuity < RawDiscontinuity ~ "improved",
      HarmonisedDiscontinuity > RawDiscontinuity ~ "worsened",
      TRUE ~ "no_change"
    ),
    Label = unname(VARIABLE_LABELS[Variable])
  ) %>%
  arrange(
    ESM,
    desc(ReductionPct)
  )

write_csv(
  transition_metrics,
  file.path(OUT_ROOT, "05_transition_metrics.csv")
)


# ==============================================================================
# 8. ANNUAL MEANS FOR TRANSITION FIGURES
# ==============================================================================

annual_means <- monthly_spatial_means %>%
  group_by(ESM, Source, Variable, Year) %>%
  group_modify(
    ~ tibble(
      AnnualMean = if (.y$Variable[[1]] %in% CIRCULAR_VARIABLES) {
        circular_mean_deg(.x$Value)
      } else {
        mean(.x$Value, na.rm = TRUE)
      }
    )
  ) %>%
  ungroup() %>%
  mutate(
    Label = unname(VARIABLE_LABELS[Variable]),
    Source = factor(
      Source,
      levels = c(
        "raw_past1000",
        "harmonised_past1000",
        "raw_historical"
      ),
      labels = c(
        "Raw past1000",
        "Harmonised past1000",
        "Raw historical"
      )
    )
  )

write_csv(
  annual_means,
  file.path(OUT_ROOT, "06_annual_spatial_means.csv")
)


# ==============================================================================
# 9. FIGURES — ANNUAL TRANSITION
# ==============================================================================

for (esm in ESMS) {

  p <- annual_means %>%
    filter(ESM == esm) %>%
    ggplot(
      aes(
        x = Year,
        y = AnnualMean,
        group = Source,
        linetype = Source,
        shape = Source
      )
    ) +
    geom_vline(
      xintercept = 1849.5,
      linewidth = 0.4,
      alpha = 0.7
    ) +
    geom_line(
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      size = 1.2,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~ Label,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      title = paste0(
        esm,
        ": past1000-historical transition"
      ),
      subtitle = paste0(
        "Raw vs harmonised past1000 (1840-1849) ",
        "compared with raw historical (1850-1859)"
      ),
      x = NULL,
      y = "Spatial mean",
      linetype = NULL,
      shape = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(
      OUT_ROOT,
      paste0(
        "Fig01_annual_transition_",
        esm,
        ".png"
      )
    ),
    plot = p,
    width = 13,
    height = 10,
    dpi = 300
  )
}


# ==============================================================================
# 10. FIGURE — MONTHLY CLIMATOLOGICAL JUMP BEFORE VS AFTER
# ==============================================================================

metrics_plot <- transition_metrics %>%
  filter(
    is.finite(RawDiscontinuity),
    is.finite(HarmonisedDiscontinuity)
  ) %>%
  select(
    ESM,
    Variable,
    Label,
    RawDiscontinuity,
    HarmonisedDiscontinuity
  ) %>%
  pivot_longer(
    cols = c(
      RawDiscontinuity,
      HarmonisedDiscontinuity
    ),
    names_to = "Version",
    values_to = "MeanAbsoluteJump"
  ) %>%
  mutate(
    Version = recode(
      Version,
      RawDiscontinuity = "Before harmonisation",
      HarmonisedDiscontinuity = "After harmonisation"
    )
  )

p_jump <- ggplot(
  metrics_plot,
  aes(
    x = MeanAbsoluteJump,
    y = reorder(Label, MeanAbsoluteJump),
    shape = Version
  )
) +
  geom_point(
    size = 2.4,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(
    ~ ESM,
    scales = "free_x"
  ) +
  labs(
    title = "Discontinuity at the past1000-historical boundary",
    subtitle = paste0(
      "Mean absolute difference between monthly climatologies ",
      "(1840-1849 vs 1850-1859)"
    ),
    x = "Mean absolute climatological discontinuity",
    y = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(
    OUT_ROOT,
    "Fig02_discontinuity_before_after.png"
  ),
  plot = p_jump,
  width = 11,
  height = 8,
  dpi = 300
)


# ==============================================================================
# 11. FIGURE — PERCENT REDUCTION
# ==============================================================================

p_reduction <- transition_metrics %>%
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
    linewidth = 0.4
  ) +
  geom_point(
    size = 2.5
  ) +
  facet_wrap(
    ~ ESM,
    scales = "free_y"
  ) +
  labs(
    title = "Reduction in the 1850 climatological discontinuity",
    subtitle = paste0(
      "Positive values indicate that harmonisation moved ",
      "past1000 closer to historical conditions"
    ),
    x = "Reduction in mean absolute discontinuity (%)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(
    OUT_ROOT,
    "Fig03_discontinuity_reduction_percent.png"
  ),
  plot = p_reduction,
  width = 10,
  height = 8,
  dpi = 300
)


# ==============================================================================
# 12. CONSOLE REPORT
# ==============================================================================

message("\n============================================================")
message("1850 HARMONISATION VALIDATION — SUMMARY")
message("============================================================")

for (esm in ESMS) {

  message("\n", esm)
  message(strrep("-", nchar(esm)))

  tab <- transition_metrics %>%
    filter(ESM == esm) %>%
    select(
      Variable,
      RawDiscontinuity,
      HarmonisedDiscontinuity,
      ReductionPct,
      Status
    )

  print(tab, n = Inf)
}

status_summary <- transition_metrics %>%
  count(ESM, Status)

message("\nStatus summary:")
print(status_summary, n = Inf)

message("\nOutputs written to:")
message(OUT_ROOT)

message("\nKey file:")
message(
  file.path(
    OUT_ROOT,
    "05_transition_metrics.csv"
  )
)

message("\nKey figures:")
message(
  file.path(
    OUT_ROOT,
    "Fig02_discontinuity_before_after.png"
  )
)
message(
  file.path(
    OUT_ROOT,
    "Fig03_discontinuity_reduction_percent.png"
  )
)

message("\n============================================================")
message("VALIDATION COMPLETE")
message("============================================================")
