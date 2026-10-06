# ==============================================================================
# 00_query_cmip6_srw.R
# ==============================================================================
# Southern right whale (SRW) project
#
# PURPOSE
#   Query the ESGF / CMIP6 catalogue for the environmental data used in the SRW
#   project. This script ONLY performs the catalogue query and saves the raw
#   query results. It does not select the final ensemble members and it does not
#   download NetCDF files.
#
# OUTPUTS
#   00_query_PMIP_past1000_raw.csv
#   00_query_CMIP_historical_SSP_raw.csv
#
# NEXT SCRIPT
#   01_prepare_cmip6_download_manifest.R
# ============================================================================== 

# ---- 1. Packages --------------------------------------------------------------

library(epwshiftr)
library(dplyr)

# ---- 2. Configuration ---------------------------------------------------------

QUERY_DIR <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/00enviro/queryCMIP6"

dir.create(QUERY_DIR, recursive = TRUE, showWarnings = FALSE)

# Final ESMs retained for the SRW analysis.
ESMS <- c("MIROC-ES2L", "MRI-ESM2-0")

# Ocean variables used to build the environmental stacks.
# mlotst is available for MRI-ESM2-0 but not for MIROC-ES2L; this is expected.
OCEAN_VARIABLES <- c(
  "tos",      # sea surface temperature
  "sos",      # sea surface salinity
  "so",       # 3-D salinity
  "mlotst",   # mixed layer depth
  "zos",      # sea surface height
  "uo",       # eastward ocean current
  "vo",       # northward ocean current
  "tauuo",    # eastward wind stress
  "tauvo",    # northward wind stress
  "wo",       # vertical ocean velocity
  "thetao"    # 3-D potential temperature
)

# Sea-ice concentration was queried separately in the original workflow.
SEA_ICE_VARIABLES <- "siconc"

# Candidate native resolutions used during the original query stage.
RESOLUTIONS <- c("100 km", "1x1 degree", "250 km")

# ---- 3. Query helper ----------------------------------------------------------

run_query <- function(activity, experiment, variables) {
  epwshiftr::init_cmip6_index(
    activity = activity,
    variable = variables,
    source = ESMS,
    frequency = "mon",
    experiment = experiment,
    variant = NULL,
    replica = FALSE,
    latest = TRUE,
    resolution = RESOLUTIONS
  )
}

# ---- 4. PMIP / past1000 -------------------------------------------------------

cat("\n============================================================\n")
cat("QUERY 1/2: PMIP / past1000\n")
cat("============================================================\n")

past_ocean <- run_query(
  activity = "PMIP",
  experiment = "past1000",
  variables = OCEAN_VARIABLES
) %>%
  mutate(query_group = "ocean")

past_ice <- run_query(
  activity = "PMIP",
  experiment = "past1000",
  variables = SEA_ICE_VARIABLES
) %>%
  mutate(query_group = "sea_ice")

query_past1000 <- bind_rows(past_ocean, past_ice) %>%
  distinct(file_id, .keep_all = TRUE)

past_file <- file.path(QUERY_DIR, "00_query_PMIP_past1000_raw.csv")
write.csv(query_past1000, past_file, row.names = FALSE)

cat("Rows returned:", nrow(query_past1000), "\n")
cat("Saved:", past_file, "\n")

# ---- 5. CMIP + ScenarioMIP / historical + SSPs -------------------------------

cat("\n============================================================\n")
cat("QUERY 2/2: CMIP + ScenarioMIP / historical + SSP126 + SSP585\n")
cat("============================================================\n")

future_ocean <- run_query(
  activity = c("CMIP", "ScenarioMIP"),
  experiment = c("historical", "ssp126", "ssp585"),
  variables = OCEAN_VARIABLES
) %>%
  mutate(query_group = "ocean")

future_ice <- run_query(
  activity = c("CMIP", "ScenarioMIP"),
  experiment = c("historical", "ssp126", "ssp585"),
  variables = SEA_ICE_VARIABLES
) %>%
  mutate(query_group = "sea_ice")

query_historical_ssp <- bind_rows(future_ocean, future_ice) %>%
  distinct(file_id, .keep_all = TRUE)

future_file <- file.path(QUERY_DIR, "00_query_CMIP_historical_SSP_raw.csv")
write.csv(query_historical_ssp, future_file, row.names = FALSE)

cat("Rows returned:", nrow(query_historical_ssp), "\n")
cat("Saved:", future_file, "\n")

# ---- 6. Quick summary ---------------------------------------------------------

query_all <- bind_rows(query_past1000, query_historical_ssp)

summary_query <- query_all %>%
  count(source_id, experiment_id, member_id, variable_id, nominal_resolution, name = "n_files") %>%
  arrange(source_id, experiment_id, variable_id, member_id)

cat("\nQuery summary:\n")
print(summary_query, n = Inf)

cat("\n00_query_cmip6_srw.R finished successfully.\n")
cat("Next: run 01_prepare_cmip6_download_manifest.R\n")
