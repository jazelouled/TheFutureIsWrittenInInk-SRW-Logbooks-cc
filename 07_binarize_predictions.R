library(raster)
library(sf)
library(dplyr)
library(stringr)
library(lubridate)

# ====== CONFIG ======
base_bin_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
ESMs         <- c("MIROC-ES2L", "MRI-ES2-0", "MRI-ES2-0", "MRI-ESM2-0") %>% unique()  # per si varies algun nom
ESMs         <- c("MIROC-ES2L", "MRI-ESM2-0")
algorithms   <- c("GBM", "RF", "SVM")   # ajusta si cal
scenarios    <- c("past1000","historical","ssp126","ssp585")

regions_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/ker95_by_region"

# Sortida (si vols)
# out_total_csv   <- "binary_habitat_area_total.csv"
# out_regions_csv <- "binary_habitat_area_by_region.csv"

# ====== LLEGIR REGIONS I CORREGIR NOMS ======
shp_files <- list.files(regions_dir, pattern = "\\.shp$", full.names = TRUE)
stopifnot(length(shp_files) > 0)

regions_list <- lapply(shp_files, function(f) {
  # extreu nom de regió del fitxer, treient prefixos i extensió
  base <- basename(f)
  region_guess <- base %>%
    str_remove("\\.shp$") %>%
    str_remove("^ker95_") %>%
    str_replace_all("_", " ") %>%
    str_trim() %>%
    str_to_title()
  
  sf_obj <- st_read(f, quiet = TRUE)
  sf_obj$Region_raw <- region_guess
  sf_obj
})

regions_sf <- do.call(rbind, regions_list) %>% st_make_valid()

# Correcció: intercanvi South Africa <-> New Zealand
regions_sf$Region <- regions_sf$Region_raw
regions_sf$Region[regions_sf$Region == "South Africa"] <- "TMP__New Zealand"
regions_sf$Region[regions_sf$Region == "New Zealand"]  <- "South Africa"
regions_sf$Region[regions_sf$Region == "TMP__New Zealand"] <- "New Zealand"

# (Opcional) neteja més noms si cal:
regions_sf$Region <- regions_sf$Region %>% str_squish()

# ====== HELPERS ======
# Donar àrea de cel·la (km²). Si lon/lat, usa raster::area(). Si no, projecta a equal-area.
cell_area_km2 <- function(r) {
  if (isLonLat(r)) {
    raster::area(r) # km²
  } else {
    # Projecta a Equal-Area global (EPSG:6933 ~ Cylindrical Equal Area, metres)
    cea <- "+proj=cea +lon_0=0 +lat_ts=30 +datum=WGS84 +units=m +no_defs"
    r_m <- projectRaster(r, crs = cea, method = "ngb")
    # Àrea per cel·la = res_x * res_y (m²) -> km²
    resm <- res(r_m) # metres
    a_km2 <- (resm[1] * resm[2]) / 1e6
    # Tornem a la malla original amb un valor constant per cel·la
    r_out <- setValues(r, a_km2)
    r_out
  }
}

# Parse de metadades (data, ESM, alg, escenari) des del path i nom
parse_meta <- function(file_path) {
  b <- basename(file_path)
  # format: yyyymmdd_ESM_ALG_prediction_binary.tif
  date_str <- str_sub(b, 1, 8)
  date     <- suppressWarnings(ymd(date_str))
  parts <- strsplit(b, "_")[[1]]
  esm  <- parts[2]
  alg  <- parts[3]
  # escenari ve del directori
  path_parts <- strsplit(dirname(file_path), .Platform$file.sep)[[1]]
  scen <- tail(path_parts, 1)
  list(date = date, esm = esm, alg = alg, scenario = scen)
}




library(raster)
library(sf)
library(parallel)
library(stringr)
library(dplyr)

n_cores <- parallel::detectCores() - 3

# --- Funció per processar una combinació completa ---
process_combination <- function(esm, alg, scen, base_bin_dir, regions_sf) {
  folder <- file.path(base_bin_dir, esm, alg, scen)
  if (!dir.exists(folder)) return(NULL)
  
  message("\n=== ", esm, " | ", alg, " | ", scen, " ===")
  tif_files <- list.files(folder, pattern = "_prediction_binary\\.tif$", full.names = TRUE)
  if (!length(tif_files)) {
    message("   ⚠️  Cap fitxer trobat.")
    return(NULL)
  }
  
  total_files <- length(tif_files)
  message("   ", total_files, " fitxers trobats.")
  
  # Inicialitzem resultats
  ts_total <- list()
  ts_by_region <- list()
  
  for (i in seq_along(tif_files)) {
    f <- tif_files[i]
    meta <- parse_meta(f)
    if (is.na(meta$date)) next
    
    pct <- round((i / total_files) * 100, 1)
    message(sprintf("   → [%s | %s | %s] %s (%d/%d, %.1f%%)",
                    esm, alg, scen, format(meta$date, "%Y-%m-%d"), i, total_files, pct))
    
    r <- try(raster(f), silent = TRUE)
    if (inherits(r, "try-error")) {
      message("     ⚠️  Error llegint fitxer: ", basename(f))
      next
    }
    
    # Binari + àrea
    r <- calc(r, fun = function(x) as.numeric(x >= 1))
    a <- cell_area_km2(r)
    hab_km2_r <- r * a
    
    # --- Total global ---
    total_km2 <- cellStats(hab_km2_r, sum, na.rm = TRUE)
    ts_total[[length(ts_total) + 1]] <- data.frame(
      date       = meta$date,
      esm        = meta$esm,
      algorithm  = meta$alg,
      scenario   = meta$scenario,
      habitat_km2 = total_km2,
      stringsAsFactors = FALSE
    )
    
    # --- Per regió ---
    reg_crs <- try(st_transform(regions_sf, crs(r)), silent = TRUE)
    if (inherits(reg_crs, "try-error")) {
      reg_crs <- st_transform(regions_sf, st_crs(r) %||% 4326)
    }
    reg_sp <- as(reg_crs, "Spatial")
    
    reg_names <- unique(reg_crs$Region)
    for (reg_name in reg_names) {
      reg_one <- reg_sp[which(reg_sp$Region == reg_name), ]
      if (is.null(reg_one) || nrow(reg_one) == 0) next
      
      hab_reg <- try(mask(hab_km2_r, reg_one), silent = TRUE)
      if (inherits(hab_reg, "try-error")) next
      
      reg_km2 <- cellStats(hab_reg, sum, na.rm = TRUE)
      ts_by_region[[length(ts_by_region) + 1]] <- data.frame(
        date        = meta$date,
        esm         = meta$esm,
        algorithm   = meta$alg,
        scenario    = meta$scenario,
        region      = reg_name,
        habitat_km2 = reg_km2,
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(list(ts_total = do.call(rbind, ts_total),
              ts_by_region = do.call(rbind, ts_by_region)))
}

# --- Combinacions a processar ---
combinations <- expand.grid(esm = ESMs, alg = algorithms, scen = scenarios, stringsAsFactors = FALSE)

# --- Executar en paral·lel ---
cat("\n🚀 Processant prediccions binàries en paral·lel amb", n_cores, "nuclis...\n")

results_list <- mclapply(seq_len(nrow(combinations)), function(i) {
  row <- combinations[i, ]
  process_combination(row$esm, row$alg, row$scen, base_bin_dir, regions_sf)
}, mc.cores = n_cores)

cat("\n✅ Totes les combinacions processades.\n")

# --- Consolidar resultats ---
ts_total_all     <- do.call(rbind, lapply(results_list, `[[`, "ts_total"))
ts_by_region_all <- do.call(rbind, lapply(results_list, `[[`, "ts_by_region"))




library(ggplot2)
library(dplyr)
library(patchwork)

# --- Prepara dades ---
ts_total_all$date <- as.Date(ts_total_all$date)
ts_by_region_all$date <- as.Date(ts_by_region_all$date)

# --- Gràfic total (global) ---
p_total <- ts_total_all %>%
  group_by(date, esm, algorithm, scenario) %>%
  summarise(habitat_km2 = mean(habitat_km2, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = date, y = habitat_km2, color = scenario, group = interaction(esm, algorithm, scenario))) +
  geom_line(alpha = 0.6) +
  labs(title = "Total habitat area over time", x = "Year", y = "Habitat area (km²)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# --- Gràfic per regió ---
p_region <- ts_by_region_all %>%
  group_by(date, region, scenario) %>%
  summarise(habitat_km2 = mean(habitat_km2, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = date, y = habitat_km2, color = scenario, group = scenario)) +
  geom_line() +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Habitat area by region", x = "Year", y = "Habitat area (km²)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# --- Mostra-ho tot en un grid ---
p_total + p_region + plot_layout(ncol = 1, heights = c(1, 2))





# --- Gràfic total amb LOESS ---
p_total <- ts_total_all %>%
  group_by(date, esm, algorithm, scenario) %>%
  summarise(habitat_km2 = mean(habitat_km2, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = as.numeric(format(date, "%Y")), y = habitat_km2, color = scenario)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 0.8) +
  labs(title = "Total habitat area (LOESS smoothed)", x = "Year", y = "Habitat area (km²)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# --- Gràfic per regió amb LOESS ---
p_region <- ts_by_region_all %>%
  group_by(date, region, scenario) %>%
  summarise(habitat_km2 = mean(habitat_km2, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = as.numeric(format(date, "%Y")), y = habitat_km2, color = scenario)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 0.8) +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Habitat area by region (LOESS smoothed)", x = "Year", y = "Habitat area (km²)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# --- Mostra-ho tot en un grid ---
p_total + p_region + plot_layout(ncol = 1, heights = c(1, 2))




library(dplyr)
library(ggplot2)
library(patchwork)

# --- Calcular anomalies per cada escenari ---
ts_total_all_anomalia <- ts_total_all %>%
  group_by(esm, algorithm, scenario) %>%
  mutate(habitat_anom = habitat_km2 - mean(habitat_km2, na.rm = TRUE)) %>%
  ungroup()

ts_by_region_all_anomalia <- ts_by_region_all %>%
  group_by(esm, algorithm, scenario, region) %>%
  mutate(habitat_anom = habitat_km2 - mean(habitat_km2, na.rm = TRUE)) %>%
  ungroup()


p_total_anom <- ggplot(ts_total_all_anomalia, aes(x = as.numeric(format(date, "%Y")),
                                                  y = habitat_anom,
                                                  color = scenario)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 1) +
  scale_color_manual(values = c(
    "past1000" = "black",
    "historical" = "#1f78b4",
    "ssp126" = "#33a02c",
    "ssp585" = "#e31a1c"
  )) +
  labs(title = "Total habitat area anomaly (LOESS smoothed)",
       x = "Year", y = expression("Habitat anomaly (km"^2*")")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

p_regions_anom <- ggplot(ts_by_region_all_anomalia, aes(x = as.numeric(format(date, "%Y")),
                                                        y = habitat_anom,
                                                        color = scenario)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 0.9) +
  facet_wrap(~ region, scales = "free_y") +
  scale_color_manual(values = c(
    "past1000" = "black",
    "historical" = "#1f78b4",
    "ssp126" = "#33a02c",
    "ssp585" = "#e31a1c"
  )) +
  labs(title = "Habitat area anomaly by region (LOESS smoothed)",
       x = "Year", y = expression("Habitat anomaly (km"^2*")")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")







library(dplyr)
library(ggplot2)

# 1️⃣ Seleccionem la baseline (past1000)
baseline_df <- ts_by_region_all %>%
  filter(scenario == "past1000") %>%
  group_by(region) %>%
  summarise(baseline_mean = mean(habitat_km2, na.rm = TRUE))

# 2️⃣ Calculem l'anomalia respecte la baseline comuna (past1000)
ts_anom_global <- ts_total_all %>%
  mutate(reference = mean(habitat_km2[scenario == "past1000"], na.rm = TRUE),
         habitat_anom = habitat_km2 - reference)

ts_anom_region <- ts_by_region_all %>%
  left_join(baseline_df, by = "region") %>%
  mutate(habitat_anom = habitat_km2 - baseline_mean)

# 3️⃣ Representem (LOESS + CI reals)
cols <- c("past1000" = "black",
          "historical" = "#1f78b4",
          "ssp126" = "#33a02c",
          "ssp585" = "#e31a1c")

p_total <- ggplot(ts_anom_global, aes(x = date, y = habitat_anom, color = scenario)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.2, linewidth = 0.8) +
  scale_color_manual(values = cols) +
  labs(title = "Total habitat area anomaly (relative to past1000)",
       x = "Year", y = expression("Habitat anomaly (km"^2*")")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

p_region <- ggplot(ts_anom_region, aes(x = date, y = habitat_anom, color = scenario)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 0.8) +
  facet_wrap(~region, scales = "free_y") +
  scale_color_manual(values = cols) +
  labs(title = "Habitat area anomaly by region (relative to past1000)",
       x = "Year", y = expression("Habitat anomaly (km"^2*")")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

p_total
p_region

library(dplyr)
library(ggplot2)

ts_total_all <- ts_total_all %>%
  mutate(year = as.numeric(format(as.Date(date), "%Y"))) %>%
  filter(year >= 1600, year <= 2100)

ts_by_region_all <- ts_by_region_all %>%
  mutate(year = as.numeric(format(as.Date(date), "%Y"))) %>%
  filter(year >= 1600, year <= 2100)


cols_scen <- c(
  past1000  = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)


p_total_models <- ggplot(ts_total_all,
                         aes(x = year, y = habitat_km2,
                             color = scenario, group = interaction(esm, algorithm, scenario))) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 0.8, alpha = 0.2) +
  facet_grid(esm ~ algorithm, scales = "free_y") +
  scale_color_manual(values = cols_scen) +
  labs(
    x = "Year",
    y = expression("Habitat area (km"^2*")"),
    title = "Total habitat area by ESM and algorithm (LOESS smoothed)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_total_models


library(raster)
library(dplyr)
library(ggplot2)
library(stringr)
library(zoo)

# --- Config ---
base_dir <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/01envStacks_Raw_0s/MIROC-ES2L"
varname  <- "SeaSurfaceTemperature"

# --- Funció per extreure any i mes dels noms de fitxer ---
extract_date <- function(filename) {
  d <- str_extract(basename(filename), "\\d{8}")
  if (is.na(d)) return(NA)
  as.Date(d, format = "%Y%m%d")
}

# --- Llista de fitxers ---
past_dir <- file.path(base_dir, "past1000")
hist_dir <- file.path(base_dir, "historical")

past_files <- list.files(past_dir, pattern = "\\.grd$", full.names = TRUE)
hist_files <- list.files(hist_dir, pattern = "\\.grd$", full.names = TRUE)

# --- Extreure només els últims anys de past1000 i primers d'historical ---
past_files <- past_files[str_detect(past_files, "18(3[5-9]|4[0-9])")]     # 1835–1849
hist_files <- hist_files[str_detect(hist_files, "18(4[8-9]|5[0-9]|6[0-1])")] # 1848–1861

# --- Funció per llegir i calcular la mitjana espacial ---
get_mean_sst <- function(files, scenario) {
  res <- lapply(files, function(f) {
    s <- try(stack(f), silent = TRUE)
    if (inherits(s, "try-error") || !(varname %in% names(s))) return(NULL)
    mean_val <- cellStats(s[[varname]], mean, na.rm = TRUE)
    date_val <- extract_date(f)
    data.frame(date = date_val, scenario = scenario, mean_sst = mean_val)
  })
  do.call(rbind, res)
}

# --- Calcular mitjanes per escenari ---
df_past <- get_mean_sst(past_files, "past1000")
df_hist <- get_mean_sst(hist_files, "historical")
df_all <- bind_rows(df_past, df_hist) %>% arrange(date)

# --- Suavitzar amb LOESS per veure tendència ---
ggplot(df_all, aes(x = date, y = mean_sst, color = scenario)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("past1000" = "black", "historical" = "#1f78b4")) +
  labs(
    title = "Sea Surface Temperature (MIROC-ES2L)",
    subtitle = "Últims anys de past1000 vs primers d'historical",
    x = "Any",
    y = "Temperatura mitjana (°C)",
    color = "Escenari"
  ) +
  theme_minimal(base_size = 13)



base_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks"

rds_files <- list.files(base_dir, pattern = "Cutoff", recursive = TRUE, full.names = TRUE)

for (f in rds_files) {
  cat("\n=== ", f, " ===\n")
  print(readRDS(f))
}







library(dplyr)
library(ggplot2)

# --- Prepara dades ---
ts_total_all <- ts_total_all %>%
  mutate(year = as.numeric(format(as.Date(date), "%Y"))) %>%
  filter(year >= 1600, year <= 2100)

cols_scen <- c(
  past1000  = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)

# --- Gràfic ---
p_total_models <- ggplot(
  ts_total_all,
  aes(x = year, y = habitat_km2,
      color = scenario,
      group = interaction(esm, algorithm, scenario))
) +
  # Punts per any
  geom_point(alpha = 0.4, size = 1) +
  # Línia LOESS suau amb IC
  geom_smooth(method = "loess", se = TRUE, span = 0.25,
              linewidth = 0.8, alpha = 0.2) +
  facet_grid(esm ~ algorithm, scales = "free_y") +
  scale_color_manual(values = cols_scen) +
  labs(
    x = "Year",
    y = expression("Habitat area (km"^2*")"),
    title = "Total habitat area by ESM and algorithm (LOESS smoothed + yearly points)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

p_total_models




library(ggplot2)
library(dplyr)

# Prepara dades
ts_total_all <- ts_total_all %>%
  mutate(
    year = as.numeric(format(as.Date(date), "%Y")),
    scenario = factor(scenario, levels = c("past1000", "historical", "ssp126", "ssp585"))
  ) %>%
  filter(year >= 1600, year <= 2100)

# Colors pels escenaris
cols_scen <- c(
  past1000  = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)

# --- Plot: boxplots comparatius entre MIROC i MRI ---
p_box_compare <- ggplot(ts_total_all,
                        aes(x = scenario, y = habitat_km2,
                            fill = esm)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7, alpha = 0.8, outlier.size = 0.6) +
  facet_wrap(~ algorithm, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("MIROC-ES2L" = "#fdae61",
                               "MRI-ESM2-0" = "#2c7bb6")) +
  labs(
    x = "Scenario",
    y = expression("Habitat area (km"^2*")"),
    title = "Comparison of total habitat area by ESM and algorithm",
    fill = "ESM"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_box_compare







library(raster)
library(dplyr)
library(ggplot2)
library(stringr)
library(parallel)

# === CONFIGURACIÓ ===
base_dir <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/01envStacks_Raw_0s"
ESMs <- c("MIROC-ES2L", "MRI-ESM2-0")
scenarios <- c("past1000", "historical", "ssp126", "ssp585")

# Nombre de nuclis per utilitzar
n_cores <- max(1, detectCores() - 1)

cat("\n===== INICI DE LECTURA DE SST EN PARA·LEL =====\n")
cat("💻 Utilitzant", n_cores, "nuclis\n\n")

# --- Helper: processa un fitxer concret ---
process_file <- function(f, esm, scen) {
  date_str <- str_extract(basename(f), "\\d{8}")
  date <- as.Date(date_str, format = "%Y%m%d")
  
  cat(sprintf("     📂 %s | %s | Fitxer: %s\n", esm, scen, basename(f)))
  
  r <- try(stack(f), silent = TRUE)
  if (inherits(r, "try-error")) {
    cat("     ⚠️  Error llegint fitxer!\n")
    return(NULL)
  }
  
  if (!"SeaSurfaceTemperature" %in% names(r)) {
    cat("     ⚠️  Cap 'SeaSurfaceTemperature' al fitxer.\n")
    return(NULL)
  }
  
  sst <- r[["SeaSurfaceTemperature"]]
  
  mean_sst <- cellStats(sst, mean, na.rm = TRUE)
  min_sst  <- cellStats(sst, min, na.rm = TRUE)
  max_sst  <- cellStats(sst, max, na.rm = TRUE)
  
  cat(sprintf("     🌡️  SST: mean = %.2f°C | min = %.2f°C | max = %.2f°C\n", mean_sst, min_sst, max_sst))
  
  data.frame(
    date = date,
    esm = esm,
    scenario = scen,
    mean_sst = mean_sst,
    min_sst = min_sst,
    max_sst = max_sst
  )
}

# --- Crear totes les combinacions ---
combos <- expand.grid(esm = ESMs, scen = scenarios, stringsAsFactors = FALSE)

# --- Bucle principal per ESM i escenari ---
sst_all <- list()

for (i in seq_len(nrow(combos))) {
  esm  <- combos$esm[i]
  scen <- combos$scen[i]
  
  cat("\n🟦 Processant:", esm, "| Escenari:", scen, "\n")
  
  stack_dir <- file.path(base_dir, esm, scen)
  if (!dir.exists(stack_dir)) {
    cat("  ⚠️  No existeix:", stack_dir, "\n")
    next
  }
  
  files <- list.files(stack_dir, pattern = paste0("monthly_stack_", esm, "_", scen, "_.*\\.grd$"), full.names = TRUE)
  if (length(files) == 0) {
    cat("  ⚠️  Cap fitxer trobat a", stack_dir, "\n")
    next
  }
  
  cat("  🔹 Fitxers trobats:", length(files), "\n")
  
  # --- Processar en paral·lel ---
  results <- mclapply(files, function(f) process_file(f, esm, scen), mc.cores = n_cores)
  
  sst_all[[paste(esm, scen, sep = "_")]] <- bind_rows(results)
}

cat("\n===== FI DE LECTURA EN PARA·LEL =====\n")

# --- Consolidar resultats ---
sst_df <- bind_rows(sst_all) %>%
  mutate(year = as.numeric(format(date, "%Y")))

# --- Plot per comparar evolució temporal ---
cols <- c(past1000 = "black", historical = "#1f78b4", ssp126 = "#33a02c", ssp585 = "#e31a1c")

ggplot(sst_df, aes(x = year, y = mean_sst, color = scenario)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25) +
  facet_wrap(~ esm, ncol = 1, scales = "free_y") +
  scale_color_manual(values = cols) +
  labs(
    x = "Year",
    y = "Global mean Sea Surface Temperature (°C)",
    title = "Temporal evolution of global SST per ESM and scenario"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )







library(ggplot2)

cols <- c(
  past1000 = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)

ggplot(sst_df, aes(x = year, y = mean_sst, color = scenario, linetype = esm)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 1) +
  scale_color_manual(values = cols) +
  labs(
    x = "Year",
    y = expression("Global mean SST (°C)"),
    title = "Comparison of SST temporal evolution — MIROC-ES2L vs MRI-ESM2-0",
    color = "Scenario",
    linetype = "ESM"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )




ggplot(sst_df, aes(x = scenario, y = mean_sst, fill = esm)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.7) +
  scale_fill_manual(values = c("MIROC-ES2L" = "#fdae61", "MRI-ESM2-0" = "#2b83ba")) +
  labs(
    x = "Scenario",
    y = "Global mean SST (°C)",
    title = "Distribution of mean global SST per scenario and ESM",
    fill = "ESM"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



library(ggplot2)
library(dplyr)

# --- Preparem dades ---
ts_total_all <- ts_total_all %>%
  mutate(year = as.numeric(format(as.Date(date), "%Y"))) %>%
  filter(year >= 1600, year <= 2100)

# Colors per escenaris
cols_scen <- c(
  past1000  = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)

# Colors per models
cols_esm <- c(
  "MIROC-ES2L" = "#ffb74d",
  "MRI-ESM2-0" = "#64b5f6"
)

# --- Plot comparatiu ---
ggplot(ts_total_all, aes(x = year, y = habitat_km2,
                         color = scenario,
                         linetype = esm)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 1) +
  scale_color_manual(values = cols_scen) +
  scale_linetype_manual(values = c("MIROC-ES2L" = "solid", "MRI-ESM2-0" = "dashed")) +
  facet_wrap(~ algorithm, ncol = 1, scales = "free_y") +
  labs(
    x = "Year",
    y = expression("Habitat area (km"^2*")"),
    title = "Comparison of total habitat area evolution — MIROC-ES2L vs MRI-ESM2-0",
    color = "Scenario",
    linetype = "ESM"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )





library(dplyr)
library(ggplot2)
library(readr)

# ==== RUTES ====
path_miroc <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw/MIROC-ES2L"
path_mri   <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw/MRI-ESM2-0"

# ==== LLEGIR FITXERS ====
files_miroc <- list.files(path_miroc, pattern = "\\.csv$", full.names = TRUE)
files_mri   <- list.files(path_mri, pattern = "\\.csv$", full.names = TRUE)

# Llegeix i combina
df_miroc <- files_miroc %>%
  lapply(read_csv, show_col_types = FALSE) %>%
  bind_rows() %>%
  filter(PresAbs == 1) %>%
  mutate(ESM = "MIROC-ES2L")

df_mri <- files_mri %>%
  lapply(read_csv, show_col_types = FALSE) %>%
  bind_rows() %>%
  filter(PresAbs == 1) %>%
  mutate(ESM = "MRI-ESM2-0")

# ==== COMBINAR ====
df_all <- bind_rows(df_miroc, df_mri) %>%
  select(ESM, TemperatureSurface) %>%
  filter(!is.na(TemperatureSurface), TemperatureSurface < 40, TemperatureSurface > -2)

# ==== GRAFIC ====
cols <- c("MIROC-ES2L" = "#1f78b4", "MRI-ESM2-0" = "#33a02c")

ggplot(df_all, aes(x = TemperatureSurface, fill = ESM, color = ESM)) +
  geom_density(alpha = 0.35, linewidth = 1.2) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(
    x = "Sea Surface Temperature (°C)",
    y = "Density",
    title = "Distribution of Sea Surface Temperature at presence points",
    subtitle = "Comparison between MIROC-ES2L and MRI-ESM2-0 presences"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )







library(raster)
library(fs)
library(dplyr)
library(stringr)
library(lubridate)
library(parallel)

# 📂 Paths base
base_dir <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/01envStacks_Raw_0s"
ESMs     <- c("MIROC-ES2L", "MRI-ESM2-0")
scenarios <- c("past1000", "historical", "ssp126", "ssp585")

# 🧭 Llegeix punts de presència
pres_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw"
miroc_df <- read.csv(file.path(pres_dir, "MIROC-ES2L/SRW_extractRaw_150km_MIROC-ES2L.csv"))
mri_df   <- read.csv(file.path(pres_dir, "MRI-ESM2-0/SRW_extractRaw_150km_MRI-ESM2-0.csv"))

coords_all <- bind_rows(
  miroc_df %>% select(Lon, Lat),
  mri_df   %>% select(Lon, Lat)
) %>%
  distinct() %>%
  mutate(
    Lon = as.numeric(Lon),
    Lat = as.numeric(Lat)
  ) %>%
  filter(!is.na(Lon), !is.na(Lat))

pts <- SpatialPoints(coords_all, proj4string = CRS("+proj=longlat +datum=WGS84"))
bbox_region <- extent(pts)

# --- Extracció de valors SST
n_cores <- detectCores() - 2
sst_summary <- list()

for (esm in ESMs) {
  cat("\n=== ", esm, " ===\n")
  
  for (scen in scenarios) {
    folder <- file.path(base_dir, esm, scen)
    grd_files <- list.files(folder, pattern = "\\.grd$", full.names = TRUE)
    if (!length(grd_files)) next
    
    # extreure any
    years <- as.numeric(str_extract(basename(grd_files), "(?<=_)\\d{8}(?=\\.grd)") %>% substr(1, 4))
    df <- data.frame(file = grd_files, year = years)
    
    # treure primers 10 anys de l’històric
    if (scen == "historical") {
      min_year <- min(df$year, na.rm = TRUE)
      df <- df %>% filter(year > min_year + 10)
    }
    
    # mostra de 100 fitxers
    df <- df %>% slice_sample(n = min(200, nrow(.)))
    
    total_files <- nrow(df)
    message("  -> ", scen, ": ", total_files, " fitxers.")
    
    vals <- mclapply(seq_len(total_files), function(i) {
      f <- df$file[i]
      r <- try(raster::stack(f)$SeaSurfaceTemperature, silent = TRUE)
      if (inherits(r, "try-error")) return(NULL)
      
      # retalla al bounding box de les dades
      r_crop <- crop(r, bbox_region)
      mean(r_crop[], na.rm = TRUE)
    }, mc.cores = n_cores)
    
    vals <- unlist(vals)
    df$mean_sst <- vals
    
    sst_summary[[paste(esm, scen, sep = "_")]] <- df
  }
}

sst_df <- bind_rows(sst_summary, .id = "source") %>%
  mutate(
    esm = str_split_fixed(source, "_", 2)[,1],
    scenario = str_split_fixed(source, "_", 2)[,2]
  )

# --- Plot comparatiu
library(ggplot2)
cols <- c("past1000" = "black", "historical" = "#1f78b4", "ssp126" = "#33a02c", "ssp585" = "#e31a1c")

ggplot(sst_df, aes(x = year, y = mean_sst, color = scenario)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 1) +
  facet_wrap(~ esm, ncol = 1, scales = "free_y") +
  scale_color_manual(values = cols) +
  labs(
    x = "Year",
    y = "Regional mean SST (°C)",
    title = "Temporal evolution of SST within the region with data"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )






library(ggplot2)
library(dplyr)

# --- Prepara les dades ---
ts_total_all <- ts_total_all %>%
  mutate(
    year = as.numeric(format(as.Date(date), "%Y")),
    scenario = factor(scenario, levels = c("past1000", "historical", "ssp126", "ssp585"))
  ) %>%
  filter(year >= 1600, year <= 2100)

# --- Colors coherents ---
cols_scen <- c(
  past1000  = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)

# --- Gràfic global ---
p_total <- ggplot(ts_total_all, 
                  aes(x = year, y = habitat_km2, color = scenario)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 0.9) +
  facet_grid(esm ~ algorithm, scales = "free_y") +
  scale_color_manual(values = cols_scen) +
  labs(
    x = "Year",
    y = expression("Habitat area (km"^2*")"),
    title = "Temporal evolution of total suitable habitat area",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

print(p_total)




# ==============================================================
# PDP REAL DE TEMPERATURA: USANT MODELS JA GUARDATS
# ==============================================================

library(pdp)
library(caret)
library(dplyr)
library(ggplot2)
library(fs)

# --- Configuració
ESMs <- c("MIROC-ES2L", "MRI-ESM2-0")
algorithms <- c("GBM", "RF", "SVM")
predictor <- "SeaSurfaceTemperature"

# --- Directoris (RAW)
base_model_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/02modelling/00raw"
base_data_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw"

# --- Contenidors
pdp_list <- list()

for (esm in ESMs) {
  for (alg in algorithms) {
    
    message("\n→ ", esm, " | ", alg)
    
    # --- Carregar model
    model_path <- file.path(base_model_dir, esm, "01FitModels-noDepth", alg,
                            paste0(alg, "_model.RDS"))
    if (!file.exists(model_path)) {
      message("   ⚠️ No model found: ", model_path)
      next
    }
    model_fit <- readRDS(model_path)
    
    # --- Carregar dades (RAW)
    data_path <- file.path(base_data_dir, esm,
                           paste0("SRW_extractRaw_150km_", esm, ".csv"))
    if (!file.exists(data_path)) {
      message("   ⚠️ No data found: ", data_path)
      next
    }
    df <- read.csv(data_path)
    
    # --- Comprovació i neteja
    if (!("PresAbs" %in% names(df))) {
      message("   ⚠️ Column 'PresAbs' not found in ", esm)
      next
    }
    
    df <- df %>% dplyr::select(-PresAbs) %>% na.omit()
    
    if (!(predictor %in% names(df))) {
      message("   ⚠️ Predictor not found in data: ", predictor)
      next
    }
    
    # --- Calcular PDP
    pd <- tryCatch({
      pdp::partial(model_fit,
                   pred.var = predictor,
                   grid.resolution = 100,
                   train = df,
                   progress = "text")
    }, error = function(e) {
      message("   ⚠️ Error generating PDP: ", e$message)
      return(NULL)
    })
    
    if (!is.null(pd)) {
      pd$ESM <- esm
      pd$Algorithm <- alg
      pdp_list[[paste(esm, alg, sep = "_")]] <- pd
      message("   ✅ PDP generated successfully.")
    }
  }
}

# --- Combinar
pdp_all <- bind_rows(pdp_list)

# --- Gràfic combinat
ggplot(pdp_all, aes(x = SeaSurfaceTemperature, y = yhat,
                    color = ESM, linetype = Algorithm)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 1) +
  scale_color_manual(values = c("MIROC-ES2L" = "#1b9e77", "MRI-ESM2-0" = "#d95f02")) +
  labs(
    x = "Sea Surface Temperature (°C)",
    y = "Predicted habitat suitability",
    title = "Partial Dependence of Temperature across ESMs and algorithms (RAW models)",
    color = "ESM",
    linetype = "Algorithm"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )




library(ggplot2)
library(dplyr)

# Si yhat és en escala logit, fem la inversa:
pdp_all <- pdp_all %>%
  mutate(yhat_prob = 1 / (1 + exp(-yhat)))  # desfem logit

ggplot(pdp_all, aes(x = SeaSurfaceTemperature, y = yhat_prob, color = ESM)) +
  # geom_point(alpha = 0.15, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 1) +
  scale_color_manual(values = c("MIROC-ES2L" = "#1b9e77", "MRI-ESM2-0" = "#d95f02")) +
  labs(
    x = "Sea Surface Temperature (°C)",
    y = "Predicted habitat suitability (probability)",
    title = "Temperature–suitability relationship (LOESS smoothed, back-transformed)",
    color = "ESM"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )







library(dplyr)
library(ggplot2)
library(lubridate)

# Afegim columna de temporada
ts_total_all <- ts_total_all %>%
  mutate(
    month = month(as.Date(date)),
    season = case_when(
      month %in% 5:10 ~ "breeding",
      month %in% c(11, 12, 1, 2, 3, 4) ~ "nonbreeding",
      TRUE ~ NA_character_
    )
  )

# Colors per escenari
cols_scen <- c(
  past1000  = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)



library(dplyr)
library(ggplot2)
library(lubridate)

# --- Assegurem columnes correctes ---
ts_total_all <- ts_total_all %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date),
    season = case_when(
      month %in% 5:10 ~ "breeding",
      TRUE ~ "nonbreeding"
    )
  )

# --- Resum per any, escenari i temporada (integrant els ESMs) ---
ts_summary_integrated <- ts_total_all %>%
  group_by(year, scenario, season) %>%
  summarise(mean_habitat = mean(habitat_km2, na.rm = TRUE), .groups = "drop")

# --- Colors per escenari ---
cols_scen <- c(
  past1000  = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)

# --- Gràfic integrat ---
ggplot(ts_summary_integrated,
       aes(x = year, y = mean_habitat,
           color = scenario, group = scenario)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 1.1, alpha = 0.2) +
  # facet_wrap(~ season, ncol = 1, scales = "free_y") +
  scale_color_manual(values = cols_scen) +
  labs(
    x = "Year",
    y = expression("Habitat area (km"^2*")"),
    title = "Integrated temporal evolution of total habitat area by scenario and season",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )





library(dplyr)
library(ggplot2)
library(lubridate)

# --- Preparar les dades ---
ts_total_all <- ts_total_all %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date),
    season = case_when(
      month %in% 5:10 ~ "breeding",
      TRUE ~ "nonbreeding"
    )
  )

# --- Calcular la mitjana anual de la superfície d’hàbitat per ESM i escenari ---
ts_annual_mean <- ts_total_all %>%
  group_by(esm, scenario, year) %>%
  summarise(mean_habitat = mean(habitat_km2, na.rm = TRUE), .groups = "drop")

# --- Colors per escenari ---
cols_scen <- c(
  past1000  = "black",
  historical = "#1f78b4",
  ssp126 = "#33a02c",
  ssp585 = "#e31a1c"
)

# --- Gràfic ---
ggplot(ts_annual_mean,
       aes(x = year, y = mean_habitat,
           color = scenario, group = scenario)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, linewidth = 1.1, alpha = 0.25) +
  facet_wrap(~ esm, ncol = 1, scales = "free_y") +
  scale_color_manual(values = cols_scen) +
  labs(
    x = "Year",
    y = expression("Mean annual habitat area (km"^2*")"),
    title = "Annual mean habitat area per ESM and scenario",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )









library(ggplot2)
library(dplyr)
library(readr)


# --- Paths ---
path_miroc <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw/MIROC-ES2L/SRW_extractRaw_150km_MIROC-ES2L.csv"
path_mri   <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw/MRI-ESM2-0/SRW_extractRaw_150km_MRI-ESM2-0.csv"

# --- Read data ---
df_miroc <- read_csv(path_miroc)
df_mri   <- read_csv(path_mri)

# --- Filter presences ---
pres_miroc <- df_miroc %>%
  filter(PresAbs == 1) %>%
  mutate(ESM = "MIROC-ES2L")

pres_mri <- df_mri %>%
  filter(PresAbs == 1) %>%
  mutate(ESM = "MRI-ESM2-0")

# --- Combine ---
pres_all <- bind_rows(pres_miroc, pres_mri)

# --- Plot histogram ---
ggplot(pres_all, aes(x = SeaSurfaceTemperature, fill = ESM)) +
  geom_histogram(position = "identity", alpha = 0.3, bins = 100, color = "grey20") +
  scale_fill_manual(values = c("MIROC-ES2L" = "#1b9e77", "MRI-ESM2-0" = "#d95f02")) +
  labs(
    x = "Sea Surface Temperature (°C)",
    y = "Count of presences",
    title = "Distribution of sea surface temperature at presence points",
    fill = "ESM"
  ) +
  # facet_wrap(~ Model, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )





# ============================================================
# Function to read and extract info (safe numeric extraction)
# ============================================================
read_cutoff <- function(file) {
  cutoff_value <- tryCatch({
    val <- readRDS(file)
    # Si és una llista o té camps interns
    if (is.list(val)) {
      as.numeric(val$CutoffValue %||% val$cutoff %||% val[[1]])
    } else {
      as.numeric(val)
    }
  }, error = function(e) NA_real_)
  
  esm <- stringr::str_extract(file, "(MIROC-ES2L|MRI-ESM2-0)")
  alg <- stringr::str_extract(file, "(GBM|RF|SVM)")
  
  tibble(
    File = basename(file),
    ESM = esm,
    Algorithm = alg,
    Cutoff = cutoff_value
  )
}

# ============================================================
# Read and combine all cutoffs
# ============================================================
cutoff_df <- purrr::map_dfr(rds_files, read_cutoff) %>%
  dplyr::filter(!is.na(Cutoff) & is.finite(Cutoff))

# ============================================================
# Plot bar chart (numeric-safe)
# ============================================================
ggplot(cutoff_df, aes(x = Algorithm, y = Cutoff, fill = ESM)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  geom_text(
    aes(label = sprintf("%.2f", Cutoff)),
    position = position_dodge(width = 0.7),
    vjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(values = c("MIROC-ES2L" = "#1b9e77", "MRI-ESM2-0" = "#d95f02")) +
  labs(
    title = "Cutoff values per ESM and algorithm (standardized models)",
    x = "Algorithm",
    y = "Cutoff (probability threshold)",
    fill = "ESM"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )





library(raster)
library(dplyr)
library(purrr)
library(ggplot2)
library(parallel)

# ===== CONFIG =====
ESMs <- c("MIROC-ES2L", "MRI-ESM2-0")
algorithms <- c("GBM", "RF", "SVM")
scenarios <- c("past1000", "historical")

base_bin_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
base_sst_dir <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/01envStacks_Raw_0s"

# ===== FUNCIO =====
extract_temp_from_binary <- function(bin_file, sst_file, esm, alg, scenario) {
  message(sprintf("→ %s | %s | %s | %s", esm, alg, scenario, basename(bin_file)))
  
  r_bin <- try(raster(bin_file), silent = TRUE)
  if (inherits(r_bin, "try-error") || is.null(r_bin)) return(NULL)
  
  sst_stack <- try(stack(sst_file), silent = TRUE)
  if (inherits(sst_stack, "try-error") || is.null(sst_stack)) return(NULL)
  
  layer_name <- names(sst_stack)[grepl("SeaSurfaceTemperature", names(sst_stack), ignore.case = TRUE)]
  if (length(layer_name) == 0) return(NULL)
  r_sst <- sst_stack[[layer_name[1]]]
  
  safe_compare <- try(compareRaster(r_bin, r_sst, extent = TRUE, rowcol = TRUE, crs = TRUE, stopIfNotEqual = FALSE), silent = TRUE)
  if (inherits(safe_compare, "try-error")) {
    r_sst <- try(resample(r_sst, r_bin, method = "bilinear"), silent = TRUE)
    if (inherits(r_sst, "try-error")) return(NULL)
  }
  
  mask_sst <- try(mask(r_sst, r_bin, maskvalue = 0), silent = TRUE)
  if (inherits(mask_sst, "try-error")) return(NULL)
  
  vals <- values(mask_sst)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NULL)
  
  data.frame(
    SST = vals,
    ESM = esm,
    Algorithm = alg,
    Scenario = scenario
  )
}

# ===== LOOP GLOBAL =====
results <- list()

for (esm in ESMs) {
  for (alg in algorithms) {
    message("\n===== ", esm, " | ", alg, " =====")
    
    for (scen in scenarios) {
      bin_dir <- file.path(base_bin_dir, esm, alg, scen)
      sst_dir <- file.path(base_sst_dir, esm, scen)
      
      bin_files <- list.files(bin_dir, pattern = "_binary\\.tif$", full.names = TRUE)
      sst_files <- list.files(sst_dir, pattern = "\\.grd$", full.names = TRUE)
      
      # Subsample per anar ràpid
      set.seed(1)
      bin_files <- head(bin_files, 50)
      sst_files <- head(sst_files, 50)
      
      df_temp <- mcmapply(
        extract_temp_from_binary,
        bin_files, sst_files,
        MoreArgs = list(esm = esm, alg = alg, scenario = scen),
        mc.cores = max(1, detectCores() - 2),
        SIMPLIFY = FALSE
      )
      
      df_temp <- bind_rows(df_temp)
      results[[paste(esm, alg, scen, sep = "_")]] <- df_temp
    }
  }
}

df_all <- bind_rows(results)

# ===== GRÀFIC =====
cols_scen <- c("past1000" = "#1b9e77", "historical" = "#d95f02")

ggplot(df_all, aes(x = SST, fill = Scenario)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity", bins = 50) +
  scale_fill_manual(values = cols_scen) +
  facet_grid(ESM ~ Algorithm, scales = "free_y") +
  labs(
    x = "Sea Surface Temperature (°C)",
    y = "Density",
    title = "Distribution of SST within suitable habitat (binary=1)\nacross ESMs and algorithms"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )





library(ggplot2)
library(dplyr)

# Ordenem els nivells del factor perquè 'historical' vagi al damunt
df_all$scenario <- factor(df_all$Scenario, levels = c("past1000", "historical"))

ggplot(df_all, aes(x = SST, fill = scenario)) +
  geom_histogram(
    position = "identity",
    alpha = 0.8,
    bins = 40,
    color = "grey30"
  ) +
  facet_grid(ESM ~ Algorithm) +
  scale_fill_manual(values = c("past1000" = "#66c2a5", "historical" = "#fc8d62")) +
  labs(
    x = "Sea Surface Temperature (°C)",
    y = "Number of suitable cells",
    title = "Distribution of SST within suitable habitat (binary=1)\nacross ESMs and algorithms"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )







library(ggplot2)
library(dplyr)

# Àrea aproximada per cel·la a 0.25º (~770 km²)
cell_area_km2 <- 770

df_bin_temp <- df_all

df_bin_temp$Scenario <- factor(df_bin_temp$Scenario,
                               levels = c("past1000", "historical"))

ggplot(df_bin_temp, aes(x = SST, fill = Scenario)) +
  # Past1000 (opac, sota)
  geom_histogram(
    data = subset(df_bin_temp, Scenario == "past1000"),
    aes(weight = cell_area_km2),  # converteix counts a km²
    bins = 50,
    position = "identity",
    color = "black",
    alpha = 1
  ) +
  # Historical (transparent, sobre)
  geom_histogram(
    data = subset(df_bin_temp, Scenario == "historical"),
    aes(weight = cell_area_km2),
    bins = 50,
    position = "identity",
    color = "black",
    alpha = 0.7
  ) +
  facet_grid(ESM ~ Algorithm, scales = "free_y") +
  scale_fill_manual(
    values = c("past1000" = "#66c2a5", "historical" = "#fc8d62")
  ) +
  labs(
    x = "Sea Surface Temperature (°C)",
    y = "Suitable habitat area (km²)",
    title = "Distribution of SST within suitable habitat (binary=1)\nacross ESMs and algorithms"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )




library(raster)
library(stringr)
library(dplyr)
library(ggplot2)
library(parallel)
library(lubridate)

# ===================== CONFIG =====================
base_bin_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
base_sst_dir <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/01envStacks_Raw_0s"

ESMs       <- c("MIROC-ES2L", "MRI-ESM2-0")
algorithms <- c("GBM", "RF", "SVM")
scenarios  <- c("past1000", "historical")

# Mostra només una part dels fitxers per combo
max_files_per_combo <- 80  # ajusta segons la velocitat que vulguis

ncores <- max(1, detectCores() - 2)

# ===================== FUNCIONS =====================
extract_sst_in_habitat <- function(fbin, esm, scen) {
  date_str <- str_sub(basename(fbin), 1, 8)
  dtt <- suppressWarnings(ymd(date_str))
  
  sst_stack_path <- file.path(
    base_sst_dir, esm, scen,
    sprintf("monthly_stack_%s_%s_%s.grd", esm, scen, date_str)
  )
  if (!file.exists(sst_stack_path)) return(NULL)
  
  rbin <- try(raster(fbin), silent = TRUE)
  if (inherits(rbin, "try-error")) return(NULL)
  
  sst_stk <- try(stack(sst_stack_path), silent = TRUE)
  if (inherits(sst_stk, "try-error")) return(NULL)
  if (!("SeaSurfaceTemperature" %in% names(sst_stk))) return(NULL)
  
  sst <- sst_stk$SeaSurfaceTemperature
  
  sst2 <- try(projectRaster(sst, rbin, method = "bilinear"), silent = TRUE)
  if (inherits(sst2, "try-error")) return(NULL)
  
  rbin01 <- calc(rbin, fun = function(x) as.numeric(x >= 1))
  sst_mask <- try(mask(sst2, rbin01, maskvalue = 0), silent = TRUE)
  if (inherits(sst_mask, "try-error")) return(NULL)
  
  vals <- getValues(sst_mask)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(NULL)
  
  data.frame(
    date = dtt,
    esm = esm,
    scenario = scen,
    mean_sst = mean(vals, na.rm = TRUE)
  )
}

# ===================== BUCLE PRINCIPAL =====================
results <- list()

for (esm in ESMs) {
  for (alg in algorithms) {
    for (scen in scenarios) {
      
      folder_bin <- file.path(base_bin_dir, esm, alg, scen)
      if (!dir.exists(folder_bin)) next
      
      tif_files <- list.files(folder_bin, pattern = "_prediction_binary\\.tif$", full.names = TRUE)
      if (!length(tif_files)) next
      
      # Subsample aleatori per anar més ràpid
      set.seed(42)
      tif_files <- sample(tif_files, min(length(tif_files), max_files_per_combo))
      
      message(sprintf("\n=== %s | %s | %s ===", esm, alg, scen))
      message("Fitxers processats: ", length(tif_files))
      
      combo_df <- mclapply(
        tif_files,
        extract_sst_in_habitat,
        esm = esm,
        scen = scen,
        mc.cores = ncores
      )
      
      combo_df <- do.call(rbind, combo_df)
      if (!is.null(combo_df) && nrow(combo_df) > 0) {
        combo_df$algorithm <- alg
        results[[length(results) + 1]] <- combo_df
      }
    }
  }
}

sst_in_hab <- do.call(rbind, results)
sst_in_hab$year <- year(sst_in_hab$date)

# ===================== GRÀFIC =====================
cols_scen <- c(past1000 = "black", historical = "#1f78b4")

ggplot(sst_in_hab,
       aes(x = year, y = mean_sst,
           color = scenario)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.3, linewidth = 0.9) +
  facet_grid(esm ~ algorithm, scales = "free_y") +
  scale_color_manual(values = cols_scen) +
  labs(
    x = "Any",
    y = "SST mitjana dins hàbitat (°C)",
    title = "Evolució de la SST dins de les cel·les amb hàbitat (subsample)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )




# ==============================================================================
# SST ALS PUNTS DE PRESÈNCIA (PresAbs = 1)
# ==============================================================================

library(raster)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(purrr)

# --- DIRECTORIS ---------------------------------------------------------------
base_presabs_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw"
base_sst_dir <- "/Volumes/MyPassport/2024_SouthernRightWhaleLogbooks/02envStacks"

# --- CONFIGURACIÓ -------------------------------------------------------------
ESMs <- c("MIROC-ES2L", "MRI-ESM2-0")
scenarios <- c("past1000", "historical")

# --- FUNCIÓ PER EXTREURE SST ALS PUNTS DE PRESÈNCIA --------------------------
extract_sst_at_presences <- function(presabs_csv, esm, scen, n_sample = NULL) {
  df <- read.csv(presabs_csv)
  df <- df %>% filter(PresAbs == 1)
  if (nrow(df) == 0) return(NULL)
  
  # opcionalment submostreig
  if (!is.null(n_sample) && nrow(df) > n_sample) {
    df <- df %>% sample_n(n_sample)
  }
  
  coords <- df %>% dplyr::select(Lon, Lat)
  pts <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # obtenir data del fitxer
  date_str <- str_sub(basename(presabs_csv), 1, 8)
  dtt <- suppressWarnings(ymd(date_str))
  
  # trobar stack mensual corresponent
  sst_stack_path <- file.path(
    base_sst_dir, esm, scen,
    sprintf("monthly_stack_%s_%s_%s.grd", esm, scen, date_str)
  )
  
  if (!file.exists(sst_stack_path)) return(NULL)
  
  sst_stk <- try(stack(sst_stack_path), silent = TRUE)
  if (inherits(sst_stk, "try-error")) return(NULL)
  if (!("SeaSurfaceTemperature" %in% names(sst_stk))) return(NULL)
  
  sst <- sst_stk$SeaSurfaceTemperature
  vals <- raster::extract(sst, pts)
  
  data.frame(
    date = dtt,
    esm = esm,
    scenario = scen,
    mean_sst_pres = mean(vals, na.rm = TRUE),
    n_points = length(vals[is.finite(vals)])
  )
}

# --- BUCLE PER RECOLLIR DADES ------------------------------------------------
sst_pres_all <- list()

for (esm in ESMs) {
  for (scen in scenarios) {
    message("\n→ Processant ", esm, " — ", scen)
    
    presabs_dir <- file.path(base_presabs_dir, esm)
    presabs_files <- list.files(presabs_dir, pattern = sprintf("%s.*\\.csv$", esm), full.names = TRUE)
    
    if (length(presabs_files) == 0) next
    
    # mostreig d’uns pocs per velocitat
    presabs_files <- sample(presabs_files, min(20, length(presabs_files)))
    
    df_list <- map(presabs_files, ~ extract_sst_at_presences(.x, esm, scen, n_sample = 1000))
    df_combined <- bind_rows(df_list)
    
    if (!is.null(df_combined)) sst_pres_all[[paste(esm, scen, sep = "_")]] <- df_combined
  }
}

sst_pres_all <- bind_rows(sst_pres_all) %>% na.omit()
sst_pres_all$year <- year(sst_pres_all$date)

# --- GRÀFIC -------------------------------------------------------------------
ggplot(sst_pres_all, aes(x = year, y = mean_sst_pres, color = scenario)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25, size = 1.1) +
  facet_wrap(~ esm, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("past1000" = "#66c2a5", "historical" = "#fc8d62")) +
  labs(
    title = "Evolució de la SST als punts de presència",
    subtitle = "Comparació entre escenaris past1000 i historical",
    x = "Any",
    y = "SST mitjana a presències (°C)",
    color = "Escenari"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )


# ANÀLISI ESPACIAL


library(raster)
library(dplyr)
library(ggplot2)

# === CONFIG ===
base_sst_dir <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/01envStacks_Raw_0s"
ESMs <- c("MIROC-ES2L", "MRI-ESM2-0")

# --- Funció per obtenir els 10 primers / últims fitxers ---
get_subset_files <- function(dir_path, n = 10) {
  files <- list.files(dir_path, pattern = "monthly_stack_.*\\.grd$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  dates <- as.Date(sub(".*_(\\d{8})\\.grd$", "\\1", basename(files)), format = "%Y%m%d")
  df <- data.frame(file = files, date = dates)
  df <- df[order(df$date), ]
  list(first = head(df$file, n), last = tail(df$file, n))
}

# --- Bucle principal ---
diff_list <- list()

for (esm in ESMs) {
  message("\n→ Processing ", esm)
  
  dir_past <- file.path(base_sst_dir, esm, "past1000")
  dir_hist <- file.path(base_sst_dir, esm, "historical")
  
  if (!dir.exists(dir_past) | !dir.exists(dir_hist)) {
    message("   ⚠️ Missing directories for ", esm)
    next
  }
  
  subset_past <- get_subset_files(dir_past, n = 10)
  subset_hist <- get_subset_files(dir_hist, n = 10)
  
  if (is.null(subset_past) | is.null(subset_hist)) next
  
  message("   Averaging last 10 past1000 and first 10 historical files...")
  
  # --- Past1000 ---
  past_stack <- stack(subset_past$last)
  vars_past <- grep("SeaSurfaceTemperature", names(past_stack), value = TRUE)
  if (length(vars_past) == 0) {
    message("   ⚠️ No SeaSurfaceTemperature found in past1000 for ", esm)
    next
  }
  past_mean <- mean(past_stack[[vars_past]], na.rm = TRUE)
  
  # --- Historical ---
  hist_stack <- stack(subset_hist$first)
  vars_hist <- grep("SeaSurfaceTemperature", names(hist_stack), value = TRUE)
  if (length(vars_hist) == 0) {
    message("   ⚠️ No SeaSurfaceTemperature found in historical for ", esm)
    next
  }
  hist_mean <- mean(hist_stack[[vars_hist]], na.rm = TRUE)
  
  # Resample històric a malla del past1000
  hist_mean_rs <- try(resample(hist_mean, past_mean, method = "bilinear"), silent = TRUE)
  if (inherits(hist_mean_rs, "try-error")) next
  
  # Diferència: històric - passat
  diff_r <- hist_mean_rs - past_mean
  diff_list[[esm]] <- diff_r
  
  # Imprimir el promig global del canvi
  delta_mean <- cellStats(diff_r, mean, na.rm = TRUE)
  message(sprintf("   Mean ΔSST (%s): %.3f °C", esm, delta_mean))
}

# --- Passar a dataframe per ggplot ---
plot_list <- list()
for (esm in names(diff_list)) {
  df <- as.data.frame(rasterToPoints(diff_list[[esm]]))
  colnames(df) <- c("lon", "lat", "diff_sst")
  df$ESM <- esm
  plot_list[[esm]] <- df
}
df_all <- bind_rows(plot_list)

# --- Gràfic ---
ggplot(df_all, aes(lon, lat, fill = diff_sst)) +
  geom_raster() +
  facet_wrap(~ ESM, ncol = 1) +
  scale_fill_gradient2(
    name = expression(Delta * "SST (°C)"),
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-2, 2), oob = scales::squish
  ) +
  coord_quickmap() +
  labs(
    title = "ΔSST between late past1000 (last 10 years) and early historical (first 10 years)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )


