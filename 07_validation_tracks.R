# ============================================================
# VALIDACIÓ BINÀRIA AMB DADES DE TRACKING (2014–2024)
# ============================================================

library(raster)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(ggplot2)
library(fs)

# === CONFIGURACIÓ ===
ESMs <- c("MIROC-ES2L", "MRI-ESM2-0")
algorithms <- c("GBM", "RF", "SVM")
scenario <- "ssp126"

base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
tracking_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"

# === 1. LLEGIR DADES DE TRACKING ===
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>% bind_rows() %>%
  mutate(
    date = suppressWarnings(lubridate::dmy_hm(date)),
    year = year(date),
    month = month(date),
    season = ifelse(month %in% 5:10, "breeding", "nonbreeding"),
    lon_round = round(lon / 0.25) * 0.25,
    lat_round = round(lat / 0.25) * 0.25
  ) %>%
  filter(!is.na(lon_round) & !is.na(lat_round))

# Llista d’anys reals en les dades
years_tracking <- sort(unique(tracking_all$year))
message("Tracking disponible entre: ", min(years_tracking), "–", max(years_tracking))

# === 2. FUNCIÓ AUXILIAR: rasteritzar tracking ===
create_tracking_raster <- function(df, template_raster) {
  coords <- df %>% dplyr::select(lon_round, lat_round)
  if (nrow(coords) == 0) return(NULL)
  sp_pts <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
  r <- rasterize(sp_pts, template_raster, field = 1, fun = "max", background = 0)
  return(r)
}

# === 3. CÀLCUL DE VALIDACIÓ ===
results <- list()

for (esm in ESMs) {
  for (alg in algorithms) {
    message("\n→ Processant ", esm, " — ", alg)
    
    pred_dir <- file.path(base_pred_dir, esm, alg, scenario)
    pred_files <- list.files(pred_dir, pattern = "binary\\.tif$", full.names = TRUE)
    
    if (length(pred_files) == 0) {
      message("   ⚠️ No s'han trobat fitxers per ", esm, "-", alg)
      next
    }
    
    # --- Iterem per any i època ---
    yearly_results <- data.frame()
    
    for (yy in years_tracking) {
      for (ss in c("breeding", "nonbreeding")) {
        
        # Subset del tracking per any i època
        df_sub <- tracking_all %>%
          filter(year == yy, season == ss)
        
        if (nrow(df_sub) == 0) next
        
        # Seleccionem predicció més propera temporalment
        pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
        if (all(is.na(pred_dates))) next
        
        target_date <- median(df_sub$date)
        # Converteix la diferència de dates a dies (numèric)
        date_diffs <- abs(as.numeric(difftime(pred_dates, target_date, units = "days")))
        
        # Agafa la predicció temporalment més propera
        idx_close <- which.min(date_diffs)       
        if (length(idx_close) == 0) next
        
        f_pred <- pred_files[idx_close]
        r_pred <- try(raster(f_pred), silent = TRUE)
        if (inherits(r_pred, "try-error")) next
        
        # Raster de presència tracking sobre la mateixa malla
        r_track <- create_tracking_raster(df_sub, r_pred)
        if (is.null(r_track)) next
        
        # Extreure valors
        vals_pred <- getValues(r_pred)
        vals_track <- getValues(r_track)
        
        # Evitar NA
        valid <- which(!is.na(vals_pred) & !is.na(vals_track))
        vals_pred <- vals_pred[valid]
        vals_track <- vals_track[valid]
        
        # Calcular mètriques
        tp <- sum(vals_pred == 1 & vals_track == 1)
        fn <- sum(vals_pred == 0 & vals_track == 1)
        fp <- sum(vals_pred == 1 & vals_track == 0)
        tn <- sum(vals_pred == 0 & vals_track == 0)
        
        sensitivity <- ifelse((tp + fn) > 0, tp / (tp + que ), NA)
        precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
        
        yearly_results <- bind_rows(yearly_results, data.frame(
          ESM = esm,
          Algorithm = alg,
          Year = yy,
          Season = ss,
          TP = tp,
          FN = fn,
          FP = fp,
          TN = tn,
          Sensitivity = sensitivity,
          Precision = precision
        ))
      }
    }
    
    # Guardar resultats per combinació
    if (nrow(yearly_results) > 0) {
      yearly_results <- yearly_results %>%
        group_by(ESM, Algorithm, Season) %>%
        summarise(across(c(Sensitivity, Precision), mean, na.rm = TRUE)) %>%
        mutate(Scenario = scenario)
      results[[paste(esm, alg, sep = "_")]] <- yearly_results
    }
  }
}

results_all <- bind_rows(results)

# === 4. GRÀFICS ===

ggplot(results_all, aes(x = Algorithm, y = Sensitivity, fill = ESM)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Season) +
  scale_fill_manual(values = c("MIROC-ES2L" = "#1b9e77", "MRI-ESM2-0" = "#d95f02")) +
  labs(
    title = "Model validation using tracking presence data",
    subtitle = "Mean sensitivity by ESM, algorithm, and season (breeding vs nonbreeding)",
    y = "Sensitivity", x = "Algorithm", fill = "ESM"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# --- També pots veure precisió ---
ggplot(results_all, aes(x = Algorithm, y = Precision, fill = ESM)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Season) +
  scale_fill_manual(values = c("MIROC-ES2L" = "#1b9e77", "MRI-ESM2-0" = "#d95f02")) +
  labs(
    title = "Precision of model predictions",
    subtitle = "Proportion of predicted suitable cells with tracking presences",
    y = "Precision", x = "Algorithm", fill = "ESM"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))









library(raster)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(sf)

# === CONFIGURACIÓ ===
esm <- "MIROC-ES2L"
alg <- "GBM"
scenario <- "ssp126"
year_target <- 2015
months_breeding <- 5:10  # maig a octubre

base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
tracking_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"

# === LLEGIR FITXERS BINARIS ===
pred_dir <- file.path(base_pred_dir, esm, alg, scenario)
pred_files <- list.files(pred_dir, pattern = "binary\\.tif$", full.names = TRUE)

# Extreure dates
pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
sel_idx <- which(year(pred_dates) == year_target & month(pred_dates) %in% months_breeding)
pred_files <- pred_files[sel_idx]

if (length(pred_files) == 0) stop("No s'han trobat prediccions per la temporada breeding de ", year_target)

# Carregar i combinar prediccions
r_stack <- stack(pred_files)

# Cel·les bones si han estat 1 en almenys un mes
r_any <- calc(r_stack, fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))

# === LLEGIR DADES DE TRACKING ===
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>%
  bind_rows() %>%
  mutate(
    date = suppressWarnings(dmy_hm(date)),
    lon_round = round(lon / 0.25) * 0.25,
    lat_round = round(lat / 0.25) * 0.25
  ) %>%
  filter(year(date) == year_target & month(date) %in% months_breeding)

pts_sf <- st_as_sf(tracking_all, coords = c("lon_round", "lat_round"), crs = 4326)

# === CONVERTIR RASTER A DATAFRAME ===
r_df <- as.data.frame(r_any, xy = TRUE)
colnames(r_df) <- c("x", "y", "pred")
r_df$pred <- as.factor(ifelse(r_df$pred == 1, "Suitable", "Unsuitable"))

# === MAPA ===
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = pred)) +
  scale_fill_manual(
    values = c("Unsuitable" = "grey85", "Suitable" = "#1b9e77"),
    name = "Predicció"
  ) +
  geom_sf(data = pts_sf, aes(geometry = geometry), inherit.aes = FALSE,
          color = "red", fill = "red", shape = 21, size = 0.05, alpha = 0.8)  +
  coord_sf(
    xlim = range(r_df$x, na.rm = TRUE),
    ylim = range(r_df$y, na.rm = TRUE),
    expand = FALSE
  ) +
  labs(
    title = paste("Predicció d'hàbitat (1 = adequat en algun mes)"),
    subtitle = paste(esm, "-", alg, "-", scenario, " | Breeding season", year_target, "(maig–octubre)"),
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )







library(sf)
library(raster)
library(dplyr)

# --- BUFFER RADIUS (en graus) ---
buffer_radius <- 0.25  # ~25 km a l'equador

# --- Converteix el raster binari en polígons d’hàbitat adequat ---
r_poly <- rasterToPolygons(r_any, fun = function(x) x == 1, dissolve = TRUE)
r_poly_sf <- st_as_sf(r_poly)

# --- Fes buffer al voltant dels punts de tracking ---
pts_buffer <- st_buffer(pts_sf, dist = buffer_radius)

# --- Comprova quants buffers intersecten amb zones adequades ---
intersections <- st_intersects(pts_buffer, r_poly_sf, sparse = FALSE)
pts_sf$hit <- apply(intersections, 1, any)

# --- Càlcul de sensibilitat ---
total_points <- nrow(pts_sf)
suitable_points <- sum(pts_sf$hit, na.rm = TRUE)
sensitivity <- suitable_points / total_points

cat("Sensibilitat (amb buffer de", buffer_radius, "°): ",
    round(sensitivity, 3), "\n")
cat("Total punts:", total_points, "| Punts dins hàbitat adequat:", suitable_points, "\n")

# --- Mapa visual per comprovar ---
ggplot() +
  geom_raster(data = as.data.frame(r_any, xy = TRUE),
              aes(x = x, y = y, fill = as.factor(layer))) +
  scale_fill_manual(values = c("0" = "grey85", "1" = "#1b9e77"),
                    name = "Predicció") +
  geom_sf(data = r_poly_sf, fill = NA, color = "darkgreen", size = 0.3) +
  geom_sf(data = pts_sf, aes(color = hit), size = 0.8) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue"),
                     name = "Coincideix?") +
  # # coord_sf(xlim = range(r_poly_sf$geometry[[1]][[1]][,1]),
  # #          ylim = range(r_poly_sf$geometry[[1]][[1]][,2]),
  #          expand = FALSE) +
  labs(
    title = paste("Validació per buffer (", esm, "-", alg, "-", scenario, ")", sep = ""),
    subtitle = paste("Temporada breeding", year_target, " — buffer =", buffer_radius, "°"),
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )





library(sf)
library(raster)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# --- BUFFER RADIUS ---
buffer_radius <- 0.25  # en graus

# --- Convertir raster binari en polígons adequats ---
r_poly <- rasterToPolygons(r_any, fun = function(x) x == 1, dissolve = TRUE)
r_poly_sf <- st_as_sf(r_poly)

# --- Fer buffer al voltant dels punts de tracking ---
pts_buffer <- st_buffer(pts_sf, dist = buffer_radius)

# --- Comprovar interseccions ---
intersections <- st_intersects(pts_buffer, r_poly_sf, sparse = FALSE)
pts_sf$hit <- apply(intersections, 1, any)

# --- Càlcul de sensibilitat ---
total_points <- nrow(pts_sf)
suitable_points <- sum(pts_sf$hit, na.rm = TRUE)
sensitivity <- suitable_points / total_points

cat("Sensibilitat (amb buffer de", buffer_radius, "°): ",
    round(sensitivity, 3), "\n")
cat("Total punts:", total_points, "| Punts dins hàbitat adequat:", suitable_points, "\n")

# --- MAPA DEL MÓN ---
world <- ne_countries(scale = "medium", returnclass = "sf")

# --- DEFINIR LÍMITS PER A L'ARGENTINA ---
zoom_extent <- list(
  xlim = c(-70, -45),
  ylim = c(-60, -30)
)

# --- MAPA ---
ggplot() +
  # Base del món
  geom_sf(data = world, fill = "grey95", color = "grey60", size = 0.3) +
  
  # Hàbitat adequat predit
  geom_raster(data = as.data.frame(r_any, xy = TRUE),
              aes(x = x, y = y, fill = as.factor(layer)), alpha = 0.8) +
  scale_fill_manual(values = c("0" = "grey85", "1" = "#1b9e77"),
                    name = "Predicció") +
  
  # Polígons d’hàbitat
  geom_sf(data = r_poly_sf, fill = NA, color = "darkgreen", size = 0.2) +
  
  # Punts de tracking (blau = encert, vermell = error)
  geom_sf(data = pts_sf, aes(color = hit), size = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue"),
                     name = "Coincideix amb hàbitat") +
  
  # Zoom sobre l’Argentina
  coord_sf(
    xlim = zoom_extent$xlim,
    ylim = zoom_extent$ylim,
    expand = FALSE
  ) +
  
  labs(
    title = paste("Validació visual de predicció d'hàbitat"),
    subtitle = paste(esm, "-", alg, "-", scenario, "| Breeding season", year_target, "| Buffer =", buffer_radius, "°"),
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal(base_size = 13) 






library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)

# --- CONFIGURACIÓ ---
esm <- "MRI-ESM2-0"
alg <- "SVM"
scenario <- "ssp126"
year_target <- 2015
buffer_deg <- 0.25

base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
tracking_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"

# --- LLEGIR RASTER BINARI ---
pred_dir <- file.path(base_pred_dir, esm, alg, scenario)
pred_files <- list.files(pred_dir, pattern = "2015.*binary\\.tif$", full.names = TRUE)

if (length(pred_files) == 0) stop("No binary files for 2015 in ", pred_dir)

r_bin <- raster(pred_files[1])  # agafem un fitxer d'exemple
r_df <- as.data.frame(r_bin, xy = TRUE)
colnames(r_df) <- c("x", "y", "pred")

# --- LLEGIR TRACKING ---
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>%
  bind_rows() %>%
  mutate(date = suppressWarnings(lubridate::dmy_hm(date))) %>%
  filter(year(date) == year_target & month(date) %in% 5:10)  # breeding

# --- CONVERTIR A SF ---
pts_sf <- st_as_sf(tracking_all, coords = c("lon", "lat"), crs = 4326)

# --- COMPROVAR VALOR BINARI EN CADA PUNT ---
vals <- raster::extract(r_bin, pts_sf)
pts_sf$pred_val <- vals
pts_sf$hit <- pts_sf$pred_val == 1

# --- COMPTAR QUANTS PUNTS CAUEN EN CEL·LES NUL·LES ---
na_ratio <- sum(is.na(pts_sf$pred_val)) / nrow(pts_sf)
zero_ratio <- sum(pts_sf$pred_val == 0, na.rm = TRUE) / nrow(pts_sf)
one_ratio <- sum(pts_sf$pred_val == 1, na.rm = TRUE) / nrow(pts_sf)

cat("Proporció de punts sobre cel·les NA:", round(na_ratio * 100, 2), "%\n")
cat("Proporció de punts sobre cel·les amb 0:", round(zero_ratio * 100, 2), "%\n")
cat("Proporció de punts sobre cel·les amb 1:", round(one_ratio * 100, 2), "%\n")

# --- MAPA: zoom a Argentina ---
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = as.factor(pred))) +
  scale_fill_manual(
    values = c("0" = "grey80", "1" = "#1b9e77"),
    name = "Predicció binària"
  ) +
  geom_sf(data = world, fill = NA, color = "black", size = 0.2) +
  geom_sf(data = pts_sf, aes(color = hit), size = 1.3, alpha = 0.9) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue"),
                     name = "Coincideix amb hàbitat") +
  coord_sf(xlim = c(-70, -45), ylim = c(-55, -30), expand = FALSE) +
  labs(
    title = paste(esm, "-", alg, "-", scenario),
    subtitle = paste("Breeding season", year_target, "| Buffer:", buffer_deg, "°"),
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )





# --- ELIMINEM PUNTS EN CEL·LES NA (fora del domini del model) ---
pts_valid <- pts_sf %>% filter(!is.na(pred_val))

# --- RECALCULEM PROPORCIONS SOBRE CEL·LES MARINES ---
na_ratio_new <- 0  # ja no hi ha NA
zero_ratio_new <- sum(pts_valid$pred_val == 0, na.rm = TRUE) / nrow(pts_valid)
one_ratio_new <- sum(pts_valid$pred_val == 1, na.rm = TRUE) / nrow(pts_valid)

cat("Recalculat només per punts dins del domini del model (oceà):\n")
cat("  Cel·les 0:", round(zero_ratio_new * 100, 2), "%\n")
cat("  Cel·les 1:", round(one_ratio_new * 100, 2), "%\n")

# --- NOVA SENSIBILITAT ---
sensitivity_new <- one_ratio_new
cat("\nNova sensibilitat (limitada a àrea marina):", round(sensitivity_new * 100, 2), "%\n")



library(sf)
library(raster)
library(dplyr)
library(ggplot2)

# === BUFFER DE 0.5° ===
buffer_deg <- 0.5  # radi del buffer (graus)

# Convertim punts a sf si no ho són
if (!inherits(pts_sf, "sf")) {
  pts_sf <- st_as_sf(tracking_all, coords = c("lon_round", "lat_round"), crs = 4326)
}

# Fem el buffer (en graus; CRS 4326 -> metres si projectes)
pts_buf <- st_buffer(pts_sf, dist = buffer_deg)

# Convertim el raster a polígon per superposar
pred_poly <- rasterToPolygons(r_any, dissolve = FALSE) %>% st_as_sf()

# Ens quedem només amb cel·les "1"
pred_suitable <- pred_poly %>% filter(layer == 1)

# Intersecció buffer–hàbitat
intersect_sf <- st_intersects(pts_buf, pred_suitable, sparse = FALSE)
hits <- rowSums(intersect_sf) > 0

# Afegim columna de coincidència
pts_sf$coincideix <- hits

# === RECALCUL SENSIBILITAT (només punts dins del domini marí) ===
valid_pts <- pts_sf %>% filter(!is.na(coincideix))
sensitivity_buf <- mean(valid_pts$coincideix, na.rm = TRUE)

cat("\nSensibilitat amb buffer de", buffer_deg, "°:", round(sensitivity_buf * 100, 2), "%\n")

# === MAPA DE RESULTAT ===
ggplot() +
  geom_sf(data = pred_suitable, fill = "#66c2a5", color = NA, alpha = 0.4) +
  geom_sf(data = pts_buf, fill = NA, color = "black", size = 0.1, alpha = 0.3) +
  geom_sf(data = pts_sf, aes(color = coincideix), size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue")) +
  coord_sf(xlim = c(-70, -55), ylim = c(-55, -35), expand = FALSE) +
  labs(
    title = paste(esm, "-", alg, "-", scenario),
    subtitle = paste("Breeding season", year_target, "| Buffer =", buffer_deg, "°"),
    x = "Longitud", y = "Latitud", color = "Coincideix"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))









library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(rnaturalearth)

# ===== CONFIG =====
esm       <- "MIROC-ES2L"
esm       <- "MRI-ESM2-0"
alg       <- "GBM"
scenario  <- "ssp126"              # 2015 és ja SSP
year_tgt  <- 2015
months_br <- c(1,2,3,4,11,12)                  # breeding: maig–oct

# Argentina approx (ajusta si cal)
bbox_lon <- c(-75, -45)
bbox_lat <- c(-60, -30)

base_env_dir  <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/02envStacks_Standardized_0s"
base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
tracking_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"

# ===== 1) TRACKING: filtrar Argentina + Breeding 2015 =====
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>% bind_rows()

tracking_all <- tracking_all %>%
  mutate(date = suppressWarnings(lubridate::dmy_hm(date))) %>%
  filter(!is.na(date),
         lubridate::year(date) == year_tgt,
         lubridate::month(date) %in% months_br,
         lon >= bbox_lon[1], lon <= bbox_lon[2],
         lat >= bbox_lat[1], lat <= bbox_lat[2])

pts_sf <- st_as_sf(tracking_all, coords = c("lon","lat"), crs = 4326)

# ===== 2) PREDICCIONS BINÀRIES: agregat (almenys un mes = 1) =====
pred_dir   <- file.path(base_pred_dir, esm, alg, scenario)
pred_files <- list.files(pred_dir, pattern = "binary\\.tif$", full.names = TRUE)

pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
sel_idx <- which(year(pred_dates) == year_tgt & month(pred_dates) %in% months_br)
stopifnot(length(sel_idx) > 0)

r_stack <- stack(pred_files[sel_idx])
r_any   <- calc(r_stack, fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))

# Retalla a la regió d'interès
bbox_sp <- as(extent(bbox_lon[1], bbox_lon[2], bbox_lat[1], bbox_lat[2]), "SpatialPolygons")
proj4string(bbox_sp) <- CRS("+proj=longlat +datum=WGS84")
r_any_roi <- crop(r_any, bbox_sp)

# ===== 3) MÀSCARA DE RESOLUCIÓ (qualsevol .grd vàlid del mateix ESM/escenari) =====
env_any <- list.files(file.path(base_env_dir, esm, scenario), pattern="\\.grd$", full.names = TRUE)[1]
env_stack <- stack(env_any)
env_mask  <- env_stack[[1]]; env_mask[!is.na(env_mask)] <- 1
env_mask  <- crop(env_mask, bbox_sp)
# Reprojecta la màscara si cal per coincidir amb la malla de predicció
if (!compareRaster(env_mask, r_any_roi, crs=TRUE, res=TRUE, stopIfNotEqual=FALSE)) {
  env_mask <- projectRaster(env_mask, r_any_roi, method="ngb")
}

# ===== 4) CLASSIFICAR PUNTS (inside / coastline_gap / outside) =====
vals_pred <- raster::extract(r_any_roi, st_coordinates(pts_sf))
vals_mask <- raster::extract(env_mask,  st_coordinates(pts_sf))

pts_sf$cat <- "outside"
pts_sf$cat[!is.na(vals_pred) & vals_pred == 1] <- "inside"
pts_sf$cat[is.na(vals_mask)] <- "coastline_gap"

# ===== 5) PREPARAR DFs per ggplot =====
r_df   <- as.data.frame(r_any_roi, xy = TRUE); colnames(r_df) <- c("x","y","pred")
env_df <- as.data.frame(env_mask,   xy = TRUE); colnames(env_df) <- c("x","y","mask")

world  <- ne_countries(scale="medium", returnclass="sf")

# ===== 6) MAPA (ZOOM ARGENTINA, breeding 2015) =====
ggplot() +
  # Màscara de resolució (oceà resolt)
  geom_raster(data = subset(env_df, !is.na(mask)), aes(x = x, y = y), fill = "grey90") +
  # Hàbitat binari (almenys un mes = 1)
  geom_raster(data = subset(r_df, pred == 1), aes(x = x, y = y), fill = "#1b9e77", alpha = 0.65) +
  # Costes
  geom_sf(data = world, fill = "grey70", color = "white", linewidth = 0.2) +
  # Punts classificats
  geom_sf(data = subset(pts_sf, cat == "inside"),        color = "blue3",   size = 0.3, alpha = 0.4) +
  geom_sf(data = subset(pts_sf, cat == "coastline_gap"), color = "orange",  size = 0.3, alpha = 0.4) +
  geom_sf(data = subset(pts_sf, cat == "outside"),       color = "red3",    size = 0.3, alpha = 0.4) +
  coord_sf(xlim = bbox_lon, ylim = bbox_lat, expand = FALSE) +
  labs(
    title = "Breeding season 2015 — Argentina",
    subtitle = paste(esm, "|", alg, "|", scenario, "— prediccions (≥1 mes=1) vs tracking"),
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    plot.subtitle= element_text(hjust = 0.5),
    panel.grid   = element_blank()
  )
















# ======================================================================
# VALIDACIÓ ESPACIAL PER TOTS ELS ESMs I ALGORISMES (Argentina, 2015)
# ======================================================================
library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(rnaturalearth)

# ===== CONFIG =====
ESMs       <- c("MIROC-ES2L", "MRI-ESM2-0")
ALGORITHMS <- c("GBM", "RF", "SVM")
scenario   <- "ssp126"
year_tgt   <- 2015
months_br  <- c(5,6,7,8,9,10)  # breeding
months_nb  <- c(1,2,3,4,11,12) # non-breeding

# Argentina approx (ajusta si cal)
bbox_lon <- c(-75, -45)
bbox_lat <- c(-60, -30)

base_env_dir  <- "/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/02envStacks_Standardized_0s"
base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
tracking_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"

# ===== TRACKING: Argentina + breeding 2015 =====
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>% bind_rows() %>%
  mutate(date = suppressWarnings(lubridate::dmy_hm(date))) %>%
  filter(!is.na(date),
         lubridate::year(date) == year_tgt,
         lubridate::month(date) %in% months_br,
         lon >= bbox_lon[1], lon <= bbox_lon[2],
         lat >= bbox_lat[1], lat <= bbox_lat[2])
pts_sf <- st_as_sf(tracking_all, coords = c("lon","lat"), crs = 4326)

# ===== PREP ENV MASK (del primer ESM com a referència) =====
env_any <- list.files(file.path(base_env_dir, ESMs[1], scenario), pattern="\\.grd$", full.names = TRUE)[1]
env_stack <- stack(env_any)
env_mask  <- env_stack[[1]]; env_mask[!is.na(env_mask)] <- 1
bbox_sp <- as(extent(bbox_lon[1], bbox_lon[2], bbox_lat[1], bbox_lat[2]), "SpatialPolygons")
proj4string(bbox_sp) <- CRS("+proj=longlat +datum=WGS84")
env_mask  <- crop(env_mask, bbox_sp)

# ===== BUCLE PRINCIPAL =====
maps_list <- list()

for (esm in ESMs) {
  for (alg in ALGORITHMS) {
    message("\n=== Processant ", esm, " — ", alg, " ===")
    
    pred_dir   <- file.path(base_pred_dir, esm, alg, scenario)
    pred_files <- list.files(pred_dir, pattern = "binary\\.tif$", full.names = TRUE)
    
    if (length(pred_files) == 0) {
      message("⚠️ No prediccions trobades per ", esm, "-", alg)
      next
    }
    
    pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
    sel_idx <- which(year(pred_dates) == year_tgt & month(pred_dates) %in% months_br)
    if (length(sel_idx) == 0) {
      message("⚠️ No fitxers per breeding 2015.")
      next
    }
    
    # --- Stack & calc ANY ---
    r_stack <- stack(pred_files[sel_idx])
    r_any   <- calc(r_stack, fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))
    r_any_roi <- try(crop(r_any, bbox_sp), silent = TRUE)
    if (inherits(r_any_roi, "try-error") || is.null(r_any_roi)) next
    
    # --- Ajustar màscara de resolució ---
    env_mask_tmp <- try(env_mask, silent = TRUE)
    if (inherits(env_mask_tmp, "try-error") || is.null(env_mask_tmp)) next
    
    # Només si tots dos són raster
    if (inherits(env_mask_tmp, "Raster") && inherits(r_any_roi, "Raster")) {
      same_res <- try(compareRaster(env_mask_tmp, r_any_roi, crs=TRUE, res=TRUE, stopIfNotEqual=FALSE), silent=TRUE)
      if (inherits(same_res, "try-error") || is.logical(same_res)) {
        env_mask_tmp <- try(projectRaster(env_mask_tmp, r_any_roi, method="ngb"), silent=TRUE)
      }
    } else {
      next
    }
    
    # --- Ajustar màscara de resolució ---
    env_mask_tmp <- env_mask
    if (!compareRaster(env_mask_tmp, r_any_roi, crs=TRUE, res=TRUE, stopIfNotEqual=FALSE)) {
      env_mask_tmp <- projectRaster(env_mask_tmp, r_any_roi, method="ngb")
    }
    
    # --- Classificar punts ---
    vals_pred <- raster::extract(r_any_roi, st_coordinates(pts_sf))
    vals_mask <- raster::extract(env_mask_tmp,  st_coordinates(pts_sf))
    pts_sf$cat <- "outside"
    pts_sf$cat[!is.na(vals_pred) & vals_pred == 1] <- "inside"
    pts_sf$cat[is.na(vals_mask)] <- "coastline_gap"
    
    # --- DataFrames per ggplot ---
    r_df   <- as.data.frame(r_any_roi, xy = TRUE); colnames(r_df) <- c("x","y","pred")
    env_df <- as.data.frame(env_mask_tmp,   xy = TRUE); colnames(env_df) <- c("x","y","mask")
    world  <- ne_countries(scale="medium", returnclass="sf")
    
    # --- MAPA ---
    g <- ggplot() +
      geom_raster(data = subset(env_df, !is.na(mask)), aes(x = x, y = y), fill = "grey90") +
      geom_raster(data = subset(r_df, pred == 1), aes(x = x, y = y), fill = "#1b9e77", alpha = 0.65) +
      geom_sf(data = world, fill = "grey70", color = "white", linewidth = 0.2) +
      geom_sf(data = subset(pts_sf, cat == "inside"),        color = "blue3",   size = 0.3, alpha = 0.4) +
      geom_sf(data = subset(pts_sf, cat == "coastline_gap"), color = "orange",  size = 0.3, alpha = 0.4) +
      geom_sf(data = subset(pts_sf, cat == "outside"),       color = "red3",    size = 0.3, alpha = 0.4) +
      coord_sf(xlim = bbox_lon, ylim = bbox_lat, expand = FALSE) +
      labs(
        title = "Breeding season 2015 — Argentina",
        subtitle = paste(esm, "|", alg, "|", scenario, "— prediccions (≥1 mes=1) vs tracking"),
        x = "Longitud", y = "Latitud"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title   = element_text(face = "bold", hjust = 0.5),
        plot.subtitle= element_text(hjust = 0.5),
        panel.grid   = element_blank()
      )
    
    maps_list[[paste(esm, alg, sep = "_")]] <- g
  }
}





sum_binaries_for <- function(reg_sf, season_tag, bin_index) {
  stopifnot(season_tag %in% names(season_months))
  months_ok <- season_months[[season_tag]]
  
  idx <- bin_index %>%
    mutate(date = ymd(str_extract(basename(file), "\\d{8}")),
           year = year(date),
           month = month(date)) %>%
    filter(!is.na(date),
           year %in% val_years,
           month %in% months_ok)
  
  if (nrow(idx) == 0) return(NULL)
  
  # Intentem agafar un raster vàlid per fer servir de plantilla
  first_valid <- NULL
  for (ff in idx$file) {
    r_try <- try(raster(ff), silent = TRUE)
    if (!inherits(r_try, "try-error")) { first_valid <- r_try; break }
  }
  if (is.null(first_valid)) {
    message("⚠️ Cap raster vàlid per crear plantilla.")
    return(NULL)
  }
  
  if (is.na(crs(first_valid))) crs(first_valid) <- CRS("+proj=longlat +datum=WGS84")
  if (st_crs(reg_sf)$epsg %in% c(NA, 0)) reg_sf <- st_set_crs(reg_sf, 4326)
  
  reg2 <- if (st_crs(reg_sf) != st_crs(first_valid)) st_transform(reg_sf, st_crs(first_valid)) else reg_sf
  reg_sp <- as(reg2, "Spatial")
  
  tmpl_reg <- try(mask(crop(first_valid, reg_sp), reg_sp), silent = TRUE)
  if (inherits(tmpl_reg, "try-error") || is.null(tmpl_reg)) {
    message("⚠️ No s'ha pogut crear plantilla a la regió.")
    return(NULL)
  }
  
  sum_r <- setValues(tmpl_reg, 0)
  n_added <- 0L
  
  for (f in idx$file) {
    r <- try(raster(f), silent = TRUE)
    if (inherits(r, "try-error")) next
    
    # Retalla a la regió i comprova que no torni un objecte lògic
    r_reg <- try(mask(crop(r, reg_sp), reg_sp), silent = TRUE)
    if (inherits(r_reg, "try-error") || is.logical(r_reg) || is.null(r_reg)) next
    
    if (!compareRaster(r_reg, tmpl_reg, crs = TRUE, res = TRUE, extent = TRUE, stopIfNotEqual = FALSE)) {
      r_reg <- try(projectRaster(r_reg, tmpl_reg, method = "ngb"), silent = TRUE)
      if (inherits(r_reg, "try-error") || is.null(r_reg)) next
    }
    
    r_bin <- try(calc(r_reg, function(x) as.numeric(x >= 1)), silent = TRUE)
    if (inherits(r_bin, "try-error") || is.null(r_bin)) next
    
    sum_r <- try(sum_r + r_bin, silent = TRUE)
    if (inherits(sum_r, "try-error")) next
    
    n_added <- n_added + 1L
  }
  
  if (n_added == 0L) {
    message("⚠️ Cap capa afegida correctament per ", season_tag)
    return(NULL)
  }
  
  return(list(sum_r = sum_r, n_layers = n_added))
}

















library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(rnaturalearth)

# ===== CONFIGURACIÓ GENERAL =====
ESMs       <- c("MIROC-ES2L", "MRI-ESM2-0")
ALGORITHMS <- c("GBM", "RF", "SVM")
scenario   <- "ssp126"
years_val  <- 2015:2024
months_breeding    <- c(5,6,7,8,9,10)
months_nonbreeding <- c(11,12,1,2,3,4)

# Paths
base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
kernel_dir    <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/ker95_by_region"
tracking_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"

# Llegim shapefiles dels kernels (95%)
kernel_files <- list.files(kernel_dir, pattern = "\\.shp$", full.names = TRUE)
regions_sf <- do.call(rbind, lapply(kernel_files, st_read, quiet = TRUE))
regions_sf$Region <- basename(tools::file_path_sans_ext(kernel_files))

# Llegim dades de tracking (totes les regions)
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>% bind_rows() %>%
  mutate(date = suppressWarnings(dmy_hm(date))) %>%
  filter(year(date) %in% years_val) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

# ===== BUCLE PER REGIÓ =====
for (reg_name in regions_sf$Region) {
  message("\n=== REGIÓ: ", reg_name, " ===")
  
  reg <- regions_sf %>% filter(Region == reg_name)
  
  # 1️⃣ — Si el shapefile està en UTM, reprojectem
  if (max(abs(st_bbox(reg)[c("xmin","xmax","ymin","ymax")])) > 180) {
    message("   Reprojectant ", reg_name, " a WGS84...")
    reg <- st_transform(reg, 4326)
  }
  
  # Converteix a Spatial
  reg_sp <- as(reg, "Spatial")
  
  # Filtra punts de tracking dins la regió
  pts_reg <- tracking_all[st_intersects(tracking_all, reg, sparse = FALSE), ]
  if (nrow(pts_reg) == 0) {
    message("   ⚠️ Sense punts de tracking dins la regió.")
    next
  }
  
  for (season_tag in c("breeding", "nonbreeding")) {
    message("  - Temporada: ", season_tag)
    
    months_target <- if (season_tag == "breeding") months_breeding else months_nonbreeding
    
    # 2️⃣ — Agreguem binaris per cada ESM + algoritme dins el període
    rasters_all <- list()
    for (esm in ESMs) {
      for (alg in ALGORITHMS) {
        pred_dir <- file.path(base_pred_dir, esm, alg, scenario)
        pred_files <- list.files(pred_dir, pattern = "binary\\.tif$", full.names = TRUE)
        if (length(pred_files) == 0) next
        
        # Filtra per dates dins el període i mesos de la temporada
        pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
        sel <- which(year(pred_dates) %in% years_val & month(pred_dates) %in% months_target)
        if (length(sel) == 0) next
        
        pred_files <- pred_files[sel]
        message("     · ", esm, " — ", alg, " (", length(pred_files), " arxius)")
        
        # Carrega i suma (si almenys un mes és 1 → 1)
        stack_temp <- try(stack(pred_files), silent = TRUE)
        if (inherits(stack_temp, "try-error")) next
        r_any <- calc(stack_temp, fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))
        
        rasters_all[[paste(esm, alg, sep="_")]] <- r_any
      }
    }
    
    if (length(rasters_all) == 0) {
      message("   ⚠️ Cap raster binari carregat per ", reg_name, " ", season_tag)
      next
    }
    
    # 3️⃣ — Uniformitzar resolucions i extents
    # --- Comparació segura entre rasters ---
    same_grid <- function(r1, r2) {
      tryCatch({
        same_res <- all(res(r1) == res(r2))
        same_crs <- identical(crs(r1), crs(r2))
        same_ext <- all(abs(extent(r1)@xmin - extent(r2)@xmin) < 1e-6,
                        abs(extent(r1)@xmax - extent(r2)@xmax) < 1e-6,
                        abs(extent(r1)@ymin - extent(r2)@ymin) < 1e-6,
                        abs(extent(r1)@ymax - extent(r2)@ymax) < 1e-6)
        all(same_res, same_crs, same_ext)
      }, error = function(e) FALSE)
    }
    
    # --- Reprojecció robusta ---
    r_ref <- rasters_all[[1]]
    
    for (i in seq_along(rasters_all)) {
      r_cur <- rasters_all[[i]]
      
      if (!inherits(r_cur, "RasterLayer")) {
        message("⚠️ ", names(rasters_all)[i], " no és un raster vàlid.")
        next
      }
      
      if (!same_grid(r_cur, r_ref)) {
        message("Reprojectant ", names(rasters_all)[i], " per igualar resolució/extent/CRS...")
        r_proj <- try(projectRaster(r_cur, r_ref, method = "ngb"), silent = TRUE)
        if (!inherits(r_proj, "try-error")) {
          rasters_all[[i]] <- r_proj
        } else {
          message("⚠️ No s'ha pogut reprojectar ", names(rasters_all)[i], ", omès.")
          rasters_all[[i]] <- NULL
        }
      }
    }
    
    # Neteja
    rasters_all <- rasters_all[!sapply(rasters_all, is.null)]
    
    # 4️⃣ — Suma total: píxel = 1 si qualsevol model l’identifica com a adequat
    r_union <- calc(stack(rasters_all), fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))
    
    # 5️⃣ — Crop i màscara a la regió
    r_union_reg <- try({
      reg_expanded <- extent(reg_sp) + 0.1
      mask(crop(r_union, reg_expanded), reg_sp)
    }, silent = TRUE)
    if (inherits(r_union_reg, "try-error") || is.null(r_union_reg)) {
      message("   ⚠️ Error fent crop/mask a ", reg_name)
      next
    }
    
    # 6️⃣ — Classificar punts de tracking dins / fora / costa
    vals <- raster::extract(r_union_reg, st_coordinates(pts_reg))
    pts_reg$cat <- ifelse(vals == 1, "inside", "outside")
    vals_mask <- !is.na(r_union_reg[])
    vals_mask_pts <- raster::extract(r_union_reg, st_coordinates(pts_reg))
    pts_reg$cat[is.na(vals_mask_pts)] <- "coastline_gap"
    
    # 7️⃣ — Càlcul de sensibilitat excloent buits costaners
    pts_eval <- pts_reg %>% filter(cat != "coastline_gap")
    sensitivity <- sum(pts_eval$cat == "inside") / nrow(pts_eval)
    message("     ➤ Sensibilitat (sense buits costaners): ", round(sensitivity, 3))
    
    # 8️⃣ — Mapa
    r_df <- as.data.frame(r_union_reg, xy = TRUE)
    colnames(r_df) <- c("x", "y", "pred")
    r_df$pred <- as.factor(ifelse(r_df$pred == 1, "Suitable", NA))
    
    world <- ne_countries(scale = "medium", returnclass = "sf")
    
    gg <- ggplot() +
      geom_raster(data = subset(r_df, !is.na(pred)), aes(x = x, y = y), fill = "#1b9e77", alpha = 0.6) +
      geom_sf(data = world, fill = "grey80", color = "white", linewidth = 0.2) +
      geom_sf(data = subset(pts_reg, cat == "inside"), color = "blue", size = 0.5, alpha = 0.7) +
      geom_sf(data = subset(pts_reg, cat == "outside"), color = "red", size = 0.5, alpha = 0.7) +
      geom_sf(data = subset(pts_reg, cat == "coastline_gap"), color = "orange", size = 0.5, alpha = 0.7) +
      coord_sf(xlim = range(r_df$x, na.rm = TRUE), ylim = range(r_df$y, na.rm = TRUE), expand = FALSE) +
      labs(
        title = paste0(reg_name, " — ", season_tag, " (", years_val[1], "-", tail(years_val, 1), ")"),
        subtitle = paste("Predicció agregada (", scenario, ") | Sensibilitat:", round(sensitivity, 3)),
        x = "Longitud", y = "Latitud"
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "none")
    
    print(gg)
  }
}












library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(rnaturalearth)

# ===== CONFIGURACIÓ GENERAL =====
ESMs       <- c("MIROC-ES2L", "MRI-ESM2-0")
ALGORITHMS <- c("GBM", "RF", "SVM")
scenario   <- "ssp126"
years_val  <- 2015:2024
months_breeding    <- c(5,6,7,8,9,10)
months_nonbreeding <- c(11,12,1,2,3,4)

# Paths
base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
kernel_dir    <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/ker95_by_region"
tracking_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"

# Llegim shapefiles dels kernels (95%)
kernel_files <- list.files(kernel_dir, pattern = "\\.shp$", full.names = TRUE)
regions_sf <- do.call(rbind, lapply(kernel_files, st_read, quiet = TRUE))
regions_sf$Region <- basename(tools::file_path_sans_ext(kernel_files))

# Llegim dades de tracking
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>% bind_rows() %>%
  mutate(date = suppressWarnings(dmy_hm(date))) %>%
  filter(year(date) %in% years_val) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

# ===== BUCLE PER REGIÓ =====
for (reg_name in regions_sf$Region) {
  message("\n=== REGIÓ: ", reg_name, " ===")
  
  reg <- regions_sf %>% filter(Region == reg_name)
  
  # Si el shapefile està en UTM, reprojecta
  if (max(abs(st_bbox(reg)[c("xmin","xmax","ymin","ymax")])) > 180) {
    message("   Reprojectant ", reg_name, " a WGS84...")
    reg <- st_transform(reg, 4326)
  }
  
  # Filtra punts de tracking dins la regió
  pts_reg <- tracking_all[st_intersects(tracking_all, reg, sparse = FALSE), ]
  if (nrow(pts_reg) == 0) {
    message("   ⚠️ Sense punts de tracking dins la regió.")
    next
  }
  
  for (season_tag in c("breeding", "nonbreeding")) {
    message("  - Temporada: ", season_tag)
    
    months_target <- if (season_tag == "breeding") months_breeding else months_nonbreeding
    
    # --- Agreguem binaris per cada ESM + algoritme dins el període ---
    rasters_all <- list()
    for (esm in ESMs) {
      for (alg in ALGORITHMS) {
        pred_dir <- file.path(base_pred_dir, esm, alg, scenario)
        pred_files <- list.files(pred_dir, pattern = "binary\\.tif$", full.names = TRUE)
        if (length(pred_files) == 0) next
        
        pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
        sel <- which(year(pred_dates) %in% years_val & month(pred_dates) %in% months_target)
        if (length(sel) == 0) next
        
        pred_files <- pred_files[sel]
        message("     · ", esm, " — ", alg, " (", length(pred_files), " arxius)")
        
        stack_temp <- try(stack(pred_files), silent = TRUE)
        if (inherits(stack_temp, "try-error")) next
        
        r_any <- calc(stack_temp, fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))
        rasters_all[[paste(esm, alg, sep="_")]] <- r_any
      }
    }
    
    if (length(rasters_all) == 0) next
    
    # --- Unificar resolucions i CRS ---
    same_grid <- function(r1, r2) {
      tryCatch({
        same_res <- all(res(r1) == res(r2))
        same_crs <- identical(crs(r1), crs(r2))
        all(same_res, same_crs)
      }, error = function(e) FALSE)
    }
    
    r_ref <- rasters_all[[1]]
    for (i in seq_along(rasters_all)) {
      r_cur <- rasters_all[[i]]
      if (!same_grid(r_cur, r_ref)) {
        r_proj <- try(projectRaster(r_cur, r_ref, method = "ngb"), silent = TRUE)
        if (!inherits(r_proj, "try-error")) rasters_all[[i]] <- r_proj
      }
    }
    
    # --- Suma total: píxel = 1 si qualsevol model l’identifica com a adequat ---
    r_union <- calc(stack(rasters_all), fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))
    
    # --- Classificar punts ---
    vals <- raster::extract(r_union, st_coordinates(pts_reg))
    pts_reg$cat <- ifelse(vals == 1, "inside", "outside")
    
    # --- Calcular sensibilitat
    sensitivity <- sum(pts_reg$cat == "inside") / nrow(pts_reg)
    message("     ➤ Sensibilitat (total): ", round(sensitivity, 3))
    
    # --- Bounding box per fer zoom segons el kernel ---
    bb <- st_bbox(reg)
    xlim <- c(bb["xmin"] - 2, bb["xmax"] + 2)
    ylim <- c(bb["ymin"] - 2, bb["ymax"] + 2)
    
    # --- Plot ---
    r_df <- as.data.frame(r_union, xy = TRUE)
    colnames(r_df) <- c("x", "y", "pred")
    r_df$pred <- as.factor(ifelse(r_df$pred == 1, "Suitable", NA))
    world <- ne_countries(scale = "medium", returnclass = "sf")
    
    gg <- ggplot() +
      geom_raster(data = subset(r_df, !is.na(pred)), aes(x = x, y = y), fill = "#1b9e77", alpha = 0.5) +
      geom_sf(data = world, fill = "grey80", color = "white", linewidth = 0.2) +
      geom_sf(data = reg, fill = NA, color = "black", linewidth = 0.7, linetype = "dashed") +
      geom_sf(data = subset(pts_reg, cat == "inside"), color = "blue", size = 0.6, alpha = 0.8) +
      geom_sf(data = subset(pts_reg, cat == "outside"), color = "red", size = 0.6, alpha = 0.8) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      labs(
        title = paste0(reg_name, " — ", season_tag, " (", years_val[1], "-", tail(years_val, 1), ")"),
        subtitle = paste("Predicció agregada (", scenario, ") | Sensibilitat:", round(sensitivity, 3)),
        x = "Longitud", y = "Latitud"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none"
      )
    
    print(gg)
  }
}







library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(rnaturalearth)

# ===== CONFIGURACIÓ GENERAL =====
ESMs       <- c("MIROC-ES2L", "MRI-ESM2-0")
ALGORITHMS <- c("GBM", "RF", "SVM")
scenario   <- "ssp126"
years_val  <- 2015:2024
months_breeding    <- c(5,6,7,8,9,10)
months_nonbreeding <- c(11,12,1,2,3,4)

# --- Flexibilitat del consens ---
#  threshold_fraction = 0.5  → majoria simple (≥50%)
#  threshold_fraction = 1.0  → tots els models han de coincidir
#  threshold_fraction = 1/length(ESMs)/length(ALGORITHMS) → almenys un model (union)
threshold_fraction <- 0.5   # canvia-ho fàcilment segons vulguis

# Paths
base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
kernel_dir    <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/ker95_by_region"
tracking_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"

# Llegim shapefiles dels kernels (95%)
kernel_files <- list.files(kernel_dir, pattern = "\\.shp$", full.names = TRUE)
regions_sf <- do.call(rbind, lapply(kernel_files, st_read, quiet = TRUE))
regions_sf$Region <- basename(tools::file_path_sans_ext(kernel_files))

# Llegim dades de tracking
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>% bind_rows() %>%
  mutate(date = suppressWarnings(dmy_hm(date))) %>%
  filter(year(date) %in% years_val) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

# ===== BUCLE PER REGIÓ =====
for (reg_name in regions_sf$Region) {
  message("\n=== REGIÓ: ", reg_name, " ===")
  
  reg <- regions_sf %>% filter(Region == reg_name)
  
  # Si el shapefile està en UTM, reprojecta
  if (max(abs(st_bbox(reg)[c("xmin","xmax","ymin","ymax")])) > 180) {
    message("   Reprojectant ", reg_name, " a WGS84...")
    reg <- st_transform(reg, 4326)
  }
  
  # Filtra punts de tracking dins la regió
  pts_reg <- tracking_all[st_intersects(tracking_all, reg, sparse = FALSE), ]
  if (nrow(pts_reg) == 0) {
    message("   ⚠️ Sense punts de tracking dins la regió.")
    next
  }
  
  for (season_tag in c("breeding", "nonbreeding")) {
    message("  - Temporada: ", season_tag)
    
    months_target <- if (season_tag == "breeding") months_breeding else months_nonbreeding
    
    # --- Agreguem binaris per cada ESM + algoritme dins el període ---
    rasters_all <- list()
    for (esm in ESMs) {
      for (alg in ALGORITHMS) {
        pred_dir <- file.path(base_pred_dir, esm, alg, scenario)
        pred_files <- list.files(pred_dir, pattern = "binary\\.tif$", full.names = TRUE)
        if (length(pred_files) == 0) next
        
        pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
        sel <- which(year(pred_dates) %in% years_val & month(pred_dates) %in% months_target)
        if (length(sel) == 0) next
        
        pred_files <- pred_files[sel]
        message("     · ", esm, " — ", alg, " (", length(pred_files), " arxius)")
        
        stack_temp <- try(stack(pred_files), silent = TRUE)
        if (inherits(stack_temp, "try-error")) next
        
        # Hàbitat = adequat si almenys un mes és 1
        r_any <- calc(stack_temp, fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))
        rasters_all[[paste(esm, alg, sep="_")]] <- r_any
      }
    }
    
    if (length(rasters_all) == 0) next
    
    # --- Uniformitzar resolucions i CRS ---
    same_grid <- function(r1, r2) {
      tryCatch({
        same_res <- all(res(r1) == res(r2))
        same_crs <- identical(crs(r1), crs(r2))
        all(same_res, same_crs)
      }, error = function(e) FALSE)
    }
    
    r_ref <- rasters_all[[1]]
    for (i in seq_along(rasters_all)) {
      r_cur <- rasters_all[[i]]
      if (!same_grid(r_cur, r_ref)) {
        r_proj <- try(projectRaster(r_cur, r_ref, method = "ngb"), silent = TRUE)
        if (!inherits(r_proj, "try-error")) rasters_all[[i]] <- r_proj
      }
    }
    
    # --- Agregació flexible segons llindar ---
    n_models <- length(rasters_all)
    r_agregat <- calc(stack(rasters_all), fun = function(x) {
      frac <- sum(x == 1, na.rm = TRUE) / sum(!is.na(x))
      as.numeric(frac >= threshold_fraction)
    })
    
    # --- Classificar punts ---
    vals <- raster::extract(r_agregat, st_coordinates(pts_reg))
    pts_reg$cat <- ifelse(vals == 1, "inside", "outside")
    
    # --- Calcular sensibilitat
    sensitivity <- sum(pts_reg$cat == "inside") / nrow(pts_reg)
    message("     ➤ Sensibilitat (llindar=", threshold_fraction, "): ", round(sensitivity, 3))
    
    # --- Bounding box per fer zoom segons el kernel ---
    bb <- st_bbox(reg)
    xlim <- c(bb["xmin"] - 2, bb["xmax"] + 2)
    ylim <- c(bb["ymin"] - 2, bb["ymax"] + 2)
    
    # --- Plot ---
    r_df <- as.data.frame(r_agregat, xy = TRUE)
    colnames(r_df) <- c("x", "y", "pred")
    r_df$pred <- as.factor(ifelse(r_df$pred == 1, "Suitable", NA))
    world <- ne_countries(scale = "medium", returnclass = "sf")
    
    gg <- ggplot() +
      geom_raster(data = subset(r_df, !is.na(pred)), aes(x = x, y = y), fill = "#1b9e77", alpha = 0.5) +
      geom_sf(data = world, fill = "grey80", color = "white", linewidth = 0.2) +
      geom_sf(data = reg, fill = NA, color = "black", linewidth = 0.7, linetype = "dashed") +
      geom_sf(data = subset(pts_reg, cat == "inside"), color = "blue", size = 0.6, alpha = 0.8) +
      geom_sf(data = subset(pts_reg, cat == "outside"), color = "red", size = 0.6, alpha = 0.8) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      labs(
        title = paste0(reg_name, " — ", season_tag, " (", years_val[1], "-", tail(years_val, 1), ")"),
        subtitle = paste0("Consens ≥", round(threshold_fraction*100), "% | Sensibilitat: ", round(sensitivity, 3)),
        x = "Longitud", y = "Latitud"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none"
      )
    
    print(gg)
  }
}




library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(rnaturalearth)

# ===== CONFIGURACIÓ =====
ESMs       <- c("MIROC-ES2L", "MRI-ESM2-0")
# ESMs       <- c("MRI-ESM2-0")
ALGORITHMS <- c("GBM", "RF", "SVM")
scenario   <- "ssp126"
years_val  <- 2015:2024
months_breeding    <- c(5,6,7,8,9,10)
months_nonbreeding <- c(11,12,1,2,3,4)

# ===== DIRECTORIS =====
base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
kernel_dir    <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/ker95_by_region"
tracking_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"
example_res_layer <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/03predict/01standardized/MIROC-ES2L/GBM/historical/18500116_MIROC-ES2L_GBM_prediction.tif"

# ===== LLEGIR KERNELS =====
kernel_files <- list.files(kernel_dir, pattern = "\\.shp$", full.names = TRUE)
regions_sf <- do.call(rbind, lapply(kernel_files, st_read, quiet = TRUE))
regions_sf$Region <- basename(tools::file_path_sans_ext(kernel_files))

# Correcció noms intercanviats
regions_sf$Region <- gsub("South_Africa", "tmp", regions_sf$Region)
regions_sf$Region <- gsub("NZ", "South_Africa", regions_sf$Region)
regions_sf$Region <- gsub("tmp", "NZ", regions_sf$Region)

# ===== LLEGIR TRACKING =====
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>% bind_rows() %>%
  mutate(date = suppressWarnings(dmy_hm(date))) %>%
  filter(year(date) %in% years_val) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

# ===== LAYER DE RESOLUCIÓ =====
res_layer <- raster(example_res_layer)
res_layer_mask <- !is.na(res_layer) * 1

# ===== BUCLE PRINCIPAL =====
for (reg_name in regions_sf$Region) {
  cat("\n=== REGIÓ:", reg_name, "===\n")
  reg <- regions_sf %>% filter(Region == reg_name)
  
  # Reprojectar si cal
  if (max(abs(st_bbox(reg)[c("xmin","xmax","ymin","ymax")])) > 180) {
    reg <- st_transform(reg, 4326)
  }
  
  # 🧩 Eliminar parts que creuen l'antimeridià
  reg <- st_make_valid(reg)
  lon_mean <- mean(st_coordinates(reg)[,1], na.rm = TRUE)
  if (lon_mean > 100) {
    reg <- suppressWarnings(st_crop(reg, xmin = 0, xmax = 180, ymin = -90, ymax = 90))
  } else if (lon_mean < -100) {
    reg <- suppressWarnings(st_crop(reg, xmin = -180, xmax = 0, ymin = -90, ymax = 90))
  }
  reg <- st_make_valid(reg)
  reg_sp <- as(reg, "Spatial")
  
  # Punts dins la regió
  pts_reg <- tracking_all[st_intersects(tracking_all, reg, sparse = FALSE), ]
  if (nrow(pts_reg) == 0) {
    message("   ⚠️ Sense punts de tracking dins la regió.")
    next
  }
  
  for (season_tag in c("breeding", "nonbreeding")) {
    message("  - Temporada:", season_tag)
    months_target <- if (season_tag == "breeding") months_breeding else months_nonbreeding
    
    # --- Agregació binaris ---
    rasters_all <- list()
    for (esm in ESMs) {
      for (alg in ALGORITHMS) {
        pred_dir <- file.path(base_pred_dir, esm, alg, scenario)
        pred_files <- list.files(pred_dir, pattern = "binary\\.tif$", full.names = TRUE)
        if (length(pred_files) == 0) next
        
        pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
        sel <- which(year(pred_dates) %in% years_val & month(pred_dates) %in% months_target)
        if (length(sel) == 0) next
        
        stack_temp <- try(stack(pred_files[sel]), silent = TRUE)
        if (inherits(stack_temp, "try-error")) next
        r_any <- calc(stack_temp, fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))
        rasters_all[[paste(esm, alg, sep="_")]] <- r_any
      }
    }
    
    if (length(rasters_all) == 0) next
    
    # Uniformitzar resolucions
    r_ref <- rasters_all[[1]]
    same_grid <- function(r1, r2) {
      tryCatch({
        all(res(r1) == res(r2)) && identical(crs(r1), crs(r2))
      }, error = function(e) FALSE)
    }
    for (i in seq_along(rasters_all)) {
      r_cur <- rasters_all[[i]]
      if (!same_grid(r_cur, r_ref)) {
        r_proj <- try(projectRaster(r_cur, r_ref, method = "ngb"), silent = TRUE)
        if (!inherits(r_proj, "try-error")) rasters_all[[i]] <- r_proj
      }
    }
    
    # --- Consens total (≥1 model = 1) ---
    r_union <- calc(stack(rasters_all), fun = function(x) as.numeric(any(x == 1, na.rm = TRUE)))
    
    # --- Ampliar l'extensió per a la capa de resolució (donar "cancha") ---
    ext_big <- extent(reg_sp) + 20  # marge en graus
    env_mask_wide <- crop(res_layer_mask, ext_big)
    
    # --- Classificar punts ---
    vals_pred <- raster::extract(r_union, st_coordinates(pts_reg))
    vals_mask <- raster::extract(env_mask_wide, st_coordinates(pts_reg))
    
    pts_reg$cat <- "outside"
    pts_reg$cat[vals_pred == 1] <- "inside"
    pts_reg$cat[is.na(vals_mask) | vals_mask == 0] <- "coastline_gap"
    
    pts_eval <- pts_reg %>% filter(cat != "coastline_gap")
    sensitivity <- sum(pts_eval$cat == "inside") / nrow(pts_eval)
    message("     ➤ Sensibilitat (sense buits): ", round(sensitivity, 3))
    
    # --- Dades per al mapa ---
    r_df <- as.data.frame(r_union, xy = TRUE)
    colnames(r_df) <- c("x", "y", "pred")
    r_df$pred <- as.factor(ifelse(r_df$pred == 1, "Suitable", NA))
    env_df <- as.data.frame(env_mask_wide, xy = TRUE)
    colnames(env_df) <- c("x", "y", "mask")
    world <- ne_countries(scale = "medium", returnclass = "sf")
    
    bb <- st_bbox(reg)
    xlim <- c(bb["xmin"] - 5, bb["xmax"] + 5)
    ylim <- c(bb["ymin"] - 5, bb["ymax"] + 5)
    
    gg <- ggplot() +
      geom_raster(data = subset(env_df, mask == 1), aes(x, y), fill = "lightblue", alpha=0.5) +
      geom_raster(data = subset(r_df, !is.na(pred)), aes(x = x, y = y), fill = "#1b9e77", alpha = 0.65) +
      geom_sf(data = world, fill = "grey70", color = "white", linewidth = 0.2) +
      geom_sf(data = reg, fill = NA, color = "black", linetype = "dotted", linewidth = 0.5) +
      geom_sf(data = subset(pts_reg, cat == "inside"), color = "blue", size = 0.3, alpha = 0.4) +
      geom_sf(data = subset(pts_reg, cat == "outside"), color = "red", size = 0.3, alpha = 0.4) +
      geom_sf(data = subset(pts_reg, cat == "coastline_gap"), color = "orange", size = 0.3, alpha = 0.4) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      labs(
        title = paste0(reg_name, " — ", season_tag, " (", years_val[1], "-", tail(years_val, 1), ")"),
        subtitle = paste("Consens total (≥1 model) | Sensibilitat:", round(sensitivity, 2)),
        x = "Longitud", y = "Latitud"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
      )
    
    print(gg)
  }
}












library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(patchwork)

# ===== CONFIG =====
ESMs       <- c("MIROC-ES2L", "MRI-ESM2-0")
ALGORITHMS <- c("GBM", "RF", "SVM")
scenario   <- "ssp126"
years_val  <- 2015:2024
months_breeding    <- c(5,6,7,8,9,10)
months_nonbreeding <- c(11,12,1,2,3,4)

# ===== PATHS =====
base_pred_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/01standardizedBinaryPred"
kernel_dir    <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/04binaryPredictions/ker95_by_region"
tracking_dir  <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00input/01trackingDataValidation"
example_res_layer <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/03predict/01standardized/MIROC-ES2L/GBM/historical/18500116_MIROC-ES2L_GBM_prediction.tif"

# ===== KERNELS =====
kernel_files <- list.files(kernel_dir, pattern = "\\.shp$", full.names = TRUE)
regions_raw <- lapply(kernel_files, function(f){
  sf <- st_read(f, quiet = TRUE)
  sf$RegionFile <- tools::file_path_sans_ext(basename(f))
  sf
})
regions_sf <- do.call(rbind, regions_raw)
regions_sf$Region <- regions_sf$RegionFile %>%
  gsub("^ver95_|^ker95_", "", .)

# intercanvi SA ↔︎ NZ
regions_sf$Region <- gsub("^South_Africa$", "tmp", regions_sf$Region)
regions_sf$Region <- gsub("^New_Zealand$", "South_Africa", regions_sf$Region)
regions_sf$Region <- gsub("^tmp$", "New_Zealand", regions_sf$Region)

region_order <- c("Argentina","South_Africa","Australia","New_Zealand")
regions_sf <- regions_sf %>% filter(Region %in% region_order)
regions_sf$Region <- factor(regions_sf$Region, levels = region_order)

# ===== TRACKS =====
track_files <- list.files(tracking_dir, pattern = "fit_mpm_.*\\.csv$", full.names = TRUE)
tracking_all <- lapply(track_files, read.csv) %>%
  bind_rows() %>%
  mutate(date = suppressWarnings(dmy_hm(date))) %>%
  filter(year(date) %in% years_val) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

# ===== RESOLUTION MASK =====
res_layer <- raster(example_res_layer)
res_mask  <- !is.na(res_layer) * 1

# ===== FUNCIONS =====
same_grid <- function(r1, r2){
  tryCatch({
    all(res(r1) == res(r2)) && identical(crs(r1), crs(r2))
  }, error = function(e) FALSE)
}

make_map <- function(reg_name, season_tag){
  reg <- regions_sf %>% filter(Region == reg_name)
  if (nrow(reg) == 0) return(ggplot() + theme_void())
  
  if (max(abs(st_bbox(reg)[c("xmin","xmax","ymin","ymax")])) > 180) reg <- st_transform(reg, 4326)
  reg <- st_make_valid(reg)
  lon_mean <- mean(st_coordinates(reg)[,1], na.rm = TRUE)
  if (is.finite(lon_mean)) {
    if (lon_mean > 100) reg <- suppressWarnings(st_crop(reg, xmin = 0, xmax = 180, ymin = -90, ymax = 90))
    if (lon_mean < -100) reg <- suppressWarnings(st_crop(reg, xmin = -180, xmax = 0, ymin = -90, ymax = 90))
  }
  
  pts_reg_all <- tracking_all[st_intersects(tracking_all, reg, sparse = FALSE), ]
  months_target <- if (season_tag == "breeding") months_breeding else months_nonbreeding
  pts_reg <- pts_reg_all %>% filter(month(date) %in% months_target)
  if (nrow(pts_reg) == 0) {
    world <- ne_countries(scale="medium", returnclass="sf")
    bb <- st_bbox(reg); xlim <- c(bb["xmin"]-5, bb["xmax"]+5); ylim <- c(bb["ymin"]-5, bb["ymax"]+5)
    return(ggplot() +
             geom_sf(data=world, fill="grey80", color="white", linewidth=0.2) +
             geom_sf(data=reg, fill=NA, color="black", linetype="dotted", linewidth=0.5) +
             coord_sf(xlim=xlim, ylim=ylim, expand=FALSE) +
             theme_minimal() +
             theme(panel.border=element_rect(colour="black", fill=NA, linewidth=0.6),
                   axis.title=element_blank()))
  }
  
  # Agregació binària (consens total ≥1 model)
  rasters_all <- list()
  for (esm in ESMs) {
    for (alg in ALGORITHMS) {
      pred_dir <- file.path(base_pred_dir, esm, alg, scenario)
      pred_files <- list.files(pred_dir, pattern="binary\\.tif$", full.names=TRUE)
      if (!length(pred_files)) next
      pred_dates <- ymd(str_extract(basename(pred_files), "\\d{8}"))
      sel <- which(year(pred_dates) %in% years_val & month(pred_dates) %in% months_target)
      if (!length(sel)) next
      
      stk <- try(stack(pred_files[sel]), silent=TRUE)
      if (inherits(stk, "try-error")) next
      r_any <- calc(stk, fun=function(x) as.numeric(any(x==1, na.rm=TRUE)))
      rasters_all[[paste(esm, alg, sep="_")]] <- r_any
    }
  }
  if (!length(rasters_all)) return(ggplot() + theme_void())
  
  r_ref <- rasters_all[[1]]
  for (i in seq_along(rasters_all)) {
    r_cur <- rasters_all[[i]]
    if (!same_grid(r_cur, r_ref)) {
      r_proj <- try(projectRaster(r_cur, r_ref, method="ngb"), silent=TRUE)
      if (!inherits(r_proj, "try-error")) rasters_all[[i]] <- r_proj
    }
  }
  
  r_union <- calc(stack(rasters_all), fun=function(x) as.numeric(any(x==1, na.rm=TRUE)))
  ext_big <- extent(as(reg,"Spatial")) + 20
  env_mask_wide <- crop(res_mask, ext_big)
  
  vals_pred <- raster::extract(r_union, st_coordinates(pts_reg))
  vals_mask <- raster::extract(env_mask_wide, st_coordinates(pts_reg))
  pts_reg$cat <- "outside"
  pts_reg$cat[vals_pred == 1] <- "inside"
  pts_reg$cat[is.na(vals_mask) | vals_mask == 0] <- "coastline_gap"
  
  pts_eval <- pts_reg %>% filter(cat != "coastline_gap")
  sensitivity <- if (nrow(pts_eval)>0) sum(pts_eval$cat=="inside")/nrow(pts_eval) else NA
  
  r_df <- as.data.frame(r_union, xy=TRUE); colnames(r_df)<-c("x","y","pred")
  r_df$pred <- ifelse(r_df$pred==1,"Suitable",NA)
  env_df <- as.data.frame(env_mask_wide, xy=TRUE); colnames(env_df)<-c("x","y","mask")
  world <- ne_countries(scale="medium", returnclass="sf")
  bb <- st_bbox(reg); xlim <- c(bb["xmin"]-5, bb["xmax"]+5); ylim <- c(bb["ymin"]-5, bb["ymax"]+5)
  
  ggplot() +
    geom_raster(data=subset(env_df, mask==1), aes(x,y), fill="lightblue", alpha=0.5) +
    geom_raster(data=subset(r_df, pred=="Suitable"), aes(x,y), fill="#1b9e77", alpha=0.65) +
    geom_sf(data=world, fill="grey70", color="white", linewidth=0.2) +
    geom_sf(data=reg, fill=NA, color="black", linetype="dotted", linewidth=0.5) +
    geom_sf(data=subset(pts_reg, cat=="inside"), color="blue", size=0.05, alpha=0.3) +
    geom_sf(data=subset(pts_reg, cat=="outside"), color="red", size=0.05, alpha=0.3) +
    geom_sf(data=subset(pts_reg, cat=="coastline_gap"), color="orange", size=0.1, alpha=0.3) +
    coord_sf(xlim=xlim, ylim=ylim, expand=FALSE) +
    theme_minimal(base_size=11) +
    theme(
      axis.title=element_blank(),
      panel.border=element_rect(colour="black", fill=NA, linewidth=0.6),
      plot.margin=margin(2,2,2,2)
    )
}

# ===== MAPS 2×4 =====
regs <- c("Argentina","South_Africa","Australia","New_Zealand")
p_list <- list(
  make_map("Argentina","breeding"),
  make_map("South_Africa","breeding"),
  make_map("Australia","breeding"),
  make_map("New_Zealand","breeding"),
  make_map("Argentina","nonbreeding"),
  make_map("South_Africa","nonbreeding"),
  make_map("Australia","nonbreeding"),
  make_map("New_Zealand","nonbreeding")
)

# Fila superior = breeding, inferior = nonbreeding
# Només latitud al primer mapa (columna 1), longitud només als de sota (fila 2)
left_cols <- c(1,5)
bottom_row <- 5:8

mods <- lapply(seq_along(p_list), function(i){
  p_list[[i]] +
    theme(
      axis.text.y = if (i %in% left_cols) element_text(size=8) else element_blank(),
      axis.ticks.y = if (i %in% left_cols) element_line() else element_blank(),
      axis.text.x = if (i %in% bottom_row) element_text(size=8) else element_blank(),
      axis.ticks.x = if (i %in% bottom_row) element_line() else element_blank()
    )
})

final_plot <- (mods[[1]] | mods[[2]] | mods[[3]] | mods[[4]]) /
  (mods[[5]] | mods[[6]] | mods[[7]] | mods[[8]])

final_plot
