#==============================================================================
# FIT MODELS: RF, GBM, SVM WITH 5-FOLD CV
#==============================================================================

# 1. Load libraries
library(caret)
library(dplyr)
library(randomForest)
library(gbm)
library(pdp)
library(e1071)
library(tidyr)

# 1.0. MIROC
# # 1.1. RAW
# data_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/01tracking/00extract/00raw/MIROC-ES2L/SRW_extractRaw_150km_MIROC-ES2L.csv"
# df <- read.csv(data_path)

# data_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw/MIROC-ES2L/SRW_extractRaw_150km_MIROC-ES2L.csv"
# df <- read.csv(data_path)


# 1.2. STANDARDIZED
# data_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/01tracking/00extract/01standardized/MIROC-ES2L/SRW_extractStandardized_150km_MIROC-ES2L.csv"
# df <- read.csv(data_path)

data_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/01standardized/MIROC-ES2L/SRW_extractStandardized_150km_MIROC-ES2L.csv"
df <- read.csv(data_path)


# 2.0 MRI
# # 2.1. RAW
# data_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/01tracking/00extract/00raw/MRI-ESM2-0/SRW_extractRaw_150km_MRI-ESM2-0.csv"
# df <- read.csv(data_path)

# data_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/00raw/MRI-ESM2-0/SRW_extractRaw_150km_MRI-ESM2-0.csv"
# df <- read.csv(data_path)


# 2.2. STANDARDIZED
# data_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/01tracking/00extract/01standardized/MRI-ESM2-0/SRW_extractStandardized_150km_MRI-ESM2-0.csv"
# df <- read.csv(data_path)

# data_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/01presAbsExtract/01standardized/MRI-ESM2-0/SRW_extractStandardized_150km_MRI-ESM2-0.csv"
# df <- read.csv(data_path)


# Balance the dataset
presences <- df %>% filter(PresAbs == 1)
absences <- df %>% filter(PresAbs == 0)
set.seed(1234)
balanced_absences <- absences %>% sample_n(nrow(presences))
balanced_df <- bind_rows(presences, balanced_absences)

# Select predictors and response
predictor_columns <- names(df)[22:(ncol(df)-1)]  # Adjust as needed
response_column <- "PresAbs"

balanced_df <- balanced_df %>%
  dplyr::select(all_of(c(response_column, predictor_columns))) %>%
  na.omit()

# Convert response to factor
balanced_df$PresAbs <- ifelse(balanced_df$PresAbs == 1, "pres", "abs")
balanced_df$PresAbs <- factor(balanced_df$PresAbs, levels = c("pres", "abs"))

# Optional: drop variable if needed
# MIROC
# balanced_df_ <- balanced_df %>% dplyr::select(-SeaSurfaceTemperature)
# balanced_df_ <- balanced_df_ %>% dplyr::select(-TempDifferenceSurfaceDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-DistanceToIceEdge)
balanced_df_ <- balanced_df_ %>% dplyr::select(-SalinityDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-TemperatureSurface)
balanced_df_ <- balanced_df_ %>% dplyr::select(-SalinitySurface )
balanced_df_ <- balanced_df_ %>% dplyr::select(-WCurrentSurface)
balanced_df_ <- balanced_df_ %>% dplyr::select(-TemperatureDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-TemperatureGradientDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-SalinityGradientDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-SeaIceConcentration)
balanced_df_ <- balanced_df_ %>% dplyr::select(-UCurrentDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-VCurrentDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-WCurrentDepth)
colnames(balanced_df_)


# MRI
balanced_df_ <- balanced_df %>% dplyr::select(-DistanceToIceEdge)
balanced_df_ <- balanced_df_ %>% dplyr::select(-SalinityDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-TemperatureSurface)
balanced_df_ <- balanced_df_ %>% dplyr::select(-SalinitySurface )
balanced_df_ <- balanced_df_ %>% dplyr::select(-WCurrentSurface)
balanced_df_ <- balanced_df_ %>% dplyr::select(-TemperatureDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-TemperatureGradientDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-SalinityGradientDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-SeaIceConcentration)
balanced_df_ <- balanced_df_ %>% dplyr::select(-UCurrentDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-VCurrentDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-WCurrentDepth)
balanced_df_ <- balanced_df_ %>% dplyr::select(-MixedLayerDepth)




# 3. Define trainControl
tc <- trainControl(
  method = "repeatedcv",
  number = 5,             # 5-fold
  repeats = 1,           # Repeat 10 times
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = "final",
  verboseIter = TRUE      # Optional: see progress
)



# 4. Define tuning grids
rfGrid <- expand.grid(
  mtry = c(2, 4, 6),
  splitrule = "gini",
  min.node.size = c(5, 10, 15)
)

gbmGrid <- expand.grid(
  interaction.depth = c(1, 3, 5),
  n.trees = (1:5)*100,
  shrinkage = c(0.1, 0.01),
  n.minobsinnode = 10
)

svmGrid <- expand.grid(C = c(0.25, 0.5, 1, 2, 4, 8))




# Set reproducible global seed
set.seed(123)

# Calculate number of tuning combinations
n_rf     <- nrow(rfGrid)   # 3 mtry × 3 min.node.size = 9
n_gbm    <- nrow(gbmGrid)  # 3 depths × 5 trees × 2 shrinkage = 30
n_svm    <- nrow(svmGrid)  # 6 values of C

# Set number of resampling iterations
folds <- 5
repeats <- 10
n_iter <- folds * repeats  # = 50

# Function to generate seed list for a given number of tuning combinations
make_seed_list <- function(n_tune, n_iter) {
  seeds <- vector("list", n_iter + 1)
  for (i in 1:n_iter) {
    set.seed(111)
    seeds[[i]] <- sample.int(1000000, n_tune)
  }
  set.seed(111)
  seeds[[n_iter + 1]] <- sample.int(1000000, 1)  # Final model
  return(seeds)
}

# Generate seed lists
rf_seeds  <- make_seed_list(n_rf, n_iter)
gbm_seeds <- make_seed_list(n_gbm, n_iter)
svm_seeds <- make_seed_list(n_svm, n_iter)


# # 5. Create output folder
# MIROC
# outdir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/02modelling/00raw/MIROC-ES2L/01FitModels-noDepth/"
# if (!dir.exists(outdir)) dir.create(outdir)
# 
# # 5. Create output folder
# # MRI
# outdir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/02modelling/00raw/MRI-ESM2-0/01FitModels-noDepth/"
# if (!dir.exists(outdir)) dir.create(outdir)


# 5. Create output folder
# MIROC
# outdir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/02modelling/01standardized/MIROC-ES2L/01FitModels-noDepth/"
# if (!dir.exists(outdir)) dir.create(outdir)
# 
# # MRI
# outdir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output/02modelling/01standardized/MRI-ESM2-0/01FitModels-noDepth/"
# if (!dir.exists(outdir)) dir.create(outdir)


# 5. Create output folder
# MIROC
# outdir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/02modelling/00raw/MIROC-ES2L/01FitModels-noDepth/"
# if (!dir.exists(outdir)) dir.create(outdir)

# MIROC
outdir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/02modelling/01standardized/MIROC-ES2L/01FitModels-noDepth-noSST/"
if (!dir.exists(outdir)) dir.create(outdir)





# 5. Create output folder
# MRI
# outdir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/02modelling/00raw/MRI-ESM2-0/01FitModels-noDepth/"
# if (!dir.exists(outdir)) dir.create(outdir)


# 5. Create output folder
# MRI
outdir <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/02modelling/01standardized/MRI-ESM2-0/01FitModels-noDepth/"
if (!dir.exists(outdir)) dir.create(outdir)






#==============================================================================
# 6. Train and Save Models
#==============================================================================

# --- RANDOM FOREST ---
tc_rf <- trainControl(
  method = "repeatedcv",
  number = 5,             # 5-fold
  repeats = 10,           # Repeat 10 times
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = "final",
  verboseIter = TRUE,      # Optional: see progress
  seeds = rf_seeds
)


set.seed(123)
time_rf <- system.time({
  rf_model <- train(PresAbs ~ ., data = balanced_df_,
                    method = "ranger",
                    trControl = tc_rf,
                    tuneGrid = rfGrid,
                    importance = "impurity",
                    metric = "ROC")
})

outdirRF <- paste0(outdir, "/RF")
if (!dir.exists(outdirRF)) dir.create(outdirRF)
saveRDS(rf_model, file = file.path(outdirRF, "RF_model.rds"))
write.csv(rf_model$results, file = file.path(outdirRF, "RF_results.csv"), row.names = FALSE)
write.csv(varImp(rf_model)$importance, file = file.path(outdirRF, "RF_variableImportance.csv"))


# --
# Plot performance RF
# Step 1: Combine tuning parameters into a single label (adjust based on model type)
resample_data <- rf_model$resample %>%
  mutate(config = "rf_default")  # or build from mtry/splitrule/min.node.size if tuning

# Step 2: Convert data to long format for faceting
resample_long <- resample_data %>%
  pivot_longer(cols = c(ROC, Sens, Spec), names_to = "Metric", values_to = "Value")

# Step 3: Plot
plotPerformanceRF <- ggplot(resample_long, aes(x = config, y = Value)) +
  geom_violin(fill = "lightblue", alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  facet_wrap(~ Metric, scales = "fixed") +
  ylim(0.5,1)+
  labs(title = "Model Performance Across Repeats",
       x = "Model Configuration",
       y = "Metric Value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(outdirRF,"/performaceViolins_RF.pdf"), plotPerformanceRF)



# --- GBM ---
tc_gbm <- trainControl(
  method = "repeatedcv",
  number = 5,             # 5-fold
  repeats = 10,           # Repeat 10 times
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = "final",
  verboseIter = TRUE,      # Optional: see progress
  seeds = gbm_seeds
)


set.seed(123)
time_gbm <- system.time({
  gbm_model <- train(PresAbs ~ ., data = balanced_df_,
                     method = "gbm",
                     trControl = tc_gbm,
                     tuneGrid = gbmGrid,
                     verbose = FALSE,
                     metric = "ROC")
})

outdirGBM <- paste0(outdir, "/GBM2")
if (!dir.exists(outdirGBM)) dir.create(outdirGBM)
saveRDS(gbm_model, file = file.path(outdirGBM, "GBM_model.rds"))
write.csv(gbm_model$results, file = file.path(outdirGBM, "GBM_results.csv"), row.names = FALSE)
write.csv(varImp(gbm_model)$importance, file = file.path(outdirGBM, "GBM_variableImportance.csv"))

# --
# Plot performance GBM

# Step 1: Combine tuning parameters into a single label (adjust based on model type)
resample_data_GBM <- gbm_model$resample %>%
  mutate(config = "gbm_default")  # or build from mtry/splitrule/min.node.size if tuning

# Step 2: Convert data to long format for faceting
resample_long_GBM <- resample_data_GBM %>%
  pivot_longer(cols = c(ROC, Sens, Spec), names_to = "Metric", values_to = "Value")

# Step 3: Plot
plotPerformanceGBM <- ggplot(resample_long_GBM, aes(x = config, y = Value)) +
  geom_violin(fill = "lightblue", alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  facet_wrap(~ Metric, scales = "fixed") +
  ylim(0.5,1)+
  labs(title = "Model Performance Across Repeats",
       x = "Model Configuration",
       y = "Metric Value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(outdirGBM,"/performaceViolins_GBM.pdf"), plotPerformanceGBM)



# --- SVM ---
tc_svm <- trainControl(
  method = "repeatedcv",
  number = 5,             # 5-fold
  repeats = 10,           # Repeat 10 times
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = "final",
  verboseIter = TRUE,      # Optional: see progress
  seeds = svm_seeds
)

set.seed(123)
time_svm <- system.time({
  svm_model <- train(PresAbs ~ ., data = balanced_df_,
                     method = "svmRadialCost",
                     trControl = tc_svm,
                     tuneGrid = svmGrid,
                     metric = "ROC")
})

outdirSVM <- paste0(outdir, "/SVM")
if (!dir.exists(outdirSVM)) dir.create(outdirSVM)
saveRDS(svm_model, file = file.path(outdirSVM, "SVM_model.rds"))
write.csv(svm_model$results, file = file.path(outdirSVM, "SVM_results.csv"), row.names = FALSE)
write.csv(varImp(svm_model)$importance, file = file.path(outdirSVM, "SVM_variableImportance.csv"))

# --
# Plot performance SVM

# Step 1: Combine tuning parameters into a single label (adjust based on model type)
resample_data_SVM <- svm_model$resample %>%
  mutate(config = "svm_default")  # or build from mtry/splitrule/min.node.size if tuning

# Step 2: Convert data to long format for faceting
resample_long_SVM <- resample_data_SVM %>%
  pivot_longer(cols = c(ROC, Sens, Spec), names_to = "Metric", values_to = "Value")

# Step 3: Plot
plotPerformanceSVM <- ggplot(resample_long_SVM, aes(x = config, y = Value)) +
  geom_violin(fill = "lightblue", alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  facet_wrap(~ Metric, scales = "fixed") +
  ylim(0.5,1)+
  labs(title = "Model Performance Across Repeats",
       x = "Model Configuration",
       y = "Metric Value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(outdirSVM,"/performaceViolins_SVM.pdf"), plotPerformanceSVM)






#==============================================================================
# 7. Print training time
#==============================================================================

print(paste("RF time:", time_rf["elapsed"]))
print(paste("GBM time:", time_gbm["elapsed"]))
print(paste("SVM time:", time_svm["elapsed"]))



s <- raster::stack("/Volumes/23327-15324-35628/2024_SouthernRightWhaleLogbooks/StackRestoration/02envStacks_Standardized_0s/MIROC-ES2L/historical/monthly_stack_stand_MIROC-ES2L_historical_18500116.grd")
r <- raster::raster("/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/03predict/01standardized/MIROC-ES2L/GBM/historical/18500116_MIROC-ES2L_GBM_prediction.tif")
plot(r)

predictors <- c("Bathymetry",
                "SeaFloorSlope", "DistanceToShore", 
                "TemperatureGradientSurface", "TempDifferenceSurfaceDepth",
                "SeaSurfaceSalinity", "SalinityGradientSurface", 
                "SeaSurfaceHeight", "UCurrentSurface", "VCurrentSurface",
                "EddyKineticEnergy", 
                "UWindStress", "VWindStress", "WindStressMagnitude",
                "WindStressDirection"
)



model <- gbm_model
available_vars <- names(s)

common_vars <- predictors[predictors %in% available_vars]

if (length(common_vars) == 0) {
  message("No matching variables found in ", fname)
  next
}

s <- subset(s, common_vars)

# Predict probability of presence ("pres")
p_2 <- raster::predict(s, model, type = "prob", index = 1, progress = "text")
r <- raster::raster("/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/03predict/01standardized/MIROC-ES2L/GBM/historical/18500116_MIROC-ES2L_GBM_prediction.tif")
r_cc <- raster::raster("/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz/InputOutput/00output_RestoredStacks/03predict/01standardized/MIROC-ES2L/GBM/ssp585/20990116_MIROC-ES2L_GBM_prediction.tif")

plot(r_cc)

library(raster)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(viridis)
library(patchwork)

# ---- 1) Convert all rasters to data frames ----
p_df   <- as.data.frame(rasterToPoints(p))
p2_df  <- as.data.frame(rasterToPoints(p_2))
r_df   <- as.data.frame(rasterToPoints(r))
colnames(p_df)  <- c("x", "y", "value")
colnames(p2_df) <- c("x", "y", "value")
colnames(r_df)  <- c("x", "y", "value")

# ---- 2) Keep only Southern Hemisphere ----
p_df  <- subset(p_df,  y <= 0)
p2_df <- subset(p2_df, y <= 0)
r_df  <- subset(r_df,  y <= 0)

# ---- 3) Load and crop world shapefile ----
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_make_valid(world)
south_bbox <- st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = st_crs(4326))
world_south <- st_crop(world, south_bbox)

# ---- 4) Determine a shared fill scale range ----
fill_limits <- range(c(p_df$value, p2_df$value, r_df$value), na.rm = TRUE)

# ---- 5) Create individual plots ----
p_plot <- ggplot() +
  geom_raster(data = p_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = world_south, color = "grey30", fill = NA, linewidth = 0.2) +
  scale_fill_viridis(name = "Suitability", option = "viridis", limits = fill_limits) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 0), expand = FALSE) +
  labs(title = "Predicted (p)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.title = element_blank())

p2_plot <- ggplot() +
  geom_raster(data = p2_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = world_south, color = "grey30", fill = NA, linewidth = 0.2) +
  scale_fill_viridis(name = "Suitability", option = "viridis", limits = fill_limits) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 0), expand = FALSE) +
  labs(title = "Predicted (p_2)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.title = element_blank())

r_plot <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = world_south, color = "grey30", fill = NA, linewidth = 0.2) +
  scale_fill_viridis(name = "Suitability", option = "viridis", limits = fill_limits) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 0), expand = FALSE) +
  labs(title = "Saved raster (r)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.title = element_blank())

# ---- 6) Combine side by side ----
(p_plot | p2_plot | r_plot)






library(raster)
library(ggplot2)
library(rnaturalearth)
library(sf)

# ---- 1) Binary reclassification (choose your threshold)
thresh <- 0.1
p_bin  <- calc(p,  function(x) ifelse(x >= thresh, 1, 0))
p2_bin <- calc(p_2,function(x) ifelse(x >= thresh, 1, 0))
r_bin  <- calc(r,  function(x) ifelse(x >= thresh, 1, 0))

# ---- 2) Encode model combinations as unique values
# Each model gets a binary weight: p=1, p2=2, r=4 → combinations up to 7
combo <- (p_bin * 1) + (p2_bin * 2) + (r_bin * 4)

# ---- 3) Convert to dataframe for ggplot
df <- as.data.frame(rasterToPoints(combo))
colnames(df) <- c("x", "y", "code")

# Only Southern Hemisphere
df <- subset(df, y <= 0)

# ---- 4) Decode combinations into labels
df$label <- factor(df$code,
                   levels = c(0, 1, 2, 3, 4, 5, 6, 7),
                   labels = c(
                     "Unsuitable in all",
                     "Only NO SST",
                     "Only NO TEMPDIFF AND SST",
                     "NO SST + NO TEMPDIFF AND SST",
                     "Only ALL VARS",
                     "NO SST + ALL VARS",
                     "NO TEMPDIFF + ALL VARS",
                     "All three"
                   )
)

library(rnaturalearth)
library(sf)
library(ggplot2)

# ---- 1) Load and fix the world shapefile ----
world <- ne_countries(scale = "medium", returnclass = "sf")

# Ensure geometries are valid (fix polygons that break ggplot)
world <- suppressWarnings(st_make_valid(world))

# ---- 2) Crop manually using bounding box geometry (avoids st_crop issues) ----
south_box <- st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = st_crs(4326)))
world_south <- st_intersection(world, south_box)

# ---- 3) Plot with your already-prepared df ----
ggplot(df) +
  geom_raster(aes(x = x, y = y, fill = label)) +
  geom_sf(data = world, color = "grey30", fill = NA, linewidth = 0.1) +
  scale_fill_manual(
    name = "Models predicting suitability",
    values = c(
      "Unsuitable in all" = "white",
      "Only NO SST" = "#fde725",
      "Only NO TEMPDIFF AND SST" = "#5ec962",
      "NO SST + NO TEMPDIFF AND SST" = "#21918c",
      "Only ALL VARS" = "#3b528b",
      "NO SST + ALL VARS" = "#440154",
      "NO TEMPDIFF + ALL VARS" = "#482878",
      "All three" = "black"
    )
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 0), expand = FALSE) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    legend.position = "right") +
  labs(title = "Combinations of models predicting suitable habitat")
