# ==============================================================================
# CREATE / UPDATE README.md AND PUSH TO GITHUB
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. PROJECT DIRECTORY
# ------------------------------------------------------------------------------

PROJECT_DIR <- "/Users/jazelouled-cheikhbonan/Dropbox/2024_SouthernRightWhaleLogbooks_MRuiz"

setwd(PROJECT_DIR)

cat("Working directory:\n", getwd(), "\n\n")


# ------------------------------------------------------------------------------
# 2. README CONTENT
# ------------------------------------------------------------------------------

readme <- c(
  
  "# Southern Right Whale Habitat Suitability through Time",
  
  "",
  
  "![R](https://img.shields.io/badge/R-4.x-blue)",
  "![Status](https://img.shields.io/badge/status-in%20development-orange)",
  
  "",
  
  "## Overview",
  
  "",
  
  "This repository contains the analytical workflow used to reconstruct and model long-term changes in habitat suitability for the Southern right whale (*Eubalaena australis*) across the Southern Hemisphere.",
  
  "",
  
  "The project combines historical whale occurrence data with environmental reconstructions from two Earth System Models (ESMs):",
  
  "",
  
  "- **MIROC-ES2L**",
  "- **MRI-ESM2-0**",
  
  "",
  
  "The workflow includes environmental stack construction, harmonisation of climate experiments, predictor selection, species distribution modelling, historical reconstruction and future climate projections.",
  
  "",
  
  "A major methodological component is the harmonisation of the CMIP6 `past1000` and `historical` experiments around the 1850 transition, where discontinuities in several environmental variables were detected.",
  
  "",
  
  "---",
  
  "",
  
  "## Workflow",
  
  "",
  
  "```text",
  "Historical whale occurrence data",
  "              │",
  "              ▼",
  "Environmental reconstruction",
  "              │",
  "              ▼",
  "ESM harmonisation",
  "              │",
  "              ▼",
  "Predictor selection",
  "              │",
  "              ▼",
  "       RF · GBM · SVM",
  "              │",
  "              ▼",
  "1700 ───────────────────── 2100",
  "       Habitat suitability",
  "```",
  
  "",
  
  "---",
  
  "",
  
  "## Repository structure",
  
  "",
  
  "```text",
  ".",
  "├── README.md",
  "├── scripts/",
  "│   ├── 01_build_environmental_stacks.R",
  "│   ├── 02_calculate_ESM_deltas.R",
  "│   ├── 03_harmonise_past1000_raw_rebuildDerived.R",
  "│   ├── 04_validate_harmonisation.R",
  "│   ├── 05_fit_final_SDM_models.R",
  "│   └── ...",
  "├── data/",
  "├── figures/",
  "└── outputs/",
  "```",
  
  "",
  
  "Large environmental rasters and intermediate products are not stored directly in this GitHub repository.",
  
  "",
  
  "---",
  
  "",
  
  "## 1. Environmental reconstruction",
  
  "",
  
  "Monthly environmental stacks were generated independently for:",
  
  "",
  
  "- MIROC-ES2L",
  "- MRI-ESM2-0",
  
  "",
  
  "and across the following CMIP6 experiments:",
  
  "",
  
  "- `past1000`",
  "- `historical`",
  "- `ssp126`",
  "- `ssp585`",
  
  "",
  
  "Environmental information includes temperature, salinity, currents, sea ice, sea-surface height, wind stress, bathymetry and several derived predictors.",
  
  "",
  
  "### Static variables",
  
  "",
  
  "The following predictors are treated as time invariant:",
  
  "",
  
  "- `Bathymetry`",
  "- `SeaFloorSlope`",
  "- `DistanceToShore`",
  
  "",
  
  "### Derived environmental variables",
  
  "",
  
  "Derived predictors include:",
  
  "",
  
  "```text",
  "EddyKineticEnergy =",
  "    (UCurrentSurface² + VCurrentSurface²) / 2",
  "",
  "TempDifferenceSurfaceDepth =",
  "    TemperatureSurface - TemperatureDepth",
  "",
  "WindStressMagnitude =",
  "    sqrt(UWindStress² + VWindStress²)",
  "",
  "WindStressDirection =",
  "    atan2(VWindStress, UWindStress)",
  "```",
  
  "",
  
  "Temperature and salinity gradients are derived spatially, while `DistanceToIceEdge` is calculated from the 15% sea-ice concentration threshold.",
  
  "",
  
  "---",
  
  "",
  
  "## 2. Harmonisation of `past1000` and `historical`",
  
  "",
  
  "Inspection of the environmental time series revealed discontinuities between the end of the `past1000` experiment and the beginning of `historical`:",
  
  "",
  
  "```text",
  "past1000        historical",
  "   1849    →       1850",
  "```",
  
  "",
  
  "These discontinuities differed between ESMs and were particularly evident for variables such as sea-surface temperature and sea-surface height.",
  
  "",
  
  "To minimise artificial experiment-boundary effects, monthly spatial correction fields are calculated independently for each ESM.",
  
  "",
  
  "### Reference periods",
  
  "",
  
  "```text",
  "past1000:     1840–1849",
  "historical:   1850–1859",
  "```",
  
  "",
  
  "For each ESM, variable, grid cell and calendar month:",
  
  "",
  
  "```text",
  "Delta(x, y, month) =",
  "    historical climatology",
  "    -",
  "    past1000 climatology",
  "```",
  
  "",
  
  "The correction is then applied to the full `past1000` experiment.",
  
  "",
  
  "---",
  
  "",
  
  "## 3. Data integrity",
  
  "",
  
  "> **Original RAW environmental stacks are treated as immutable source data.**",
  
  "",
  
  "The original RAW directory is never modified or overwritten.",
  
  "",
  
  "All harmonised products are generated as new files in separate versioned directories.",
  
  "",
  
  "```text",
  "01envStacks_Raw_0s/",
  "        │",
  "        │ READ ONLY",
  "        ▼",
  "02_ESM_harmonisation/",
  "        │",
  "        ▼",
  "03envStacks_Raw_Harmonised_v1/",
  "```",
  
  "",
  
  "Static variables remain unchanged during harmonisation.",
  
  "",
  
  "Dynamic base variables are corrected using monthly spatial deltas, while derived predictors are recalculated from the corrected base variables to preserve mathematical consistency.",
  
  "",
  
  "---",
  
  "",
  
  "## 4. Harmonisation validation",
  
  "",
  
  "Harmonised stacks are validated by comparing:",
  
  "",
  
  "- original RAW `past1000`",
  "- harmonised `past1000`",
  "- original RAW `historical`",
  
  "",
  
  "Validation includes:",
  
  "",
  
  "- raster geometry and layer consistency",
  "- RAW vs harmonised spatial comparisons",
  "- maps of `harmonised - RAW`",
  "- verification that static layers remain unchanged",
  "- verification of derived-variable equations",
  "- inspection of environmental continuity across the 1849–1850 transition",
  
  "",
  
  "---",
  
  "",
  
  "## 5. Predictor selection",
  
  "",
  
  "Predictor selection is performed after environmental harmonisation.",
  
  "",
  
  "Candidate predictors are assessed for collinearity and ecological relevance before fitting the final species distribution models.",
  
  "",
  
  "---",
  
  "",
  
  "## 6. Species distribution modelling",
  
  "",
  
  "Three modelling algorithms are used:",
  
  "",
  
  "- Random Forest (`ranger`)",
  "- Gradient Boosting Machine (`gbm`)",
  "- Support Vector Machine with radial kernel",
  
  "",
  
  "Models are trained independently for each ESM.",
  
  "",
  
  "```text",
  "MIROC-ES2L",
  "├── RF",
  "├── GBM",
  "└── SVM",
  "",
  "MRI-ESM2-0",
  "├── RF",
  "├── GBM",
  "└── SVM",
  "```",
  
  "",
  
  "Model training uses repeated cross-validation:",
  
  "",
  
  "```text",
  "5 folds × 10 repetitions",
  "```",
  
  "",
  
  "Model selection is based on ROC performance, and presence/absence samples are balanced prior to training.",
  
  "",
  
  "---",
  
  "",
  
  "## 7. Model interpretation",
  
  "",
  
  "Model interpretation includes:",
  
  "",
  
  "- variable importance",
  "- partial dependence plots",
  "- comparison among algorithms",
  "- comparison between ESMs",
  
  "",
  
  "Variable importance is normalised within each algorithm so that the total importance sums to 100%.",
  
  "",
  
  "---",
  
  "",
  
  "## 8. Historical and future projections",
  
  "",
  
  "Final models are projected across environmental conditions covering:",
  
  "",
  
  "```text",
  "past1000",
  "historical",
  "SSP1-2.6",
  "SSP5-8.5",
  "```",
  
  "",
  
  "This enables reconstruction of long-term changes in potentially suitable Southern right whale habitat from historical conditions to future climate scenarios.",
  
  "",
  
  "---",
  
  "",
  
  "## Main software",
  
  "",
  
  "Core R packages include:",
  
  "",
  
  "```text",
  "terra",
  "raster",
  "ncdf4",
  "caret",
  "ranger",
  "gbm",
  "e1071",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "sf",
  "```",
  
  "",
  
  "---",
  
  "",
  
  "## Reproducibility",
  
  "",
  
  "Analysis scripts are organised sequentially according to workflow order:",
  
  "",
  
  "```text",
  "01_ environmental stack construction",
  "02_ ESM transition diagnostics and delta calculation",
  "03_ past1000 harmonisation and reconstruction of derived variables",
  "04_ harmonisation validation",
  "05_ final SDM fitting",
  "...",
  "```",
  
  "",
  
  "Exact software versions, environmental data sources and final workflow documentation will be added before public release.",
  
  "",
  
  "---",
  
  "",
  
  "## Authors",
  
  "",
  
  "**Jazel Ouled-Cheikh Bonan**  ",
  "Universitat de València",
  
  "",
  
  "**Marc Ruiz-Sagalés**",
  
  "",
  
  "Additional contributors and affiliations will be added before release.",
  
  "",
  
  "---",
  
  "",
  
  "## Citation",
  
  "",
  
  "Citation information will be added upon publication.",
  
  "",
  
  "---",
  
  "",
  
  "## License",
  
  "",
  
  "License information will be added before public release."
  
)


# ------------------------------------------------------------------------------
# 3. WRITE README.md
# ------------------------------------------------------------------------------

README_FILE <- file.path(
  PROJECT_DIR,
  "README.md"
)

writeLines(
  readme,
  README_FILE,
  useBytes = TRUE
)

cat("\nREADME created:\n")
cat(README_FILE, "\n")


# ------------------------------------------------------------------------------
# 4. CHECK THAT THIS IS A GIT REPOSITORY
# ------------------------------------------------------------------------------

if (!dir.exists(file.path(PROJECT_DIR, ".git"))) {
  
  stop(
    paste0(
      "\nThis directory is not currently a Git repository:\n",
      PROJECT_DIR,
      "\n\nRun git init / connect it to GitHub first."
    )
  )
}


# ------------------------------------------------------------------------------
# 5. CHECK REMOTE
# ------------------------------------------------------------------------------

remotes <- system(
  "git remote",
  intern = TRUE
)

if (!"origin" %in% remotes) {
  
  stop(
    paste0(
      "\nNo GitHub remote called 'origin' is configured.\n",
      "README.md has been created locally, but nothing has been pushed."
    )
  )
}


cat("\nGit remote:\n")

print(
  system(
    "git remote -v",
    intern = TRUE
  )
)


# ------------------------------------------------------------------------------
# 6. GIT ADD
# ------------------------------------------------------------------------------

status_add <- system(
  "git add README.md"
)

if (status_add != 0) {
  stop("git add failed.")
}


# ------------------------------------------------------------------------------
# 7. GIT COMMIT
# ------------------------------------------------------------------------------

commit_message <- "Add project README"

status_commit <- system2(
  "git",
  args = c(
    "commit",
    "-m",
    shQuote(commit_message)
  )
)


# A non-zero commit status can simply mean there was nothing new to commit.
if (status_commit != 0) {
  
  cat(
    "\nNo new commit created. README may already be up to date.\n"
  )
}


# ------------------------------------------------------------------------------
# 8. DETECT CURRENT BRANCH
# ------------------------------------------------------------------------------

branch <- system(
  "git branch --show-current",
  intern = TRUE
)

branch <- trimws(
  branch[1]
)

if (
  length(branch) == 0 ||
  is.na(branch) ||
  branch == ""
) {
  
  stop(
    "Could not determine current Git branch."
  )
}


cat(
  "\nCurrent branch:",
  branch,
  "\n"
)


# ------------------------------------------------------------------------------
# 9. PUSH TO GITHUB
# ------------------------------------------------------------------------------

cat(
  "\nPushing README to GitHub...\n"
)


push_status <- system2(
  "git",
  args = c(
    "push",
    "origin",
    branch
  )
)


if (push_status != 0) {
  
  stop(
    paste0(
      "\nGit push failed.\n",
      "The README is still safely stored locally at:\n",
      README_FILE
    )
  )
}


# ------------------------------------------------------------------------------
# 10. DONE
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("README SUCCESSFULLY PUSHED TO GITHUB\n")
cat("============================================================\n")