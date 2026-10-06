#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Southern Right Whale CMIP6 preprocessing
# Ultimate end-to-end regridding script
#
# Consolidates the logic learned across the historical SRW shell scripts:
# - process MIROC-ES2L and MRI-ESM2-0
# - process past1000, historical, ssp126, ssp585
# - handle single-layer variables directly
# - handle 3-D variables by extracting:
#     * the shallowest available level
#     * the available level closest to 200 m
# - process siconc separately but through the same generic pipeline
# - remap to 1 degree using CDO remapnn
# - crop to Southern Hemisphere (-180,180,-90,0)
# - split into one NetCDF per timestep/month
# - skip existing outputs by default
# - write logs and a summary TSV
#
# IMPORTANT:
# - Reads RAW downloaded NetCDFs only.
# - Writes only under OUTPUT_ROOT.
# - Does not modify original RAW files.
#
# Requirements:
# - bash
# - CDO
# - grid_1deg.txt
# ==============================================================================


# ==============================================================================
# 0. USER CONFIGURATION
# ==============================================================================

INPUT_ROOT="/Volumes/tichodroma/2024_SouthernRightWhaleLogbooks/00rawData/UsedESM"
OUTPUT_ROOT="/Volumes/MyPassport/2024_SouthernRightWhaleLogbooks/01regriddedData"

GRID_FILE="grid_1deg.txt"
CROP_EXTENT="-180,180,-90,0"

ESMS=("MIROC-ES2L" "MRI-ESM2-0")
SCENARIOS=("past1000" "historical" "ssp126" "ssp585")

SINGLE_LAYER_VARIABLES=(
  "tos"
  "sos"
  "zos"
  "tauuo"
  "tauvo"
  "mlotst"
  "siconc"
)

MULTI_LAYER_VARIABLES=(
  "so"
  "thetao"
  "uo"
  "vo"
  "wo"
)

OVERWRITE=false
DRY_RUN=false


# ==============================================================================
# 1. INITIAL CHECKS
# ==============================================================================

command -v cdo >/dev/null 2>&1 || {
  echo "ERROR: CDO is not available on PATH."
  exit 1
}

[[ -f "$GRID_FILE" ]] || {
  echo "ERROR: Grid file not found: $GRID_FILE"
  exit 1
}

[[ -d "$INPUT_ROOT" ]] || {
  echo "ERROR: INPUT_ROOT does not exist: $INPUT_ROOT"
  exit 1
}

mkdir -p "$OUTPUT_ROOT"

LOG_DIR="$OUTPUT_ROOT/_regrid_logs"
mkdir -p "$LOG_DIR"

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$LOG_DIR/regrid_${RUN_STAMP}.log"
ERROR_FILE="$LOG_DIR/regrid_errors_${RUN_STAMP}.log"
SUMMARY_FILE="$LOG_DIR/regrid_summary_${RUN_STAMP}.tsv"

printf "status\tesm\tscenario\tvariable\tinput_file\tlevel\tdate\toutput_file\tmessage\n" > "$SUMMARY_FILE"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "============================================================"
echo "SRW CMIP6 REGRIDDING"
echo "Started: $(date)"
echo "Input : $INPUT_ROOT"
echo "Output: $OUTPUT_ROOT"
echo "Grid  : $GRID_FILE"
echo "Crop  : $CROP_EXTENT"
echo "Dry run: $DRY_RUN"
echo "Overwrite: $OVERWRITE"
echo "============================================================"


# ==============================================================================
# 2. HELPERS
# ==============================================================================

in_array() {
  local needle="$1"
  shift
  local x
  for x in "$@"; do
    [[ "$x" == "$needle" ]] && return 0
  done
  return 1
}

log_summary() {
  local status="$1"
  local esm="$2"
  local scenario="$3"
  local variable="$4"
  local input_file="$5"
  local level="$6"
  local date="$7"
  local output_file="$8"
  local message="$9"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$status" "$esm" "$scenario" "$variable" "$input_file" \
    "$level" "$date" "$output_file" "$message" >> "$SUMMARY_FILE"
}

run_cmd() {
  echo "+ $*"
  if [[ "$DRY_RUN" == false ]]; then
    "$@"
  fi
}

safe_rm() {
  if [[ "$DRY_RUN" == false ]]; then
    rm -f "$@"
  fi
}

extract_variable_name() {
  local input_file="$1"
  local variable

  variable="$(cdo -s showname "$input_file" 2>/dev/null | awk '{print $1}')"

  if [[ -z "$variable" ]]; then
    variable="$(basename "$input_file" | sed -n 's/.*\.Omon\.\([^.]*\)\..*/\1/p')"
  fi

  if [[ -z "$variable" ]]; then
    variable="$(basename "$input_file" | sed -n 's/.*\.SImon\.\([^.]*\)\..*/\1/p')"
  fi

  echo "$variable"
}

get_levels() {
  local input_file="$1"
  cdo -s showlevel "$input_file" 2>/dev/null \
    | tr ' ' '\n' \
    | sed '/^[[:space:]]*$/d' \
    | awk '$0 ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ { print $0 }'
}

choose_surface_level() {
  awk '
    BEGIN { found=0 }
    {
      x=$1+0
      a=(x<0 ? -x : x)
      if (!found || a<best_abs) {
        best=x
        best_abs=a
        found=1
      }
    }
    END {
      if (found) print best
    }
  '
}

choose_deep_level() {
  awk '
    BEGIN { found=0; target=200 }
    {
      x=$1+0
      d=x-target
      if (d<0) d=-d
      if (!found || d<best_d) {
        best=x
        best_d=d
        found=1
      }
    }
    END {
      if (found) print best
    }
  '
}

format_level_for_filename() {
  local level="$1"
  awk -v x="$level" 'BEGIN {
    if (x == int(x)) printf "%d", x;
    else printf "%g", x;
  }'
}

make_output_dir() {
  local esm="$1"
  local scenario="$2"
  local variable="$3"

  local outdir="$OUTPUT_ROOT/$esm/$scenario/1degree/$variable"
  mkdir -p "$outdir"
  echo "$outdir"
}

should_skip_existing() {
  local outfile="$1"

  if [[ -f "$outfile" && "$OVERWRITE" == false ]]; then
    return 0
  fi

  return 1
}

detect_scenario_from_path() {
  local path="$1"
  local scenario

  for scenario in "${SCENARIOS[@]}"; do
    if [[ "$path" == *"/$scenario/"* ]]; then
      echo "$scenario"
      return 0
    fi
  done

  echo "unknown"
}

detect_esm_from_path() {
  local path="$1"
  local esm

  for esm in "${ESMS[@]}"; do
    if [[ "$path" == *"/$esm/"* ]]; then
      echo "$esm"
      return 0
    fi
  done

  echo "unknown"
}


# ==============================================================================
# 3. PROCESS ONE 2-D FIELD
# ==============================================================================

process_2d_field() {
  local source_file="$1"
  local esm="$2"
  local scenario="$3"
  local variable="$4"
  local output_stem="$5"
  local level_label="$6"

  local outdir
  outdir="$(make_output_dir "$esm" "$scenario" "$variable")"

  local ntime
  ntime="$(cdo -s ntime "$source_file" 2>/dev/null || true)"

  if [[ -z "$ntime" || "$ntime" -lt 1 ]]; then
    echo "ERROR: Cannot determine timesteps: $source_file" | tee -a "$ERROR_FILE"
    log_summary "ERROR" "$esm" "$scenario" "$variable" "$source_file" "$level_label" "" "" "Cannot determine timesteps"
    return 1
  fi

  local tmp_regrid
  tmp_regrid="$(mktemp "/tmp/srw_regrid_${variable}_XXXXXX.nc")"

  echo "Regridding: $variable | $esm | $scenario | level=$level_label"

  if ! run_cmd cdo remapnn,"$GRID_FILE" "$source_file" "$tmp_regrid"; then
    echo "ERROR: remapnn failed: $source_file" | tee -a "$ERROR_FILE"
    log_summary "ERROR" "$esm" "$scenario" "$variable" "$source_file" "$level_label" "" "" "remapnn failed"
    safe_rm "$tmp_regrid"
    return 1
  fi

  local i
  for ((i=1; i<=ntime; i++)); do
    local date
    date="$(cdo -s showdate -seltimestep,"$i" "$tmp_regrid" 2>/dev/null \
      | grep -o '[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}' \
      | head -n 1 \
      | tr -d '-')"

    if [[ -z "$date" ]]; then
      echo "WARNING: Could not determine date for timestep $i in $source_file" | tee -a "$ERROR_FILE"
      log_summary "ERROR" "$esm" "$scenario" "$variable" "$source_file" "$level_label" "" "" "Missing date at timestep $i"
      continue
    fi

    local outfile="$outdir/${output_stem}_${date}.nc"

    if should_skip_existing "$outfile"; then
      echo "SKIP existing: $outfile"
      log_summary "SKIP" "$esm" "$scenario" "$variable" "$source_file" "$level_label" "$date" "$outfile" "Already exists"
      continue
    fi

    local tmp_crop
    tmp_crop="$(mktemp "/tmp/srw_crop_${variable}_XXXXXX.nc")"

    if run_cmd cdo sellonlatbox,"$CROP_EXTENT" -seltimestep,"$i" "$tmp_regrid" "$tmp_crop"; then
      if [[ "$DRY_RUN" == false ]]; then
        mv "$tmp_crop" "$outfile"
      else
        safe_rm "$tmp_crop"
      fi

      echo "SAVED: $outfile"
      log_summary "OK" "$esm" "$scenario" "$variable" "$source_file" "$level_label" "$date" "$outfile" ""
    else
      echo "ERROR: crop failed: $source_file timestep $i" | tee -a "$ERROR_FILE"
      log_summary "ERROR" "$esm" "$scenario" "$variable" "$source_file" "$level_label" "$date" "$outfile" "crop failed"
      safe_rm "$tmp_crop"
    fi
  done

  safe_rm "$tmp_regrid"
}


# ==============================================================================
# 4. SINGLE-LAYER VARIABLES
# ==============================================================================

process_single_layer_file() {
  local input_file="$1"
  local esm="$2"
  local scenario="$3"
  local variable="$4"

  echo
  echo "------------------------------------------------------------"
  echo "SINGLE-LAYER: $variable"
  echo "$input_file"
  echo "------------------------------------------------------------"

  local stem
  stem="$(basename "$input_file" .nc)"

  process_2d_field "$input_file" "$esm" "$scenario" "$variable" "$stem" "single"
}


# ==============================================================================
# 5. MULTI-LAYER VARIABLES
# ==============================================================================

process_multi_layer_file() {
  local input_file="$1"
  local esm="$2"
  local scenario="$3"
  local variable="$4"

  echo
  echo "------------------------------------------------------------"
  echo "MULTI-LAYER: $variable"
  echo "$input_file"
  echo "------------------------------------------------------------"

  mapfile -t levels < <(get_levels "$input_file")

  if [[ "${#levels[@]}" -eq 0 ]]; then
    echo "ERROR: No depth levels detected: $input_file" | tee -a "$ERROR_FILE"
    log_summary "ERROR" "$esm" "$scenario" "$variable" "$input_file" "" "" "" "No depth levels detected"
    return 1
  fi

  echo "Available levels: ${levels[*]}"

  local surface_level
  local deep_level

  surface_level="$(printf "%s\n" "${levels[@]}" | choose_surface_level)"
  deep_level="$(printf "%s\n" "${levels[@]}" | choose_deep_level)"

  if [[ -z "$surface_level" ]]; then
    echo "ERROR: Cannot choose surface level: $input_file" | tee -a "$ERROR_FILE"
    log_summary "ERROR" "$esm" "$scenario" "$variable" "$input_file" "" "" "" "Cannot choose surface level"
    return 1
  fi

  if [[ -z "$deep_level" ]]; then
    echo "ERROR: Cannot choose ~200 m level: $input_file" | tee -a "$ERROR_FILE"
    log_summary "ERROR" "$esm" "$scenario" "$variable" "$input_file" "" "" "" "Cannot choose deep level"
    return 1
  fi

  if [[ "$surface_level" == "$deep_level" ]]; then
    echo "ERROR: Surface and deep level resolve to the same value ($surface_level): $input_file" | tee -a "$ERROR_FILE"
    log_summary "ERROR" "$esm" "$scenario" "$variable" "$input_file" "$surface_level" "" "" "Surface and deep levels identical"
    return 1
  fi

  local surface_label
  local deep_label
  surface_label="$(format_level_for_filename "$surface_level")"
  deep_label="$(format_level_for_filename "$deep_level")"

  echo "Selected surface level: $surface_level"
  echo "Selected deep level   : $deep_level"

  local stem
  stem="$(basename "$input_file" .nc)"

  local tmp_surface
  local tmp_deep
  tmp_surface="$(mktemp "/tmp/srw_surface_${variable}_XXXXXX.nc")"
  tmp_deep="$(mktemp "/tmp/srw_deep_${variable}_XXXXXX.nc")"

  if ! run_cmd cdo sellevel,"$surface_level" "$input_file" "$tmp_surface"; then
    echo "ERROR: Cannot extract surface level $surface_level from $input_file" | tee -a "$ERROR_FILE"
    log_summary "ERROR" "$esm" "$scenario" "$variable" "$input_file" "$surface_level" "" "" "Surface extraction failed"
    safe_rm "$tmp_surface" "$tmp_deep"
    return 1
  fi

  if ! run_cmd cdo sellevel,"$deep_level" "$input_file" "$tmp_deep"; then
    echo "ERROR: Cannot extract deep level $deep_level from $input_file" | tee -a "$ERROR_FILE"
    log_summary "ERROR" "$esm" "$scenario" "$variable" "$input_file" "$deep_level" "" "" "Deep extraction failed"
    safe_rm "$tmp_surface" "$tmp_deep"
    return 1
  fi

  process_2d_field \
    "$tmp_surface" "$esm" "$scenario" "$variable" \
    "${stem}_layer${surface_label}" "$surface_label"

  process_2d_field \
    "$tmp_deep" "$esm" "$scenario" "$variable" \
    "${stem}_layer${deep_label}" "$deep_label"

  safe_rm "$tmp_surface" "$tmp_deep"
}


# ==============================================================================
# 6. PROCESS ONE RAW NETCDF
# ==============================================================================

process_file() {
  local input_file="$1"

  local esm
  local scenario
  local variable

  esm="$(detect_esm_from_path "$input_file")"
  scenario="$(detect_scenario_from_path "$input_file")"
  variable="$(extract_variable_name "$input_file")"

  if [[ "$esm" == "unknown" ]]; then
    echo "SKIP: Unknown ESM: $input_file"
    return 0
  fi

  if [[ "$scenario" == "unknown" ]]; then
    echo "SKIP: Unknown scenario: $input_file"
    return 0
  fi

  if [[ -z "$variable" ]]; then
    echo "SKIP: Could not identify variable: $input_file"
    log_summary "SKIP" "$esm" "$scenario" "" "$input_file" "" "" "" "Unknown variable"
    return 0
  fi

  if in_array "$variable" "${SINGLE_LAYER_VARIABLES[@]}"; then
    process_single_layer_file "$input_file" "$esm" "$scenario" "$variable"
    return $?
  fi

  if in_array "$variable" "${MULTI_LAYER_VARIABLES[@]}"; then
    local nlevels
    nlevels="$(cdo -s nlevel "$input_file" 2>/dev/null | awk '{print $1}' || true)"

    if [[ -n "$nlevels" && "$nlevels" -le 1 ]]; then
      echo "NOTE: $variable has only one level in this file; processing as 2-D."
      process_single_layer_file "$input_file" "$esm" "$scenario" "$variable"
    else
      process_multi_layer_file "$input_file" "$esm" "$scenario" "$variable"
    fi

    return $?
  fi

  echo "SKIP unsupported variable '$variable': $input_file"
  log_summary "SKIP" "$esm" "$scenario" "$variable" "$input_file" "" "" "" "Variable not part of SRW processing set"
}


# ==============================================================================
# 7. MAIN LOOP
# ==============================================================================

for esm in "${ESMS[@]}"; do
  for scenario in "${SCENARIOS[@]}"; do

    scenario_dir="$INPUT_ROOT/$esm/$scenario"

    if [[ ! -d "$scenario_dir" ]]; then
      echo
      echo "NOTE: Missing directory, skipping: $scenario_dir"
      continue
    fi

    echo
    echo "============================================================"
    echo "PROCESSING $esm | $scenario"
    echo "============================================================"

    while IFS= read -r -d '' input_file; do
      process_file "$input_file" || true
    done < <(
      find "$scenario_dir" \
        -type f \
        -name "*.nc" \
        ! -name "._*" \
        ! -name "*.gr.*" \
        -print0
    )

  done
done


# ==============================================================================
# 8. FINISH
# ==============================================================================

echo
echo "============================================================"
echo "Finished: $(date)"
echo "Main log   : $LOG_FILE"
echo "Error log  : $ERROR_FILE"
echo "Summary TSV: $SUMMARY_FILE"
echo "============================================================"
