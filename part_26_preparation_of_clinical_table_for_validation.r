# ======================================================================
# V6a — BUILD CLINICAL/SURVIVAL TABLE FROM GEO CHARACTERISTICS
#
#
# INPUT (from V2a)
# -----------------
# /content/GSE72094_raw/clean/GSE72094_sample_characteristics_long.csv
#
# OUTPUT
# ------
# /content/GSE72094_validation/V6a_GSE72094_clinical_wide.csv
# ======================================================================

GSE72094_CLINICAL_OVERRIDE_FILE <- NA
# ^ EDIT THIS only if GEO's own characteristics (printed below) do not
#   contain what you need. Otherwise leave as NA.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

output_dir <- "/content/GSE72094_validation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

characteristics_file <- "/content/GSE72094_raw/clean/GSE72094_sample_characteristics_long.csv"

if (!is.na(GSE72094_CLINICAL_OVERRIDE_FILE)) {

  cat("Using manually supplied clinical file:", GSE72094_CLINICAL_OVERRIDE_FILE, "\n")
  if (!file.exists(GSE72094_CLINICAL_OVERRIDE_FILE)) {
    stop("ERROR: ", GSE72094_CLINICAL_OVERRIDE_FILE, " not found.")
  }
  clinical_wide <- read_csv(GSE72094_CLINICAL_OVERRIDE_FILE, show_col_types = FALSE)

} else {

  if (!file.exists(characteristics_file)) {
    stop("ERROR: ", characteristics_file, " not found. Run V2a first, or ",
         "set GSE72094_CLINICAL_OVERRIDE_FILE above.")
  }

  characteristics_long <- read_csv(characteristics_file, show_col_types = FALSE)

  cat("Characteristic fields available from GEO for GSE72094:\n")
  print(sort(unique(characteristics_long$characteristic_name)))
  cat("\n")

  # ---- pivot to one row per sample, one column per characteristic ----
  clinical_wide <- characteristics_long %>%
    dplyr::distinct(geo_sample_id, characteristic_name, .keep_all = TRUE) %>%
    tidyr::pivot_wider(
      id_cols = geo_sample_id,
      names_from = characteristic_name,
      values_from = characteristic_value
    )
}

clinical_wide$sample_id <- trimws(as.character(clinical_wide$geo_sample_id))

write_csv(clinical_wide, file.path(output_dir, "V6a_GSE72094_clinical_wide.csv"))

cat("=====================================================================\n")
cat("CLINICAL/SURVIVAL TABLE BUILT\n")
cat("=====================================================================\n")
cat("Columns available (pass the right ones into V6 below if the\n")
cat("auto-detected time/status columns look wrong):\n")
print(colnames(clinical_wide))
cat("\nSaved to:", file.path(output_dir, "V6a_GSE72094_clinical_wide.csv"), "\n")
cat("Preview:\n")
print(utils::head(clinical_wide))
