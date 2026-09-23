# ======================================================================
# V1 — LOAD AND VERIFY THE FROZEN TCGA CLASSIFIER
#
# PURPOSE
# -------
# Load the classifier artifacts produced by FREEZE_TCGA_classifier.R
# and confirm they are structurally valid before using them on
# an independent cohort. This step performs NO fitting - it only reads and checks.
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

frozen_dir <- "/content/TCGA_frozen_classifier"

required_files <- c(
  file.path(frozen_dir, "FROZEN_module_stats.csv"),
  file.path(frozen_dir, "FROZEN_state_centroids.csv"),
  file.path(frozen_dir, "FROZEN_classifier_metadata.csv")
)

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  cat("Missing frozen classifier files:\n")
  print(missing)
  stop("Upload the TCGA_frozen_classifier/ folder (produced by ",
       "FREEZE_TCGA_classifier.R) into /content/ before continuing.")
}

frozen_stats     <- read_csv(required_files[1], show_col_types = FALSE)
frozen_centroids <- read_csv(required_files[2], show_col_types = FALSE)
frozen_metadata  <- read_csv(required_files[3], show_col_types = FALSE)

cat("=====================================================================\n")
cat("FROZEN TCGA CLASSIFIER LOADED\n")
cat("=====================================================================\n")
print(frozen_metadata)

cat("\nFrozen modules (n =", nrow(frozen_stats), "):\n")
print(frozen_stats$module)

cat("\nFrozen states:\n")
print(frozen_centroids$ecosystem_state)

# ---- structural checks ----
if (!all(c("module", "tcga_mean", "tcga_sd") %in% colnames(frozen_stats))) {
  stop("FROZEN_module_stats.csv is missing expected columns.")
}
if (any(!is.finite(frozen_stats$tcga_sd)) || any(frozen_stats$tcga_sd <= 0)) {
  stop("FROZEN_module_stats.csv contains a zero or non-finite tcga_sd - ",
       "cannot standardize a new cohort with this classifier.")
}
if (nrow(frozen_centroids) < 2) {
  stop("Fewer than 2 frozen state centroids found - check FREEZE_TCGA_classifier.R output.")
}

cat("\nStructural checks passed. Classifier is ready to apply to GSE72094 in V5.\n")
