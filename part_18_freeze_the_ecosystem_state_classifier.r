# ======================================================================
# M3.2-LOCK — FREEZE THE TCGA ECOSYSTEM-STATE CLASSIFIER
#
# PURPOSE
# -------
# Produce a portable, frozen version of the M3.2 ecosystem-state
# classifier so it can be APPLIED (never re-fit) to an independent
# cohort, following the "locked classifier" design:
#
#   1. Freeze the TCGA module list, the standardization parameters
#      (mean/sd per module), and the State_1/State_2 centroids.
#   2. In the new cohort: standardize using these FROZEN mean/sd
#      values (never recompute mean/sd on the new cohort) and
#      assign each patient to the nearest FROZEN centroid.
#   3. Test survival separation under that locked assignment.
#
#
# INPUTS (already produced by M3.2 — unchanged)
# ----------------------------------------------
# /content/M3_LIONESS_entropy/M3_LIONESS_module_entropy.csv
# /content/M3_LIONESS_entropy/M3.2_module_axes/M3.2_patient_ecosystem_states.csv
# /content/M3_LIONESS_entropy/M3.2_module_axes/M3.2_standardized_module_entropy.csv
#
# OUTPUT
# ------
# /content/TCGA_frozen_classifier/
#     FROZEN_module_stats.csv        (module, tcga_mean, tcga_sd)
#     FROZEN_state_centroids.csv     (ecosystem_state x module, standardized)
#     FROZEN_classifier_metadata.csv
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

input_dir  <- "/content/M3_LIONESS_entropy"
axes_dir   <- file.path(input_dir, "M3.2_module_axes")
output_dir <- "/content/TCGA_frozen_classifier"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

module_entropy_file <- file.path(input_dir, "M3_LIONESS_module_entropy.csv")
states_file         <- file.path(axes_dir, "M3.2_patient_ecosystem_states.csv")
saved_Z_file        <- file.path(axes_dir, "M3.2_standardized_module_entropy.csv")

if (!file.exists(module_entropy_file)) {
  stop("ERROR: raw module entropy file not found:\n", module_entropy_file,
       "\nRun M3 (network + entropy) first.")
}
if (!file.exists(states_file)) {
  stop("ERROR: M3.2 state-assignment file not found:\n", states_file,
       "\nRun M3.2 first.")
}

cat("=====================================================================\n")
cat("FREEZING TCGA ECOSYSTEM-STATE CLASSIFIER\n")
cat("=====================================================================\n\n")

module_entropy <- read_csv(module_entropy_file, show_col_types = FALSE)
patient_states <- read_csv(states_file, show_col_types = FALSE)

# ---- reproduce sample_id exactly as M3.2 did ----
sample_candidates <- c("sample_id", "Sample_ID", "sample", "patient", "Patient")
existing_candidates <- sample_candidates[sample_candidates %in% colnames(module_entropy)]
sample_col <- if (length(existing_candidates) > 0) existing_candidates[1] else colnames(module_entropy)[1]

module_entropy$sample_id <- trimws(as.character(module_entropy[[sample_col]]))
module_entropy <- module_entropy[!duplicated(module_entropy$sample_id), , drop = FALSE]

# ---- reproduce module QC exactly as M3.2 did ----
module_columns <- grep("^S_module_", colnames(module_entropy), value = TRUE)
if (length(module_columns) < 3) {
  stop("ERROR: fewer than 3 S_module_ features detected in the raw entropy file.")
}

for (m in module_columns) {
  module_entropy[[m]] <- suppressWarnings(as.numeric(module_entropy[[m]]))
}

module_qc <- data.frame(
  module = module_columns,
  n  = sapply(module_columns, function(x) sum(is.finite(module_entropy[[x]]))),
  sd = sapply(module_columns, function(x) {
    y <- module_entropy[[x]]
    y <- y[is.finite(y)]
    if (length(y) > 1) sd(y) else NA_real_
  }),
  stringsAsFactors = FALSE
)

usable_modules <- module_qc %>%
  filter(n >= 0.90 * nrow(module_entropy), is.finite(sd), sd > 0) %>%
  pull(module)

cat("Usable modules recovered from raw entropy file:", length(usable_modules), "\n")

if (length(usable_modules) < 3) {
  stop("ERROR: fewer than 3 usable modules after QC.")
}

# ---- reproduce patient x module matrix + median imputation exactly as M3.2 did ----
X <- module_entropy[, c("sample_id", usable_modules), drop = FALSE]
X_matrix <- as.matrix(X[, usable_modules, drop = FALSE])
rownames(X_matrix) <- X$sample_id

for (j in seq_len(ncol(X_matrix))) {
  bad <- !is.finite(X_matrix[, j])
  if (any(bad)) {
    med <- median(X_matrix[, j], na.rm = TRUE)
    if (!is.finite(med)) {
      stop("ERROR: module has no valid values: ", colnames(X_matrix)[j])
    }
    X_matrix[bad, j] <- med
  }
}

# ---- reproduce standardization exactly as M3.2 did, and SAVE the parameters ----
tcga_mean <- apply(X_matrix, 2, mean)
tcga_sd   <- apply(X_matrix, 2, sd)

Z <- scale(X_matrix, center = tcga_mean, scale = tcga_sd)
Z <- as.matrix(Z)

if (any(!is.finite(Z))) {
  stop("ERROR: non-finite values after standardization - cannot freeze classifier.")
}

# ---- sanity check against M3.2's own saved standardized matrix ----
if (file.exists(saved_Z_file)) {

  saved_Z <- read.csv(saved_Z_file, row.names = 1, check.names = FALSE)

  common_ids <- intersect(rownames(Z), rownames(saved_Z))
  common_mod <- intersect(colnames(Z), colnames(saved_Z))

  max_diff <- max(
    abs(Z[common_ids, common_mod] - as.matrix(saved_Z[common_ids, common_mod])),
    na.rm = TRUE
  )

  cat("\nVerification against M3.2's own saved standardized matrix:\n")
  cat("  Common patients:", length(common_ids), "\n")
  cat("  Common modules :", length(common_mod), "\n")
  cat("  Max abs difference:", max_diff, "\n")

  if (!is.finite(max_diff) || max_diff > 1e-6) {
    warning("Recomputed matrix does NOT exactly match M3.2's saved output. ",
            "Inspect before trusting the frozen classifier.")
  } else {
    cat("  PASSED — recomputation matches M3.2 exactly.\n")
  }

} else {
  cat("\nNOTE: M3.2_standardized_module_entropy.csv not found — skipping the\n")
  cat("cross-check (freezing will still proceed from the raw entropy file).\n")
}

# ---- attach the already-saved M3.2 state labels ----
patient_states$sample_id <- trimws(as.character(patient_states$sample_id))
state_lookup <- patient_states %>% select(sample_id, ecosystem_state)

Z_df <- as.data.frame(Z)
Z_df$sample_id <- rownames(Z)
Z_df <- Z_df %>% inner_join(state_lookup, by = "sample_id")

if (nrow(Z_df) != nrow(Z)) {
  warning(nrow(Z) - nrow(Z_df),
          " patients had no matching M3.2 state label and were dropped",
          " from the centroid calculation.")
}

# ---- centroids = mean of standardized values per state (== kmeans centroids) ----
frozen_centroids <- Z_df %>%
  group_by(ecosystem_state) %>%
  summarise(across(all_of(usable_modules), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

cat("\nFrozen state centroids (n patients per state):\n")
print(Z_df %>% count(ecosystem_state))

# ---- save frozen artifacts ----
frozen_module_stats <- data.frame(
  module    = usable_modules,
  tcga_mean = as.numeric(tcga_mean[usable_modules]),
  tcga_sd   = as.numeric(tcga_sd[usable_modules]),
  stringsAsFactors = FALSE
)

write_csv(frozen_module_stats, file.path(output_dir, "FROZEN_module_stats.csv"))
write_csv(frozen_centroids,    file.path(output_dir, "FROZEN_state_centroids.csv"))

frozen_metadata <- data.frame(
  n_tcga_patients = nrow(Z),
  n_modules       = length(usable_modules),
  n_states        = nrow(frozen_centroids),
  states          = paste(frozen_centroids$ecosystem_state, collapse = ", "),
  frozen_on       = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_csv(frozen_metadata, file.path(output_dir, "FROZEN_classifier_metadata.csv"))

cat("\n=====================================================================\n")
cat("TCGA ECOSYSTEM-STATE CLASSIFIER FROZEN\n")
cat("=====================================================================\n")
print(frozen_metadata)

cat("\nFrozen files written to:", output_dir, "\n")
cat("  - FROZEN_module_stats.csv      (module, tcga_mean, tcga_sd)\n")
cat("  - FROZEN_state_centroids.csv   (per-state centroid, standardized space)\n")
cat("  - FROZEN_classifier_metadata.csv\n\n")

cat("IMPORTANT:\n")
cat("Download this entire TCGA_frozen_classifier/ folder (or copy it to\n")
cat("Google Drive) and upload it into the GSE31210 validation notebook.\n")
cat("Do NOT regenerate these files from a re-clustered cohort — they are\n")
cat("only valid as long as they come from this exact TCGA M3.2 run.\n")
