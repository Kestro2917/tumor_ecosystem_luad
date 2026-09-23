# ======================================================================
# V5 — APPLY THE FROZEN TCGA CLASSIFIER TO GSE72094
#
# INPUTS
# ------
# /content/TCGA_frozen_classifier/FROZEN_module_stats.csv       (from V1)
# /content/TCGA_frozen_classifier/FROZEN_state_centroids.csv    (from V1)
# /content/GSE72094_M3_LIONESS_entropy/M3_LIONESS_module_entropy.csv (from V4)
#
# OUTPUT
# ------
# /content/GSE72094_validation/GSE72094_ecosystem_states.csv
# /content/GSE72094_validation/V5_module_coverage_check.csv
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

gse_entropy_file <- "/content/GSE72094_M3_LIONESS_entropy/M3_LIONESS_module_entropy.csv"
output_dir <- "/content/GSE72094_validation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(gse_entropy_file)) {
  stop("ERROR: GSE72094 module entropy file not found:\n", gse_entropy_file,
       "\nRun V3 and V4 first.")
}

# frozen_stats / frozen_centroids are expected already in memory from V1.
# If you're running this cell fresh, uncomment the two lines below:
# frozen_stats     <- read_csv("/content/TCGA_frozen_classifier/FROZEN_module_stats.csv", show_col_types = FALSE)
# frozen_centroids <- read_csv("/content/TCGA_frozen_classifier/FROZEN_state_centroids.csv", show_col_types = FALSE)

gse_entropy <- read_csv(gse_entropy_file, show_col_types = FALSE)

# ---- identify sample id column (same convention used throughout the pipeline) ----
sample_candidates <- c("sample_id", "Sample_ID", "sample", "patient", "Patient")
existing_candidates <- sample_candidates[sample_candidates %in% colnames(gse_entropy)]
sample_col <- if (length(existing_candidates) > 0) existing_candidates[1] else colnames(gse_entropy)[1]

gse_entropy$sample_id <- trimws(as.character(gse_entropy[[sample_col]]))

# The reused M3-entropy script derives sample_id from each LIONESS network
# filename ("{id}_LIONESS_network.csv"), and its filename-parsing regex
# leaves a stray "_LIONESS" suffix attached (e.g. "GSM772610_LIONESS").
# Your original TCGA M3.6 script silently corrects for this same artifact
# via sub("_LIONESS$", "", x) before its own clinical merge - applying the
# identical fix here, at the point the ID is first read, so every
# downstream step (this file's output, V6a, V6b, V7) sees the clean ID.
gse_entropy$sample_id <- sub("_LIONESS$", "", gse_entropy$sample_id)

gse_entropy <- gse_entropy[!duplicated(gse_entropy$sample_id), , drop = FALSE]

frozen_modules <- frozen_stats$module

# ---- module coverage check ----
missing_modules <- setdiff(frozen_modules, colnames(gse_entropy))

if (length(missing_modules) > 0) {
  cat("WARNING:", length(missing_modules),
      "frozen TCGA modules are missing entirely from GSE72094's output:\n")
  print(missing_modules)
  cat("\nThis usually means no genes from that MSigDB category survived\n")
  cat("GSE72094's own network/QC steps. Missing modules are set to the\n")
  cat("FROZEN TCGA mean for every GSE72094 patient below (i.e. treated as\n")
  cat("'average' rather than assigning them an arbitrary value) — this\n")
  cat("keeps the projection purely frozen rather than inventing new\n")
  cat("information from GSE72094 itself.\n\n")
  for (m in missing_modules) gse_entropy[[m]] <- NA_real_
}

for (m in frozen_modules) {
  gse_entropy[[m]] <- suppressWarnings(as.numeric(gse_entropy[[m]]))
}

coverage <- data.frame(
  module = frozen_modules,
  present_in_gse = !(frozen_modules %in% missing_modules),
  pct_missing_in_gse = sapply(frozen_modules, function(m) mean(!is.finite(gse_entropy[[m]]))),
  stringsAsFactors = FALSE
)
write_csv(coverage, file.path(output_dir, "V5_module_coverage_check.csv"))

cat("Frozen modules used:", length(frozen_modules), "\n")
cat("Modules with any missing values in GSE72094:",
    sum(coverage$pct_missing_in_gse > 0), "\n\n")

# ---- build the GSE72094 patient x module matrix, in the SAME module order ----
X_new <- as.matrix(gse_entropy[, frozen_modules, drop = FALSE])
rownames(X_new) <- gse_entropy$sample_id

tcga_mean <- setNames(frozen_stats$tcga_mean, frozen_stats$module)
tcga_sd   <- setNames(frozen_stats$tcga_sd,   frozen_stats$module)

# ---- median-impute any missing values using the FROZEN TCGA mean
#      (never a GSE72094-specific value — keeps the projection frozen) ----
for (j in seq_len(ncol(X_new))) {
  bad <- !is.finite(X_new[, j])
  if (any(bad)) X_new[bad, j] <- tcga_mean[colnames(X_new)[j]]
}

# ---- standardize using the FROZEN TCGA mean/sd — NEVER GSE72094's own ----
Z_new <- sweep(X_new, 2, tcga_mean[colnames(X_new)], "-")
Z_new <- sweep(Z_new, 2, tcga_sd[colnames(X_new)],   "/")

if (any(!is.finite(Z_new))) {
  stop("ERROR: non-finite values after standardization — check ",
       "FROZEN_module_stats.csv for a zero-sd module.")
}

# ---- nearest-centroid assignment (Euclidean distance, standardized space) ----
centroid_matrix <- as.matrix(frozen_centroids[, frozen_modules, drop = FALSE])
rownames(centroid_matrix) <- frozen_centroids$ecosystem_state

assign_nearest_centroid <- function(x, centroids) {
  d <- apply(centroids, 1, function(c) sqrt(sum((x - c)^2)))
  names(which.min(d))
}

assigned_state <- apply(Z_new, 1, assign_nearest_centroid, centroids = centroid_matrix)

# ---- assignment-confidence margin: relative gap between nearest and
#      second-nearest centroid (0 = ambiguous/borderline, 1 = confident) ----
distance_margin <- apply(Z_new, 1, function(x) {
  d <- sort(apply(centroid_matrix, 1, function(c) sqrt(sum((x - c)^2))))
  if (length(d) < 2 || d[2] == 0) return(NA_real_)
  (d[2] - d[1]) / d[2]
})

gse_states <- data.frame(
  sample_id = rownames(Z_new),
  ecosystem_state = assigned_state,
  assignment_margin = distance_margin,
  stringsAsFactors = FALSE
)

write_csv(gse_states, file.path(output_dir, "GSE72094_ecosystem_states.csv"))

cat("=====================================================================\n")
cat("GSE72094 ECOSYSTEM-STATE ASSIGNMENT (FROZEN TCGA CLASSIFIER)\n")
cat("=====================================================================\n")
print(table(gse_states$ecosystem_state))

cat("\nMean assignment margin:", round(mean(gse_states$assignment_margin, na.rm = TRUE), 3), "\n")
cat("(margin near 0 = the patient sits almost equidistant between the two\n")
cat(" frozen centroids, i.e. a borderline assignment; margin near 1 = a\n")
cat(" confident assignment)\n\n")

cat("Saved to:", file.path(output_dir, "GSE72094_ecosystem_states.csv"), "\n")
cat("Next: run V6 to test survival separation under this locked assignment.\n")
