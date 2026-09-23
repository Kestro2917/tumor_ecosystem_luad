# ======================================================================
# V2c — CONVERT HARMONIZED MATRIX INTO M1-COMPATIBLE FORMAT
#
# PURPOSE
# -------
# The reused M3 LIONESS network script (V3, unmodified from your TCGA
# pipeline) expects a RAW/LINEAR-scale expression matrix and applies
# its own log2(x+1) transform internally — matching how the TCGA TPM
# matrix was processed. V2a/V2b (correctly) produced a LOG2-scale
# matrix. This step un-logs it (2^x - 1) so V3's own log2(x+1) is
# applied exactly once, avoiding a silent double-log-transform.
#
# INPUT (from V2b)
# -----------------
# /content/GSE72094_raw/clean/GSE72094_gene_expression_log2_harmonized.csv
#
# OUTPUT (matches the exact contract V3 expects)
# ------------------------------------------------
# /content/GSE72094/M1_expression_matrix.csv
#     columns: matrix_gene_id, <sample_1>, <sample_2>, ...
# /content/GSE72094/02_ProteinCoding_GeneAnnotation.csv
#     columns: matrix_gene_id, gene_symbol
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

harmonized_file <- "/content/GSE72094_raw/clean/GSE72094_gene_expression_log2_harmonized.csv"
output_dir <- "/content/GSE72094"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(harmonized_file)) {
  stop("ERROR: ", harmonized_file, " not found. Run V2a and V2b first.")
}

harmonized <- read_csv(harmonized_file, show_col_types = FALSE)

gene_col <- colnames(harmonized)[1L]
sample_cols <- setdiff(colnames(harmonized), gene_col)

gene_matrix <- as.matrix(harmonized[, sample_cols, drop = FALSE])
storage.mode(gene_matrix) <- "double"

cat("Harmonized matrix (log2 scale): ", nrow(gene_matrix), "genes x", ncol(gene_matrix), "samples\n")
cat("Value range (log2 scale):", round(min(gene_matrix, na.rm = TRUE), 3), "to",
    round(max(gene_matrix, na.rm = TRUE), 3), "\n")

# ---- un-log back to linear scale, matching TCGA's TPM (raw, non-negative) convention ----
linear_matrix <- 2^gene_matrix - 1
linear_matrix[linear_matrix < 0] <- 0

cat("Converted to linear scale for V3's own log2(x + 1) step.\n")
cat("Value range (linear scale):", round(min(linear_matrix, na.rm = TRUE), 3), "to",
    round(max(linear_matrix, na.rm = TRUE), 3), "\n\n")

# ---- write M1-compatible expression matrix ----
expr_out <- as.data.frame(linear_matrix)
expr_out <- cbind(matrix_gene_id = harmonized[[gene_col]], expr_out)

write_csv(expr_out, file.path(output_dir, "M1_expression_matrix.csv"))

# ---- write M1-compatible annotation (identity mapping: rows are already gene symbols) ----
annotation_out <- data.frame(
  matrix_gene_id = harmonized[[gene_col]],
  gene_symbol    = toupper(trimws(harmonized[[gene_col]])),
  stringsAsFactors = FALSE
)
write_csv(annotation_out, file.path(output_dir, "02_ProteinCoding_GeneAnnotation.csv"))

cat("=====================================================================\n")
cat("GSE72094 READY IN M1-COMPATIBLE FORMAT\n")
cat("=====================================================================\n")
cat("Expression matrix:", file.path(output_dir, "M1_expression_matrix.csv"), "\n")
cat("  Genes:", nrow(expr_out), " Samples:", length(sample_cols), "\n")
cat("Annotation file:", file.path(output_dir, "02_ProteinCoding_GeneAnnotation.csv"), "\n\n")
cat("Next: run V3 (M3 network, reused unmodified) on these two files.\n")
