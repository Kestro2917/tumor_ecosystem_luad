# ======================================================================
# V2b — PROBE-TO-GENE HARMONIZATION
#
# -----------------
# /content/GSE72094_raw/clean/GSE72094_expression_log2_matrix.csv
# /content/GSE72094_raw/clean/GSE72094_probe_annotation_raw.csv
#
# OUTPUT
# ------
# /content/GSE72094_raw/clean/GSE72094_gene_expression_log2_harmonized.csv
# /content/GSE72094_raw/tables/GSE72094_expression_harmonization_QC.csv
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
})

DATA_DIR  <- "/content"
OUT_DIR   <- file.path(DATA_DIR, "GSE72094_raw")
CLEAN_DIR <- file.path(OUT_DIR, "clean")
TABLE_DIR <- file.path(OUT_DIR, "tables")

expression_csv_file <- file.path(CLEAN_DIR, "GSE72094_expression_log2_matrix.csv")
annotation_csv_file <- file.path(CLEAN_DIR, "GSE72094_probe_annotation_raw.csv")

if (!file.exists(expression_csv_file) || !file.exists(annotation_csv_file)) {
  stop("V2a output not found. Run V2a first.")
}

# ---- helper functions (unchanged from your tested notebook) ----
detect_column <- function(data, candidates, required = TRUE) {
  exact_match <- intersect(candidates, colnames(data))
  if (length(exact_match) > 0L) return(exact_match[1L])
  lower_names <- tolower(colnames(data))
  lower_candidates <- tolower(candidates)
  matched_index <- match(lower_candidates, lower_names)
  matched_index <- matched_index[!is.na(matched_index)]
  if (length(matched_index) > 0L) return(colnames(data)[matched_index[1L]])
  if (required) {
    stop("None of the expected columns were found: ", paste(candidates, collapse = ", "),
         "\nAvailable columns:\n", paste(colnames(data), collapse = ", "))
  }
  NA_character_
}

clean_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x == "---" | x == "NA" | x == "N/A" | x == "null"] <- NA_character_
  x <- gsub('^"|"$', "", x)
  trimws(x)
}

is_ambiguous_symbol <- function(x) {
  x <- as.character(x)
  is.na(x) | !nzchar(trimws(x)) | grepl("///|//|;|\\||,", x)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ---- load V2a output ----
expression_table <- readr::read_csv(expression_csv_file, show_col_types = FALSE)
probe_column <- colnames(expression_table)[1L]
expression_log2 <- as.matrix(expression_table[, -1L, drop = FALSE])
storage.mode(expression_log2) <- "double"
rownames(expression_log2) <- as.character(expression_table[[probe_column]])

feature_annotation <- readr::read_csv(annotation_csv_file, show_col_types = FALSE)

# ---- identify probe / gene-symbol columns in the platform annotation ----
probe_id_col <- detect_column(
  feature_annotation,
  c("probe_id", "ID", "ID_REF", "ProbeID", "probe", "SPOT_ID")
)
gene_symbol_col <- detect_column(
  feature_annotation,
  c("GeneSymbol", "Gene Symbol", "Gene symbol", "GENE_SYMBOL", "gene_symbol", "Symbol",
    "SYMBOL", "Gene.Symbol", "gene_assignment", "Gene Assignment")
)
message("Probe identifier column: ", probe_id_col)
message("Gene-symbol column: ", gene_symbol_col)

# ---- build clean probe -> gene mapping ----
probe_annotation_clean <- feature_annotation %>%
  dplyr::transmute(
    probe_id = trimws(as.character(.data[[probe_id_col]])),
    original_gene_symbol = clean_gene_symbol(.data[[gene_symbol_col]])
  ) %>%
  dplyr::mutate(
    ambiguous_mapping = is_ambiguous_symbol(original_gene_symbol),
    gene_symbol = dplyr::if_else(ambiguous_mapping, NA_character_, original_gene_symbol),
    gene_symbol = trimws(gene_symbol)
  ) %>%
  dplyr::filter(!is.na(probe_id), nzchar(probe_id)) %>%
  dplyr::distinct(probe_id, .keep_all = TRUE)

# ---- probe overlap check ----
expression_probe_ids <- rownames(expression_log2)
common_probes <- intersect(expression_probe_ids, probe_annotation_clean$probe_id)
if (length(common_probes) == 0L) {
  stop("No expression probes matched the probe-annotation table.")
}
probe_match_rate <- length(common_probes) / length(expression_probe_ids)
message("Matched probes: ", length(common_probes), " of ", length(expression_probe_ids),
        " (", round(100 * probe_match_rate, 2), "%)")

# ---- merge, drop unmapped/ambiguous, collapse duplicates by median ----
expression_probe_table <- as.data.frame(expression_log2) %>%
  tibble::rownames_to_column("probe_id") %>%
  dplyr::mutate(probe_id = as.character(probe_id))

probe_expression_annotated <- expression_probe_table %>%
  dplyr::inner_join(probe_annotation_clean, by = "probe_id")

sample_columns <- colnames(expression_log2)

probe_expression_mapped <- probe_expression_annotated %>%
  dplyr::filter(!is.na(gene_symbol), nzchar(gene_symbol))

if (nrow(probe_expression_mapped) == 0L) {
  stop("No probes remained after gene-symbol mapping.")
}

# Median collapsing: microarray probe intensities are not additive.
gene_expression_table <- probe_expression_mapped %>%
  dplyr::select(gene_symbol, dplyr::all_of(sample_columns)) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sample_columns), ~ stats::median(safe_numeric(.x), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(sample_columns), ~ dplyr::if_else(is.nan(.x), NA_real_, .x)))

gene_expression_log2 <- gene_expression_table %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()
storage.mode(gene_expression_log2) <- "double"

# ---- drop genes with no usable expression at all ----
usable_gene <- apply(gene_expression_log2, 1L, function(x) sum(is.finite(x)) > 0L)
gene_expression_log2 <- gene_expression_log2[usable_gene, , drop = FALSE]

if (nrow(gene_expression_log2) == 0L) {
  stop("No genes remained after expression quality control.")
}
if (anyDuplicated(rownames(gene_expression_log2))) {
  stop("Duplicated gene symbols remain after probe collapsing.")
}

# ---- save harmonized gene-level matrix ----
saveRDS(gene_expression_log2, file.path(CLEAN_DIR, "GSE72094_gene_expression_log2_harmonized.rds"))
readr::write_csv(
  as.data.frame(gene_expression_log2) %>% tibble::rownames_to_column("gene_symbol"),
  file.path(CLEAN_DIR, "GSE72094_gene_expression_log2_harmonized.csv")
)

probe_count_per_gene <- probe_expression_mapped %>%
  dplyr::count(gene_symbol, name = "number_of_probes") %>%
  dplyr::arrange(dplyr::desc(number_of_probes), gene_symbol)

readr::write_csv(probe_annotation_clean, file.path(TABLE_DIR, "GSE72094_probe_to_gene_mapping_clean.csv"))
readr::write_csv(probe_count_per_gene, file.path(TABLE_DIR, "GSE72094_probe_count_per_gene.csv"))

harmonization_qc <- tibble::tibble(
  item = c("Expression samples", "Expression probes before mapping", "Probes matched to annotation",
           "Probe-to-annotation match rate", "Probes with unambiguous gene symbols",
           "Unique genes after median collapsing"),
  value = c(ncol(expression_log2), nrow(expression_log2), length(common_probes),
            round(probe_match_rate, 6), nrow(probe_expression_mapped), nrow(gene_expression_log2))
)
readr::write_csv(harmonization_qc, file.path(TABLE_DIR, "GSE72094_expression_harmonization_QC.csv"))

cat("\n=====================================================================\n")
cat("GSE72094 EXPRESSION HARMONIZATION SUMMARY\n")
cat("=====================================================================\n")
print(harmonization_qc, n = Inf, width = Inf)

message("\nV2b complete. Harmonized matrix: ", nrow(gene_expression_log2), " genes x ",
        ncol(gene_expression_log2), " samples.")
