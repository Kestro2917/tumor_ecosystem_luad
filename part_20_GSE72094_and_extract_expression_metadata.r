# ======================================================================
# V2a — DOWNLOAD GSE72094 AND EXTRACT EXPRESSION / METADATA
#
# PURPOSE
# -------
# Downloads the GSE72094 Series Matrix directly from GEO (no manually
# prepared file needed) and extracts the raw expression matrix, sample
# metadata, and probe annotation. This is adapted from the data-
# processing approach in your own External_Validation_GSE72094.ipynb
# (Section 9, Cell 1) — same GEOquery workflow, same scale-check logic
# — reused here rather than re-invented, since it was already tested
# on this exact GEO series.
#
# OUTPUT
# ------
# /content/GSE72094_raw/clean/GSE72094_expression_log2_matrix.csv
# /content/GSE72094_raw/clean/GSE72094_probe_annotation_raw.csv
# /content/GSE72094_raw/clean/GSE72094_sample_characteristics_long.csv
# /content/GSE72094_raw/tables/GSE72094_input_qc_summary.csv
# ======================================================================

options(stringsAsFactors = FALSE)

cran_packages <- c("dplyr", "tidyr", "tibble", "readr", "stringr", "purrr")
missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioconductor_packages <- c("GEOquery", "Biobase")
missing_bioc <- setdiff(bioconductor_packages, rownames(installed.packages()))
if (length(missing_bioc) > 0L) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
})

GEO_ACCESSION <- "GSE72094"
DATA_DIR <- "/content"
OUT_DIR   <- file.path(DATA_DIR, "GSE72094_raw")
RAW_DIR   <- file.path(OUT_DIR, "raw")
CLEAN_DIR <- file.path(OUT_DIR, "clean")
TABLE_DIR <- file.path(OUT_DIR, "tables")

for (d in c(OUT_DIR, RAW_DIR, CLEAN_DIR, TABLE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("Downloading ", GEO_ACCESSION, " Series Matrix from GEO...")

gse_list <- GEOquery::getGEO(GEO_ACCESSION, GSEMatrix = TRUE, getGPL = TRUE, destdir = RAW_DIR)

if (!is.list(gse_list) || length(gse_list) == 0L) {
  stop("getGEO() did not return any ExpressionSet objects for ", GEO_ACCESSION, ".")
}

message("Number of platform-specific ExpressionSet objects: ", length(gse_list))

# ---- select the ExpressionSet with the most samples, in case of >1 platform ----
platform_summary <- purrr::map_dfr(seq_along(gse_list), function(i) {
  eset_i <- gse_list[[i]]
  tibble::tibble(
    list_index = i,
    platform = Biobase::annotation(eset_i),
    n_features = nrow(Biobase::exprs(eset_i)),
    n_samples = ncol(Biobase::exprs(eset_i))
  )
})
print(platform_summary)

selected_index <- platform_summary$list_index[which.max(platform_summary$n_samples)]
gse_eset <- gse_list[[selected_index]]
message("Selected platform: ", Biobase::annotation(gse_eset))

# ---- extract expression matrix ----
expression_matrix <- as.matrix(Biobase::exprs(gse_eset))
storage.mode(expression_matrix) <- "double"

if (nrow(expression_matrix) == 0L || ncol(expression_matrix) == 0L) {
  stop("The downloaded expression matrix is empty.")
}
if (anyDuplicated(colnames(expression_matrix))) {
  stop("Duplicated GEO sample identifiers were detected.")
}

# ---- extract sample metadata ----
sample_metadata <- Biobase::pData(gse_eset) %>%
  tibble::rownames_to_column("geo_sample_id") %>%
  tibble::as_tibble()

if (!all(colnames(expression_matrix) %in% sample_metadata$geo_sample_id)) {
  stop("Expression sample identifiers do not match metadata.")
}
sample_metadata <- sample_metadata[
  match(colnames(expression_matrix), sample_metadata$geo_sample_id), , drop = FALSE
]

# ---- extract probe/feature annotation ----
feature_annotation <- Biobase::fData(gse_eset) %>%
  tibble::rownames_to_column("probe_id") %>%
  tibble::as_tibble()

if (nrow(feature_annotation) == 0L) {
  warning("No feature annotation was present in the Series Matrix.")
}

# ---- SCALE CHECK: only log2-transform if values look linear, never twice ----
maximum_expression <- max(expression_matrix, na.rm = TRUE)
minimum_expression <- min(expression_matrix, na.rm = TRUE)

cat("\nExpression quantiles:\n")
print(stats::quantile(expression_matrix, probs = c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

appears_log2 <- (maximum_expression < 50 && minimum_expression > -20)

if (appears_log2) {
  expression_log2 <- expression_matrix
  transformation_status <- "Series Matrix values retained; data appeared log2-scaled."
} else {
  if (any(expression_matrix < 0, na.rm = TRUE)) {
    stop("Expression values do not appear log2-scaled, but negative values ",
         "prevent log2 transformation.")
  }
  expression_log2 <- log2(expression_matrix + 1)
  transformation_status <- "Series Matrix values transformed using log2(x + 1)."
}
message(transformation_status)

# ---- parse characteristics_ch1 columns into a long clinical/phenotype table ----
characteristic_columns <- grep("^characteristics_ch1", colnames(sample_metadata), value = TRUE)

if (length(characteristic_columns) > 0L) {
  characteristics_long <- sample_metadata %>%
    dplyr::select(geo_sample_id, dplyr::all_of(characteristic_columns)) %>%
    tidyr::pivot_longer(cols = -geo_sample_id, names_to = "source_column", values_to = "characteristic") %>%
    dplyr::filter(!is.na(characteristic), nzchar(trimws(as.character(characteristic)))) %>%
    dplyr::mutate(
      characteristic = trimws(as.character(characteristic)),
      characteristic_name = dplyr::if_else(
        stringr::str_detect(characteristic, ":"),
        stringr::str_trim(stringr::str_extract(characteristic, "^[^:]+")),
        source_column
      ),
      characteristic_value = dplyr::if_else(
        stringr::str_detect(characteristic, ":"),
        stringr::str_trim(stringr::str_replace(characteristic, "^[^:]+:", "")),
        characteristic
      )
    )
} else {
  characteristics_long <- tibble::tibble(
    geo_sample_id = character(), source_column = character(),
    characteristic = character(), characteristic_name = character(),
    characteristic_value = character()
  )
  warning("No characteristics_ch1 columns were found.")
}

# ---- save everything ----
saveRDS(gse_eset, file.path(RAW_DIR, "GSE72094_selected_ExpressionSet.rds"))
saveRDS(expression_log2, file.path(CLEAN_DIR, "GSE72094_expression_log2_matrix.rds"))

readr::write_csv(
  as.data.frame(expression_log2) %>% tibble::rownames_to_column("probe_id"),
  file.path(CLEAN_DIR, "GSE72094_expression_log2_matrix.csv")
)
readr::write_csv(sample_metadata, file.path(CLEAN_DIR, "GSE72094_sample_metadata_raw.csv"))
readr::write_csv(feature_annotation, file.path(CLEAN_DIR, "GSE72094_probe_annotation_raw.csv"))
readr::write_csv(characteristics_long, file.path(CLEAN_DIR, "GSE72094_sample_characteristics_long.csv"))
readr::write_csv(platform_summary, file.path(TABLE_DIR, "GSE72094_platform_summary.csv"))

input_qc <- tibble::tibble(
  item = c("GEO accession", "Selected platform", "Number of probes", "Number of samples",
           "Minimum expression value", "Maximum expression value", "Expression transformation",
           "Number of metadata columns", "Number of characteristic columns"),
  value = c(GEO_ACCESSION, Biobase::annotation(gse_eset), nrow(expression_log2), ncol(expression_log2),
            minimum_expression, maximum_expression, transformation_status,
            ncol(sample_metadata), length(characteristic_columns))
)
readr::write_csv(input_qc, file.path(TABLE_DIR, "GSE72094_input_qc_summary.csv"))

cat("\n=====================================================================\n")
cat("GSE72094 DOWNLOAD AND EXTRACTION SUMMARY\n")
cat("=====================================================================\n")
print(input_qc)

cat("\nUnique characteristic names available for clinical/survival lookup (used in V6):\n")
if (nrow(characteristics_long) > 0L) {
  print(sort(unique(characteristics_long$characteristic_name)))
}

message("\nV2a complete. Output directory: ", OUT_DIR)
