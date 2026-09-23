# ======================================================================
# M2b Gene-signature library construction and enrichment scoring
# Expression:
#   M1_expression_matrix.csv
# Clinical:
#   M1_clinical.csv
# Gene annotation:
#   02_ProteinCoding_GeneAnnotation.csv
# MSigDB signature libraries (from the M2 MSigDB Signature Builder):
#   M2_MSigDB_LUAD/M2_CellType_MSigDB.csv
#   M2_MSigDB_LUAD/M2_CellState_MSigDB.csv
#   M2_MSigDB_LUAD/M2_TumorStroma_MSigDB.csv
#
# Outputs:
#   M2b_CellType_EnrichmentScores.csv
#   M2b_CellState_EnrichmentScores.csv
#   M2b_TumorStroma_EnrichmentScores.csv
#   M2b_CellType_signature_overlap.csv
#   M2b_CellState_signature_overlap.csv
#   M2b_TumorStroma_signature_overlap.csv
#   M2b_sample_matching.csv
#   M2b_run_summary.csv
#
# ======================================================================


# ======================================================================
# 0. FILE SETTINGS
# ======================================================================

EXPRESSION_FILE <- "/content/M1_expression_matrix.csv"
CLINICAL_FILE <- "/content/M1_clinical.csv"
ANNOTATION_FILE <- "/content/02_ProteinCoding_GeneAnnotation.csv"

MSIGDB_DIR <- "/content/M2_MSigDB_LUAD"

LIBRARY_FILES <- list(
  CellType    = file.path(MSIGDB_DIR, "M2_CellType_MSigDB.csv"),
  CellState   = file.path(MSIGDB_DIR, "M2_CellState_MSigDB.csv"),
  TumorStroma = file.path(MSIGDB_DIR, "M2_TumorStroma_MSigDB.csv")
)

MATCHING_FILE <- "/content/M2b_sample_matching.csv"
SUMMARY_FILE <- "/content/M2b_run_summary.csv"

# Minimum number of marker genes required, per signature, for that
# signature to be scored (same convention as the original M2 script).
MIN_MARKERS <- 3


# ======================================================================
# 1. START
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("M2b - MSigDB-BASED ssGSEA ENRICHMENT SCORING\n")
cat("============================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# 2. INSTALL / LOAD PACKAGES
# ======================================================================

cat("Checking required packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GSVA", quietly = TRUE)) {
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
}

cran_packages <- c("data.table", "dplyr")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GSVA)
})

cat("Packages ready.\n")


# ======================================================================
# 3. CHECK INPUT FILES
# ======================================================================

cat("\nChecking input files...\n")

core_required_files <- c(EXPRESSION_FILE, CLINICAL_FILE, ANNOTATION_FILE)
missing_core <- core_required_files[!file.exists(core_required_files)]

if (length(missing_core) > 0) {
  cat("\nMissing files:\n")
  print(missing_core)
  stop("One or more required core files are missing.")
}

missing_libraries <- names(LIBRARY_FILES)[!file.exists(unlist(LIBRARY_FILES))]

if (length(missing_libraries) > 0) {
  cat("\nMissing MSigDB library files:\n")
  print(unlist(LIBRARY_FILES[missing_libraries]))
  stop("Run the M2 MSigDB Signature Builder script first - one or more library CSVs are missing.")
}

cat("All required files found.\n")


# ======================================================================
# 4. READ + PREPARE EXPRESSION MATRIX
#    (same logic as the original M2 ssGSEA script - unchanged)
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("READING EXPRESSION MATRIX\n")
cat("============================================================\n")

expr_df <- fread(EXPRESSION_FILE, data.table = FALSE, check.names = FALSE)

cat("Rows:", nrow(expr_df), "\n")
cat("Columns:", ncol(expr_df), "\n")

gene_column <- "matrix_gene_id"

if (!(gene_column %in% colnames(expr_df))) {
  stop("matrix_gene_id column was not found in expression matrix.")
}

expression_samples <- setdiff(colnames(expr_df), gene_column)

cat("Number of expression samples:", length(expression_samples), "\n")

for (sample_name in expression_samples) {
  expr_df[, sample_name] <- suppressWarnings(as.numeric(expr_df[, sample_name]))
}

cat("\n")
cat("============================================================\n")
cat("READING GENE ANNOTATION\n")
cat("============================================================\n")

annotation <- fread(ANNOTATION_FILE, data.table = FALSE, check.names = FALSE)

if (!("matrix_gene_id" %in% colnames(annotation))) {
  stop("matrix_gene_id is missing from annotation file.")
}
if (!("gene_symbol" %in% colnames(annotation))) {
  stop("gene_symbol is missing from annotation file.")
}

annotation$matrix_gene_id <- trimws(as.character(annotation$matrix_gene_id))
annotation$gene_symbol <- trimws(as.character(annotation$gene_symbol))

annotation_index <- match(as.character(expr_df$matrix_gene_id), annotation$matrix_gene_id)
gene_symbols <- annotation$gene_symbol[annotation_index]

expr_df$gene_symbol <- gene_symbols
expr_df <- expr_df[!is.na(expr_df$gene_symbol) & expr_df$gene_symbol != "", ]

cat("Genes remaining after annotation:", nrow(expr_df), "\n")

duplicate_symbols <- duplicated(expr_df$gene_symbol)
cat("Duplicate gene symbols:", sum(duplicate_symbols), "\n")

expr_df <- expr_df[!duplicate_symbols, ]
cat("Unique gene symbols:", nrow(expr_df), "\n")

expr_matrix <- as.matrix(expr_df[, expression_samples, drop = FALSE])
rownames(expr_matrix) <- toupper(expr_df$gene_symbol)
storage.mode(expr_matrix) <- "numeric"

expr_matrix[is.na(expr_matrix)] <- 0
expr_matrix[is.infinite(expr_matrix)] <- 0

cat("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")


# ======================================================================
# 5. READ CLINICAL DATA + SAMPLE MATCHING
#    (same logic as the original M2 ssGSEA script - unchanged)
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("SAMPLE MATCHING\n")
cat("============================================================\n")

clinical <- fread(CLINICAL_FILE, data.table = FALSE, check.names = FALSE)

if (!("sample_barcode" %in% colnames(clinical))) {
  stop("sample_barcode column was not found in clinical file.")
}

clinical_samples <- trimws(as.character(clinical$sample_barcode))
expression_samples <- colnames(expr_matrix)

common_samples <- intersect(expression_samples, clinical_samples)

cat("Expression samples:", length(expression_samples), "\n")
cat("Clinical samples:", length(clinical_samples), "\n")
cat("Matched samples:", length(common_samples), "\n")

if (length(common_samples) < 10) {
  stop("Fewer than 10 samples matched.")
}

matching_table <- data.frame(
  expression_sample_id = expression_samples,
  matched = expression_samples %in% common_samples,
  stringsAsFactors = FALSE
)
matching_table$clinical_sample_id <- clinical_samples[
  match(matching_table$expression_sample_id, clinical_samples)
]

write.csv(matching_table, MATCHING_FILE, row.names = FALSE, quote = FALSE)
cat("Saved:", MATCHING_FILE, "\n")

matched_samples <- expression_samples[expression_samples %in% common_samples]
expr_matrix <- expr_matrix[, matched_samples, drop = FALSE]

cat("Final expression samples:", ncol(expr_matrix), "\n")
cat("Sample matching completed successfully.\n")


# ======================================================================
# 6. PREPARE EXPRESSION FOR ssGSEA (same log2 heuristic as original)
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("EXPRESSION SCALE CHECK\n")
cat("============================================================\n")

minimum_value <- min(expr_matrix, na.rm = TRUE)
maximum_value <- max(expr_matrix, na.rm = TRUE)

cat("Minimum:", minimum_value, "\n")
cat("Maximum:", maximum_value, "\n")

if (maximum_value > 50) {
  cat("\nLarge expression values detected. Applying log2(x + 1)...\n")
  expr_ssgsea <- log2(expr_matrix + 1)
} else {
  cat("\nExpression values used as supplied.\n")
  expr_ssgsea <- expr_matrix
}


# ======================================================================
# 7. HELPER — RUN ssGSEA (modern GSVA API, same pattern as original)
# ======================================================================

run_ssgsea <- function(expr, gene_sets, min_size) {
  result <- NULL

  tryCatch(
    {
      cat("\nTrying modern GSVA API...\n")

      ssgsea_parameters <- GSVA::ssgseaParam(
        exprData = expr,
        geneSets = gene_sets,
        minSize = min_size,
        maxSize = Inf,
        normalize = TRUE
      )

      result <- GSVA::gsva(ssgsea_parameters, verbose = TRUE)
    },
    error = function(e) {
      cat("\nModern GSVA API failed.\n")
      cat("Error:", conditionMessage(e), "\n")
    }
  )

  result
}


# ======================================================================
# 8. HELPER — SCORE ONE MSigDB LIBRARY END TO END
#
# For a given library CSV (signature, category, gene, source,
# collection):
#   1. Build one gene set per SIGNATURE (not per category), filtered
#      to genes present in the expression matrix.
#   2. Drop signatures with fewer than MIN_MARKERS overlapping genes,
#      reporting the same kind of overlap table the original script
#      produced per cell type.
#   3. Run ssGSEA -> raw signature-level scores.
#   4. Z-score each signature across samples.
#   5. Aggregate to CATEGORY-level scores by averaging member
#      signature z-scores (mirrors the "domain score" convention used
#      later in this pipeline for combining z-scored components).
#   6. Save one wide, sample-level CSV: sample_id, sig_<name>_raw,
#      sig_<name>_z, score_<category>_z.
# ======================================================================

score_library <- function(library_name, library_path, expr_ssgsea,
                           expression_genes, min_markers, outdir = "/content") {

  cat("\n")
  cat("============================================================\n")
  cat("SCORING LIBRARY:", library_name, "\n")
  cat("============================================================\n")

  lib <- fread(library_path, data.table = FALSE, check.names = FALSE)

  required_cols <- c("signature", "category", "gene")
  missing_cols <- setdiff(required_cols, colnames(lib))
  if (length(missing_cols) > 0) {
    stop(paste0(
      "Library file ", library_path,
      " is missing required column(s): ", paste(missing_cols, collapse = ", ")
    ))
  }

  lib$signature <- trimws(as.character(lib$signature))
  lib$category <- trimws(as.character(lib$category))
  lib$gene <- toupper(trimws(as.character(lib$gene)))

  cat("Signatures in library:", n_distinct(lib$signature), "\n")
  cat("Categories in library:", n_distinct(lib$category), "\n")

  # ---- Build one gene set per signature, restricted to expression genes
  signature_names <- unique(lib$signature)
  gene_sets <- list()
  overlap_results <- data.frame(
    signature = character(0), category = character(0),
    genes_original = integer(0), genes_in_expression = integer(0),
    overlap_fraction = numeric(0), stringsAsFactors = FALSE
  )

  for (sig in signature_names) {
    sig_genes <- unique(lib$gene[lib$signature == sig])
    sig_genes <- sig_genes[!is.na(sig_genes) & sig_genes != ""]

    genes_present <- intersect(sig_genes, expression_genes)
    gene_sets[[sig]] <- genes_present

    overlap_fraction <- if (length(sig_genes) > 0) {
      length(genes_present) / length(sig_genes)
    } else {
      0
    }

    overlap_results <- rbind(overlap_results, data.frame(
      signature = sig,
      category = lib$category[lib$signature == sig][1],
      genes_original = length(sig_genes),
      genes_in_expression = length(genes_present),
      overlap_fraction = overlap_fraction,
      stringsAsFactors = FALSE
    ))
  }

  overlap_file <- file.path(outdir, paste0("M2b_", library_name, "_signature_overlap.csv"))
  write.csv(overlap_results, overlap_file, row.names = FALSE, quote = FALSE)
  cat("Saved:", overlap_file, "\n")

  # ---- Drop signatures below the minimum marker threshold
  usable_signatures <- names(gene_sets)[lengths(gene_sets) >= min_markers]
  dropped_signatures <- setdiff(names(gene_sets), usable_signatures)

  if (length(dropped_signatures) > 0) {
    cat("\nDropping", length(dropped_signatures),
        "signature(s) with fewer than", min_markers, "genes in expression data.\n")
  }

  gene_sets <- gene_sets[usable_signatures]

  if (length(gene_sets) == 0) {
    stop(paste0("No usable signatures remained for library: ", library_name))
  }

  cat("Usable signatures for ssGSEA:", length(gene_sets), "\n")

  # ---- Run ssGSEA
  ssgsea_result <- run_ssgsea(expr_ssgsea, gene_sets, min_markers)

  if (is.null(ssgsea_result)) {
    stop(paste0("ssGSEA failed for library: ", library_name))
  }

  score_matrix <- as.matrix(ssgsea_result)
  cat("ssGSEA result:", nrow(score_matrix), "signatures x", ncol(score_matrix), "samples\n")

  # ---- Transpose to sample-level table, keep raw + z per signature
  score_table <- as.data.frame(t(score_matrix), check.names = FALSE, stringsAsFactors = FALSE)
  score_table$sample_id <- rownames(score_table)
  rownames(score_table) <- NULL

  signature_cols_raw <- character(0)
  signature_cols_z <- character(0)
  sig_to_category <- setNames(overlap_results$category, overlap_results$signature)

  for (sig in rownames(score_matrix)) {
    raw_col <- paste0("sig_", sig, "_raw")
    z_col <- paste0("sig_", sig, "_z")

    values <- as.numeric(score_table[[sig]])
    values_sd <- sd(values, na.rm = TRUE)

    score_table[[raw_col]] <- values
    score_table[[z_col]] <- if (is.na(values_sd) || values_sd == 0) {
      0
    } else {
      as.numeric(scale(values))
    }

    score_table[[sig]] <- NULL
    signature_cols_raw <- c(signature_cols_raw, raw_col)
    signature_cols_z <- c(signature_cols_z, z_col)
  }

  # ---- Category-level aggregate score = mean of member signature z-scores
  categories_present <- unique(sig_to_category[rownames(score_matrix)])

  for (cat_name in categories_present) {
    member_signatures <- names(sig_to_category)[sig_to_category == cat_name]
    member_signatures <- intersect(member_signatures, rownames(score_matrix))
    member_z_cols <- paste0("sig_", member_signatures, "_z")
    member_z_cols <- intersect(member_z_cols, colnames(score_table))

    category_col <- paste0("score_", cat_name, "_z")

    if (length(member_z_cols) > 0) {
      score_table[[category_col]] <- rowMeans(
        as.matrix(score_table[, member_z_cols, drop = FALSE]), na.rm = TRUE
      )
    } else {
      score_table[[category_col]] <- NA_real_
    }
  }

  category_cols <- paste0("score_", categories_present, "_z")

  final_cols <- c("sample_id", signature_cols_raw, signature_cols_z, category_cols)
  score_table <- score_table[, final_cols, drop = FALSE]

  output_file <- file.path(outdir, paste0("M2b_", library_name, "_EnrichmentScores.csv"))
  write.csv(score_table, output_file, row.names = FALSE, quote = FALSE)
  cat("Saved:", output_file, "\n")

  list(
    library_name = library_name,
    n_signatures_scored = length(gene_sets),
    n_signatures_dropped = length(dropped_signatures),
    n_categories = length(categories_present),
    n_samples = nrow(score_table),
    output_file = output_file,
    overlap_file = overlap_file
  )
}


# ======================================================================
# 9. RUN SCORING FOR EACH MSigDB LIBRARY
# ======================================================================

expression_genes <- rownames(expr_ssgsea)

run_results <- list()

for (library_name in names(LIBRARY_FILES)) {
  run_results[[library_name]] <- score_library(
    library_name = library_name,
    library_path = LIBRARY_FILES[[library_name]],
    expr_ssgsea = expr_ssgsea,
    expression_genes = expression_genes,
    min_markers = MIN_MARKERS
  )
}


# ======================================================================
# 10. RUN SUMMARY
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("M2b RUN SUMMARY\n")
cat("============================================================\n")

summary_table <- do.call(rbind, lapply(run_results, function(r) {
  data.frame(
    library = r$library_name,
    signatures_scored = r$n_signatures_scored,
    signatures_dropped = r$n_signatures_dropped,
    categories = r$n_categories,
    samples = r$n_samples,
    output_file = r$output_file,
    stringsAsFactors = FALSE
  )
}))
rownames(summary_table) <- NULL

print(summary_table, row.names = FALSE)

write.csv(summary_table, SUMMARY_FILE, row.names = FALSE, quote = FALSE)
cat("\nSaved:", SUMMARY_FILE, "\n")


# ======================================================================
# 11. FINISHED
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("M2b FINISHED SUCCESSFULLY\n")
cat("============================================================\n")
cat("Finished:", as.character(Sys.time()), "\n")
cat("\n")
cat("NOTE: this enrichment layer is independent of the original\n")
cat("marker-based M2_celltype_scores.csv. Downstream scripts that\n")
cat("expect score_ClMAC_z / score_AltMAC_z / etc. should keep using\n")
cat("that file; this one adds broader MSigDB-standard signature and\n")
cat("category scores alongside it, distinguished by the 'sig_' and\n")
cat("category-level 'score_<Category>_z' naming used here.\n")
cat("============================================================\n")