# ======================================================================
# M4 — IMMUNE INFILTRATION ANALYSIS (quanTIseq)
#
# GOOGLE COLAB-READY R SCRIPT
#
# PURPOSE
# -------
# Deconvolve bulk expression into immune cell-type fractions using
# quanTIseq, then test whether those fractions differ between the
# FROZEN M3.2 ecosystem states.
#
# IMPORTANT
# ---------
#   - M3.2 ecosystem states are FROZEN. This script does not recluster
#     patients and does not redefine State_1 / State_2.
#   - quanTIseq is run once on the full cohort expression matrix;
#     ecosystem state is used only as a downstream grouping variable
#     for comparison, never as an input to deconvolution.
#
# INPUTS
# -------
# /content/M1_expression_matrix.csv
# /content/M3_LIONESS_entropy/M3.2_module_axes/
#     M3.2_patient_ecosystem_states.csv
#
# OUTPUTS
# --------
# /content/M4_immune_infiltration/
#     M4_quanTIseq_fractions.csv
#     M4_immune_state_comparison.csv
#     M4_immune_boxplots.pdf
#     M4_immune_heatmap_significant.pdf
#
# ======================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

cat("\n")
cat("=====================================================================\n")
cat("M4 — IMMUNE INFILTRATION ANALYSIS (quanTIseq)\n")
cat("=====================================================================\n\n")

cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# 1. PACKAGES
# ======================================================================

cat("Checking required R packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# quantiseqr implements the quanTIseq deconvolution method and is
# distributed via Bioconductor.
if (!requireNamespace("quantiseqr", quietly = TRUE)) {
  cat("Installing: quantiseqr\n")
  BiocManager::install("quantiseqr", ask = FALSE, update = FALSE)
}

required_cran <- c("dplyr", "tidyr", "readr", "tibble", "ggplot2", "pheatmap")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(quantiseqr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

cat("Packages ready.\n\n")


# ======================================================================
# 2. FILE SETTINGS
# ======================================================================

EXPRESSION_FILE <- "/content/M1_expression_matrix.csv"

ANNOTATION_FILE <- "/content/02_ProteinCoding_GeneAnnotation.csv"

STATE_FILE <- "/content/M3_LIONESS_entropy/M3.2_module_axes/M3.2_patient_ecosystem_states.csv"

OUTPUT_DIR <- "/content/M4_immune_infiltration"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

cat("Output directory:", normalizePath(OUTPUT_DIR), "\n\n")


# ======================================================================
# 3. LOAD INPUTS AND MAP GENES TO SYMBOLS
# ======================================================================

if (!file.exists(EXPRESSION_FILE)) {
  stop(paste0("Expression file not found: ", EXPRESSION_FILE))
}

if (!file.exists(ANNOTATION_FILE)) {
  stop(paste0("Gene annotation file not found: ", ANNOTATION_FILE))
}

if (!file.exists(STATE_FILE)) {
  stop(paste0(
    "FROZEN M3.2 ecosystem-state file not found: ", STATE_FILE,
    "\nRun M3.2 before M4."
  ))
}

cat("Loading expression matrix...\n")
expr_df <- read_csv(EXPRESSION_FILE, show_col_types = FALSE)

id_col <- colnames(expr_df)[1]
colnames(expr_df)[1] <- "matrix_gene_id"

expression_samples <- setdiff(colnames(expr_df), "matrix_gene_id")

cat("Loading protein-coding gene annotation...\n")
annotation <- read_csv(ANNOTATION_FILE, show_col_types = FALSE)

if (!("matrix_gene_id" %in% colnames(annotation))) {
  stop("Annotation file does not contain 'matrix_gene_id'.")
}
if (!("gene_symbol" %in% colnames(annotation))) {
  stop("Annotation file does not contain 'gene_symbol'.")
}

annotation$matrix_gene_id <- trimws(as.character(annotation$matrix_gene_id))
annotation$gene_symbol    <- toupper(trimws(as.character(annotation$gene_symbol)))

# ---- Map expression identifiers to gene symbols (mirrors M3) ----
cat("Mapping expression genes to gene symbols...\n")

expr_df$matrix_gene_id <- trimws(as.character(expr_df$matrix_gene_id))

annotation_index <- match(expr_df$matrix_gene_id, annotation$matrix_gene_id)
mapped_symbols   <- annotation$gene_symbol[annotation_index]

matched_count <- sum(!is.na(mapped_symbols) & mapped_symbols != "")
cat("Expression genes:", nrow(expr_df), "\n")
cat("Mapped genes:", matched_count, "\n")

expr_df$gene_symbol <- mapped_symbols

# ---- Remove unmapped genes ----
expr_df <- expr_df[!is.na(expr_df$gene_symbol) & expr_df$gene_symbol != "", , drop = FALSE]

# ---- Remove duplicate gene symbols ----
cat("Removing duplicated gene symbols...\n")
expr_df$gene_symbol <- toupper(trimws(expr_df$gene_symbol))
duplicate_symbols   <- duplicated(expr_df$gene_symbol)
cat("Duplicate genes:", sum(duplicate_symbols), "\n")

expr_df <- expr_df[!duplicate_symbols, , drop = FALSE]
cat("Unique genes:", nrow(expr_df), "\n\n")

# ---- Build final expression matrix, rownames = gene symbols ----
expr_matrix <- as.matrix(expr_df[, expression_samples, drop = FALSE])
rownames(expr_matrix) <- expr_df$gene_symbol
storage.mode(expr_matrix) <- "numeric"

cat("Expression matrix:", nrow(expr_matrix), "genes x",
    ncol(expr_matrix), "samples\n\n")

cat("Loading FROZEN ecosystem-state assignments...\n")
ecosystem_states <- read_csv(STATE_FILE, show_col_types = FALSE)

if (!"sample_id" %in% colnames(ecosystem_states) ||
    !"ecosystem_state" %in% colnames(ecosystem_states)) {
  stop("M3.2 state file must contain 'sample_id' and 'ecosystem_state' columns.")
}

kv <- function(label, value) cat(label, ":", value, "\n")

kv("Patients with frozen ecosystem state", nrow(ecosystem_states))
print(table(ecosystem_states$ecosystem_state))
cat("\n")

# ======================================================================
# 4. RUN quanTIseq DECONVOLUTION
# ======================================================================

cat("Running quanTIseq deconvolution (TIL10 signature)...\n")

quantiseq_result <- quantiseqr::run_quantiseq(
  expression_data = expr_matrix,
  signature_matrix = "TIL10",
  is_arraydata     = FALSE,
  is_tumordata     = TRUE,
  scale_mRNA       = TRUE
)

# run_quantiseq() returns a data frame with one row per sample and one
# column per cell type (plus a "Sample" identifier column).
quantiseq_result <- as.data.frame(quantiseq_result)

sample_col <- intersect(c("Sample", "sample", "sample_id"),
                         colnames(quantiseq_result))[1]

if (is.na(sample_col)) {
  stop("Could not identify the sample-identifier column returned by run_quantiseq().")
}

quantiseq_result <- quantiseq_result %>%
  rename(sample_id = all_of(sample_col))

cell_type_cols <- setdiff(colnames(quantiseq_result), "sample_id")

kv("quanTIseq cell types scored", length(cell_type_cols))
cat(paste(cell_type_cols, collapse = ", "), "\n\n")

write_csv(
  quantiseq_result,
  file.path(OUTPUT_DIR, "M4_quanTIseq_fractions.csv")
)


# ======================================================================
# 5. MERGE WITH FROZEN ECOSYSTEM STATE (ID-FORMAT NORMALIZED)
# ======================================================================

# TCGA/GEO sample identifiers frequently pick up cosmetic differences
# between files (e.g. "-" vs "." from R's automatic column-name
# mangling, or trailing suffixes). Join on a normalized key rather
# than the raw string, then keep the ORIGINAL ecosystem_states
# sample_id going forward.


normalize_id <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("[.\\-]", "_", x)      # collapse "-" and "." to a common "_"
  x <- sub("_LIONESS$", "", x)     # strip the "_LIONESS" suffix added to state file IDs
  x
}

quantiseq_result <- quantiseq_result %>%
  mutate(.join_key = normalize_id(sample_id))

ecosystem_states_keyed <- ecosystem_states %>%
  mutate(.join_key = normalize_id(sample_id))

merged <- ecosystem_states_keyed %>%
  select(.join_key, sample_id, ecosystem_state) %>%
  inner_join(
    quantiseq_result %>% select(-sample_id),
    by = ".join_key"
  ) %>%
  select(-.join_key)

kv("Patients with both quanTIseq fractions and frozen state", nrow(merged))

if (nrow(merged) == 0) {
  cat("\nStill 0 matches after normalization. Example IDs:\n")
  cat("  quanTIseq (normalized):",
      paste(head(normalize_id(quantiseq_result$sample_id), 5), collapse = " | "), "\n")
  cat("  states (normalized):   ",
      paste(head(normalize_id(ecosystem_states$sample_id), 5), collapse = " | "), "\n")
  stop("No overlapping sample IDs between quanTIseq output and frozen ecosystem states. ",
       "Compare the printed ID formats above and adjust normalize_id() accordingly ",
       "(e.g. truncate TCGA barcodes to 15/12 characters if one file has extra suffix fields).")
}

if (n_distinct(merged$ecosystem_state) != 2) {
  stop("Expected exactly two frozen ecosystem states (State_1, State_2).")
}

merged$ecosystem_state <- factor(merged$ecosystem_state)

state_levels <- levels(merged$ecosystem_state)


# ======================================================================
# 6. DIFFERENTIAL IMMUNE-FRACTION TESTING (STATE_1 vs STATE_2)
# ======================================================================

cat("\nTesting immune-cell fractions by frozen ecosystem state...\n")

rank_biserial <- function(x, group) {
  # Rank-biserial correlation effect size for a two-group Wilcoxon test.
  w <- suppressWarnings(wilcox.test(x ~ group)$statistic)
  n1 <- sum(group == levels(group)[1])
  n2 <- sum(group == levels(group)[2])
  1 - (2 * w) / (n1 * n2)
}

immune_results <- lapply(cell_type_cols, function(ct) {

  values <- merged[[ct]]
  group  <- merged$ecosystem_state

  test <- tryCatch(
    wilcox.test(values ~ group),
    error = function(e) NULL
  )

  medians <- tapply(values, group, median, na.rm = TRUE)

  data.frame(
    cell_type    = ct,
    median_state1 = medians[state_levels[1]],
    median_state2 = medians[state_levels[2]],
    p_value      = if (is.null(test)) NA_real_ else test$p.value,
    effect_size_rb = if (is.null(test)) NA_real_ else rank_biserial(values, group),
    stringsAsFactors = FALSE
  )
})

immune_results_df <- bind_rows(immune_results)

immune_results_df$FDR <- p.adjust(immune_results_df$p_value, method = "BH")

immune_results_df <- immune_results_df %>%
  arrange(FDR, p_value)

print(immune_results_df)

write_csv(
  immune_results_df,
  file.path(OUTPUT_DIR, "M4_immune_state_comparison.csv")
)

n_significant <- sum(immune_results_df$FDR < 0.05, na.rm = TRUE)
kv("Cell types significant at FDR < 0.05", n_significant)


# ======================================================================
# 7. FIGURE — BOXPLOTS OF ALL IMMUNE-CELL FRACTIONS BY STATE
# ======================================================================

cat("\nGenerating immune-fraction boxplots...\n")

long_immune <- merged %>%
  select(sample_id, ecosystem_state, all_of(cell_type_cols)) %>%
  pivot_longer(
    cols = all_of(cell_type_cols),
    names_to = "cell_type",
    values_to = "fraction"
  )

boxplot_immune <- ggplot(
  long_immune,
  aes(x = ecosystem_state, y = fraction, fill = ecosystem_state)
) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(
    title = "quanTIseq Immune-Cell Fractions by Frozen Ecosystem State",
    x = "Ecosystem state",
    y = "Estimated fraction"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  file.path(OUTPUT_DIR, "M4_immune_boxplots.pdf"),
  boxplot_immune,
  width = 12, height = 9
)


# ======================================================================
# 8. FIGURE — HEATMAP OF FDR-SIGNIFICANT IMMUNE POPULATIONS
# ======================================================================

significant_cell_types <- immune_results_df %>%
  filter(FDR < 0.05) %>%
  pull(cell_type)

if (length(significant_cell_types) >= 2) {

  cat("Generating heatmap of FDR-significant immune populations...\n")

  heatmap_matrix <- merged %>%
    select(sample_id, all_of(significant_cell_types)) %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>%
    scale()

  annotation_df <- merged %>%
    select(sample_id, ecosystem_state) %>%
    column_to_rownames("sample_id")

  pdf(file.path(OUTPUT_DIR, "M4_immune_heatmap_significant.pdf"),
      width = 10, height = 8)

  pheatmap(
    t(heatmap_matrix),
    annotation_col = annotation_df,
    show_colnames = FALSE,
    main = "FDR-Significant Immune Populations by Ecosystem State"
  )

  dev.off()

} else {

  cat("Fewer than two FDR-significant immune populations; ",
      "heatmap skipped.\n")
}

cat("\nM4 immune infiltration analysis complete.\n")
cat("Finished:", as.character(Sys.time()), "\n")
