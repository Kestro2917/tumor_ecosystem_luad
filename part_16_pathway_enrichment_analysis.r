# ======================================================================
# M7 — PATHWAY ENRICHMENT ANALYSIS (GSVA), STATE_1 vs STATE_2
# PURPOSE
# -------
# Score canonical pathways per patient using GSVA and test which
# pathways differ between the FROZEN M3.2 ecosystem states.
#
# GENE-SET LIBRARY
# ------------------
# Two sources are combined:
#   1. The full MSigDB Hallmark collection (H, 50 gene sets), covering
#      broad cancer-relevant programs (interferon signaling, EMT,
#      hypoxia, KRAS signaling, oxidative phosphorylation, cell-cycle/
#      G2M and E2F targets, inflammatory response, etc.).
#   2. A small set of curated antigen-presentation gene sets from
#      MSigDB C2/C5 (keyword: "ANTIGEN_PROCESSING"), since antigen
#      presentation is not represented as its own Hallmark set.
#
#
# METHOD
# ------
# 1. Load expression matrix and map matrix_gene_id -> gene_symbol.
# 2. Build the combined gene-set list via msigdbr (version-tolerant).
# 3. Run GSVA on log2(x + 1)-transformed expression to obtain a
#    pathway x patient enrichment-score matrix.
# 4. Merge with the FROZEN M3.2 ecosystem-state assignment.
# 5. Wilcoxon rank-sum test per pathway (State_1 vs State_2), BH-FDR
#    across all pathways tested.
# 6. Figures: heatmap of top FDR-significant pathways ordered by
#    state, and a dot plot summarizing effect direction/magnitude and
#    significance for the top pathways.
#
# INPUTS
# -------
# /content/M1_expression_matrix.csv
# /content/02_ProteinCoding_GeneAnnotation.csv
# /content/M3_LIONESS_entropy/M3.2_module_axes/
#     M3.2_patient_ecosystem_states.csv
#
# OUTPUTS
# --------
# /content/M7_pathway_enrichment/
#     M7_GSVA_pathway_scores.csv
#     M7_pathway_state_comparison.csv
#     M7_pathway_enrichment_heatmap.pdf
#     M7_pathway_dotplot.pdf
#
# ======================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

cat("\n")
cat("=====================================================================\n")
cat("M7 — PATHWAY ENRICHMENT ANALYSIS (GSVA), STATE_1 vs STATE_2\n")
cat("=====================================================================\n\n")

cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# 1. PACKAGES
# ======================================================================

cat("Checking required R packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GSVA", quietly = TRUE)) {
  cat("Installing: GSVA\n")
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
}

required_cran <- c("msigdbr", "dplyr", "tidyr", "readr", "tibble",
                    "ggplot2", "pheatmap")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
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

OUTPUT_DIR <- "/content/M7_pathway_enrichment"

FDR_THRESHOLD <- 0.05
TOP_N_PATHWAYS_FOR_FIGURES <- 30

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

cat("Output directory:", normalizePath(OUTPUT_DIR), "\n\n")

kv <- function(label, value) cat(label, ":", value, "\n")

# Version-tolerant msigdbr wrapper: tries the modern `collection`
# argument first (msigdbr >= 10.0.0), and falls back to the older
# `category` argument if that fails (older msigdbr releases).
get_msigdb_sets <- function(species = "Homo sapiens", set, subset = NULL) {
  result <- tryCatch(
    {
      if (is.null(subset)) {
        msigdbr(species = species, collection = set)
      } else {
        msigdbr(species = species, collection = set, subcollection = subset)
      }
    },
    error = function(e) NULL
  )

  if (is.null(result) || nrow(result) == 0) {
    result <- tryCatch(
      {
        if (is.null(subset)) {
          msigdbr(species = species, category = set)
        } else {
          msigdbr(species = species, category = set, subcategory = subset)
        }
      },
      error = function(e) NULL
    )
  }

  if (is.null(result)) {
    result <- tibble()
  }

  result
}


# ======================================================================
# 3. LOAD EXPRESSION MATRIX AND MAP GENES TO SYMBOLS
# ======================================================================

if (!file.exists(EXPRESSION_FILE)) {
  stop(paste0("Expression file not found: ", EXPRESSION_FILE))
}

if (!file.exists(ANNOTATION_FILE)) {
  stop(paste0("Gene annotation file not found: ", ANNOTATION_FILE))
}

if (!file.exists(STATE_FILE)) {
  stop(paste0("FROZEN M3.2 ecosystem-state file not found: ", STATE_FILE,
              "\nRun M3.2 before M7."))
}

cat("Loading expression matrix...\n")
expr_df <- read_csv(EXPRESSION_FILE, show_col_types = FALSE)
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

# ---- Map expression identifiers to gene symbols ----
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

# Log2 transform if the matrix appears to be on a linear (TPM-like) scale.
if (max(expr_matrix, na.rm = TRUE) > 50) {
  cat("Applying log2(x + 1) transformation prior to GSVA.\n")
  expr_matrix <- log2(expr_matrix + 1)
}

kv("Expression matrix (genes x samples)",
   paste(nrow(expr_matrix), "x", ncol(expr_matrix)))

ecosystem_states <- read_csv(STATE_FILE, show_col_types = FALSE)

if (!"sample_id" %in% colnames(ecosystem_states)) {
  colnames(ecosystem_states)[1] <- "sample_id"
}

ecosystem_states <- ecosystem_states %>%
  select(sample_id, ecosystem_state)

kv("Patients with frozen ecosystem state", nrow(ecosystem_states))


# ======================================================================
# 4. BUILD GENE-SET LIBRARY (HALLMARK + ANTIGEN PROCESSING)
# ======================================================================

cat("\nBuilding pathway gene-set library from MSigDB...\n")

hallmark_sets <- get_msigdb_sets(set = "H")

antigen_sets <- bind_rows(
  get_msigdb_sets(set = "C2"),
  get_msigdb_sets(set = "C5")
) %>%
  filter(grepl("ANTIGEN_PROCESSING", gs_name, ignore.case = TRUE))

kv("Hallmark gene-set rows retrieved", nrow(hallmark_sets))
kv("Antigen-processing gene-set rows retrieved", nrow(antigen_sets))

if (nrow(hallmark_sets) == 0) {
  stop(paste0(
    "msigdbr returned 0 rows for the Hallmark collection under both the ",
    "'collection' and 'category' argument names. Check the installed ",
    "msigdbr version with packageVersion('msigdbr') and confirm valid ",
    "collection codes with msigdbr_collections()."
  ))
}

pathway_table <- bind_rows(hallmark_sets, antigen_sets) %>%
  distinct(gs_name, gene_symbol) %>%
  mutate(gene_symbol = toupper(trimws(gene_symbol)))

gene_set_list <- split(pathway_table$gene_symbol, pathway_table$gs_name)

# Restrict each gene set to genes actually present in the expression matrix.
gene_set_list <- lapply(gene_set_list, function(genes) {
  intersect(unique(genes), rownames(expr_matrix))
})

gene_set_list <- gene_set_list[lengths(gene_set_list) >= 5]

kv("Pathways included (Hallmark + antigen processing)", length(gene_set_list))

if (length(gene_set_list) == 0) {
  stop(paste0(
    "No pathways retained after intersecting MSigDB gene sets with the ",
    "expression matrix's gene symbols. This usually means the gene-symbol ",
    "mapping in Section 3 did not work as expected — check that ",
    "02_ProteinCoding_GeneAnnotation.csv's 'gene_symbol' column uses the ",
    "same HGNC symbols MSigDB expects (e.g. compare head(rownames(expr_matrix)) ",
    "against head(pathway_table$gene_symbol))."
  ))
}


# ======================================================================
# 5. RUN GSVA
# ======================================================================

cat("\nRunning GSVA (this may take several minutes)...\n")

gsva_param <- GSVA::gsvaParam(
  exprData = expr_matrix,
  geneSets = gene_set_list,
  kcdf     = "Gaussian"
)

gsva_scores <- GSVA::gsva(gsva_param, verbose = FALSE)

kv("GSVA output (pathways x samples)",
   paste(nrow(gsva_scores), "x", ncol(gsva_scores)))

gsva_scores_df <- as.data.frame(gsva_scores) %>%
  rownames_to_column("pathway")

write_csv(
  gsva_scores_df,
  file.path(OUTPUT_DIR, "M7_GSVA_pathway_scores.csv")
)


# ======================================================================
# 6. MERGE WITH FROZEN ECOSYSTEM STATE
# ======================================================================

normalize_id <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("[.\\-]", "_", x)
  x <- sub("_LIONESS$", "", x)
  x
}

gsva_long_raw <- gsva_scores_df %>%
  pivot_longer(cols = -pathway, names_to = "sample_id", values_to = "score") %>%
  mutate(.join_key = normalize_id(sample_id))

ecosystem_states_keyed <- ecosystem_states %>%
  mutate(.join_key = normalize_id(sample_id))

gsva_long <- gsva_long_raw %>%
  select(-sample_id) %>%
  inner_join(
    ecosystem_states_keyed %>% select(.join_key, sample_id, ecosystem_state),
    by = ".join_key"
  ) %>%
  select(-.join_key)

kv("Patient-pathway rows with matched ecosystem state", nrow(gsva_long))

if (nrow(gsva_long) == 0) {
  cat("\nExample IDs for troubleshooting:\n")
  cat("  GSVA columns (normalized):",
      paste(head(unique(normalize_id(gsva_long_raw$sample_id)), 5), collapse = " | "), "\n")
  cat("  states (normalized):      ",
      paste(head(unique(ecosystem_states_keyed$.join_key), 5), collapse = " | "), "\n")
  stop("No overlapping sample IDs between GSVA output and frozen ecosystem states.")
}

if (n_distinct(gsva_long$ecosystem_state) != 2) {
  stop("Expected exactly two frozen ecosystem states (State_1, State_2).")
}

gsva_long$ecosystem_state <- factor(gsva_long$ecosystem_state)
state_levels <- levels(gsva_long$ecosystem_state)


# ======================================================================
# 7. DIFFERENTIAL PATHWAY TESTING (STATE_1 vs STATE_2)
# ======================================================================

cat("\nTesting pathway activity by frozen ecosystem state...\n")

pathway_results <- gsva_long %>%
  group_by(pathway) %>%
  summarise(
    median_state1 = median(score[ecosystem_state == state_levels[1]], na.rm = TRUE),
    median_state2 = median(score[ecosystem_state == state_levels[2]], na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(score ~ ecosystem_state)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    delta_median = median_state1 - median_state2,
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(FDR, p_value)

print(head(pathway_results, 15))

write_csv(
  pathway_results,
  file.path(OUTPUT_DIR, "M7_pathway_state_comparison.csv")
)

n_significant <- sum(pathway_results$FDR < FDR_THRESHOLD, na.rm = TRUE)
kv(paste0("Pathways significant at FDR < ", FDR_THRESHOLD), n_significant)


# ======================================================================
# 8. FIGURE — HEATMAP OF TOP SIGNIFICANT PATHWAYS
# ======================================================================

top_pathways <- pathway_results %>%
  filter(!is.na(FDR)) %>%
  arrange(FDR) %>%
  slice_head(n = TOP_N_PATHWAYS_FOR_FIGURES) %>%
  pull(pathway)

if (length(top_pathways) >= 2) {

  cat("\nGenerating pathway-enrichment heatmap...\n")

  heatmap_matrix <- gsva_scores_df %>%
    filter(pathway %in% top_pathways) %>%
    column_to_rownames("pathway") %>%
    as.matrix()

  # Map heatmap columns (raw GSVA sample IDs) to states via the same
  # normalized-key join used above, then keep the original column names.
  col_lookup <- tibble(sample_id = colnames(heatmap_matrix)) %>%
    mutate(.join_key = normalize_id(sample_id)) %>%
    inner_join(
      ecosystem_states_keyed %>% select(.join_key, ecosystem_state),
      by = ".join_key"
    )

  sample_order <- col_lookup %>%
    arrange(ecosystem_state) %>%
    pull(sample_id)

  heatmap_matrix <- heatmap_matrix[, sample_order, drop = FALSE]

  annotation_col <- col_lookup %>%
    filter(sample_id %in% sample_order) %>%
    distinct(sample_id, ecosystem_state) %>%
    column_to_rownames("sample_id")

  pdf(file.path(OUTPUT_DIR, "M7_pathway_enrichment_heatmap.pdf"),
      width = 12, height = 10)

  pheatmap(
    heatmap_matrix,
    scale = "row",
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    show_colnames = FALSE,
    fontsize_row = 6,
    main = "Top Differential Pathways by Frozen Ecosystem State (GSVA)"
  )

  dev.off()

} else {
  cat("Fewer than two significant pathways; heatmap skipped.\n")
}


# ======================================================================
# 9. FIGURE — DOT PLOT OF TOP PATHWAYS
# ======================================================================

cat("Generating pathway dot plot...\n")

dotplot_data <- pathway_results %>%
  filter(!is.na(FDR)) %>%
  arrange(FDR) %>%
  slice_head(n = TOP_N_PATHWAYS_FOR_FIGURES) %>%
  mutate(
    pathway = gsub("^HALLMARK_|^GOBP_", "", pathway),
    pathway = reorder(pathway, -log10(FDR)),
    direction = ifelse(delta_median > 0,
                        paste0("Higher in ", state_levels[1]),
                        paste0("Higher in ", state_levels[2]))
  )

dot_plot <- ggplot(
  dotplot_data,
  aes(x = -log10(FDR), y = pathway, size = abs(delta_median), color = direction)
) +
  geom_point() +
  labs(
    title = "Top Differential Pathways by Frozen Ecosystem State",
    x = expression(-log[10]("FDR")),
    y = NULL,
    size = "|median difference|",
    color = "Direction"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7))

ggsave(
  file.path(OUTPUT_DIR, "M7_pathway_dotplot.pdf"),
  dot_plot, width = 10, height = 9
)

cat("\nM7 pathway enrichment analysis complete.\n")
cat("Finished:", as.character(Sys.time()), "\n")