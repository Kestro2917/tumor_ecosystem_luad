# ============================================================
# M2 - Gene-signature library construction and enrichment scoring
#
# Creates:
#   1. M2_CellType_MSigDB.csv
#   2. M2_CellState_MSigDB.csv
#   3. M2_TumorStroma_MSigDB.csv
#   4. M2_MSigDB_signature_summary.csv
#
# Purpose:
#   Build curated MSigDB-derived gene-signature libraries for
#   downstream TCGA-LUAD molecular ecosystem analysis.
#
# ============================================================

options(stringsAsFactors = FALSE)
options(timeout = 600)

cat("\n")
cat("============================================================\n")
cat("M2 - MSigDB SIGNATURE BUILDER FOR TCGA-LUAD\n")
cat("============================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")


# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

required_packages <- c("msigdbr", "dplyr", "stringr", "readr", "tibble")

cat("Checking required R packages...\n")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(msigdbr)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)  # provides tribble(), used for the pattern tables below
})

cat("All required packages loaded successfully.\n\n")


# ============================================================
# 2. SMALL HELPER — LABELED VALUE PRINTING
# ============================================================

kv <- function(label, value) {
  cat(label, ":", value, "\n")
}


# ============================================================
# 3. CREATE OUTPUT DIRECTORY
# ============================================================

outdir <- "M2_MSigDB_LUAD"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

cat("Output directory:", normalizePath(outdir), "\n\n")


# ============================================================
# 4. DOWNLOAD MSigDB HUMAN GENE SETS
# ============================================================

cat("Downloading MSigDB human gene sets...\n")

msig <- tryCatch(
  msigdbr(species = "Homo sapiens"),
  error = function(e) {
    stop(paste0("\nERROR while downloading MSigDB:\n", conditionMessage(e), "\n"))
  }
)

if (nrow(msig) == 0) {
  stop("MSigDB returned zero rows. Cannot continue.")
}

kv("Total MSigDB rows", format(nrow(msig), big.mark = ","))
cat("\n")


# ============================================================
# 5. VERIFY REQUIRED COLUMNS
# ============================================================

cat("Available MSigDB columns:\n")
print(colnames(msig))
cat("\n")

required_columns <- c("gs_name", "gene_symbol", "gs_collection")
missing_columns <- setdiff(required_columns, colnames(msig))

if (length(missing_columns) > 0) {
  stop(paste0(
    "The following required MSigDB columns are missing: ",
    paste(missing_columns, collapse = ", "), "\n"
  ))
}


# ============================================================
# 6. STANDARDIZE GENE NAMES
# ============================================================

cat("Cleaning MSigDB gene-set data...\n")

msig <- msig %>%
  mutate(
    gs_name = as.character(gs_name),
    gene_symbol = toupper(as.character(gene_symbol)),
    gs_collection = as.character(gs_collection)
  ) %>%
  filter(
    !is.na(gs_name), !is.na(gene_symbol), !is.na(gs_collection),
    gs_name != "", gene_symbol != "", gs_collection != ""
  ) %>%
  distinct(gs_name, gene_symbol, .keep_all = TRUE)

kv("Clean MSigDB rows", format(nrow(msig), big.mark = ","))
kv("Unique gene sets", format(n_distinct(msig$gs_name), big.mark = ","))
cat("\n")


# ============================================================
# 7. GENE-SET SIZE + FILTER
# ============================================================

MIN_SIZE <- 10
MAX_SIZE <- 300

set_sizes <- msig %>%
  group_by(gs_name) %>%
  summarise(n_genes = n_distinct(gene_symbol), .groups = "drop")

msig <- msig %>% left_join(set_sizes, by = "gs_name")

cat("Applying gene-set size filter:", MIN_SIZE, "to", MAX_SIZE, "genes\n")

msig_filtered <- msig %>% filter(n_genes >= MIN_SIZE, n_genes <= MAX_SIZE)

kv("Gene sets remaining after size filtering",
   format(n_distinct(msig_filtered$gs_name), big.mark = ","))
cat("\n")


# ============================================================
# 8. HELPERS: PATTERN MATCHING, CATEGORY ASSIGNMENT, JACCARD
# ============================================================

make_pattern <- function(terms) paste(terms, collapse = "|")

# Assigns a category to each gene-set name using an ORDERED table of
# (category, pattern) pairs - first matching pattern wins, exactly
# reproducing case_when()'s top-to-bottom semantics. Rows with no
# match are dropped, matching the original filter(!is.na(category)).
assign_category <- function(data, pattern_table) {
  upper_name <- toupper(data$gs_name)
  category <- rep(NA_character_, nrow(data))

  for (i in seq_len(nrow(pattern_table))) {
    still_unassigned <- is.na(category)
    hit <- still_unassigned & str_detect(upper_name, pattern_table$pattern[i])
    category[hit] <- pattern_table$category[i]
  }

  data$category <- category
  data %>% filter(!is.na(category))
}

jaccard <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  union_genes <- union(a, b)
  if (length(union_genes) == 0) return(0)
  length(intersect(a, b)) / length(union_genes)
}

# Two signatures are considered redundant when Jaccard similarity
# >= threshold. Within the same category, the larger signature wins.
reduce_redundancy <- function(data, threshold = 0.80) {
  if (nrow(data) == 0) return(data)

  sets <- data %>%
    group_by(category, gs_name) %>%
    summarise(genes = list(unique(gene_symbol)), n_genes = n_distinct(gene_symbol),
              .groups = "drop") %>%
    arrange(category, desc(n_genes), gs_name)

  keep <- rep(TRUE, nrow(sets))

  if (nrow(sets) > 1) {
    for (i in seq_len(nrow(sets) - 1)) {
      if (!keep[i]) next

      for (j in seq.int(i + 1, nrow(sets))) {
        if (!keep[j] || sets$category[i] != sets$category[j]) next

        similarity <- jaccard(sets$genes[[i]], sets$genes[[j]])
        if (is.finite(similarity) && similarity >= threshold) keep[j] <- FALSE
      }
    }
  }

  retained <- sets %>% filter(keep) %>% select(category, gs_name)
  data %>% inner_join(retained, by = c("category", "gs_name"))
}

format_output <- function(data, collection_name) {
  if (nrow(data) == 0) {
    return(data.frame(
      signature = character(0), category = character(0), gene = character(0),
      source = character(0), collection = character(0), stringsAsFactors = FALSE
    ))
  }

  data %>%
    select(signature = gs_name, category, gene = gene_symbol) %>%
    distinct() %>%
    arrange(category, signature, gene) %>%
    mutate(source = "MSigDB", collection = collection_name) %>%
    select(signature, category, gene, source, collection)
}

check_sizes <- function(data, label) {
  if (nrow(data) == 0) {
    cat(label, ": no signatures found.\n")
    return(data.frame(signature = character(0), category = character(0),
                       n_genes = integer(0), stringsAsFactors = FALSE))
  }

  sizes <- data %>%
    group_by(signature, category) %>%
    summarise(n_genes = n_distinct(gene), .groups = "drop")

  bad <- sizes %>% filter(n_genes < MIN_SIZE)

  kv(paste(label, "signatures"), nrow(sizes))
  kv(paste(label, "gene entries"), nrow(data))

  if (nrow(bad) > 0) {
    warning(paste(label, "contains signatures below the minimum size of", MIN_SIZE))
  }

  sizes
}

# One end-to-end library builder: filter to collection(s), keyword
# match against gs_name, assign category, then reduce redundancy.
# Replaces the three duplicated filter -> case_when -> filter blocks.
build_signature_library <- function(msig_filtered, collections, keyword_terms,
                                     pattern_table, label, threshold = 0.80) {

  candidates <- msig_filtered %>%
    filter(gs_collection %in% collections) %>%
    filter(str_detect(toupper(gs_name), make_pattern(keyword_terms)))

  kv(paste(label, "gene sets identified"), n_distinct(candidates$gs_name))

  categorized <- assign_category(candidates, pattern_table)
  reduced <- reduce_redundancy(categorized, threshold = threshold)

  kv(paste("Remaining", tolower(label), "signatures"), n_distinct(reduced$gs_name))

  reduced
}


# ============================================================
# 9. CELL-TYPE SIGNATURES (MSigDB C8)
# ============================================================

cat("------------------------------------------------------------\n")
cat("9. Building cell-type signatures\n")
cat("------------------------------------------------------------\n")

celltype_terms <- c(
  "MACROPHAGE", "MONOCYTE", "DENDRITIC", "NEUTROPHIL", "MAST",
  "NATURAL_KILLER", "NK_CELL", "NK_CELLS", "CD8", "CD4", "TREG",
  "REGULATORY_T", "T_CELL", "B_CELL", "B_CELLS", "PLASMA_CELL",
  "FIBROBLAST", "ENDOTHELIAL", "PERICYTE"
)

celltype_patterns <- tribble(
  ~category,        ~pattern,
  "Macrophage",      "MACROPHAGE",
  "Monocyte",        "MONOCYTE",
  "Dendritic_Cell",  "DENDRITIC",
  "Neutrophil",      "NEUTROPHIL",
  "Mast_Cell",       "MAST",
  "NK_Cell",         "NATURAL_KILLER|NK_CELL|NK_CELLS",
  "Treg",            "TREG|REGULATORY_T",
  "CD8_T_Cell",      "CD8",
  "CD4_T_Cell",      "CD4",
  "T_Cell",          "T_CELL",
  "Plasma_Cell",     "PLASMA_CELL",
  "B_Cell",          "B_CELL|B_CELLS",
  "Fibroblast",      "FIBROBLAST",
  "Endothelial",     "ENDOTHELIAL",
  "Pericyte",        "PERICYTE"
)

celltype_final <- build_signature_library(
  msig_filtered,
  collections = "C8",
  keyword_terms = celltype_terms,
  pattern_table = celltype_patterns,
  label = "Cell-type"
)


# ============================================================
# 10. CELL-STATE SIGNATURES (MSigDB C7, C5)
# ============================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("10. Building cell-state signatures\n")
cat("------------------------------------------------------------\n")

state_terms <- c(
  "EXHAUST", "DYSFUNCTION", "CYTOTOX", "CYTOLYTIC", "INTERFERON", "IFN",
  "ANTIGEN", "MHC", "IMMUNE", "INFLAMMAT", "TNF", "NF_KB",
  "T_CELL_ACTIV", "T_CELL_RESPONSE", "B_CELL_ACTIV", "MACROPHAGE_ACTIV",
  "MYELOID", "SUPPRESS", "ACTIVATION"
)

cellstate_patterns <- tribble(
  ~category,               ~pattern,
  "T_Cell_Exhaustion",      "EXHAUST|DYSFUNCTION",
  "Cytotoxicity",           "CYTOTOX|CYTOLYTIC",
  "Interferon",             "INTERFERON|IFN",
  "Antigen_Presentation",   "ANTIGEN|MHC",
  "Myeloid_Activation",     "MYELOID",
  "Macrophage_Activation",  "MACROPHAGE",
  "Immune_Suppression",     "SUPPRESS",
  "Inflammation",           "INFLAMMAT|TNF|NF_KB",
  "Immune_Activation",      "T_CELL_ACTIV|T_CELL_RESPONSE|ACTIVATION",
  "B_Cell_Activation",      "B_CELL_ACTIV"
)

cellstate_final <- build_signature_library(
  msig_filtered,
  collections = c("C7", "C5"),
  keyword_terms = state_terms,
  pattern_table = cellstate_patterns,
  label = "Cell-state"
)


# ============================================================
# 11. TUMOR / STROMA SIGNATURES (MSigDB C2, C4, C5)
# ============================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("11. Building tumor/stroma signatures\n")
cat("------------------------------------------------------------\n")

tumor_terms <- c(
  "EPITHELIAL", "MALIGNANT", "CANCER", "TUMOR", "EMT", "MESENCHYM",
  "HYPOXIA", "FIBROBLAST", "CAF", "ENDOTHELIAL", "PERICYTE",
  "EXTRACELLULAR_MATRIX", "ECM", "ANGIOGEN", "PROLIFER", "CELL_CYCLE"
)

tumorstroma_patterns <- tribble(
  ~category,        ~pattern,
  "EMT",             "EMT|EPITHELIAL_MESENCHYM",
  "Hypoxia",         "HYPOXIA|HYPOXIC",
  "CAF",             "FIBROBLAST|CAF",
  "Endothelial",     "ENDOTHELIAL|ANGIOGEN",
  "Pericyte",        "PERICYTE",
  "ECM",             "EXTRACELLULAR_MATRIX|ECM",
  "Proliferation",   "PROLIFER|CELL_CYCLE",
  "Malignant",       "MALIGNANT|CANCER|TUMOR",
  "Epithelial",      "EPITHELIAL"
)

tumorstroma_final <- build_signature_library(
  msig_filtered,
  collections = c("C2", "C4", "C5"),
  keyword_terms = tumor_terms,
  pattern_table = tumorstroma_patterns,
  label = "Tumor/stroma"
)


# ============================================================
# 12. FORMAT FINAL OUTPUT TABLES
# ============================================================

celltype_final <- format_output(celltype_final, "C8")
cellstate_final <- format_output(cellstate_final, "C7_C5")
tumorstroma_final <- format_output(tumorstroma_final, "C2_C4_C5")


# ============================================================
# 13. FINAL SIGNATURE SIZE CHECKS
# ============================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("13. Final signature-size checks\n")
cat("------------------------------------------------------------\n")

celltype_sizes <- check_sizes(celltype_final, "Cell type")
cellstate_sizes <- check_sizes(cellstate_final, "Cell state")
tumorstroma_sizes <- check_sizes(tumorstroma_final, "Tumor/Stroma")


# ============================================================
# 14. SAVE FINAL CSV FILES
# ============================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("14. Writing CSV files\n")
cat("------------------------------------------------------------\n")

celltype_file <- file.path(outdir, "M2_CellType_MSigDB.csv")
cellstate_file <- file.path(outdir, "M2_CellState_MSigDB.csv")
tumorstroma_file <- file.path(outdir, "M2_TumorStroma_MSigDB.csv")
summary_file <- file.path(outdir, "M2_MSigDB_signature_summary.csv")

write.csv(celltype_final, celltype_file, row.names = FALSE, quote = FALSE)
write.csv(cellstate_final, cellstate_file, row.names = FALSE, quote = FALSE)
write.csv(tumorstroma_final, tumorstroma_file, row.names = FALSE, quote = FALSE)

cat("Created:", celltype_file, "\n")
cat("Created:", cellstate_file, "\n")
cat("Created:", tumorstroma_file, "\n")


# ============================================================
# 15. SUMMARY TABLE
# ============================================================

summary_table <- bind_rows(
  celltype_sizes %>% mutate(library = "CellType"),
  cellstate_sizes %>% mutate(library = "CellState"),
  tumorstroma_sizes %>% mutate(library = "TumorStroma")
) %>%
  select(library, category, signature, n_genes) %>%
  arrange(library, category, signature)

write.csv(summary_table, summary_file, row.names = FALSE, quote = FALSE)

cat("Created:", summary_file, "\n")


# ============================================================
# 16. FINAL RESULTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("M2 COMPLETE\n")
cat("============================================================\n")

kv("Cell-type signatures", n_distinct(celltype_final$signature))
kv("Cell-state signatures", n_distinct(cellstate_final$signature))
kv("Tumor/stroma signatures", n_distinct(tumorstroma_final$signature))

print_categories <- function(data, label) {
  cat("\n", label, " categories:\n", sep = "")
  if (nrow(data) > 0) {
    print(table(data$category))
  } else {
    cat("No", tolower(label), "signatures found.\n")
  }
}

print_categories(celltype_final, "Cell-type")
print_categories(cellstate_final, "Cell-state")
print_categories(tumorstroma_final, "Tumor/stroma")


# ============================================================
# 17. OUTPUT FILE LOCATIONS
# ============================================================

cat("\nFiles created:\n")
cat("1.", normalizePath(celltype_file), "\n")
cat("2.", normalizePath(cellstate_file), "\n")
cat("3.", normalizePath(tumorstroma_file), "\n")
cat("4.", normalizePath(summary_file), "\n")

cat("\nOutput directory contents:\n")
print(list.files(outdir, full.names = TRUE))


# ============================================================
# END
# ============================================================

cat("\n")
cat("M2 finished successfully.\n")
cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")