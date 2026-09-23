# ======================================================================
# M6 — KEAP1 / STK11 MUTATION ENRICHMENT ANALYSIS
# PURPOSE
# -------
# Test whether KEAP1 and/or STK11 somatic mutation status is enriched
# in one FROZEN M3.2 ecosystem state versus the other, and whether
# mutation status associates with myeloid network entropy (S_myeloid,
# from M5).
#
#
# EXPECTED MUTATION INPUT
# ------------------------
# This script accepts EITHER of the following, and will auto-detect
# which is supplied. If NEITHER is present, it will attempt to fetch
# TCGA-LUAD masked somatic mutation calls directly from GDC via
# TCGAbiolinks and write them to MUTATION_MAF_FILE automatically.
#
#   (A) A binary patient x gene mutation matrix, one row per patient:
#         M6_mutation_matrix.csv
#         columns: sample_id, <gene_1>, <gene_2>, ... (0 = wild-type,
#         1 = mutant; nonsynonymous/functional calls only)
#
#   (B) A standard MAF-style long table:
#         M6_mutations.maf.csv
#         required columns: Tumor_Sample_Barcode, Hugo_Symbol,
#         Variant_Classification
#         (silent / intronic / UTR variants are excluded; all other
#         classifications are treated as mutation-positive)
#
# Only KEAP1 and STK11 are extracted for this analysis; the input file
# may contain additional genes, which are ignored here.
#
# METHOD
# ------
# 1. Load mutation calls (or auto-fetch from GDC); derive a binary
#    KEAP1 and STK11 status per patient, plus a combined
#    "KEAP1_or_STK11" co-mutation flag.
# 2. Merge with the FROZEN M3.2 ecosystem-state assignment and, where
#    available, the M5 S_myeloid score, joined at the 12-character
#    TCGA patient-barcode level.
# 3. Fisher's exact test: mutation status x ecosystem state, for
#    KEAP1, STK11, and the combined flag; BH-FDR across the 3 tests.
# 4. Wilcoxon rank-sum test: S_myeloid by mutation status, for the
#    same 3 groupings; BH-FDR across the 3 tests.
# 5. Figures: mutation-frequency barplot by state, a lightweight
#    oncoprint (tile plot) of KEAP1/STK11 status ordered by ecosystem
#    state, and boxplots of S_myeloid stratified by mutation status
#    and colored by ecosystem state.
#
# INPUTS
# -------
# /content/M3_LIONESS_entropy/M3.2_module_axes/
#     M3.2_patient_ecosystem_states.csv
# /content/M5_myeloid_entropy/M5_patient_S_myeloid.csv   (optional)
# /content/M6_mutation_matrix.csv   OR   /content/M6_mutations.maf.csv
#     (auto-fetched from GDC into the .maf.csv path if neither exists)
#
# OUTPUTS
# --------
# /content/M6_mutation_enrichment/
#     M6_mutation_status_per_patient.csv
#     M6_mutation_state_enrichment.csv
#     M6_mutation_S_myeloid_association.csv
#     M6_mutation_frequency_barplot.pdf
#     M6_oncoprint.pdf
#     M6_S_myeloid_by_mutation_boxplots.pdf
#
# ======================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

cat("\n")
cat("=====================================================================\n")
cat("M6 — KEAP1 / STK11 MUTATION ENRICHMENT ANALYSIS\n")
cat("=====================================================================\n\n")

cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# 1. PACKAGES
# ======================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "tibble")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

cat("Packages ready.\n\n")


# ======================================================================
# 2. FILE SETTINGS
# ======================================================================

STATE_FILE <- "/content/M3_LIONESS_entropy/M3.2_module_axes/M3.2_patient_ecosystem_states.csv"

S_MYELOID_FILE <- "/content/M5_myeloid_entropy/M5_patient_S_myeloid.csv"

MUTATION_MATRIX_FILE <- "/content/M6_mutation_matrix.csv"
MUTATION_MAF_FILE     <- "/content/M6_mutations.maf.csv"

GENES_OF_INTEREST <- c("KEAP1", "STK11")

EXCLUDED_VARIANT_CLASSES <- c(
  "Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
  "IGR", "RNA", "lincRNA"
)

OUTPUT_DIR <- "/content/M6_mutation_enrichment"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

cat("Output directory:", normalizePath(OUTPUT_DIR), "\n\n")

kv <- function(label, value) cat(label, ":", value, "\n")

# Normalize + collapse any sample/barcode identifier to the 12-character
# TCGA PATIENT barcode (TCGA-XX-XXXX), which is the only level at which
# MAF (DNA aliquot) and expression/entropy (RNA aliquot) barcodes for
# the same patient are guaranteed to agree.
patient_key <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("[.\\-]", "_", x)     # unify "-" and "." to "_"
  x <- sub("_LIONESS$", "", x)    # strip pipeline-added suffix, if present
  substr(x, 1, 12)                # TCGA_XX_XXXX patient-level barcode
}


# ======================================================================
# 3. LOAD FROZEN ECOSYSTEM STATE (REQUIRED) AND S_myeloid (OPTIONAL)
# ======================================================================

if (!file.exists(STATE_FILE)) {
  stop(paste0("FROZEN M3.2 ecosystem-state file not found: ", STATE_FILE,
              "\nRun M3.2 before M6."))
}

ecosystem_states <- read_csv(STATE_FILE, show_col_types = FALSE)

if (!"sample_id" %in% colnames(ecosystem_states)) {
  colnames(ecosystem_states)[1] <- "sample_id"
}

ecosystem_states <- ecosystem_states %>%
  select(sample_id, ecosystem_state)

kv("Patients with frozen ecosystem state", nrow(ecosystem_states))

has_S_myeloid <- file.exists(S_MYELOID_FILE)

if (has_S_myeloid) {
  S_myeloid_df <- read_csv(S_MYELOID_FILE, show_col_types = FALSE)
  if (!"sample_id" %in% colnames(S_myeloid_df)) {
    colnames(S_myeloid_df)[1] <- "sample_id"
  }
  cat("M5 S_myeloid scores found and will be included in this analysis.\n")
} else {
  S_myeloid_df <- NULL
  cat("M5 S_myeloid scores not found; mutation-vs-entropy association ",
      "step will be skipped (run M5 first to enable it).\n", sep = "")
}


# ======================================================================
# 3.5. AUTO-FETCH TCGA-LUAD MUTATION CALLS IF NO LOCAL FILE IS PRESENT
# ======================================================================

if (!file.exists(MUTATION_MATRIX_FILE) && !file.exists(MUTATION_MAF_FILE)) {

  cat("No local mutation file found at either expected path.\n")
  cat("Attempting to fetch TCGA-LUAD masked somatic mutation calls from GDC...\n")

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    cat("Installing: TCGAbiolinks (this can take several minutes)\n")
    BiocManager::install("TCGAbiolinks", ask = FALSE, update = FALSE)
  }

  suppressPackageStartupMessages(library(TCGAbiolinks))

  maf_query <- tryCatch(
    GDCquery(
      project = "TCGA-LUAD",
      data.category = "Simple Nucleotide Variation",
      data.type = "Masked Somatic Mutation",
      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
    ),
    error = function(e) {
      cat("GDCquery failed with workflow.type filter:", conditionMessage(e), "\n")
      cat("Retrying without a workflow.type filter...\n")
      tryCatch(
        GDCquery(
          project = "TCGA-LUAD",
          data.category = "Simple Nucleotide Variation",
          data.type = "Masked Somatic Mutation"
        ),
        error = function(e2) {
          cat("GDCquery failed again:", conditionMessage(e2), "\n")
          NULL
        }
      )
    }
  )

  if (!is.null(maf_query)) {

    GDCdownload(maf_query)
    maf_data <- GDCprepare(maf_query)

    maf_export <- as.data.frame(maf_data) %>%
      select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
      distinct()

    write_csv(maf_export, MUTATION_MAF_FILE)

    cat("Fetched", nrow(maf_export), "mutation records for",
        n_distinct(maf_export$Tumor_Sample_Barcode), "samples.\n")
    cat("Saved to:", MUTATION_MAF_FILE, "\n\n")

  } else {

    stop(paste0(
      "Automatic GDC fetch failed (e.g. no internet from this session, ",
      "or the GDC data.type/workflow.type naming has changed). Download a ",
      "TCGA-LUAD MAF manually (GDC Data Portal or cBioPortal) and save it as:\n  ",
      MUTATION_MAF_FILE,
      "\nwith at least the columns Tumor_Sample_Barcode, Hugo_Symbol, ",
      "Variant_Classification."
    ))
  }
}


# ======================================================================
# 4. LOAD MUTATION CALLS (AUTO-DETECT FORMAT) AND DERIVE STATUS
# ======================================================================

cat("\nLoading mutation calls...\n")

mutation_status <- NULL

if (file.exists(MUTATION_MATRIX_FILE)) {

  cat("Detected binary mutation matrix:", MUTATION_MATRIX_FILE, "\n")

  mut_matrix <- read_csv(MUTATION_MATRIX_FILE, show_col_types = FALSE)

  if (!"sample_id" %in% colnames(mut_matrix)) {
    colnames(mut_matrix)[1] <- "sample_id"
  }

  genes_present <- intersect(GENES_OF_INTEREST, colnames(mut_matrix))

  if (length(genes_present) == 0) {
    stop("Neither KEAP1 nor STK11 found as columns in the mutation matrix.")
  }

  mutation_status <- mut_matrix %>%
    select(sample_id, all_of(genes_present)) %>%
    mutate(across(all_of(genes_present), ~ as.integer(. > 0)))

} else if (file.exists(MUTATION_MAF_FILE)) {

  cat("Detected MAF-style mutation table:", MUTATION_MAF_FILE, "\n")

  maf <- read_csv(MUTATION_MAF_FILE, show_col_types = FALSE)

  required_maf_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol",
                          "Variant_Classification")

  if (!all(required_maf_cols %in% colnames(maf))) {
    stop(paste0(
      "MAF-style mutation file must contain columns: ",
      paste(required_maf_cols, collapse = ", ")
    ))
  }

  # MAF (DNA) aliquot barcodes end in a different portion/analyte code
  # than RNA-seq aliquot barcodes for the same patient, so mutation
  # status is derived at the 12-character PATIENT barcode level.
  maf_filtered <- maf %>%
    filter(
      Hugo_Symbol %in% GENES_OF_INTEREST,
      !(Variant_Classification %in% EXCLUDED_VARIANT_CLASSES)
    ) %>%
    mutate(patient_barcode = patient_key(Tumor_Sample_Barcode))

  all_samples <- unique(patient_key(maf$Tumor_Sample_Barcode))

  mutation_status <- tibble(sample_id = all_samples)

  for (gene in GENES_OF_INTEREST) {
    mutant_samples <- maf_filtered %>%
      filter(Hugo_Symbol == gene) %>%
      pull(patient_barcode) %>%
      unique()

    mutation_status[[gene]] <- as.integer(
      mutation_status$sample_id %in% mutant_samples
    )
  }

} else {

  stop(paste0(
    "No mutation input file found. Provide either:\n  ",
    MUTATION_MATRIX_FILE, "\nor\n  ", MUTATION_MAF_FILE
  ))
}

genes_available <- intersect(GENES_OF_INTEREST, colnames(mutation_status))

kv("Genes available for enrichment testing", paste(genes_available, collapse = ", "))

if ("KEAP1" %in% genes_available && "STK11" %in% genes_available) {
  mutation_status <- mutation_status %>%
    mutate(KEAP1_or_STK11 = as.integer(KEAP1 == 1 | STK11 == 1))
  mutation_groupings <- c(genes_available, "KEAP1_or_STK11")
} else {
  mutation_groupings <- genes_available
}

write_csv(
  mutation_status,
  file.path(OUTPUT_DIR, "M6_mutation_status_per_patient.csv")
)


# ======================================================================
# 5. MERGE MUTATION STATUS WITH FROZEN ECOSYSTEM STATE (+ S_myeloid)
#    Joined at the 12-character TCGA PATIENT barcode level, since MAF
#    (DNA) and expression/entropy (RNA) aliquot suffixes never match.
# ======================================================================

mutation_status_keyed <- mutation_status %>%
  mutate(.join_key = patient_key(sample_id))

ecosystem_states_keyed <- ecosystem_states %>%
  mutate(.join_key = patient_key(sample_id))

merged <- ecosystem_states_keyed %>%
  inner_join(
    mutation_status_keyed %>% select(-sample_id),
    by = ".join_key"
  ) %>%
  select(-.join_key)

if (has_S_myeloid) {
  S_myeloid_keyed <- S_myeloid_df %>%
    mutate(.join_key = patient_key(sample_id)) %>%
    select(-sample_id)

  merged <- merged %>%
    mutate(.join_key = patient_key(sample_id)) %>%
    left_join(S_myeloid_keyed, by = ".join_key") %>%
    select(-.join_key)
}

kv("Patients with both mutation status and frozen state", nrow(merged))

if (nrow(merged) == 0) {
  cat("\nNo overlapping patient barcodes. Example IDs:\n")
  cat("  ecosystem_states (patient-level):",
      paste(head(patient_key(ecosystem_states$sample_id), 5), collapse = " | "), "\n")
  cat("  mutation_status (patient-level): ",
      paste(head(patient_key(mutation_status$sample_id), 5), collapse = " | "), "\n")
  stop("No overlapping sample IDs between mutation data and frozen ecosystem states. ",
       "Compare the printed patient-barcode examples above.")
}

if (n_distinct(merged$ecosystem_state) != 2) {
  stop("Expected exactly two frozen ecosystem states (State_1, State_2).")
}

merged$ecosystem_state <- factor(merged$ecosystem_state)
state_levels <- levels(merged$ecosystem_state)


# ======================================================================
# 6. FISHER'S EXACT TEST: MUTATION STATUS x ECOSYSTEM STATE
# ======================================================================

cat("\nTesting mutation-status enrichment by frozen ecosystem state...\n")

enrichment_results <- lapply(mutation_groupings, function(gene) {

  status <- factor(merged[[gene]], levels = c(0, 1),
                    labels = c("WT", "Mutant"))

  tab <- table(status, merged$ecosystem_state)

  test <- fisher.test(tab)

  freq_state1 <- mean(merged[[gene]][merged$ecosystem_state == state_levels[1]] == 1, na.rm = TRUE)
  freq_state2 <- mean(merged[[gene]][merged$ecosystem_state == state_levels[2]] == 1, na.rm = TRUE)

  data.frame(
    gene              = gene,
    mutation_freq_state1 = freq_state1,
    mutation_freq_state2 = freq_state2,
    odds_ratio        = unname(test$estimate),
    p_value           = test$p.value,
    stringsAsFactors  = FALSE
  )
})

enrichment_df <- bind_rows(enrichment_results)
enrichment_df$FDR <- p.adjust(enrichment_df$p_value, method = "BH")
enrichment_df <- enrichment_df %>% arrange(FDR, p_value)

print(enrichment_df)

write_csv(
  enrichment_df,
  file.path(OUTPUT_DIR, "M6_mutation_state_enrichment.csv")
)


# ======================================================================
# 7. MUTATION STATUS x S_myeloid (WILCOXON), IF AVAILABLE
# ======================================================================

if (has_S_myeloid && "S_myeloid" %in% colnames(merged)) {

  cat("\nTesting S_myeloid by mutation status...\n")

  myeloid_assoc <- lapply(mutation_groupings, function(gene) {

    df <- merged %>%
      filter(!is.na(S_myeloid)) %>%
      mutate(status = factor(.data[[gene]], levels = c(0, 1),
                              labels = c("WT", "Mutant")))

    if (n_distinct(df$status) < 2) {
      return(data.frame(
        gene = gene, median_WT = NA, median_Mutant = NA,
        p_value = NA_real_, stringsAsFactors = FALSE
      ))
    }

    test <- wilcox.test(S_myeloid ~ status, data = df)
    medians <- tapply(df$S_myeloid, df$status, median, na.rm = TRUE)

    data.frame(
      gene          = gene,
      median_WT     = medians["WT"],
      median_Mutant = medians["Mutant"],
      p_value       = test$p.value,
      stringsAsFactors = FALSE
    )
  })

  myeloid_assoc_df <- bind_rows(myeloid_assoc)
  myeloid_assoc_df$FDR <- p.adjust(myeloid_assoc_df$p_value, method = "BH")
  myeloid_assoc_df <- myeloid_assoc_df %>% arrange(FDR, p_value)

  print(myeloid_assoc_df)

  write_csv(
    myeloid_assoc_df,
    file.path(OUTPUT_DIR, "M6_mutation_S_myeloid_association.csv")
  )

} else {
  cat("\nSkipping mutation-vs-S_myeloid association (M5 output not found).\n")
}


# ======================================================================
# 8. FIGURE — MUTATION FREQUENCY BARPLOT BY STATE
# ======================================================================

cat("\nGenerating mutation-frequency barplot...\n")

freq_long <- enrichment_df %>%
  select(gene, mutation_freq_state1, mutation_freq_state2) %>%
  pivot_longer(
    cols = c(mutation_freq_state1, mutation_freq_state2),
    names_to = "state", values_to = "frequency"
  ) %>%
  mutate(
    state = recode(state,
                    mutation_freq_state1 = state_levels[1],
                    mutation_freq_state2 = state_levels[2])
  )

barplot_freq <- ggplot(freq_long, aes(x = gene, y = frequency, fill = state)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(
    title = "KEAP1 / STK11 Mutation Frequency by Frozen Ecosystem State",
    x = "Gene", y = "Mutation frequency", fill = "Ecosystem state"
  ) +
  theme_bw()

ggsave(
  file.path(OUTPUT_DIR, "M6_mutation_frequency_barplot.pdf"),
  barplot_freq, width = 7, height = 5
)


# ======================================================================
# 9. FIGURE — LIGHTWEIGHT ONCOPRINT (TILE PLOT)
# ======================================================================

cat("Generating oncoprint-style tile plot...\n")

oncoprint_long <- merged %>%
  select(sample_id, ecosystem_state, all_of(genes_available)) %>%
  arrange(ecosystem_state) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id))) %>%
  pivot_longer(
    cols = all_of(genes_available),
    names_to = "gene", values_to = "mutant"
  ) %>%
  mutate(mutant = factor(mutant, levels = c(0, 1), labels = c("WT", "Mutant")))

oncoprint_plot <- ggplot(oncoprint_long, aes(x = sample_id, y = gene, fill = mutant)) +
  geom_tile(color = "grey90") +
  facet_grid(~ ecosystem_state, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("WT" = "grey85", "Mutant" = "firebrick")) +
  labs(title = "KEAP1 / STK11 Mutation Status by Frozen Ecosystem State",
       x = "Patient", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  file.path(OUTPUT_DIR, "M6_oncoprint.pdf"),
  oncoprint_plot, width = 12, height = 3.5
)


# ======================================================================
# 10. FIGURE — S_myeloid BY MUTATION STATUS, COLORED BY STATE
# ======================================================================

if (has_S_myeloid && "S_myeloid" %in% colnames(merged)) {

  cat("Generating S_myeloid-by-mutation-status boxplots...\n")

  box_long <- merged %>%
    filter(!is.na(S_myeloid)) %>%
    select(sample_id, ecosystem_state, S_myeloid, all_of(genes_available)) %>%
    pivot_longer(
      cols = all_of(genes_available),
      names_to = "gene", values_to = "mutant"
    ) %>%
    mutate(mutant = factor(mutant, levels = c(0, 1), labels = c("WT", "Mutant")))

  box_plot <- ggplot(box_long, aes(x = mutant, y = S_myeloid, fill = ecosystem_state)) +
    geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.75)) +
    facet_wrap(~ gene) +
    labs(
      title = "S_myeloid by Mutation Status and Frozen Ecosystem State",
      x = "Mutation status", y = "S_myeloid", fill = "Ecosystem state"
    ) +
    theme_bw()

  ggsave(
    file.path(OUTPUT_DIR, "M6_S_myeloid_by_mutation_boxplots.pdf"),
    box_plot, width = 9, height = 5
  )
}

cat("\nM6 KEAP1/STK11 mutation enrichment analysis complete.\n")
cat("Finished:", as.character(Sys.time()), "\n")