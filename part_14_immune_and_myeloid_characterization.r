# ======================================================================
# M5 — MYELOID ENTROPY DIFFERENTIAL ANALYSIS
#
# PURPOSE
# -------
# Test whether a composite myeloid network-entropy score (S_myeloid),
# built from the M3 module-level entropy features already assigned to
# myeloid-lineage MSigDB categories, differs between the FROZEN M3.2
# ecosystem states.
#
# IMPORTANT
# ---------
#   - M3.2 ecosystem states are FROZEN and are not recomputed here.
#   - S_myeloid is constructed ONLY from module-entropy features that
#     were already computed upstream (M3 / M3.2); no new network or
#     entropy calculation is performed in this script.
#
# MYELOID MODULE DEFINITION
# --------------------------
# S_myeloid is defined as the mean of the standardized (z-scored)
# module entropy across the myeloid-lineage MSigDB categories assigned
# by the M2 signature builder, referenced by their actual column names
# in the M3.2 standardized module-entropy matrix (the "S_module_"
# prefix + make.names() convention applied when M3 first wrote
# M3_LIONESS_module_entropy.csv):
#
#   S_module_Macrophage, S_module_Monocyte, S_module_Dendritic_Cell,
#   S_module_Neutrophil, S_module_Myeloid_Activation,
#   S_module_Macrophage_Activation
#
# INPUTS
# -------
# /content/M3_LIONESS_entropy/M3.2_module_axes/
#     M3.2_patient_ecosystem_states.csv
#     M3.2_standardized_module_entropy.csv
#
# OUTPUTS
# --------
# /content/M5_myeloid_entropy/
#     M5_patient_S_myeloid.csv
#     M5_myeloid_entropy_state_comparison.csv
#     M5_myeloid_entropy_violin.pdf
#
# ======================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

cat("\n")
cat("=====================================================================\n")
cat("M5 — MYELOID ENTROPY DIFFERENTIAL ANALYSIS\n")
cat("=====================================================================\n\n")

cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# 1. PACKAGES
# ======================================================================

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

STANDARDIZED_ENTROPY_FILE <- "/content/M3_LIONESS_entropy/M3.2_module_axes/M3.2_standardized_module_entropy.csv"

OUTPUT_DIR <- "/content/M5_myeloid_entropy"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

cat("Output directory:", normalizePath(OUTPUT_DIR), "\n\n")

kv <- function(label, value) cat(label, ":", value, "\n")


# ======================================================================
# 3. LOAD FROZEN INPUTS
# ======================================================================

for (f in c(STATE_FILE, STANDARDIZED_ENTROPY_FILE)) {
  if (!file.exists(f)) {
    stop(paste0("Required FROZEN M3.2 file not found: ", f,
                "\nRun M3.2 before M5."))
  }
}

ecosystem_states <- read_csv(STATE_FILE, show_col_types = FALSE)

standardized_entropy <- read_csv(STANDARDIZED_ENTROPY_FILE, show_col_types = FALSE)

# The standardized module-entropy file is sometimes written with row
# names and no header for the first column, so readr auto-names it
# "...1" (or similar) instead of "sample_id". Detect and fix that here.
if (!"sample_id" %in% colnames(standardized_entropy)) {
  first_col <- colnames(standardized_entropy)[1]
  cat("'sample_id' column not found in standardized module-entropy file; ",
      "treating first column '", first_col, "' as sample_id.\n", sep = "")
  colnames(standardized_entropy)[1] <- "sample_id"
}

if (!"sample_id" %in% colnames(ecosystem_states)) {
  first_col <- colnames(ecosystem_states)[1]
  cat("'sample_id' column not found in ecosystem-state file; ",
      "treating first column '", first_col, "' as sample_id.\n", sep = "")
  colnames(ecosystem_states)[1] <- "sample_id"
}

kv("Patients in FROZEN state file", nrow(ecosystem_states))
kv("Patients in standardized module-entropy matrix", nrow(standardized_entropy))

# ======================================================================
# 4. DEFINE MYELOID MODULE SET AND COMPUTE S_myeloid
# ======================================================================

myeloid_modules_requested <- c(
  "S_module_Macrophage",
  "S_module_Monocyte",
  "S_module_Dendritic_Cell",
  "S_module_Neutrophil",
  "S_module_Myeloid_Activation",
  "S_module_Macrophage_Activation"
)

available_modules <- setdiff(colnames(standardized_entropy), "sample_id")

myeloid_modules_used <- intersect(myeloid_modules_requested, available_modules)

kv("Myeloid modules requested", length(myeloid_modules_requested))
kv("Myeloid modules available in M3.2 matrix", length(myeloid_modules_used))

if (length(setdiff(myeloid_modules_requested, myeloid_modules_used)) > 0) {
  warning(paste0(
    "The following myeloid modules were not found in the M3.2 ",
    "standardized module-entropy matrix and were skipped: ",
    paste(setdiff(myeloid_modules_requested, myeloid_modules_used),
          collapse = ", ")
  ))
}

if (length(myeloid_modules_used) < 1) {
  stop("No myeloid-lineage modules were found in the M3.2 matrix. Cannot compute S_myeloid.")
}

cat("Modules contributing to S_myeloid:",
    paste(myeloid_modules_used, collapse = ", "), "\n\n")

S_myeloid_df <- standardized_entropy %>%
  select(sample_id, all_of(myeloid_modules_used)) %>%
  rowwise() %>%
  mutate(S_myeloid = mean(c_across(all_of(myeloid_modules_used)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(sample_id, S_myeloid)

write_csv(
  S_myeloid_df,
  file.path(OUTPUT_DIR, "M5_patient_S_myeloid.csv")
)


# ======================================================================
# 5. MERGE WITH FROZEN ECOSYSTEM STATE (ID-FORMAT NORMALIZED)
# ======================================================================

normalize_id <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("[.\\-]", "_", x)
  x <- sub("_LIONESS$", "", x)
  x
}

S_myeloid_keyed <- S_myeloid_df %>%
  mutate(.join_key = normalize_id(sample_id))

ecosystem_states_keyed <- ecosystem_states %>%
  select(sample_id, ecosystem_state) %>%
  mutate(.join_key = normalize_id(sample_id))

merged <- ecosystem_states_keyed %>%
  inner_join(
    S_myeloid_keyed %>% select(-sample_id),
    by = ".join_key"
  ) %>%
  select(-.join_key) %>%
  filter(!is.na(S_myeloid))

kv("Patients with both S_myeloid and frozen state", nrow(merged))

if (nrow(merged) == 0) {
  stop("No overlapping sample IDs between S_myeloid and frozen ecosystem states. ",
       "Check ID formats with: head(ecosystem_states$sample_id); head(S_myeloid_df$sample_id)")
}

if (n_distinct(merged$ecosystem_state) != 2) {
  stop("Expected exactly two frozen ecosystem states (State_1, State_2).")
}

merged$ecosystem_state <- factor(merged$ecosystem_state)
state_levels <- levels(merged$ecosystem_state)

# ======================================================================
# 6. DIFFERENTIAL TEST: S_myeloid BY FROZEN ECOSYSTEM STATE
# ======================================================================

cat("\nTesting S_myeloid by frozen ecosystem state...\n")

rank_biserial <- function(x, group) {
  w <- suppressWarnings(wilcox.test(x ~ group)$statistic)
  n1 <- sum(group == levels(group)[1])
  n2 <- sum(group == levels(group)[2])
  1 - (2 * w) / (n1 * n2)
}

wtest <- wilcox.test(S_myeloid ~ ecosystem_state, data = merged)

medians <- tapply(merged$S_myeloid, merged$ecosystem_state, median, na.rm = TRUE)

effect_size <- rank_biserial(merged$S_myeloid, merged$ecosystem_state)

myeloid_summary <- data.frame(
  feature          = "S_myeloid",
  n_modules_used   = length(myeloid_modules_used),
  modules_used     = paste(myeloid_modules_used, collapse = "; "),
  median_state1    = medians[state_levels[1]],
  median_state2    = medians[state_levels[2]],
  n_state1         = sum(merged$ecosystem_state == state_levels[1]),
  n_state2         = sum(merged$ecosystem_state == state_levels[2]),
  p_value          = wtest$p.value,
  effect_size_rb   = effect_size,
  stringsAsFactors = FALSE
)

print(myeloid_summary)

write_csv(
  myeloid_summary,
  file.path(OUTPUT_DIR, "M5_myeloid_entropy_state_comparison.csv")
)


# ======================================================================
# 7. FIGURE — VIOLIN PLOT OF S_myeloid BY STATE
# ======================================================================

cat("\nGenerating S_myeloid violin plot...\n")

p_label <- format.pval(wtest$p.value, digits = 3, eps = 1e-300)

violin_plot <- ggplot(merged, aes(x = ecosystem_state, y = S_myeloid,
                                   fill = ecosystem_state)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, fill = "white") +
  labs(
    title = "Myeloid Network Entropy (S_myeloid) by Frozen Ecosystem State",
    subtitle = paste0(
      "Wilcoxon p = ", p_label,
      "   |   rank-biserial r = ", round(effect_size, 3)
    ),
    x = "Ecosystem state",
    y = "S_myeloid (mean standardized module entropy)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  file.path(OUTPUT_DIR, "M5_myeloid_entropy_violin.pdf"),
  violin_plot,
  width = 6, height = 6
)

cat("\nM5 myeloid entropy differential analysis complete.\n")
cat("Finished:", as.character(Sys.time()), "\n")
