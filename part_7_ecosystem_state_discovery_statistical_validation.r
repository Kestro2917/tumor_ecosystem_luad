# ======================================================================
# M3.4 — LIONESS ECOSYSTEM STATE QUANTITATIVE VALIDATION
# PURPOSE
# -------
# M3.4 quantitatively validates the FROZEN ecosystem states generated
# in M3.2 and biologically validated in M3.3.
#
# IMPORTANT
# ---------
# This script DOES NOT recreate or redefine ecosystem states.
#
# It uses:
#   M3.2_patient_ecosystem_states.csv
#   M3.2_standardized_module_entropy.csv
#
# Main analyses:
#
#   1. Frozen state composition
#   2. Global module-state separation
#   3. PERMANOVA-style multivariate separation
#   4. Silhouette analysis
#   5. Pairwise centroid distances
#   6. Within-state / between-state distances
#   7. PCA separation
#   8. Leave-one-out classification accuracy
#   9. Bootstrap state stability
#  10. Module contribution analysis
#  11. Publication-ready figures
#
# OUTPUT DIRECTORY
# ----------------
# /content/M3_LIONESS_entropy/M3.4_quantitative_validation
#
# ======================================================================


# ======================================================================
# SECTION 1 — BASIC SETTINGS
# ======================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

cat("\n")
cat("=====================================================================\n")
cat("M3.4 — LIONESS ECOSYSTEM STATE QUANTITATIVE VALIDATION\n")
cat("=====================================================================\n\n")

cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# SECTION 2 — INSTALL / LOAD PACKAGES
# ======================================================================

cat("Checking required R packages...\n")

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "cluster",
  "vegan",
  "patchwork"
)

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    cat("Installing:", pkg, "\n")

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org",
      quiet = TRUE
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cluster)
  library(vegan)
  library(patchwork)
})

cat("Packages ready.\n\n")


# ======================================================================
# SECTION 3 — DIRECTORIES
# ======================================================================

base_dir <- "/content/M3_LIONESS_entropy"

input_dir <- file.path(
  base_dir,
  "M3.2_module_axes"
)

m33_dir <- file.path(
  base_dir,
  "M3.3_biological_validation"
)

output_dir <- file.path(
  base_dir,
  "M3.4_quantitative_validation"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

cat("Input directory:\n")
cat(input_dir, "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")


# ======================================================================
# SECTION 4 — INPUT FILES
# ======================================================================

state_file <- file.path(
  input_dir,
  "M3.2_patient_ecosystem_states.csv"
)

entropy_file <- file.path(
  input_dir,
  "M3.2_standardized_module_entropy.csv"
)

cat("Checking input files...\n\n")

if (!file.exists(state_file)) {
  stop(
    paste0(
      "ERROR: Frozen state file not found:\n",
      state_file
    )
  )
}

if (!file.exists(entropy_file)) {
  stop(
    paste0(
      "ERROR: Standardized module entropy file not found:\n",
      entropy_file
    )
  )
}

cat("FOUND:\n")
cat(state_file, "\n")
cat(entropy_file, "\n\n")

cat("All required files found.\n")


# ======================================================================
# SECTION 5 — READ FROZEN STATES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING FROZEN ECOSYSTEM STATES\n")
cat("=====================================================================\n")

states <- read_csv(
  state_file,
  show_col_types = FALSE
)

cat("Rows:", nrow(states), "\n")
cat("Columns:", ncol(states), "\n\n")

cat("Columns detected:\n")
print(names(states))

# ----------------------------------------------------------------------
# Detect sample ID
# ----------------------------------------------------------------------

possible_id_names <- c(
  "sample_id",
  "Sample_ID",
  "sample",
  "Sample",
  "patient_id",
  "Patient_ID"
)

state_id <- possible_id_names[
  possible_id_names %in% names(states)
]

if (length(state_id) == 0) {
  state_id <- names(states)[1]
} else {
  state_id <- state_id[1]
}

# ----------------------------------------------------------------------
# Detect ecosystem state
# ----------------------------------------------------------------------

possible_state_names <- c(
  "ecosystem_state",
  "Ecosystem_State",
  "state",
  "State",
  "cluster",
  "cluster_numeric"
)

state_col_candidates <- possible_state_names[
  possible_state_names %in% names(states)
]

if ("ecosystem_state" %in% names(states)) {

  state_col <- "ecosystem_state"

} else if (length(state_col_candidates) > 0) {

  state_col <- state_col_candidates[1]

} else {

  stop(
    "ERROR: Could not identify ecosystem-state column."
  )
}

cat("\nSample ID column:", state_id, "\n")
cat("Ecosystem-state column:", state_col, "\n\n")


# ----------------------------------------------------------------------
# Keep only required columns
# ----------------------------------------------------------------------

state_df <- states %>%
  select(
    sample_id_internal = all_of(state_id),
    ecosystem_state_internal = all_of(state_col)
  )

names(state_df) <- c(
  "sample_id",
  "ecosystem_state"
)

state_df$sample_id <- as.character(
  state_df$sample_id
)

state_df$ecosystem_state <- as.character(
  state_df$ecosystem_state
)

# Remove duplicated sample IDs
state_df <- state_df %>%
  distinct(sample_id, .keep_all = TRUE)

cat("Unique patients:", nrow(state_df), "\n\n")


# ======================================================================
# SECTION 6 — STATE SUMMARY
# ======================================================================

cat("Frozen ecosystem-state structure:\n")

state_summary <- state_df %>%
  count(ecosystem_state, name = "n") %>%
  mutate(
    fraction = n / sum(n),
    percentage = 100 * fraction
  ) %>%
  arrange(ecosystem_state)

print(state_summary)

write_csv(
  state_summary,
  file.path(
    output_dir,
    "M3.4_state_sizes.csv"
  )
)


# ======================================================================
# SECTION 7 — READ MODULE ENTROPY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING STANDARDIZED MODULE ENTROPY\n")
cat("=====================================================================\n")

entropy <- read_csv(
  entropy_file,
  show_col_types = FALSE
)

cat("Rows:", nrow(entropy), "\n")
cat("Columns:", ncol(entropy), "\n\n")

cat("Detected columns:\n")
print(names(entropy))


# ----------------------------------------------------------------------
# Detect sample ID
# ----------------------------------------------------------------------

entropy_id_candidates <- possible_id_names[
  possible_id_names %in% names(entropy)
]

if (length(entropy_id_candidates) > 0) {

  entropy_id <- entropy_id_candidates[1]

} else {

  entropy_id <- names(entropy)[1]
}

cat("\nEntropy sample ID column:", entropy_id, "\n")


# ----------------------------------------------------------------------
# Identify module columns
# ----------------------------------------------------------------------

module_columns <- names(entropy)[
  grepl(
    "^S_module_",
    names(entropy)
  )
]

if (length(module_columns) == 0) {

  # fallback
  module_columns <- names(entropy)[
    sapply(
      entropy,
      function(x) is.numeric(x)
    )
  ]
}

cat(
  "Module entropy features detected:",
  length(module_columns),
  "\n"
)

if (length(module_columns) < 3) {
  stop(
    "ERROR: Fewer than 3 module features detected."
  )
}


# ======================================================================
# SECTION 8 — STANDARDIZE INTERNAL DATA STRUCTURE
# ======================================================================

entropy_df <- entropy %>%
  mutate(
    sample_id = as.character(
      .data[[entropy_id]]
    )
  ) %>%
  select(
    sample_id,
    all_of(module_columns)
  )

# ----------------------------------------------------------------------
# Convert modules safely to numeric
# ----------------------------------------------------------------------

for (mod in module_columns) {

  entropy_df[[mod]] <- suppressWarnings(
    as.numeric(entropy_df[[mod]])
  )
}

# ----------------------------------------------------------------------
# Remove duplicated IDs
# ----------------------------------------------------------------------

entropy_df <- entropy_df %>%
  distinct(sample_id, .keep_all = TRUE)

cat(
  "Unique patients in entropy data:",
  nrow(entropy_df),
  "\n"
)


# ======================================================================
# SECTION 9 — MERGE FROZEN STATES WITH MODULE FEATURES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("MERGING FROZEN STATES WITH MODULE FEATURES\n")
cat("=====================================================================\n")

validation_df <- inner_join(
  state_df,
  entropy_df,
  by = "sample_id"
)

cat(
  "Patients after merge:",
  nrow(validation_df),
  "\n"
)

if (nrow(validation_df) < 10) {
  stop(
    "ERROR: Too few patients after merging states and module data."
  )
}

validation_df$ecosystem_state <- factor(
  validation_df$ecosystem_state
)

cat("\nFinal state distribution:\n")
print(
  table(validation_df$ecosystem_state)
)


# ======================================================================
# SECTION 10 — MODULE QC
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("MODULE QC\n")
cat("=====================================================================\n")

module_qc_list <- vector(
  "list",
  length(module_columns)
)

for (i in seq_along(module_columns)) {

  mod <- module_columns[i]

  x <- validation_df[[mod]]

  module_qc_list[[i]] <- data.frame(
    module = mod,
    n = sum(!is.na(x)),
    missing = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

module_qc <- bind_rows(
  module_qc_list
)

print(
  module_qc
)

write_csv(
  module_qc,
  file.path(
    output_dir,
    "M3.4_module_QC.csv"
  )
)


# ----------------------------------------------------------------------
# Remove modules with insufficient information
# ----------------------------------------------------------------------

valid_modules <- module_qc %>%
  filter(
    missing == 0,
    is.finite(sd),
    sd > 0
  ) %>%
  pull(module)

cat(
  "\nModules retained:",
  length(valid_modules),
  "\n"
)

if (length(valid_modules) < 3) {
  stop(
    "ERROR: Insufficient valid modules."
  )
}


# ======================================================================
# SECTION 11 — CREATE ANALYSIS MATRIX
# ======================================================================

analysis_matrix <- as.matrix(
  validation_df[, valid_modules, drop = FALSE]
)

rownames(analysis_matrix) <- validation_df$sample_id

# ----------------------------------------------------------------------
# Ensure numeric
# ----------------------------------------------------------------------

storage.mode(analysis_matrix) <- "numeric"

# ----------------------------------------------------------------------
# Remove any rows containing non-finite values
# ----------------------------------------------------------------------

finite_rows <- apply(
  analysis_matrix,
  1,
  function(x) all(is.finite(x))
)

if (!all(finite_rows)) {

  cat(
    "Removing",
    sum(!finite_rows),
    "patients with non-finite module values.\n"
  )

  analysis_matrix <- analysis_matrix[
    finite_rows,
    ,
    drop = FALSE
  ]

  validation_df <- validation_df[
    match(
      rownames(analysis_matrix),
      validation_df$sample_id
    ),
    ,
    drop = FALSE
  ]
}


# ======================================================================
# SECTION 12 — MODULE-WISE STATE STATISTICS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("MODULE-WISE STATE EFFECTS\n")
cat("=====================================================================\n")

module_effect_list <- vector(
  "list",
  length(valid_modules)
)

for (i in seq_along(valid_modules)) {

  mod <- valid_modules[i]

  tmp <- validation_df[
    ,
    c("ecosystem_state", mod)
  ]

  names(tmp) <- c(
    "state",
    "value"
  )

  tmp <- tmp[
    is.finite(tmp$value),
    ,
    drop = FALSE
  ]

  kw <- kruskal.test(
    value ~ state,
    data = tmp
  )

  state_means <- tapply(
    tmp$value,
    tmp$state,
    mean,
    na.rm = TRUE
  )

  state_medians <- tapply(
    tmp$value,
    tmp$state,
    median,
    na.rm = TRUE
  )

  module_effect_list[[i]] <- data.frame(
    module = mod,
    statistic = as.numeric(
      kw$statistic
    ),
    p_value = as.numeric(
      kw$p.value
    ),
    stringsAsFactors = FALSE
  )

  # Add state means and medians dynamically
  for (st in names(state_means)) {

    safe_name <- paste0(
      "mean_",
      st
    )

    module_effect_list[[i]][[safe_name]] <- (
      unname(state_means[st])
    )

    safe_median_name <- paste0(
      "median_",
      st
    )

    module_effect_list[[i]][[safe_median_name]] <- (
      unname(state_medians[st])
    )
  }
}

module_effects <- bind_rows(
  module_effect_list
)

module_effects$FDR <- p.adjust(
  module_effects$p_value,
  method = "BH"
)

module_effects <- module_effects %>%
  arrange(FDR, p_value)

cat("\nTop module effects:\n")
print(
  head(
    module_effects,
    20
  )
)

write_csv(
  module_effects,
  file.path(
    output_dir,
    "M3.4_module_state_effects.csv"
  )
)


# ======================================================================
# SECTION 13 — MULTIVARIATE PERMANOVA
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("MULTIVARIATE ECOSYSTEM-STATE SEPARATION\n")
cat("=====================================================================\n")

cat(
  "Running PERMANOVA on standardized module entropy...\n"
)

set.seed(20260830)

distance_matrix <- dist(
  analysis_matrix,
  method = "euclidean"
)

permanova_result <- adonis2(
  analysis_matrix ~ ecosystem_state,
  data = validation_df,
  method = "euclidean",
  permutations = 999
)

print(
  permanova_result
)

# Save PERMANOVA table
permanova_df <- as.data.frame(
  permanova_result
)

permanova_df$term <- rownames(
  permanova_df
)

rownames(permanova_df) <- NULL

write_csv(
  permanova_df,
  file.path(
    output_dir,
    "M3.4_PERMANOVA.csv"
  )
)


# ======================================================================
# SECTION 14 — BETADISPER / DISPERSION TEST
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("TESTING WITHIN-STATE DISPERSION\n")
cat("=====================================================================\n")

dispersion <- betadisper(
  distance_matrix,
  validation_df$ecosystem_state
)

dispersion_anova <- anova(
  dispersion
)

cat("\nDispersion ANOVA:\n")
print(
  dispersion_anova
)

dispersion_perm <- permutest(
  dispersion,
  permutations = 999
)

cat("\nPermutation test:\n")
print(
  dispersion_perm
)

capture.output(
  dispersion_anova,
  file = file.path(
    output_dir,
    "M3.4_dispersion_ANOVA.txt"
  )
)

capture.output(
  dispersion_perm,
  file = file.path(
    output_dir,
    "M3.4_dispersion_permutation_test.txt"
  )
)


# ======================================================================
# SECTION 15 — SILHOUETTE ANALYSIS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("SILHOUETTE ANALYSIS\n")
cat("=====================================================================\n")

silhouette_object <- silhouette(
  as.numeric(
    validation_df$ecosystem_state
  ),
  distance_matrix
)

silhouette_values <- as.data.frame(
  silhouette_object
)

names(silhouette_values) <- c(
  "cluster",
  "neighbor",
  "silhouette_width"
)

silhouette_values$sample_id <- rownames(
  silhouette_values
)

# Add state labels
silhouette_values$ecosystem_state <- validation_df$ecosystem_state

mean_silhouette <- mean(
  silhouette_values$silhouette_width,
  na.rm = TRUE
)

cat(
  "\nMean silhouette width:",
  round(
    mean_silhouette,
    4
  ),
  "\n"
)

silhouette_by_state <- silhouette_values %>%
  group_by(
    ecosystem_state
  ) %>%
  summarise(
    n = n(),
    mean_silhouette = mean(
      silhouette_width,
      na.rm = TRUE
    ),
    median_silhouette = median(
      silhouette_width,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(
  silhouette_by_state
)

write_csv(
  silhouette_values,
  file.path(
    output_dir,
    "M3.4_patient_silhouette_values.csv"
  )
)

write_csv(
  silhouette_by_state,
  file.path(
    output_dir,
    "M3.4_silhouette_by_state.csv"
  )
)


# ======================================================================
# SECTION 16 — SILHOUETTE PLOT
# ======================================================================

silhouette_plot_df <- silhouette_values %>%
  arrange(
    ecosystem_state,
    silhouette_width
  ) %>%
  mutate(
    rank_within_state = row_number()
  )

silhouette_plot <- ggplot(
  silhouette_plot_df,
  aes(
    x = rank_within_state,
    y = silhouette_width,
    group = ecosystem_state
  )
) +
  geom_col() +
  facet_wrap(
    ~ ecosystem_state,
    scales = "free_x"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Silhouette Width of Frozen LIONESS Ecosystem States",
    x = "Patients within ecosystem state",
    y = "Silhouette width"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "M3.4_silhouette_plot.png"
  ),
  silhouette_plot,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "M3.4_silhouette_plot.pdf"
  ),
  silhouette_plot,
  width = 10,
  height = 6
)


# ======================================================================
# SECTION 17 — CENTROID DISTANCES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("ECOSYSTEM-STATE CENTROID DISTANCES\n")
cat("=====================================================================\n")

centroid_matrix <- matrix(
  NA_real_,
  nrow = nlevels(
    validation_df$ecosystem_state
  ),
  ncol = length(valid_modules)
)

rownames(centroid_matrix) <- levels(
  validation_df$ecosystem_state
)

colnames(centroid_matrix) <- valid_modules

for (st in levels(validation_df$ecosystem_state)) {

  idx <- validation_df$ecosystem_state == st

  centroid_matrix[st, ] <- colMeans(
    analysis_matrix[idx, , drop = FALSE],
    na.rm = TRUE
  )
}

centroid_distance <- as.matrix(
  dist(
    centroid_matrix,
    method = "euclidean"
  )
)

centroid_distance_df <- as.data.frame(
  as.table(
    centroid_distance
  )
)

names(centroid_distance_df) <- c(
  "state_1",
  "state_2",
  "centroid_distance"
)

centroid_distance_df <- centroid_distance_df %>%
  filter(
    as.character(state_1) <
      as.character(state_2)
  )

print(
  centroid_distance_df
)

write_csv(
  centroid_distance_df,
  file.path(
    output_dir,
    "M3.4_state_centroid_distances.csv"
  )
)


# ======================================================================
# SECTION 18 — WITHIN / BETWEEN STATE DISTANCES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("WITHIN-STATE AND BETWEEN-STATE DISTANCES\n")
cat("=====================================================================\n")

distance_vector <- as.vector(
  distance_matrix
)

distance_matrix_full <- as.matrix(
  distance_matrix
)

state_labels <- validation_df$ecosystem_state

within_between_list <- list()

counter <- 1

for (i in seq_len(nrow(distance_matrix_full) - 1)) {

  for (j in seq(
    i + 1,
    nrow(distance_matrix_full)
  )) {

    state_i <- as.character(
      state_labels[i]
    )

    state_j <- as.character(
      state_labels[j]
    )

    relation <- ifelse(
      state_i == state_j,
      "within_state",
      "between_state"
    )

    within_between_list[[counter]] <- data.frame(
      sample_i = rownames(
        distance_matrix_full
      )[i],
      sample_j = rownames(
        distance_matrix_full
      )[j],
      state_i = state_i,
      state_j = state_j,
      relation = relation,
      distance = distance_matrix_full[i, j]
    )

    counter <- counter + 1
  }
}

within_between_df <- bind_rows(
  within_between_list
)

distance_summary <- within_between_df %>%
  group_by(
    relation
  ) %>%
  summarise(
    n = n(),
    mean_distance = mean(
      distance,
      na.rm = TRUE
    ),
    median_distance = median(
      distance,
      na.rm = TRUE
    ),
    sd_distance = sd(
      distance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(
  distance_summary
)

write_csv(
  within_between_df,
  file.path(
    output_dir,
    "M3.4_pairwise_patient_distances.csv"
  )
)

write_csv(
  distance_summary,
  file.path(
    output_dir,
    "M3.4_within_between_distance_summary.csv"
  )
)


# ======================================================================
# SECTION 19 — PCA
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("PCA OF ECOSYSTEM MODULE ENTROPY\n")
cat("=====================================================================\n")

pca_result <- prcomp(
  analysis_matrix,
  center = TRUE,
  scale. = TRUE
)

pca_variance <- data.frame(
  PC = paste0(
    "PC",
    seq_along(
      pca_result$sdev
    )
  ),
  variance_percent =
    100 *
    (
      pca_result$sdev^2 /
      sum(
        pca_result$sdev^2
      )
    )
)

print(
  head(
    pca_variance,
    10
  )
)

write_csv(
  pca_variance,
  file.path(
    output_dir,
    "M3.4_PCA_variance_explained.csv"
  )
)

pca_scores <- as.data.frame(
  pca_result$x
)

pca_scores$sample_id <- rownames(
  pca_scores
)

pca_scores$ecosystem_state <- validation_df$ecosystem_state

write_csv(
  pca_scores,
  file.path(
    output_dir,
    "M3.4_PCA_scores.csv"
  )
)


# ----------------------------------------------------------------------
# PCA plot
# ----------------------------------------------------------------------

pc1_percent <- round(
  pca_variance$variance_percent[1],
  1
)

pc2_percent <- round(
  pca_variance$variance_percent[2],
  1
)

pca_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    shape = ecosystem_state
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  stat_ellipse(
    aes(
      group = ecosystem_state
    ),
    type = "norm",
    level = 0.90,
    linewidth = 0.7
  ) +
  labs(
    title = "PCA Separation of Frozen LIONESS Ecosystem States",
    x = paste0(
      "PC1 (",
      pc1_percent,
      "%)"
    ),
    y = paste0(
      "PC2 (",
      pc2_percent,
      "%)"
    ),
    shape = "Ecosystem state"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "M3.4_PCA_ecosystem_states.png"
  ),
  pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "M3.4_PCA_ecosystem_states.pdf"
  ),
  pca_plot,
  width = 8,
  height = 6
)


# ======================================================================
# SECTION 20 — PAIRWISE STATE SEPARATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("PAIRWISE ECOSYSTEM-STATE SEPARATION\n")
cat("=====================================================================\n")

state_levels <- levels(
  validation_df$ecosystem_state
)

pairwise_results <- list()

pair_counter <- 1

if (length(state_levels) >= 2) {

  for (i in seq_len(length(state_levels) - 1)) {

    for (j in seq(
      i + 1,
      length(state_levels)
    )) {

      st1 <- state_levels[i]
      st2 <- state_levels[j]

      idx <- validation_df$ecosystem_state %in% c(
        st1,
        st2
      )

      x <- validation_df[
        idx,
        ,
        drop = FALSE
      ]

      x$pair_state <- droplevels(
        factor(
          x$ecosystem_state
        )
      )

      pair_dist <- dist(
        as.matrix(
          x[, valid_modules, drop = FALSE]
        ),
        method = "euclidean"
      )

      pair_permanova <- adonis2(
        as.matrix(
          x[, valid_modules, drop = FALSE]
        ) ~ pair_state,
        data = x,
        method = "euclidean",
        permutations = 999
      )

      pairwise_results[[pair_counter]] <- data.frame(
        state_1 = st1,
        state_2 = st2,
        F_statistic = as.numeric(
          pair_permanova$F[1]
        ),
        R2 = as.numeric(
          pair_permanova$R2[1]
        ),
        p_value = as.numeric(
          pair_permanova$`Pr(>F)`[1]
        )
      )

      pair_counter <- pair_counter + 1
    }
  }
}

pairwise_state_results <- bind_rows(
  pairwise_results
)

pairwise_state_results$FDR <- p.adjust(
  pairwise_state_results$p_value,
  method = "BH"
)

print(
  pairwise_state_results
)

write_csv(
  pairwise_state_results,
  file.path(
    output_dir,
    "M3.4_pairwise_state_PERMANOVA.csv"
  )
)


# ======================================================================
# SECTION 21 — LEAVE-ONE-OUT NEAREST-CENTROID CLASSIFICATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LEAVE-ONE-OUT NEAREST-CENTROID CLASSIFICATION\n")
cat("=====================================================================\n")

n_patients <- nrow(
  analysis_matrix
)

predicted_states <- character(
  n_patients
)

actual_states <- as.character(
  validation_df$ecosystem_state
)

for (i in seq_len(n_patients)) {

  training_idx <- seq_len(n_patients) != i

  training_matrix <- analysis_matrix[
    training_idx,
    ,
    drop = FALSE
  ]

  training_states <- actual_states[
    training_idx
  ]

  test_vector <- analysis_matrix[
    i,
    ,
    drop = FALSE
  ]

  training_centroids <- list()

  for (st in unique(training_states)) {

    training_centroids[[st]] <- colMeans(
      training_matrix[
        training_states == st,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }

  centroid_distances <- numeric(
    length(training_centroids)
  )

  names(centroid_distances) <- names(
    training_centroids
  )

  for (st in names(training_centroids)) {

    centroid_distances[st] <- sqrt(
      sum(
        (
          test_vector -
          training_centroids[[st]]
        )^2
      )
    )
  }

  predicted_states[i] <- names(
    which.min(
      centroid_distances
    )
  )
}

classification_df <- data.frame(
  sample_id = validation_df$sample_id,
  actual_state = actual_states,
  predicted_state = predicted_states,
  stringsAsFactors = FALSE
)

classification_accuracy <- mean(
  classification_df$actual_state ==
    classification_df$predicted_state
)

cat(
  "\nLeave-one-out classification accuracy:",
  round(
    classification_accuracy,
    4
  ),
  "\n"
)

confusion_matrix <- table(
  Actual = classification_df$actual_state,
  Predicted = classification_df$predicted_state
)

cat("\nConfusion matrix:\n")
print(
  confusion_matrix
)

write_csv(
  classification_df,
  file.path(
    output_dir,
    "M3.4_leave_one_out_classification.csv"
  )
)

write.csv(
  as.data.frame.matrix(
    confusion_matrix
  ),
  file.path(
    output_dir,
    "M3.4_confusion_matrix.csv"
  )
)


# ======================================================================
# SECTION 22 — PER-STATE CLASSIFICATION ACCURACY
# ======================================================================

state_accuracy <- classification_df %>%
  group_by(
    actual_state
  ) %>%
  summarise(
    n = n(),
    correct = sum(
      actual_state ==
        predicted_state
    ),
    accuracy =
      correct / n,
    .groups = "drop"
  )

print(
  state_accuracy
)

write_csv(
  state_accuracy,
  file.path(
    output_dir,
    "M3.4_state_classification_accuracy.csv"
  )
)


# ======================================================================
# SECTION 23 — BOOTSTRAP STATE STABILITY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("BOOTSTRAP STATE STABILITY\n")
cat("=====================================================================\n")

set.seed(20260830)

n_bootstrap <- 500

bootstrap_accuracy <- numeric(
  n_bootstrap
)

bootstrap_state_sizes <- matrix(
  NA_real_,
  nrow = n_bootstrap,
  ncol = length(state_levels)
)

colnames(
  bootstrap_state_sizes
) <- state_levels

cat(
  "Bootstrap replicates:",
  n_bootstrap,
  "\n"
)

for (b in seq_len(n_bootstrap)) {

  sample_indices <- sample(
    seq_len(n_patients),
    size = n_patients,
    replace = TRUE
  )

  boot_matrix <- analysis_matrix[
    sample_indices,
    ,
    drop = FALSE
  ]

  boot_states <- actual_states[
    sample_indices
  ]

  # ------------------------------------------------------------
  # Bootstrap nearest-centroid resubstitution accuracy
  # ------------------------------------------------------------

  predicted_boot <- character(
    length(boot_states)
  )

  for (i in seq_along(boot_states)) {

    centroid_distances <- numeric(
      length(state_levels)
    )

    names(centroid_distances) <- state_levels

    for (st in state_levels) {

      idx_state <- boot_states == st

      if (sum(idx_state) > 0) {

        centroid <- colMeans(
          boot_matrix[
            idx_state,
            ,
            drop = FALSE
          ],
          na.rm = TRUE
        )

        centroid_distances[st] <- sqrt(
          sum(
            (
              boot_matrix[i, ] -
              centroid
            )^2
          )
        )

      } else {

        centroid_distances[st] <- Inf
      }
    }

    predicted_boot[i] <- names(
      which.min(
        centroid_distances
      )
    )
  }

  bootstrap_accuracy[b] <- mean(
    predicted_boot ==
      boot_states
  )

  # ------------------------------------------------------------
  # State proportions
  # ------------------------------------------------------------

  boot_table <- table(
    factor(
      boot_states,
      levels = state_levels
    )
  )

  bootstrap_state_sizes[b, ] <- (
    as.numeric(
      boot_table
    ) /
      length(boot_states)
  )
}


bootstrap_summary <- data.frame(
  metric = c(
    "mean_accuracy",
    "median_accuracy",
    "sd_accuracy",
    "lower_95_CI",
    "upper_95_CI"
  ),
  value = c(
    mean(
      bootstrap_accuracy,
      na.rm = TRUE
    ),
    median(
      bootstrap_accuracy,
      na.rm = TRUE
    ),
    sd(
      bootstrap_accuracy,
      na.rm = TRUE
    ),
    quantile(
      bootstrap_accuracy,
      0.025,
      na.rm = TRUE
    ),
    quantile(
      bootstrap_accuracy,
      0.975,
      na.rm = TRUE
    )
  )
)

print(
  bootstrap_summary
)

write_csv(
  bootstrap_summary,
  file.path(
    output_dir,
    "M3.4_bootstrap_stability_summary.csv"
  )
)


# ----------------------------------------------------------------------
# Bootstrap state proportion summary
# ----------------------------------------------------------------------

bootstrap_proportion_df <- as.data.frame(
  bootstrap_state_sizes
)

bootstrap_proportion_summary <- data.frame(
  ecosystem_state = state_levels,
  mean_fraction = colMeans(
    bootstrap_state_sizes,
    na.rm = TRUE
  ),
  lower_95_CI = apply(
    bootstrap_state_sizes,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  ),
  upper_95_CI = apply(
    bootstrap_state_sizes,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )
)

print(
  bootstrap_proportion_summary
)

write_csv(
  bootstrap_proportion_summary,
  file.path(
    output_dir,
    "M3.4_bootstrap_state_proportions.csv"
  )
)


# ======================================================================
# SECTION 24 — MODULE CONTRIBUTION TO STATE SEPARATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("MODULE CONTRIBUTION TO ECOSYSTEM-STATE SEPARATION\n")
cat("=====================================================================\n")

module_contribution_list <- vector(
  "list",
  length(valid_modules)
)

for (i in seq_along(valid_modules)) {

  mod <- valid_modules[i]

  x <- validation_df[
    ,
    c("ecosystem_state", mod)
  ]

  names(x) <- c(
    "state",
    "value"
  )

  x <- x[
    is.finite(x$value),
    ,
    drop = FALSE
  ]

  overall_mean <- mean(
    x$value,
    na.rm = TRUE
  )

  state_means <- tapply(
    x$value,
    x$state,
    mean,
    na.rm = TRUE
  )

  between_variance <- mean(
    (
      state_means -
      overall_mean
    )^2,
    na.rm = TRUE
  )

  within_variance <- mean(
    tapply(
      x$value,
      x$state,
      var,
      na.rm = TRUE
    ),
    na.rm = TRUE
  )

  contribution_ratio <- (
    between_variance /
      (
        between_variance +
        within_variance +
        1e-12
      )
  )

  module_contribution_list[[i]] <- data.frame(
    module = mod,
    between_state_variance =
      between_variance,
    within_state_variance =
      within_variance,
    contribution_ratio =
      contribution_ratio
  )
}

module_contribution <- bind_rows(
  module_contribution_list
) %>%
  arrange(
    desc(
      contribution_ratio
    )
  )

print(
  head(
    module_contribution,
    20
  )
)

write_csv(
  module_contribution,
  file.path(
    output_dir,
    "M3.4_module_contribution_to_state_separation.csv"
  )
)


# ======================================================================
# SECTION 25 — TOP MODULE CONTRIBUTION PLOT
# ======================================================================

top_n_modules <- min(
  15,
  nrow(
    module_contribution
  )
)

top_modules <- module_contribution[
  seq_len(
    top_n_modules
  ),
  ,
  drop = FALSE
]

top_modules$module <- factor(
  top_modules$module,
  levels = rev(
    top_modules$module
  )
)

contribution_plot <- ggplot(
  top_modules,
  aes(
    x = module,
    y = contribution_ratio
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top Modules Contributing to Ecosystem-State Separation",
    x = "Module",
    y = "Between-state contribution ratio"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "M3.4_top_module_contributions.png"
  ),
  contribution_plot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "M3.4_top_module_contributions.pdf"
  ),
  contribution_plot,
  width = 9,
  height = 7
)


# ======================================================================
# SECTION 26 — STATE MODULE HEATMAP DATA
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("STATE MODULE CENTROID PROFILE\n")
cat("=====================================================================\n")

state_profile_list <- list()

counter <- 1

for (st in state_levels) {

  idx <- validation_df$ecosystem_state == st

  means <- colMeans(
    analysis_matrix[
      idx,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  state_profile_list[[counter]] <- data.frame(
    ecosystem_state = st,
    module = names(means),
    mean_entropy = as.numeric(
      means
    )
  )

  counter <- counter + 1
}

state_profiles <- bind_rows(
  state_profile_list
)

write_csv(
  state_profiles,
  file.path(
    output_dir,
    "M3.4_state_module_profiles.csv"
  )
)


# ======================================================================
# SECTION 27 — HEATMAP
# ======================================================================

heatmap_plot <- ggplot(
  state_profiles,
  aes(
    x = module,
    y = ecosystem_state,
    fill = mean_entropy
  )
) +
  geom_tile() +
  scale_x_discrete(
    guide = guide_axis(
      angle = 90
    )
  ) +
  labs(
    title = "Frozen LIONESS Ecosystem-State Module Profiles",
    x = "Module",
    y = "Ecosystem state",
    fill = "Mean\nentropy"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "M3.4_state_module_heatmap.png"
  ),
  heatmap_plot,
  width = 14,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "M3.4_state_module_heatmap.pdf"
  ),
  heatmap_plot,
  width = 14,
  height = 5
)


# ======================================================================
# SECTION 28 — OVERALL QUANTITATIVE VALIDATION SUMMARY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.4 QUANTITATIVE VALIDATION SUMMARY\n")
cat("=====================================================================\n")

permanova_p <- as.numeric(
  permanova_result$`Pr(>F)`[1]
)

permanova_r2 <- as.numeric(
  permanova_result$R2[1]
)

mean_silhouette_value <- mean_silhouette

n_significant_modules <- sum(
  module_effects$FDR < 0.05,
  na.rm = TRUE
)

n_significant_modules_fdr01 <- sum(
  module_effects$FDR < 0.01,
  na.rm = TRUE
)

summary_df <- data.frame(
  metric = c(
    "n_patients",
    "n_modules",
    "n_ecosystem_states",
    "PERMANOVA_R2",
    "PERMANOVA_p_value",
    "mean_silhouette",
    "LOO_classification_accuracy",
    "bootstrap_mean_accuracy",
    "bootstrap_lower_95_CI",
    "bootstrap_upper_95_CI",
    "modules_FDR_lt_0.05",
    "modules_FDR_lt_0.01"
  ),
  value = c(
    nrow(validation_df),
    length(valid_modules),
    length(state_levels),
    permanova_r2,
    permanova_p,
    mean_silhouette_value,
    classification_accuracy,
    bootstrap_summary$value[
      bootstrap_summary$metric ==
        "mean_accuracy"
    ],
    bootstrap_summary$value[
      bootstrap_summary$metric ==
        "lower_95_CI"
    ],
    bootstrap_summary$value[
      bootstrap_summary$metric ==
        "upper_95_CI"
    ],
    n_significant_modules,
    n_significant_modules_fdr01
  )
)

print(
  summary_df
)

write_csv(
  summary_df,
  file.path(
    output_dir,
    "M3.4_validation_summary.csv"
  )
)


# ======================================================================
# SECTION 29 — SAVE COMPLETE ANALYSIS DATASET
# ======================================================================

write_csv(
  validation_df,
  file.path(
    output_dir,
    "M3.4_validation_dataset.csv"
  )
)


# ======================================================================
# SECTION 30 — FINAL REPORT
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.4 COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n\n")

cat("Patients analyzed:",
    nrow(validation_df),
    "\n")

cat("Modules analyzed:",
    length(valid_modules),
    "\n")

cat("Ecosystem states:",
    length(state_levels),
    "\n\n")

cat("PERMANOVA R2:",
    round(
      permanova_r2,
      4
    ),
    "\n")

cat("PERMANOVA P:",
    format.pval(
      permanova_p,
      digits = 4
    ),
    "\n")

cat("Mean silhouette:",
    round(
      mean_silhouette_value,
      4
    ),
    "\n")

cat("Leave-one-out classification accuracy:",
    round(
      classification_accuracy,
      4
    ),
    "\n")

cat("Bootstrap mean classification accuracy:",
    round(
      bootstrap_summary$value[
        bootstrap_summary$metric ==
          "mean_accuracy"
      ],
      4
    ),
    "\n")

cat(
  "Bootstrap 95% CI:",
  round(
    bootstrap_summary$value[
      bootstrap_summary$metric ==
        "lower_95_CI"
    ],
    4
  ),
  "-",
  round(
    bootstrap_summary$value[
      bootstrap_summary$metric ==
        "upper_95_CI"
    ],
    4
  ),
  "\n"
)

cat("\nModules with FDR < 0.05:",
    n_significant_modules,
    "/",
    length(valid_modules),
    "\n")

cat("Modules with FDR < 0.01:",
    n_significant_modules_fdr01,
    "/",
    length(valid_modules),
    "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Main output files:\n\n")

output_files <- c(
  "M3.4_state_sizes.csv",
  "M3.4_module_QC.csv",
  "M3.4_module_state_effects.csv",
  "M3.4_PERMANOVA.csv",
  "M3.4_dispersion_ANOVA.txt",
  "M3.4_dispersion_permutation_test.txt",
  "M3.4_patient_silhouette_values.csv",
  "M3.4_silhouette_by_state.csv",
  "M3.4_silhouette_plot.png",
  "M3.4_silhouette_plot.pdf",
  "M3.4_state_centroid_distances.csv",
  "M3.4_pairwise_patient_distances.csv",
  "M3.4_within_between_distance_summary.csv",
  "M3.4_PCA_variance_explained.csv",
  "M3.4_PCA_scores.csv",
  "M3.4_PCA_ecosystem_states.png",
  "M3.4_PCA_ecosystem_states.pdf",
  "M3.4_pairwise_state_PERMANOVA.csv",
  "M3.4_leave_one_out_classification.csv",
  "M3.4_confusion_matrix.csv",
  "M3.4_state_classification_accuracy.csv",
  "M3.4_bootstrap_stability_summary.csv",
  "M3.4_bootstrap_state_proportions.csv",
  "M3.4_module_contribution_to_state_separation.csv",
  "M3.4_top_module_contributions.png",
  "M3.4_top_module_contributions.pdf",
  "M3.4_state_module_profiles.csv",
  "M3.4_state_module_heatmap.png",
  "M3.4_state_module_heatmap.pdf",
  "M3.4_validation_summary.csv",
  "M3.4_validation_dataset.csv"
)

for (f in output_files) {

  full_path <- file.path(
    output_dir,
    f
  )

  if (file.exists(full_path)) {
    cat(" -", f, "\n")
  }
}

cat("\n")
cat("Finished:", as.character(Sys.time()), "\n\n")

cat("=====================================================================\n")
cat("QUANTITATIVE VALIDATION COMPLETE\n")
cat("=====================================================================\n")

