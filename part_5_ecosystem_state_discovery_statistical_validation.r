# ======================================================================
# M3.2 — LIONESS MODULE-LEVEL BIOLOGICAL AXES
#          AND PATIENT ECOSYSTEM-STATE DISCOVERY
#
#
# PURPOSE
# -------
# Determine whether module-level LIONESS entropy contains reproducible
# patient-specific biological structure.
#
# IMPORTANT
# ---------
# M3.1 demonstrated that the aggregate entropy measures:
#
#   S_global
#   S_local_mean
#   S_meso_mean
#
# are highly correlated and therefore should NOT be treated as
# independent ecosystem dimensions.
#
# M3.2 therefore focuses on the 29 usable MSigDB module-level
# entropy features.
#
# ANALYSIS
# --------
# 1. Module QC
# 2. Module variation
# 3. Module correlation
# 4. PCA
# 5. Biological module axes
# 6. Patient clustering
# 7. Silhouette analysis
# 8. Bootstrap ARI stability
# 9. Ecosystem-state assignment
# 10. State module profiles
#
# ======================================================================


# ======================================================================
# SECTION 1 — START
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.2 — LIONESS MODULE-LEVEL BIOLOGICAL AXES\n")
cat("         AND PATIENT ECOSYSTEM-STATE DISCOVERY\n")
cat("=====================================================================\n\n")

cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# SECTION 2 — PACKAGES
# ======================================================================

cat("Checking required R packages...\n")

required_packages <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "pheatmap",
  "cluster"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {

  if (!pkg %in% installed) {

    cat("Installing:", pkg, "\n")

    tryCatch(
      install.packages(
        pkg,
        repos = "https://cloud.r-project.org",
        dependencies = TRUE
      ),
      error = function(e) {
        cat(
          "WARNING: Could not install ",
          pkg,
          "\n",
          sep = ""
        )
      }
    )
  }
}

suppressPackageStartupMessages({

  library(ggplot2)
  library(dplyr)
  library(readr)
  library(pheatmap)
  library(cluster)

})

cat("Packages ready.\n\n")


# ======================================================================
# SECTION 3 — DIRECTORIES
# ======================================================================

input_dir <- "/content/M3_LIONESS_entropy"

diagnostic_dir <- file.path(
  input_dir,
  "M3.1_entropy_diagnostics"
)

output_dir <- file.path(
  input_dir,
  "M3.2_module_axes"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Input directory:\n")
cat(input_dir, "\n\n")

cat("Diagnostic directory:\n")
cat(diagnostic_dir, "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")


# ======================================================================
# SECTION 4 — INPUT FILES
# ======================================================================

module_entropy_file <- file.path(
  input_dir,
  "M3_LIONESS_module_entropy.csv"
)

# IMPORTANT:
# M3.1 actually created:
#
# M3.1_module_variation.csv
#
# NOT:
# M3.1_module_entropy_variation.csv

module_variation_file <- file.path(
  diagnostic_dir,
  "M3.1_module_variation.csv"
)

cat("=====================================================================\n")
cat("CHECKING INPUT FILES\n")
cat("=====================================================================\n\n")

if (!file.exists(module_entropy_file)) {

  stop(
    paste0(
      "\nERROR: Required module entropy file was not found:\n",
      module_entropy_file,
      "\n\nPlease make sure M3.1 generated:\n",
      "M3_LIONESS_module_entropy.csv\n"
    )
  )

}

cat("FOUND:\n")
cat(module_entropy_file, "\n")

if (file.exists(module_variation_file)) {

  cat("\nFOUND M3.1 module variation file:\n")
  cat(module_variation_file, "\n")

} else {

  cat("\nNOTE:\n")
  cat(
    "M3.1_module_variation.csv was not found.\n",
    "M3.2 will calculate module variation independently.\n"
  )
}

cat("\nInput check completed.\n\n")


# ======================================================================
# SECTION 5 — READ MODULE ENTROPY
# ======================================================================

cat("=====================================================================\n")
cat("READING MODULE ENTROPY MATRIX\n")
cat("=====================================================================\n")

module_entropy <- read_csv(
  module_entropy_file,
  show_col_types = FALSE
)

cat(
  "Rows:",
  nrow(module_entropy),
  "\n"
)

cat(
  "Columns:",
  ncol(module_entropy),
  "\n\n"
)


# ======================================================================
# SECTION 6 — IDENTIFY SAMPLE ID
# ======================================================================

sample_candidates <- c(
  "sample_id",
  "Sample_ID",
  "sample",
  "patient",
  "Patient"
)

existing_candidates <- sample_candidates[
  sample_candidates %in% colnames(module_entropy)
]

if (length(existing_candidates) > 0) {

  sample_col <- existing_candidates[1]

} else {

  sample_col <- colnames(module_entropy)[1]

}

cat(
  "Sample ID column:",
  sample_col,
  "\n\n"
)

module_entropy$sample_id <- as.character(
  module_entropy[[sample_col]]
)

module_entropy$sample_id <- trimws(
  module_entropy$sample_id
)

if (any(
  is.na(module_entropy$sample_id) |
  module_entropy$sample_id == ""
)) {

  stop(
    "ERROR: Missing or empty sample IDs detected."
  )

}


# ======================================================================
# SECTION 7 — CHECK DUPLICATE SAMPLE IDS
# ======================================================================

duplicate_ids <- duplicated(
  module_entropy$sample_id
)

if (any(duplicate_ids)) {

  cat(
    "WARNING:",
    sum(duplicate_ids),
    "duplicate sample IDs detected.\n"
  )

  cat(
    "Keeping the first occurrence of each sample.\n\n"
  )

  module_entropy <- module_entropy[
    !duplicate_ids,
    ,
    drop = FALSE
  ]

}


# ======================================================================
# SECTION 8 — IDENTIFY MODULE FEATURES
# ======================================================================

module_columns <- grep(
  "^S_module_",
  colnames(module_entropy),
  value = TRUE
)

cat(
  "Module entropy features detected:",
  length(module_columns),
  "\n\n"
)

if (length(module_columns) < 3) {

  stop(
    paste0(
      "ERROR: Fewer than 3 S_module_ features were detected.\n",
      "Detected: ",
      length(module_columns)
    )
  )

}

cat("Detected module features:\n")

print(module_columns)

cat("\n")


# ======================================================================
# SECTION 9 — FORCE NUMERIC VALUES
# ======================================================================

cat("Converting module entropy values to numeric...\n")

for (m in module_columns) {

  module_entropy[[m]] <- suppressWarnings(
    as.numeric(module_entropy[[m]])
  )

}

cat("Numeric conversion completed.\n\n")


# ======================================================================
# SECTION 10 — MODULE QC
# ======================================================================

cat("=====================================================================\n")
cat("MODULE-LEVEL QC\n")
cat("=====================================================================\n")

module_qc <- data.frame(
  module = module_columns,
  n = sapply(
    module_columns,
    function(x) {
      sum(is.finite(module_entropy[[x]]))
    }
  ),
  missing = sapply(
    module_columns,
    function(x) {
      sum(!is.finite(module_entropy[[x]]))
    }
  ),
  mean = sapply(
    module_columns,
    function(x) {

      y <- module_entropy[[x]]
      y <- y[is.finite(y)]

      if (length(y) > 0) {
        mean(y)
      } else {
        NA_real_
      }

    }
  ),
  sd = sapply(
    module_columns,
    function(x) {

      y <- module_entropy[[x]]
      y <- y[is.finite(y)]

      if (length(y) > 1) {
        sd(y)
      } else {
        NA_real_
      }

    }
  ),
  stringsAsFactors = FALSE
)

module_qc$missing_fraction <-
  module_qc$missing /
  nrow(module_entropy)

print(module_qc)

write_csv(
  module_qc,
  file.path(
    output_dir,
    "M3.2_module_QC.csv"
  )
)


# ======================================================================
# SECTION 11 — REMOVE BAD MODULES
# ======================================================================

usable_modules <- module_qc %>%
  filter(
    n >= 0.90 * nrow(module_entropy),
    is.finite(sd),
    sd > 0
  ) %>%
  pull(module)

cat("\n")
cat(
  "Modules before QC:",
  length(module_columns),
  "\n"
)

cat(
  "Modules after QC:",
  length(usable_modules),
  "\n\n"
)

if (length(usable_modules) < 3) {

  stop(
    "ERROR: Fewer than 3 usable module features remain after QC."
  )

}

cat("Usable modules:\n")
print(usable_modules)
cat("\n")


# ======================================================================
# SECTION 12 — CREATE CLEAN PATIENT × MODULE MATRIX
# ======================================================================

cat("Creating patient × module matrix...\n")

X <- module_entropy[
  ,
  c("sample_id", usable_modules),
  drop = FALSE
]

X_matrix <- as.matrix(
  X[
    ,
    usable_modules,
    drop = FALSE
  ]
)

rownames(X_matrix) <- X$sample_id


# ======================================================================
# SECTION 13 — SAFE MEDIAN IMPUTATION
# ======================================================================

cat("Checking missing module values...\n")

total_missing <- sum(
  !is.finite(X_matrix)
)

cat(
  "Invalid values before imputation:",
  total_missing,
  "\n"
)

for (j in seq_len(ncol(X_matrix))) {

  bad <- !is.finite(
    X_matrix[, j]
  )

  if (any(bad)) {

    med <- median(
      X_matrix[, j],
      na.rm = TRUE
    )

    if (!is.finite(med)) {

      stop(
        paste0(
          "ERROR: Module has no valid values: ",
          colnames(X_matrix)[j]
        )
      )

    }

    X_matrix[bad, j] <- med

  }

}

cat(
  "Invalid values after imputation:",
  sum(!is.finite(X_matrix)),
  "\n\n"
)


# ======================================================================
# SECTION 14 — MODULE VARIATION
# ======================================================================

cat("=====================================================================\n")
cat("MODULE BETWEEN-PATIENT VARIATION\n")
cat("=====================================================================\n")

module_stats <- data.frame(

  module = colnames(X_matrix),

  n = apply(
    X_matrix,
    2,
    function(x) sum(is.finite(x))
  ),

  mean = apply(
    X_matrix,
    2,
    mean
  ),

  sd = apply(
    X_matrix,
    2,
    sd
  ),

  median = apply(
    X_matrix,
    2,
    median
  ),

  IQR = apply(
    X_matrix,
    2,
    IQR
  ),

  MAD = apply(
    X_matrix,
    2,
    mad
  ),

  stringsAsFactors = FALSE
)

module_stats$CV <- abs(
  module_stats$sd /
    module_stats$mean
)

module_stats <- module_stats %>%
  arrange(
    desc(sd)
  )

print(
  module_stats
)

write_csv(
  module_stats,
  file.path(
    output_dir,
    "M3.2_module_variation.csv"
  )
)


# ======================================================================
# SECTION 15 — SAVE STANDARDIZED MATRIX
# ======================================================================

cat("\n")
cat("Standardizing module entropy features...\n")

Z <- scale(
  X_matrix,
  center = TRUE,
  scale = TRUE
)

Z <- as.matrix(Z)

if (any(!is.finite(Z))) {

  stop(
    "ERROR: Non-finite values remain after standardization."
  )

}

write.csv(
  Z,
  file = file.path(
    output_dir,
    "M3.2_standardized_module_entropy.csv"
  ),
  row.names = TRUE
)

cat("Standardization completed.\n\n")


# ======================================================================
# SECTION 16 — MODULE CORRELATION
# ======================================================================

cat("=====================================================================\n")
cat("MODULE ENTROPY CORRELATION\n")
cat("=====================================================================\n")

module_cor <- cor(
  Z,
  method = "spearman",
  use = "pairwise.complete.obs"
)

module_cor[
  !is.finite(module_cor)
] <- 0

diag(module_cor) <- 1

write.csv(
  module_cor,
  file = file.path(
    output_dir,
    "M3.2_module_correlation.csv"
  ),
  row.names = TRUE
)

cat(
  "Correlation matrix generated for",
  ncol(Z),
  "modules.\n\n"
)


# ======================================================================
# SECTION 17 — MODULE CORRELATION HEATMAP
# ======================================================================

pdf(
  file.path(
    output_dir,
    "M3.2_module_correlation_heatmap.pdf"
  ),
  width = 12,
  height = 10
)

pheatmap(
  module_cor,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  main = "LIONESS Module Entropy Correlation"
)

dev.off()


# ======================================================================
# SECTION 18 — PCA
# ======================================================================

cat("=====================================================================\n")
cat("PCA OF MODULE ENTROPY\n")
cat("=====================================================================\n")

pca <- prcomp(
  Z,
  center = FALSE,
  scale. = FALSE
)

variance <- pca$sdev^2

variance_explained <- variance /
  sum(variance)

cumulative_variance <- cumsum(
  variance_explained
)

pca_variance <- data.frame(

  PC = paste0(
    "PC",
    seq_along(
      variance_explained
    )
  ),

  variance_explained =
    variance_explained,

  cumulative_variance =
    cumulative_variance

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
    "M3.2_PCA_variance.csv"
  )
)


# ======================================================================
# SECTION 19 — PCA SCORES
# ======================================================================

pca_scores <- as.data.frame(
  pca$x
)

pca_scores$sample_id <- rownames(
  pca$x
)

pca_scores <- pca_scores %>%
  select(
    sample_id,
    everything()
  )

write_csv(
  pca_scores,
  file.path(
    output_dir,
    "M3.2_PCA_scores.csv"
  )
)


# ======================================================================
# SECTION 20 — PCA LOADINGS
# ======================================================================

pca_loadings <- as.data.frame(
  pca$rotation
)

pca_loadings$module <- rownames(
  pca$rotation
)

pca_loadings <- pca_loadings %>%
  select(
    module,
    everything()
  )

write_csv(
  pca_loadings,
  file.path(
    output_dir,
    "M3.2_PCA_loadings.csv"
  )
)


# ======================================================================
# SECTION 21 — PCA PATIENT STRUCTURE
# ======================================================================

pc1_pc2 <- data.frame(

  sample_id = rownames(
    pca$x
  ),

  PC1 = pca$x[, 1],

  PC2 = pca$x[, 2],

  stringsAsFactors = FALSE

)

p_pca <- ggplot(
  pc1_pc2,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(
    size = 2.5
  ) +
  theme_bw() +
  labs(

    title =
      "Patient Structure from LIONESS Module Entropy",

    x = paste0(
      "PC1 (",
      round(
        100 *
          variance_explained[1],
        1
      ),
      "%)"
    ),

    y = paste0(
      "PC2 (",
      round(
        100 *
          variance_explained[2],
        1
      ),
      "%)"
    )

  )

ggsave(
  file.path(
    output_dir,
    "M3.2_PCA_patient_structure.pdf"
  ),
  p_pca,
  width = 8,
  height = 6
)


# ======================================================================
# SECTION 22 — BIOLOGICAL AXES / PCA LOADINGS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("PCA BIOLOGICAL AXES\n")
cat("=====================================================================\n")

n_axes <- min(
  6,
  ncol(pca$rotation)
)

axis_summary <- list()

for (pc_index in seq_len(n_axes)) {

  pc_name <- paste0(
    "PC",
    pc_index
  )

  values <- pca$rotation[
    ,
    pc_index
  ]

  ordered <- sort(
    values,
    decreasing = TRUE
  )

  top_positive <- head(
    ordered,
    10
  )

  ordered_negative <- sort(
    values,
    decreasing = FALSE
  )

  top_negative <- head(
    ordered_negative,
    10
  )

  cat("\n")
  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    pc_name,
    "variance explained:",
    round(
      100 *
        variance_explained[pc_index],
      2
    ),
    "%\n"
  )

  cat(
    "Top positive modules:\n"
  )

  print(
    data.frame(
      module = names(top_positive),
      loading = as.numeric(
        top_positive
      ),
      row.names = NULL
    )
  )

  cat(
    "Top negative modules:\n"
  )

  print(
    data.frame(
      module = names(top_negative),
      loading = as.numeric(
        top_negative
      ),
      row.names = NULL
    )
  )

  axis_summary[[pc_name]] <- data.frame(

    PC = pc_name,

    variance_explained =
      variance_explained[pc_index],

    module =
      names(c(
        top_positive,
        top_negative
      )),

    loading =
      as.numeric(
        c(
          top_positive,
          top_negative
        )
      ),

    direction = c(
      rep(
        "positive",
        length(top_positive)
      ),
      rep(
        "negative",
        length(top_negative)
      )
    ),

    stringsAsFactors = FALSE

  )

}

axis_summary_df <- bind_rows(
  axis_summary
)

write_csv(
  axis_summary_df,
  file.path(
    output_dir,
    "M3.2_biological_axes.csv"
  )
)


# ======================================================================
# SECTION 23 — MODULE HIERARCHICAL CLUSTERING
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("HIERARCHICAL CLUSTERING OF BIOLOGICAL MODULES\n")
cat("=====================================================================\n")

module_distance <- as.dist(
  1 - module_cor
)

module_distance[
  !is.finite(
    module_distance
  )
] <- 1

module_hclust <- hclust(
  module_distance,
  method = "average"
)

pdf(
  file.path(
    output_dir,
    "M3.2_module_hierarchical_clustering.pdf"
  ),
  width = 12,
  height = 8
)

plot(
  module_hclust,
  main =
    "Hierarchical Clustering of Module Entropy",
  xlab = "",
  sub = "",
  cex = 0.7
)

dev.off()


# ======================================================================
# SECTION 24 — PATIENT DISTANCE
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("PATIENT DISTANCE MATRIX\n")
cat("=====================================================================\n")

patient_distance <- dist(
  Z,
  method = "euclidean"
)

if (
  any(
    !is.finite(
      as.numeric(
        patient_distance
      )
    )
  )
) {

  stop(
    "ERROR: Patient distance matrix contains invalid values."
  )

}

cat(
  "Patient distance matrix created.\n\n"
)


# ======================================================================
# SECTION 25 — PATIENT HIERARCHICAL CLUSTERING
# ======================================================================

patient_hclust <- hclust(
  patient_distance,
  method = "ward.D2"
)

pdf(
  file.path(
    output_dir,
    "M3.2_patient_hierarchical_clustering.pdf"
  ),
  width = 12,
  height = 8
)

plot(
  patient_hclust,
  labels = FALSE,
  main =
    "Patient Clustering from Module Entropy",
  xlab = "",
  sub = "",
  hang = -1
)

dev.off()


# ======================================================================
# SECTION 26 — PATIENT MODULE HEATMAP
# ======================================================================

pdf(
  file.path(
    output_dir,
    "M3.2_patient_module_entropy_heatmap.pdf"
  ),
  width = 12,
  height = 12
)

pheatmap(
  Z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  fontsize_col = 8,
  main =
    "Standardized Patient Module Entropy"
)

dev.off()


# ======================================================================
# SECTION 27 — CLUSTER METRICS FUNCTION
# ======================================================================

calculate_cluster_metrics <- function(
  matrix_data,
  k
) {

  set.seed(
    1000 + k
  )

  km <- kmeans(
    matrix_data,
    centers = k,
    nstart = 100,
    iter.max = 1000
  )

  distance_matrix <- dist(
    matrix_data
  )

  sil <- silhouette(
    km$cluster,
    distance_matrix
  )

  mean_silhouette <- mean(
    sil[, 3]
  )

  cluster_sizes <- table(
    km$cluster
  )

  list(
    model = km,
    silhouette = mean_silhouette,
    sizes = cluster_sizes
  )

}


# ======================================================================
# SECTION 28 — TEST K = 2–6
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("PATIENT ECOSYSTEM-STATE CLUSTERING\n")
cat("=====================================================================\n")

cluster_results <- list()

cluster_summary <- data.frame(
  k = integer(),
  silhouette = numeric(),
  minimum_cluster_size = integer(),
  maximum_cluster_size = integer(),
  stringsAsFactors = FALSE
)

for (k in 2:6) {

  cat(
    "\nTesting k =",
    k,
    "\n"
  )

  result <- calculate_cluster_metrics(
    Z,
    k
  )

  cluster_results[[paste0("k", k)]] <- result

  sizes <- result$sizes

  cluster_summary <- rbind(
    cluster_summary,
    data.frame(

      k = k,

      silhouette =
        result$silhouette,

      minimum_cluster_size =
        min(sizes),

      maximum_cluster_size =
        max(sizes),

      stringsAsFactors = FALSE

    )
  )

  cat(
    "Silhouette:",
    round(
      result$silhouette,
      4
    ),
    "\n"
  )

  cat(
    "Cluster sizes:",
    paste(
      as.numeric(sizes),
      collapse = ", "
    ),
    "\n"
  )

}

print(
  cluster_summary
)

write_csv(
  cluster_summary,
  file.path(
    output_dir,
    "M3.2_cluster_quality_k2_k6.csv"
  )
)


# ======================================================================
# SECTION 29 — SELECT BEST K
# ======================================================================

best_k <- cluster_summary$k[
  which.max(
    cluster_summary$silhouette
  )
]

cat("\n")
cat(
  "Best k by mean silhouette:",
  best_k,
  "\n"
)


# ======================================================================
# SECTION 30 — FINAL CLUSTER MODEL
# ======================================================================

cat("\n")
cat("Creating final ecosystem-state model...\n")

# IMPORTANT:
# Keep the [[ ... ]] operator on ONE expression.
# This avoids the previous parser error.

final_model <- cluster_results[[paste0("k", best_k)]][["model"]]

if (is.null(final_model)) {

  stop(
    "ERROR: Final clustering model could not be retrieved."
  )

}

cat(
  "Final clustering model retrieved successfully.\n"
)


# ======================================================================
# SECTION 31 — PATIENT ECOSYSTEM STATES
# ======================================================================

patient_states <- data.frame(

  sample_id =
    rownames(Z),

  ecosystem_state =
    paste0(
      "State_",
      final_model$cluster
    ),

  cluster_numeric =
    final_model$cluster,

  stringsAsFactors = FALSE

)

write_csv(
  patient_states,
  file.path(
    output_dir,
    "M3.2_patient_ecosystem_states.csv"
  )
)

cat(
  "Patient ecosystem states saved.\n"
)


# ======================================================================
# SECTION 32 — STATE SIZES
# ======================================================================

state_sizes <- patient_states %>%
  count(
    ecosystem_state,
    name = "n"
  ) %>%
  mutate(
    fraction =
      n /
      sum(n)
  )

cat("\n")
cat("Ecosystem-state sizes:\n")
print(state_sizes)

write_csv(
  state_sizes,
  file.path(
    output_dir,
    "M3.2_state_sizes.csv"
  )
)


# ======================================================================
# SECTION 33 — STATE MODULE PROFILES
# ======================================================================

state_profile_data <- as.data.frame(
  Z
)

state_profile_data$sample_id <-
  rownames(Z)

state_profile_data <- state_profile_data %>%
  left_join(
    patient_states,
    by = "sample_id"
  )

state_module_profiles <-
  state_profile_data %>%
  group_by(
    ecosystem_state
  ) %>%
  summarise(
    across(
      all_of(
        usable_modules
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

write_csv(
  state_module_profiles,
  file.path(
    output_dir,
    "M3.2_state_module_profiles.csv"
  )
)


# ======================================================================
# SECTION 34 — STATE PROFILE HEATMAP
# ======================================================================

state_profile_matrix <- as.matrix(
  state_module_profiles[
    ,
    usable_modules,
    drop = FALSE
  ]
)

rownames(
  state_profile_matrix
) <-
  state_module_profiles$ecosystem_state

pdf(
  file.path(
    output_dir,
    "M3.2_state_module_profile_heatmap.pdf"
  ),
  width = 12,
  height = 7
)

pheatmap(
  state_profile_matrix,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  scale = "none",
  main =
    paste0(
      "Ecosystem State Module Profiles (k=",
      best_k,
      ")"
    )
)

dev.off()


# ======================================================================
# SECTION 35 — ADJUSTED RAND INDEX FUNCTION
# ======================================================================

adjusted_rand_index <- function(
  labels1,
  labels2
) {

  tab <- table(
    labels1,
    labels2
  )

  n <- sum(tab)

  if (n < 2) {
    return(NA_real_)
  }

  choose2 <- function(x) {
    x * (x - 1) / 2
  }

  sum_comb <- sum(
    choose2(tab)
  )

  row_comb <- sum(
    choose2(
      rowSums(tab)
    )
  )

  col_comb <- sum(
    choose2(
      colSums(tab)
    )
  )

  total_comb <- choose2(n)

  if (total_comb == 0) {
    return(NA_real_)
  }

  expected <-
    (
      row_comb *
      col_comb
    ) /
    total_comb

  max_index <-
    (
      row_comb +
      col_comb
    ) /
    2

  denominator <-
    max_index -
    expected

  if (
    abs(denominator) < .Machine$double.eps
  ) {

    if (
      abs(
        sum_comb -
        expected
      ) < .Machine$double.eps
    ) {

      return(1)

    } else {

      return(0)

    }

  }

  (
    sum_comb -
    expected
  ) /
  denominator

}


# ======================================================================
# SECTION 36 — BOOTSTRAP ARI
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("BOOTSTRAP CLUSTER STABILITY\n")
cat("=====================================================================\n")

set.seed(
  20260830
)

n_boot <- 500

bootstrap_ari <- data.frame(

  replicate = integer(),

  ARI = numeric(),

  stringsAsFactors = FALSE

)

original_labels <-
  final_model$cluster

for (b in seq_len(n_boot)) {

  sample_indices <- sample(
    seq_len(
      nrow(Z)
    ),
    size =
      nrow(Z),
    replace = TRUE
  )

  boot_matrix <- Z[
    sample_indices,
    ,
    drop = FALSE
  ]

  boot_model <- tryCatch(

    kmeans(
      boot_matrix,
      centers = best_k,
      nstart = 25,
      iter.max = 500
    ),

    error = function(e) {
      NULL
    }

  )

  if (!is.null(boot_model)) {

    original_boot_labels <-
      original_labels[
        sample_indices
      ]

    ari <- adjusted_rand_index(
      original_boot_labels,
      boot_model$cluster
    )

    bootstrap_ari <- rbind(

      bootstrap_ari,

      data.frame(
        replicate = b,
        ARI = ari,
        stringsAsFactors = FALSE
      )

    )

  }

}

bootstrap_ari <-
  bootstrap_ari[
    is.finite(
      bootstrap_ari$ARI
    ),
    ,
    drop = FALSE
  ]

cat(
  "\nValid bootstrap replicates:",
  nrow(bootstrap_ari),
  "\n"
)

if (nrow(bootstrap_ari) == 0) {

  stop(
    "ERROR: No valid bootstrap ARI values were produced."
  )

}

mean_ari <- mean(
  bootstrap_ari$ARI
)

median_ari <- median(
  bootstrap_ari$ARI
)

sd_ari <- sd(
  bootstrap_ari$ARI
)

ari_q05 <- quantile(
  bootstrap_ari$ARI,
  0.05
)

ari_q95 <- quantile(
  bootstrap_ari$ARI,
  0.95
)

cat(
  "Mean ARI:",
  round(
    mean_ari,
    4
  ),
  "\n"
)

cat(
  "Median ARI:",
  round(
    median_ari,
    4
  ),
  "\n"
)

cat(
  "SD ARI:",
  round(
    sd_ari,
    4
  ),
  "\n"
)

cat(
  "5th percentile:",
  round(
    ari_q05,
    4
  ),
  "\n"
)

cat(
  "95th percentile:",
  round(
    ari_q95,
    4
  ),
  "\n"
)

write_csv(
  bootstrap_ari,
  file.path(
    output_dir,
    "M3.2_bootstrap_ARI.csv"
  )
)


# ======================================================================
# SECTION 37 — BOOTSTRAP ARI DISTRIBUTION
# ======================================================================

p_ari <- ggplot(
  bootstrap_ari,
  aes(
    x = ARI
  )
) +
  geom_histogram(
    bins = 30
  ) +
  theme_bw() +
  labs(

    title =
      paste0(
        "Bootstrap Cluster Stability (k=",
        best_k,
        ")"
      ),

    x = "Adjusted Rand Index",

    y =
      "Number of bootstrap replicates"

  )

ggsave(
  file.path(
    output_dir,
    "M3.2_bootstrap_ARI_distribution.pdf"
  ),
  p_ari,
  width = 8,
  height = 6
)


# ======================================================================
# SECTION 38 — PCA WITH ECOSYSTEM STATES
# ======================================================================

pc_state_plot <- pca_scores %>%
  left_join(
    patient_states,
    by = "sample_id"
  )

p_state <- ggplot(
  pc_state_plot,
  aes(
    x = PC1,
    y = PC2,
    shape = ecosystem_state
  )
) +
  geom_point(
    size = 3
  ) +
  theme_bw() +
  labs(

    title =
      paste0(
        "LIONESS Module Entropy Ecosystem States (k=",
        best_k,
        ")"
      ),

    x = paste0(
      "PC1 (",
      round(
        100 *
          variance_explained[1],
        1
      ),
      "%)"
    ),

    y = paste0(
      "PC2 (",
      round(
        100 *
          variance_explained[2],
        1
      ),
      "%)"
    ),

    shape =
      "Ecosystem state"

  )

ggsave(
  file.path(
    output_dir,
    "M3.2_PCA_ecosystem_states.pdf"
  ),
  p_state,
  width = 9,
  height = 7
)


# ======================================================================
# SECTION 39 — TOP MODULES BY STATE
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("TOP MODULES CHARACTERIZING EACH ECOSYSTEM STATE\n")
cat("=====================================================================\n")

top_state_modules <- list()

state_names <- unique(
  patient_states$ecosystem_state
)

for (state in state_names) {

  state_data <-
    state_module_profiles %>%
    filter(
      ecosystem_state == state
    )

  state_values <- as.numeric(
    state_data[
      ,
      usable_modules,
      drop = TRUE
    ]
  )

  names(state_values) <-
    usable_modules

  ordered <- sort(
    state_values,
    decreasing = TRUE
  )

  top10 <- head(
    ordered,
    10
  )

  cat("\n")
  cat(
    state,
    ":\n"
  )

  state_table <- data.frame(

    ecosystem_state =
      state,

    module =
      names(top10),

    standardized_entropy =
      as.numeric(top10),

    stringsAsFactors = FALSE

  )

  print(
    state_table
  )

  top_state_modules[[state]] <-
    state_table

}

top_state_modules_df <-
  bind_rows(
    top_state_modules
  )

write_csv(
  top_state_modules_df,
  file.path(
    output_dir,
    "M3.2_top_modules_by_state.csv"
  )
)


# ======================================================================
# SECTION 40 — STATE SEPARATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("STATE SEPARATION\n")
cat("=====================================================================\n")

state_silhouette <- silhouette(
  final_model$cluster,
  patient_distance
)

state_silhouette_df <- data.frame(

  sample_id =
    rownames(Z),

  ecosystem_state =
    patient_states$ecosystem_state,

  silhouette_width =
    state_silhouette[, 3],

  stringsAsFactors = FALSE

)

write_csv(
  state_silhouette_df,
  file.path(
    output_dir,
    "M3.2_patient_silhouette.csv"
  )
)

mean_silhouette <- mean(
  state_silhouette_df$silhouette_width,
  na.rm = TRUE
)

cat(
  "Mean silhouette:",
  round(
    mean_silhouette,
    4
  ),
  "\n"
)


# ======================================================================
# SECTION 41 — STATE-SPECIFIC SILHOUETTE SUMMARY
# ======================================================================

state_silhouette_summary <-
  state_silhouette_df %>%
  group_by(
    ecosystem_state
  ) %>%
  summarise(

    n = n(),

    mean_silhouette =
      mean(
        silhouette_width,
        na.rm = TRUE
      ),

    median_silhouette =
      median(
        silhouette_width,
        na.rm = TRUE
      ),

    negative_fraction =
      mean(
        silhouette_width < 0,
        na.rm = TRUE
      ),

    .groups = "drop"

  )

write_csv(
  state_silhouette_summary,
  file.path(
    output_dir,
    "M3.2_state_silhouette_summary.csv"
  )
)


# ======================================================================
# SECTION 42 — COMPREHENSIVE SUMMARY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.2 FINAL SUMMARY\n")
cat("=====================================================================\n")

pc1_percent <-
  100 *
  variance_explained[1]

pc2_percent <-
  100 *
  variance_explained[2]

pc3_percent <- NA_real_

if (
  length(
    variance_explained
  ) >= 3
) {

  pc3_percent <-
    100 *
    variance_explained[3]

}

summary_table <- data.frame(

  patients =
    nrow(Z),

  modules_before_QC =
    length(module_columns),

  modules_after_QC =
    length(usable_modules),

  PC1_percent =
    pc1_percent,

  PC2_percent =
    pc2_percent,

  PC3_percent =
    pc3_percent,

  best_k =
    best_k,

  best_k_silhouette =
    cluster_summary$silhouette[
      cluster_summary$k == best_k
    ],

  mean_bootstrap_ARI =
    mean_ari,

  median_bootstrap_ARI =
    median_ari,

  bootstrap_ARI_5pct =
    as.numeric(ari_q05),

  bootstrap_ARI_95pct =
    as.numeric(ari_q95),

  mean_patient_silhouette =
    mean_silhouette,

  stringsAsFactors = FALSE

)

print(
  summary_table
)

write_csv(
  summary_table,
  file.path(
    output_dir,
    "M3.2_summary.csv"
  )
)


# ======================================================================
# SECTION 43 — INTERPRETATION GUIDE
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.2 INTERPRETATION GUIDE\n")
cat("=====================================================================\n\n")

cat(
  "M3.2 evaluates whether module-level LIONESS entropy\n",
  "contains reproducible patient-specific biological structure.\n\n"
)

cat(
  "The analysis includes:\n\n"
)

cat(
  "1. Module variability\n",
  "   Identifies modules with the greatest patient-to-patient\n",
  "   entropy variation.\n\n"
)

cat(
  "2. Module correlation\n",
  "   Identifies coordinated biological entropy programs.\n\n"
)

cat(
  "3. PCA\n",
  "   Determines whether a small number of biological axes\n",
  "   explain the module-level variation.\n\n"
)

cat(
  "4. Patient clustering\n",
  "   Tests whether patients form distinct entropy profiles.\n\n"
)

cat(
  "5. Silhouette\n",
  "   Measures separation and within-state cohesion.\n\n"
)

cat(
  "6. Bootstrap ARI\n",
  "   Tests stability of ecosystem-state assignments under\n",
  "   patient resampling.\n\n"
)

cat(
  "7. State profiles\n",
  "   Identifies modules characterizing each candidate state.\n\n"
)

cat(
  "IMPORTANT:\n",
  "The resulting states are candidate biological states.\n",
  "They should NOT be called clinically validated states until\n",
  "they are independently validated using immune, stromal,\n",
  "tumor, molecular or clinical variables.\n\n"
)


# ======================================================================
# SECTION 44 — DIAGNOSTIC THRESHOLDS
# ======================================================================

cat(
  "GENERAL DIAGNOSTIC GUIDANCE:\n\n"
)

cat(
  "Silhouette:\n",
  "  >0.50       strong separation\n",
  "  0.25–0.50   moderate structure\n",
  "  <0.25       weak/uncertain structure\n\n"
)

cat(
  "Bootstrap ARI:\n",
  "  >0.75       strong stability\n",
  "  0.50–0.75   moderate stability\n",
  "  <0.50       unstable clustering\n\n"
)

cat(
  "These are diagnostic guidelines rather than formal\n",
  "statistical significance thresholds.\n\n"
)


# ======================================================================
# SECTION 45 — OUTPUT FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.2 OUTPUT FILES\n")
cat("=====================================================================\n")

output_files <- list.files(
  output_dir,
  full.names = TRUE
)

for (i in seq_along(output_files)) {

  cat(
    i,
    ".",
    output_files[i],
    "\n"
  )

}


# ======================================================================
# SECTION 46 — FINAL CONCLUSION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.2 COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n\n")

cat(
  "Patients analyzed:",
  nrow(Z),
  "\n"
)

cat(
  "Modules before QC:",
  length(module_columns),
  "\n"
)

cat(
  "Usable biological modules:",
  length(usable_modules),
  "\n"
)

cat(
  "Best ecosystem-state k:",
  best_k,
  "\n"
)

cat(
  "PC1 variance:",
  round(
    pc1_percent,
    2
  ),
  "%\n"
)

cat(
  "PC2 variance:",
  round(
    pc2_percent,
    2
  ),
  "%\n"
)

cat(
  "Mean bootstrap ARI:",
  round(
    mean_ari,
    3
  ),
  "\n"
)

cat(
  "Mean silhouette:",
  round(
    mean_silhouette,
    3
  ),
  "\n\n"
)

cat(
  "M3.2 has completed the module-level assessment of\n",
  "patient-specific LIONESS entropy structure.\n\n"
)

cat(
  "The resulting ecosystem states are candidate states and\n",
  "require independent biological validation in subsequent\n",
  "analyses.\n"
)

cat("\n")
cat(
  "Finished:",
  as.character(Sys.time()),
  "\n"
)

cat("\n")
cat("=====================================================================\n")
cat("M3.2 FINISHED\n")
cat("=====================================================================\n")