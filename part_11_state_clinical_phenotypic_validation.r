# ======================================================================
# M3.5.2 — ORTHOGONAL BIOLOGICAL VALIDATION OF FROZEN ECOSYSTEM STATES
# PURPOSE
# -------
# Independently characterize the already-frozen M3.2 ecosystem states
# using biological variables that were NOT used to define/redefine
# those states.
#
# IMPORTANT:
#   1. M3.2 ecosystem states are FROZEN.
#   2. No clustering is performed.
#   3. No ecosystem state is redefined.
#   4. No survival analysis is performed.
#   5. No univariate hazard ratios are performed.
#   6. No multivariate hazard ratios are performed.
#   7. No Cox regression is performed.
#
# MAIN QUESTIONS
# --------------
# A. Do independent biological variables differ among frozen states?
# B. Which biological features show the strongest state separation?
# C. Are pairwise state differences present?
# D. Do the independent biological variables reproduce separation
#    of the frozen ecosystem states?
#
# INPUT
# -----
# Frozen M3.2 state file:
# /content/M3_LIONESS_entropy/M3.2_module_axes/
# M3.2_patient_ecosystem_states.csv
#
# OUTPUT
# ------
# /content/M3_LIONESS_entropy/M3.5_clinical_validation/
# M3.5.2_orthogonal_biological_validation/
#
# ======================================================================


# ======================================================================
# 0. START
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.5.2 — ORTHOGONAL BIOLOGICAL VALIDATION\n")
cat("=====================================================================\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("\n")


# ======================================================================
# 1. PACKAGE MANAGEMENT
# ======================================================================

required_packages <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "tidyr"
)

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    cat("Installing package:", pkg, "\n")

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})


# ======================================================================
# 2. PATHS
# ======================================================================

frozen_file <- paste0(
  "/content/M3_LIONESS_entropy/",
  "M3.2_module_axes/",
  "M3.2_patient_ecosystem_states.csv"
)

m3_root <- "/content/M3_LIONESS_entropy"

output_dir <- paste0(
  m3_root,
  "/M3.5_clinical_validation/",
  "M3.5.2_orthogonal_biological_validation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Frozen M3.2 file:\n")
cat(frozen_file, "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")


# ======================================================================
# 3. CHECK FROZEN STATE FILE
# ======================================================================

cat("=====================================================================\n")
cat("CHECKING FROZEN M3.2 FILE\n")
cat("=====================================================================\n")

if (!file.exists(frozen_file)) {
  stop(
    paste0(
      "\nERROR: Frozen M3.2 file was not found:\n",
      frozen_file
    )
  )
}

cat("Frozen M3.2 file found.\n\n")


# ======================================================================
# 4. READ FROZEN STATES
# ======================================================================

cat("=====================================================================\n")
cat("READING FROZEN M3.2 STATES\n")
cat("=====================================================================\n")

frozen <- read_csv(
  frozen_file,
  show_col_types = FALSE
)

cat("Rows:", nrow(frozen), "\n")
cat("Columns:", ncol(frozen), "\n\n")

print(names(frozen))


# ======================================================================
# 5. VALIDATE REQUIRED COLUMNS
# ======================================================================

required_frozen_columns <- c(
  "sample_id",
  "ecosystem_state"
)

missing_columns <- setdiff(
  required_frozen_columns,
  names(frozen)
)

if (length(missing_columns) > 0) {

  stop(
    paste0(
      "\nERROR: Required frozen-state columns are missing:\n",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# ======================================================================
# 6. STANDARDIZE FROZEN PATIENT IDS
# ======================================================================

frozen$sample_id <- as.character(frozen$sample_id)

frozen$ecosystem_state <- as.character(
  frozen$ecosystem_state
)

frozen$sample_id_clean <- trimws(
  frozen$sample_id
)

frozen <- frozen[
  !is.na(frozen$sample_id_clean) &
  frozen$sample_id_clean != "",
]

frozen <- frozen[
  !is.na(frozen$ecosystem_state) &
  frozen$ecosystem_state != "",
]

if (anyDuplicated(frozen$sample_id_clean) > 0) {

  cat("WARNING: Duplicate frozen patient IDs detected.\n")
  cat("Keeping first occurrence only.\n")

  frozen <- frozen[
    !duplicated(frozen$sample_id_clean),
  ]
}

frozen$ecosystem_state <- factor(
  frozen$ecosystem_state
)

cat("\nValid frozen patients:", nrow(frozen), "\n")
cat(
  "Frozen ecosystem states:",
  nlevels(frozen$ecosystem_state),
  "\n\n"
)


# ======================================================================
# 7. FROZEN STATE DISTRIBUTION
# ======================================================================

cat("=====================================================================\n")
cat("FROZEN ECOSYSTEM STATE DISTRIBUTION\n")
cat("=====================================================================\n")

state_table <- as.data.frame(
  table(frozen$ecosystem_state),
  stringsAsFactors = FALSE
)

names(state_table) <- c(
  "ecosystem_state",
  "n"
)

state_table$percentage <- (
  state_table$n /
    sum(state_table$n)
) * 100

print(state_table)

write_csv(
  state_table,
  file.path(
    output_dir,
    "M3.5.2_frozen_state_distribution.csv"
  )
)


# ======================================================================
# 8. SEARCH ALL M3 CSV FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("SEARCHING M3 DIRECTORY FOR BIOLOGICAL DATA\n")
cat("=====================================================================\n")

all_csv <- list.files(
  path = m3_root,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("CSV files found:", length(all_csv), "\n")


# ======================================================================
# 9. EXCLUSION RULES
#
# We do NOT use files that directly contain:
#   - frozen ecosystem state definitions
#   - clustering results
#   - state labels
#   - survival outcomes
#   - clinical validation
#   - Layer C2 outputs
#   - state-wise summaries
#   - files whose purpose is to reproduce the state labels
#
# We retain patient-level biological tables.
# ======================================================================

exclude_pattern <- paste(
  c(
    "M3\\.2_patient_ecosystem_states",
    "M3\\.2_state_",
    "M3\\.2_top_modules",
    "M3\\.2_bootstrap",
    "M3\\.2_cluster",
    "M3\\.2_module_QC",
    "M3\\.2_module_variation",
    "M3\\.2_PCA_loadings",
    "M3\\.2_PCA_variance",
    "M3\\.2_state",
    "M3\\.3_",
    "M3\\.4_",
    "M3\\.5_",
    "M3\\.5\\.1",
    "M3\\.5\\.2",
    "clinical",
    "survival",
    "hazard",
    "cox"
  ),
  collapse = "|"
)

candidate_csv <- all_csv[
  !grepl(
    exclude_pattern,
    all_csv,
    ignore.case = TRUE
  )
]

cat(
  "Files remaining after exclusion:",
  length(candidate_csv),
  "\n\n"
)


# ======================================================================
# 10. IDENTIFY PATIENT-LEVEL FILES
# ======================================================================

cat("=====================================================================\n")
cat("IDENTIFYING PATIENT-LEVEL BIOLOGICAL FILES\n")
cat("=====================================================================\n")

patient_files <- character(0)

file_diagnostic <- data.frame(
  file = character(0),
  n_rows = integer(0),
  n_columns = integer(0),
  has_sample_id = logical(0),
  stringsAsFactors = FALSE
)

for (f in candidate_csv) {

  dat <- tryCatch(
    {
      read_csv(
        f,
        n_max = 1000,
        show_col_types = FALSE
      )
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(dat)) {
    next
  }

  has_id <- (
    "sample_id" %in% names(dat)
  )

  if (!has_id) {
    next
  }

  patient_files <- c(
    patient_files,
    f
  )

  file_diagnostic <- rbind(
    file_diagnostic,
    data.frame(
      file = f,
      n_rows = nrow(dat),
      n_columns = ncol(dat),
      has_sample_id = TRUE,
      stringsAsFactors = FALSE
    )
  )
}

cat(
  "Patient-level candidate files:",
  length(patient_files),
  "\n\n"
)

if (nrow(file_diagnostic) > 0) {
  print(file_diagnostic)
}

write_csv(
  file_diagnostic,
  file.path(
    output_dir,
    "M3.5.2_patient_level_file_diagnostic.csv"
  )
)


# ======================================================================
# 11. IF NO BIOLOGICAL FILES ARE FOUND
# ======================================================================

if (length(patient_files) == 0) {

  cat("\n")
  cat("WARNING: No independent patient-level biological files found.\n")

  empty_results <- data.frame(
    variable = character(0),
    n = integer(0),
    n_states = integer(0),
    statistic = numeric(0),
    p_value = numeric(0),
    FDR = numeric(0),
    epsilon_squared = numeric(0)
  )

  write_csv(
    empty_results,
    file.path(
      output_dir,
      "M3.5.2_biological_continuous_associations.csv"
    )
  )

  summary_table <- data.frame(
    metric = c(
      "Frozen patients",
      "Frozen ecosystem states",
      "Independent biological files",
      "Biological variables tested",
      "Survival analysis",
      "Univariate hazard ratio",
      "Multivariate hazard ratio",
      "Cox regression",
      "State re-clustering"
    ),
    value = c(
      nrow(frozen),
      nlevels(frozen$ecosystem_state),
      0,
      0,
      "NOT PERFORMED",
      "NOT PERFORMED",
      "NOT PERFORMED",
      "NOT PERFORMED",
      "NOT PERFORMED"
    ),
    stringsAsFactors = FALSE
  )

  write_csv(
    summary_table,
    file.path(
      output_dir,
      "M3.5.2_overall_summary.csv"
    )
  )

  stop(
    "\nM3.5.2 stopped safely because no independent biological files were found."
  )
}


# ======================================================================
# 12. READ PATIENT-LEVEL BIOLOGICAL TABLES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING PATIENT-LEVEL BIOLOGICAL FEATURES\n")
cat("=====================================================================\n")

biological_tables <- list()

for (f in patient_files) {

  cat("\nReading:\n", f, "\n")

  dat <- tryCatch(
    {
      read_csv(
        f,
        show_col_types = FALSE
      )
    },
    error = function(e) {

      cat(
        "Could not read file. Skipping.\n"
      )

      NULL
    }
  )

  if (is.null(dat)) {
    next
  }

  if (!"sample_id" %in% names(dat)) {
    next
  }

  dat$sample_id <- as.character(
    dat$sample_id
  )

  dat$sample_id_clean <- trimws(
    dat$sample_id
  )

  dat <- dat[
    !is.na(dat$sample_id_clean) &
    dat$sample_id_clean != "",
  ]

  if (anyDuplicated(dat$sample_id_clean) > 0) {

    dat <- dat[
      !duplicated(dat$sample_id_clean),
    ]
  }

  biological_tables[[length(biological_tables) + 1]] <- dat
}

cat(
  "\nSuccessfully read:",
  length(biological_tables),
  "patient-level biological tables.\n"
)


# ======================================================================
# 13. MERGE BIOLOGICAL TABLES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("MERGING INDEPENDENT BIOLOGICAL FEATURES\n")
cat("=====================================================================\n")

merged_bio <- frozen[
  ,
  c(
    "sample_id_clean",
    "ecosystem_state"
  )
]

for (i in seq_along(biological_tables)) {

  dat <- biological_tables[[i]]

  dat <- dat[
    ,
    !duplicated(names(dat))
  ]

  keep_columns <- c(
    "sample_id_clean"
  )

  for (nm in names(dat)) {

    if (nm == "sample_id_clean") {
      next
    }

    if (nm == "sample_id") {
      next
    }

    if (nm %in% names(merged_bio)) {
      next
    }

    keep_columns <- c(
      keep_columns,
      nm
    )
  }

  dat_small <- dat[
    ,
    keep_columns,
    drop = FALSE
  ]

  merged_bio <- merge(
    merged_bio,
    dat_small,
    by = "sample_id_clean",
    all.x = TRUE,
    sort = FALSE
  )
}

cat(
  "Merged biological dataset:\n"
)

cat(
  "Rows:",
  nrow(merged_bio),
  "\n"
)

cat(
  "Columns:",
  ncol(merged_bio),
  "\n\n"
)

cat(
  "Frozen patients:",
  nrow(frozen),
  "\n"
)

cat(
  "Patients after biological merge:",
  nrow(merged_bio),
  "\n"
)

cat(
  "Patients lost:",
  nrow(frozen) -
    nrow(merged_bio),
  "\n"
)


# ======================================================================
# 14. SAVE MERGED DATASET
# ======================================================================

write_csv(
  merged_bio,
  file.path(
    output_dir,
    "M3.5.2_orthogonal_biological_validation_dataset.csv"
  )
)


# ======================================================================
# 15. IDENTIFY BIOLOGICAL VARIABLES
# ======================================================================

identifier_columns <- c(
  "sample_id",
  "sample_id_clean",
  "ecosystem_state",
  "cluster_numeric",
  "actual_state",
  "predicted_state",
  "state",
  "cluster",
  "group"
)

candidate_variables <- setdiff(
  names(merged_bio),
  identifier_columns
)

cat("\n")
cat("=====================================================================\n")
cat("IDENTIFYING BIOLOGICAL VARIABLES\n")
cat("=====================================================================\n")

cat(
  "Candidate biological variables:",
  length(candidate_variables),
  "\n"
)


# ======================================================================
# 16. CLASSIFY VARIABLES
# ======================================================================

continuous_variables <- character(0)
categorical_variables <- character(0)

for (v in candidate_variables) {

  x <- merged_bio[[v]]

  if (is.numeric(x) || is.integer(x)) {

    continuous_variables <- c(
      continuous_variables,
      v
    )

  } else {

    categorical_variables <- c(
      categorical_variables,
      v
    )
  }
}

cat(
  "Continuous variables:",
  length(continuous_variables),
  "\n"
)

cat(
  "Categorical variables:",
  length(categorical_variables),
  "\n"
)


# ======================================================================
# 17. INFORMATION FILTERING
# ======================================================================

filtered_continuous <- character(0)

for (v in continuous_variables) {

  x <- merged_bio[[v]]

  x <- x[
    is.finite(x)
  ]

  if (length(x) < 10) {
    next
  }

  if (length(unique(x)) < 3) {
    next
  }

  if (sd(x, na.rm = TRUE) == 0) {
    next
  }

  filtered_continuous <- c(
    filtered_continuous,
    v
  )
}


filtered_categorical <- character(0)

for (v in categorical_variables) {

  x <- as.character(
    merged_bio[[v]]
  )

  x <- x[
    !is.na(x) &
    x != ""
  ]

  if (length(x) < 10) {
    next
  }

  if (length(unique(x)) < 2) {
    next
  }

  if (length(unique(x)) >= nrow(merged_bio)) {
    next
  }

  filtered_categorical <- c(
    filtered_categorical,
    v
  )
}

continuous_variables <- filtered_continuous
categorical_variables <- filtered_categorical

cat("\n")
cat("After information filtering:\n")

cat(
  "Continuous:",
  length(continuous_variables),
  "\n"
)

cat(
  "Categorical:",
  length(categorical_variables),
  "\n"
)


# ======================================================================
# 18. KRUSKAL-WALLIS TEST
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C2-O — ORTHOGONAL BIOLOGICAL CONTINUOUS ASSOCIATIONS\n")
cat("=====================================================================\n")

continuous_results <- data.frame(
  variable = character(0),
  n = integer(0),
  n_states = integer(0),
  statistic = numeric(0),
  p_value = numeric(0),
  epsilon_squared = numeric(0),
  stringsAsFactors = FALSE
)

for (v in continuous_variables) {

  x <- merged_bio[[v]]
  g <- merged_bio$ecosystem_state

  keep <- (
    !is.na(x) &
    is.finite(x) &
    !is.na(g)
  )

  x2 <- x[keep]
  g2 <- g[keep]

  if (length(x2) < 10) {
    next
  }

  if (length(unique(g2)) < 2) {
    next
  }

  test <- tryCatch(
    {
      kruskal.test(
        x2 ~ g2
      )
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(test)) {
    next
  }

  H <- as.numeric(
    test$statistic
  )

  N <- length(x2)

  k <- length(
    unique(g2)
  )

  epsilon <- (
    H - k + 1
  ) / (
    N - k
  )

  if (!is.finite(epsilon)) {
    epsilon <- NA_real_
  }

  if (epsilon < 0) {
    epsilon <- 0
  }

  continuous_results <- rbind(
    continuous_results,
    data.frame(
      variable = v,
      n = N,
      n_states = k,
      statistic = H,
      p_value = test$p.value,
      epsilon_squared = epsilon,
      stringsAsFactors = FALSE
    )
  )
}

if (nrow(continuous_results) > 0) {

  continuous_results$FDR <- p.adjust(
    continuous_results$p_value,
    method = "BH"
  )

  continuous_results <- continuous_results[
    order(
      continuous_results$FDR,
      -continuous_results$epsilon_squared
    ),
  ]

} else {

  continuous_results$FDR <- numeric(0)
}

cat(
  "Continuous tests completed:",
  nrow(continuous_results),
  "\n"
)

if (nrow(continuous_results) > 0) {

  n_sig <- sum(
    continuous_results$FDR < 0.05,
    na.rm = TRUE
  )

} else {

  n_sig <- 0
}

cat(
  "FDR < 0.05:",
  n_sig,
  "\n\n"
)

if (nrow(continuous_results) > 0) {

  print(
    head(
      continuous_results,
      30
    )
  )
}


write_csv(
  continuous_results,
  file.path(
    output_dir,
    "M3.5.2_biological_continuous_associations.csv"
  )
)


# ======================================================================
# 19. PAIRWISE WILCOXON TESTS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("PAIRWISE ORTHOGONAL BIOLOGICAL COMPARISONS\n")
cat("=====================================================================\n")

pairwise_results <- data.frame(
  variable = character(0),
  state_1 = character(0),
  state_2 = character(0),
  n_1 = integer(0),
  n_2 = integer(0),
  median_1 = numeric(0),
  median_2 = numeric(0),
  p_value = numeric(0),
  stringsAsFactors = FALSE
)

states <- levels(
  merged_bio$ecosystem_state
)

if (length(states) >= 2 &&
    length(continuous_variables) > 0) {

  state_pairs <- combn(
    states,
    2,
    simplify = FALSE
  )

  for (v in continuous_variables) {

    x <- merged_bio[[v]]

    for (pair in state_pairs) {

      s1 <- pair[[1]]
      s2 <- pair[[2]]

      x1 <- x[
        merged_bio$ecosystem_state == s1
      ]

      x2 <- x[
        merged_bio$ecosystem_state == s2
      ]

      x1 <- x1[
        is.finite(x1)
      ]

      x2 <- x2[
        is.finite(x2)
      ]

      if (length(x1) < 3 ||
          length(x2) < 3) {
        next
      }

      test <- tryCatch(
        {
          wilcox.test(
            x1,
            x2,
            exact = FALSE
          )
        },
        error = function(e) {
          NULL
        }
      )

      if (is.null(test)) {
        next
      }

      pairwise_results <- rbind(
        pairwise_results,
        data.frame(
          variable = v,
          state_1 = s1,
          state_2 = s2,
          n_1 = length(x1),
          n_2 = length(x2),
          median_1 = median(
            x1,
            na.rm = TRUE
          ),
          median_2 = median(
            x2,
            na.rm = TRUE
          ),
          p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

if (nrow(pairwise_results) > 0) {

  pairwise_results$FDR <- p.adjust(
    pairwise_results$p_value,
    method = "BH"
  )
}

cat(
  "Pairwise tests completed:",
  nrow(pairwise_results),
  "\n"
)

write_csv(
  pairwise_results,
  file.path(
    output_dir,
    "M3.5.2_pairwise_biological_comparisons.csv"
  )
)


# ======================================================================
# 20. CATEGORICAL BIOLOGICAL ASSOCIATIONS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C2-O — ORTHOGONAL BIOLOGICAL CATEGORICAL ASSOCIATIONS\n")
cat("=====================================================================\n")

categorical_results <- data.frame(
  variable = character(0),
  n = integer(0),
  n_levels = integer(0),
  p_value = numeric(0),
  stringsAsFactors = FALSE
)

for (v in categorical_variables) {

  x <- as.character(
    merged_bio[[v]]
  )

  g <- as.character(
    merged_bio$ecosystem_state
  )

  keep <- (
    !is.na(x) &
    x != "" &
    !is.na(g) &
    g != ""
  )

  x2 <- x[keep]
  g2 <- g[keep]

  if (length(x2) < 10) {
    next
  }

  if (length(unique(x2)) < 2) {
    next
  }

  tab <- table(
    g2,
    x2
  )

  if (nrow(tab) < 2 ||
      ncol(tab) < 2) {
    next
  }

  test <- tryCatch(
    {
      fisher.test(
        tab,
        simulate.p.value = TRUE,
        B = 5000
      )
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(test)) {
    next
  }

  categorical_results <- rbind(
    categorical_results,
    data.frame(
      variable = v,
      n = length(x2),
      n_levels = length(unique(x2)),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  )
}

if (nrow(categorical_results) > 0) {

  categorical_results$FDR <- p.adjust(
    categorical_results$p_value,
    method = "BH"
  )

  categorical_results <- categorical_results[
    order(
      categorical_results$FDR
    ),
  ]

} else {

  categorical_results$FDR <- numeric(0)
}

cat(
  "Categorical tests completed:",
  nrow(categorical_results),
  "\n"
)

if (nrow(categorical_results) > 0) {

  cat(
    "FDR < 0.05:",
    sum(
      categorical_results$FDR < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
}

write_csv(
  categorical_results,
  file.path(
    output_dir,
    "M3.5.2_biological_categorical_associations.csv"
  )
)


# ======================================================================
# 21. TOP BIOLOGICAL FEATURES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("GENERATING TOP ORTHOGONAL BIOLOGICAL FEATURES\n")
cat("=====================================================================\n")

if (nrow(continuous_results) > 0) {

  top_n <- min(
    30,
    nrow(continuous_results)
  )

  top_features <- continuous_results[
    seq_len(top_n),
    ,
    drop = FALSE
  ]

} else {

  top_features <- data.frame(
    variable = character(0),
    n = integer(0),
    n_states = integer(0),
    statistic = numeric(0),
    p_value = numeric(0),
    epsilon_squared = numeric(0),
    FDR = numeric(0),
    stringsAsFactors = FALSE
  )
}

write_csv(
  top_features,
  file.path(
    output_dir,
    "M3.5.2_top_orthogonal_biological_features.csv"
  )
)

cat(
  "Top features:",
  nrow(top_features),
  "\n"
)


# ======================================================================
# 22. STATE-WISE CONTINUOUS PROFILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("STATE-WISE ORTHOGONAL BIOLOGICAL PROFILES\n")
cat("=====================================================================\n")

profile_results <- data.frame(
  variable = character(0),
  ecosystem_state = character(0),
  n = integer(0),
  mean = numeric(0),
  median = numeric(0),
  sd = numeric(0),
  stringsAsFactors = FALSE
)

if (length(continuous_variables) > 0) {

  for (v in continuous_variables) {

    x <- merged_bio[[v]]

    for (s in states) {

      values <- x[
        merged_bio$ecosystem_state == s
      ]

      values <- values[
        is.finite(values)
      ]

      if (length(values) == 0) {
        next
      }

      profile_results <- rbind(
        profile_results,
        data.frame(
          variable = v,
          ecosystem_state = s,
          n = length(values),
          mean = mean(
            values,
            na.rm = TRUE
          ),
          median = median(
            values,
            na.rm = TRUE
          ),
          sd = ifelse(
            length(values) > 1,
            sd(
              values,
              na.rm = TRUE
            ),
            NA_real_
          ),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

cat(
  "Profiles generated:",
  nrow(profile_results),
  "\n"
)

write_csv(
  profile_results,
  file.path(
    output_dir,
    "M3.5.2_state_wise_biological_profiles.csv"
  )
)


# ======================================================================
# 23. TOP FEATURE HEATMAP
#
# Heatmap is deliberately based only on the top independent biological
# variables. It is NOT used to redefine states.
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("GENERATING ORTHOGONAL BIOLOGICAL HEATMAP\n")
cat("=====================================================================\n")

heatmap_file <- file.path(
  output_dir,
  "M3.5.2_orthogonal_biological_heatmap.png"
)

if (nrow(top_features) >= 2) {

  heatmap_vars <- top_features$variable

  heatmap_vars <- heatmap_vars[
    heatmap_vars %in% names(merged_bio)
  ]

  if (length(heatmap_vars) >= 2) {

    heat_data <- merged_bio[
      ,
      c(
        "sample_id_clean",
        "ecosystem_state",
        heatmap_vars
      ),
      drop = FALSE
    ]

    matrix_data <- as.matrix(
      heat_data[
        ,
        heatmap_vars,
        drop = FALSE
      ]
    )

    suppressWarnings(
      storage.mode(matrix_data) <- "numeric"
    )

    keep_rows <- apply(
      matrix_data,
      1,
      function(z) {
        sum(is.finite(z)) >= 2
      }
    )

    matrix_data <- matrix_data[
      keep_rows,
      ,
      drop = FALSE
    ]

    if (nrow(matrix_data) >= 3 &&
        ncol(matrix_data) >= 2) {

      scaled_matrix <- t(
        scale(
          matrix_data
        )
      )

      scaled_matrix[
        !is.finite(scaled_matrix)
      ] <- 0

      heat_df <- as.data.frame(
        scaled_matrix
      )

      heat_df$feature <- rownames(
        heat_df
      )

      heat_long <- pivot_longer(
        heat_df,
        cols = -feature,
        names_to = "sample",
        values_to = "z"
      )

      heat_plot <- ggplot(
        heat_long,
        aes(
          x = sample,
          y = feature,
          fill = z
        )
      ) +
        geom_tile() +
        theme_minimal() +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(
            size = 7
          )
        ) +
        labs(
          title = "Top Orthogonal Biological Features",
          x = "Patients",
          y = "Biological feature",
          fill = "Z-score"
        )

      ggsave(
        heatmap_file,
        heat_plot,
        width = 12,
        height = 8,
        dpi = 300
      )

      cat(
        "Heatmap generated successfully.\n"
      )

    } else {

      cat(
        "Not enough valid data for heatmap.\n"
      )
    }

  } else {

    cat(
      "Not enough top biological variables for heatmap.\n"
    )
  }

} else {

  cat(
    "No significant biological features available for heatmap.\n"
  )
}


# ======================================================================
# 24. PCA OF ORTHOGONAL BIOLOGICAL FEATURES
#
# PCA is descriptive only.
# It does NOT create new ecosystem states.
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("ORTHOGONAL BIOLOGICAL PCA\n")
cat("=====================================================================\n")

pca_file <- file.path(
  output_dir,
  "M3.5.2_orthogonal_biological_PCA.png"
)

pca_scores_file <- file.path(
  output_dir,
  "M3.5.2_orthogonal_biological_PCA_scores.csv"
)

if (length(heatmap_vars) >= 2) {

  pca_matrix <- as.data.frame(
    merged_bio[
      ,
      heatmap_vars,
      drop = FALSE
    ]
  )

  for (j in seq_along(pca_matrix)) {

    pca_matrix[[j]] <- suppressWarnings(
      as.numeric(
        pca_matrix[[j]]
      )
    )
  }

  complete_rows <- complete.cases(
    pca_matrix
  )

  pca_matrix_complete <- pca_matrix[
    complete_rows,
    ,
    drop = FALSE
  ]

  pca_states <- merged_bio$ecosystem_state[
    complete_rows
  ]

  pca_ids <- merged_bio$sample_id_clean[
    complete_rows
  ]

  if (
    nrow(pca_matrix_complete) >= 5 &&
    ncol(pca_matrix_complete) >= 2
  ) {

    variable_sd <- apply(
      pca_matrix_complete,
      2,
      sd,
      na.rm = TRUE
    )

    variable_keep <- (
      is.finite(variable_sd) &
      variable_sd > 0
    )

    pca_matrix_complete <- pca_matrix_complete[
      ,
      variable_keep,
      drop = FALSE
    ]

    if (ncol(pca_matrix_complete) >= 2) {

      pca_fit <- tryCatch(
        {
          prcomp(
            pca_matrix_complete,
            center = TRUE,
            scale. = TRUE
          )
        },
        error = function(e) {
          NULL
        }
      )

      if (!is.null(pca_fit)) {

        pca_scores <- as.data.frame(
          pca_fit$x
        )

        pca_scores$sample_id <- pca_ids
        pca_scores$ecosystem_state <- pca_states

        write_csv(
          pca_scores,
          pca_scores_file
        )

        variance <- (
          pca_fit$sdev^2
        )

        variance <- variance /
          sum(variance)

        pc1_var <- round(
          variance[1] * 100,
          2
        )

        pc2_var <- round(
          variance[2] * 100,
          2
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
            size = 3
          ) +
          theme_minimal() +
          labs(
            title = "PCA of Orthogonal Biological Features",
            x = paste0(
              "PC1 (",
              pc1_var,
              "%)"
            ),
            y = paste0(
              "PC2 (",
              pc2_var,
              "%)"
            ),
            shape = "Frozen ecosystem state"
          )

        ggsave(
          pca_file,
          pca_plot,
          width = 8,
          height = 6,
          dpi = 300
        )

        cat(
          "PCA generated successfully.\n"
        )

      } else {

        cat(
          "PCA failed safely.\n"
        )
      }

    } else {

      cat(
        "Fewer than 2 variable dimensions remain for PCA.\n"
      )
    }

  } else {

    cat(
      "Insufficient complete observations for PCA.\n"
    )
  }

} else {

  cat(
    "Insufficient variables for PCA.\n"
  )
}


# ======================================================================
# 25. CORRELATION WITH FROZEN STATE NUMERIC CODE
#
# Descriptive only.
# NOT used to redefine states.
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("STATE-SEPARATION EFFECT SUMMARY\n")
cat("=====================================================================\n")

if (nrow(continuous_results) > 0) {

  effect_summary <- continuous_results[
    ,
    c(
      "variable",
      "statistic",
      "p_value",
      "FDR",
      "epsilon_squared"
    ),
    drop = FALSE
  ]

  effect_summary <- effect_summary[
    order(
      effect_summary$FDR,
      -effect_summary$epsilon_squared
    ),
  ]

} else {

  effect_summary <- data.frame(
    variable = character(0),
    statistic = numeric(0),
    p_value = numeric(0),
    FDR = numeric(0),
    epsilon_squared = numeric(0),
    stringsAsFactors = FALSE
  )
}

write_csv(
  effect_summary,
  file.path(
    output_dir,
    "M3.5.2_biological_effect_size_summary.csv"
  )
)


# ======================================================================
# 26. OVERALL SUMMARY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.5.2 OVERALL SUMMARY\n")
cat("=====================================================================\n")

n_continuous_sig <- 0

if (nrow(continuous_results) > 0) {

  n_continuous_sig <- sum(
    continuous_results$FDR < 0.05,
    na.rm = TRUE
  )
}

n_categorical_sig <- 0

if (nrow(categorical_results) > 0) {

  n_categorical_sig <- sum(
    categorical_results$FDR < 0.05,
    na.rm = TRUE
  )
}

n_pairwise_sig <- 0

if (nrow(pairwise_results) > 0 &&
    "FDR" %in% names(pairwise_results)) {

  n_pairwise_sig <- sum(
    pairwise_results$FDR < 0.05,
    na.rm = TRUE
  )
}

summary_table <- data.frame(
  metric = c(
    "Frozen patients",
    "Frozen ecosystem states",
    "Independent biological files",
    "Candidate biological variables",
    "Continuous variables tested",
    "Categorical variables tested",
    "Continuous FDR < 0.05",
    "Categorical FDR < 0.05",
    "Pairwise biological comparisons",
    "Pairwise FDR < 0.05",
    "State re-clustering",
    "Ecosystem states modified",
    "Survival analysis",
    "Univariate hazard ratio",
    "Multivariate hazard ratio",
    "Cox regression"
  ),
  value = c(
    nrow(frozen),
    nlevels(frozen$ecosystem_state),
    length(biological_tables),
    length(candidate_variables),
    nrow(continuous_results),
    nrow(categorical_results),
    n_continuous_sig,
    n_categorical_sig,
    nrow(pairwise_results),
    n_pairwise_sig,
    "NOT PERFORMED",
    "NO",
    "NOT PERFORMED",
    "NOT PERFORMED",
    "NOT PERFORMED",
    "NOT PERFORMED"
  ),
  stringsAsFactors = FALSE
)

print(summary_table)

write_csv(
  summary_table,
  file.path(
    output_dir,
    "M3.5.2_overall_summary.csv"
  )
)


# ======================================================================
# 27. ANALYSIS PRINCIPLES FILE
# ======================================================================

principles <- data.frame(
  analysis_component = c(
    "Frozen M3.2 ecosystem states",
    "State re-clustering",
    "State redefinition",
    "Survival analysis",
    "Univariate hazard ratios",
    "Multivariate hazard ratios",
    "Cox regression",
    "Independent biological characterization",
    "Kruskal-Wallis testing",
    "Pairwise Wilcoxon testing",
    "Categorical association testing",
    "Biological PCA",
    "Biological heatmap"
  ),
  status = c(
    "USED AS FIXED INPUT",
    "NOT PERFORMED",
    "NOT PERFORMED",
    "NOT PERFORMED",
    "NOT PERFORMED",
    "NOT PERFORMED",
    "NOT PERFORMED",
    "PERFORMED",
    "PERFORMED",
    "PERFORMED",
    "PERFORMED",
    "DESCRIPTIVE ONLY",
    "DESCRIPTIVE ONLY"
  ),
  stringsAsFactors = FALSE
)

write_csv(
  principles,
  file.path(
    output_dir,
    "M3.5.2_analysis_principles.csv"
  )
)


# ======================================================================
# 28. DIAGNOSTIC LOG
# ======================================================================

diagnostic_file <- file.path(
  output_dir,
  "M3.5.2_run_diagnostic.txt"
)

sink(
  diagnostic_file
)

cat("M3.5.2 ORTHOGONAL BIOLOGICAL VALIDATION\n")
cat("========================================\n\n")

cat(
  "Run time:",
  as.character(Sys.time()),
  "\n\n"
)

cat(
  "Frozen state file:\n",
  frozen_file,
  "\n\n"
)

cat(
  "Output directory:\n",
  output_dir,
  "\n\n"
)

cat(
  "Frozen patients:",
  nrow(frozen),
  "\n"
)

cat(
  "Frozen states:",
  nlevels(frozen$ecosystem_state),
  "\n"
)

cat(
  "Independent biological files:",
  length(biological_tables),
  "\n"
)

cat(
  "Candidate variables:",
  length(candidate_variables),
  "\n"
)

cat(
  "Continuous variables tested:",
  nrow(continuous_results),
  "\n"
)

cat(
  "Categorical variables tested:",
  nrow(categorical_results),
  "\n"
)

cat(
  "Pairwise tests:",
  nrow(pairwise_results),
  "\n\n"
)

cat(
  "Continuous FDR < 0.05:",
  n_continuous_sig,
  "\n"
)

cat(
  "Categorical FDR < 0.05:",
  n_categorical_sig,
  "\n"
)

cat(
  "Pairwise FDR < 0.05:",
  n_pairwise_sig,
  "\n\n"
)

cat("Survival analysis: NOT PERFORMED\n")
cat("Univariate hazard ratios: NOT PERFORMED\n")
cat("Multivariate hazard ratios: NOT PERFORMED\n")
cat("Cox regression: NOT PERFORMED\n")
cat("State re-clustering: NOT PERFORMED\n")
cat("State redefinition: NOT PERFORMED\n")
cat("M3.2 frozen states modified: NO\n")

sink()


# ======================================================================
# 29. FINAL MESSAGE
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.5.2 COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n")

cat(
  "Frozen patients:",
  nrow(frozen),
  "\n"
)

cat(
  "Frozen ecosystem states:",
  nlevels(frozen$ecosystem_state),
  "\n"
)

cat(
  "Independent biological files:",
  length(biological_tables),
  "\n"
)

cat(
  "Continuous biological variables tested:",
  nrow(continuous_results),
  "\n"
)

cat(
  "Continuous FDR < 0.05:",
  n_continuous_sig,
  "\n"
)

cat(
  "Categorical biological variables tested:",
  nrow(categorical_results),
  "\n"
)

cat(
  "Categorical FDR < 0.05:",
  n_categorical_sig,
  "\n"
)

cat(
  "Pairwise biological tests:",
  nrow(pairwise_results),
  "\n"
)

cat("\n")
cat("Survival analysis: NOT PERFORMED\n")
cat("Univariate hazard ratios: NOT PERFORMED\n")
cat("Multivariate hazard ratios: NOT PERFORMED\n")
cat("Cox regression: NOT PERFORMED\n")
cat("State re-clustering: NOT PERFORMED\n")
cat("M3.2 states were NOT modified.\n")

cat("\n")
cat("Output directory:\n")
cat(output_dir, "\n")

cat("\n")
cat("=====================================================================\n")
cat("M3.5.2 FINISHED\n")
cat("=====================================================================\n")