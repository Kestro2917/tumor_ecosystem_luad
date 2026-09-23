# ======================================================================
# M3.5 — FROZEN ECOSYSTEM STATE CLINICAL / PHENOTYPIC VALIDATION
# ----------------------------------------------------------------------
# LAYER A
# Baseline clinical characterization
#   - Does frozen ecosystem state associate with baseline phenotype?
#   - Continuous variables: Kruskal-Wallis
#   - Categorical variables: Fisher's exact / simulated Fisher
#   - BH-FDR correction
#   - Exploratory pairwise tests only for FDR-significant variables
#
# LAYER B
# Survival / outcome validation
#   - Uses frozen ecosystem state without modification
#   - Overall survival Kaplan-Meier
#   - Global log-rank test
#
# LAYER C
# Biological / ecosystem characterization
#   - NEVER redefines ecosystem states
#   - Uses molecular/ecosystem variables only when actually present
#   - Safely skips Layer C molecular testing if M3.2 state file contains
#     only identifiers + ecosystem state + cluster number
#
# IMPORTANT
# ----------------------------------------------------------------------
# M3.2 ecosystem states are FROZEN.
# This script does NOT recluster patients.
# This script does NOT redefine State_1 / State_2.
# This script does NOT use survival variables to define states.
# ======================================================================


# ======================================================================
# 0. INITIAL SETTINGS
# ======================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

start_time <- Sys.time()

cat("\n")
cat("=====================================================================\n")
cat("M3.5 — FROZEN ECOSYSTEM STATE CLINICAL / PHENOTYPIC VALIDATION\n")
cat("=====================================================================\n")
cat("Started:", as.character(start_time), "\n")


# ======================================================================
# 1. FILE PATHS
# ======================================================================

clinical_file <- "/content/M1_clinical.csv"

m32_state_file <- paste0(
  "/content/M3_LIONESS_entropy/",
  "M3.2_module_axes/",
  "M3.2_patient_ecosystem_states.csv"
)

output_dir <- paste0(
  "/content/M3_LIONESS_entropy/",
  "M3.5_clinical_validation"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\n")
cat("Clinical file:\n")
cat(clinical_file, "\n")

cat("\nM3.2 frozen state file:\n")
cat(m32_state_file, "\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")


# ======================================================================
# 2. INSTALL / LOAD REQUIRED PACKAGES
# ======================================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "survival",
  "survminer"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {

  if (!(pkg %in% installed)) {

    cat("Installing package:", pkg, "\n")

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org",
      dependencies = TRUE
    )
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(survminer)
})


# ======================================================================
# 3. HELPER FUNCTIONS
# ======================================================================


# ----------------------------------------------------------------------
# Clean column names
# ----------------------------------------------------------------------

clean_names_simple <- function(x) {

  x <- as.character(x)

  x <- gsub("[[:space:]]+", "_", x)

  x <- gsub("[^A-Za-z0-9_.-]", "_", x)

  x <- make.unique(x, sep = "_")

  return(x)
}


# ----------------------------------------------------------------------
# Find first existing column from candidate names
# ----------------------------------------------------------------------

find_column <- function(data, candidates) {

  available <- names(data)

  hit <- candidates[candidates %in% available]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  return(hit[1])
}


# ----------------------------------------------------------------------
# Normalize TCGA identifiers to PATIENT LEVEL
#
# Example:
#
# TCGA-05-4402
# TCGA-05-4402-01A
# TCGA-05-4402-01A-01R-1206-07
#
# all become:
#
# TCGA-05-4402
# ----------------------------------------------------------------------

normalize_tcga_patient_id <- function(x) {

  x <- as.character(x)

  x <- trimws(x)

  x <- toupper(x)

  x <- gsub("\\.", "-", x)

  x <- sub(
    "^((TCGA-[0-9A-Z]+-[0-9A-Z]+)).*$",
    "\\1",
    x
  )

  x[x %in% c("", "NA", "N/A", "NULL", "NAN")] <- NA_character_

  return(x)
}


# ----------------------------------------------------------------------
# Safe numeric conversion
# ----------------------------------------------------------------------

safe_numeric <- function(x) {

  if (is.numeric(x)) {
    return(x)
  }

  x <- as.character(x)

  x <- trimws(x)

  x[x %in% c("", "NA", "N/A", "NULL", "NaN")] <- NA_character_

  suppressWarnings(as.numeric(x))
}


# ----------------------------------------------------------------------
# Safe Fisher test
#
# Fisher's exact test can become unstable with larger tables.
# Simulated Fisher is therefore used for tables larger than 2x2.
# ----------------------------------------------------------------------

safe_fisher_test <- function(tab) {

  if (length(tab) == 0) {
    return(NA_real_)
  }

  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(NA_real_)
  }

  result <- tryCatch({

    if (nrow(tab) == 2 && ncol(tab) == 2) {

      fisher.test(tab)$p.value

    } else {

      fisher.test(
        tab,
        simulate.p.value = TRUE,
        B = 5000
      )$p.value
    }

  }, error = function(e) {

    NA_real_
  })

  return(result)
}


# ----------------------------------------------------------------------
# Safe Kruskal-Wallis
# ----------------------------------------------------------------------

safe_kruskal <- function(value, state) {

  keep <- is.finite(value) & !is.na(state)

  value <- value[keep]

  state <- droplevels(
    factor(state[keep])
  )

  if (length(value) < 3) {
    return(
      list(
        statistic = NA_real_,
        p_value = NA_real_,
        n = length(value),
        n_states = nlevels(state)
      )
    )
  }

  if (nlevels(state) < 2) {
    return(
      list(
        statistic = NA_real_,
        p_value = NA_real_,
        n = length(value),
        n_states = nlevels(state)
      )
    )
  }

  if (length(unique(value)) < 2) {
    return(
      list(
        statistic = NA_real_,
        p_value = NA_real_,
        n = length(value),
        n_states = nlevels(state)
      )
    )
  }

  result <- tryCatch({

    kt <- kruskal.test(
      value ~ state
    )

    list(
      statistic = unname(kt$statistic),
      p_value = kt$p.value,
      n = length(value),
      n_states = nlevels(state)
    )

  }, error = function(e) {

    list(
      statistic = NA_real_,
      p_value = NA_real_,
      n = length(value),
      n_states = nlevels(state)
    )
  })

  return(result)
}


# ----------------------------------------------------------------------
# Safe pairwise Wilcoxon
# ----------------------------------------------------------------------

safe_pairwise_wilcox <- function(data, variable_name) {

  result <- NULL

  if (!(variable_name %in% names(data))) {
    return(result)
  }

  temp <- data.frame(
    value = safe_numeric(data[[variable_name]]),
    ecosystem_state = as.character(data$ecosystem_state),
    stringsAsFactors = FALSE
  )

  temp <- temp[
    is.finite(temp$value) &
      !is.na(temp$ecosystem_state),
    ,
    drop = FALSE
  ]

  if (nrow(temp) < 3) {
    return(result)
  }

  temp$ecosystem_state <- factor(
    temp$ecosystem_state
  )

  states <- levels(temp$ecosystem_state)

  if (length(states) < 2) {
    return(result)
  }

  combinations <- combn(
    states,
    2,
    simplify = FALSE
  )

  output <- vector(
    "list",
    length(combinations)
  )

  counter <- 0

  for (pair in combinations) {

    x <- temp$value[
      temp$ecosystem_state == pair[1]
    ]

    y <- temp$value[
      temp$ecosystem_state == pair[2]
    ]

    if (length(x) < 2 || length(y) < 2) {
      next
    }

    test <- tryCatch({

      wt <- wilcox.test(
        x,
        y,
        exact = FALSE
      )

      wt$p.value

    }, error = function(e) {

      NA_real_
    })

    counter <- counter + 1

    output[[counter]] <- data.frame(
      variable = variable_name,
      state_1 = pair[1],
      state_2 = pair[2],
      n_state_1 = length(x),
      n_state_2 = length(y),
      median_state_1 = median(x, na.rm = TRUE),
      median_state_2 = median(y, na.rm = TRUE),
      p_value = test,
      stringsAsFactors = FALSE
    )
  }

  if (counter == 0) {
    return(NULL)
  }

  output <- output[
    seq_len(counter)
  ]

  bind_rows(output)
}


# ----------------------------------------------------------------------
# Safe median
# ----------------------------------------------------------------------

safe_median <- function(x) {

  x <- safe_numeric(x)

  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  median(x)
}


# ----------------------------------------------------------------------
# Safe minimum
# ----------------------------------------------------------------------

safe_min <- function(x) {

  x <- safe_numeric(x)

  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  min(x)
}


# ----------------------------------------------------------------------
# Safe maximum
# ----------------------------------------------------------------------

safe_max <- function(x) {

  x <- safe_numeric(x)

  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  max(x)
}


# ======================================================================
# 4. CHECK INPUT FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("CHECKING INPUT FILES\n")
cat("=====================================================================\n")

if (!file.exists(clinical_file)) {

  stop(
    paste0(
      "Clinical file does not exist: ",
      clinical_file
    )
  )
}

if (!file.exists(m32_state_file)) {

  stop(
    paste0(
      "M3.2 frozen state file does not exist: ",
      m32_state_file
    )
  )
}

cat("Clinical file found.\n")
cat("M3.2 frozen state file found.\n")


# ======================================================================
# 5. READ CLINICAL DATA
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING CLINICAL DATA\n")
cat("=====================================================================\n")

clinical <- read_csv(
  clinical_file,
  show_col_types = FALSE,
  progress = FALSE
)

names(clinical) <- clean_names_simple(
  names(clinical)
)

cat("Clinical dimensions:\n")
cat("Rows:", nrow(clinical), "\n")
cat("Columns:", ncol(clinical), "\n")


# ======================================================================
# 6. READ M3.2 FROZEN ECOSYSTEM STATES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING M3.2 FROZEN ECOSYSTEM STATES\n")
cat("=====================================================================\n")

m32 <- read_csv(
  m32_state_file,
  show_col_types = FALSE,
  progress = FALSE
)

names(m32) <- clean_names_simple(
  names(m32)
)

cat("M3.2 dimensions:\n")
cat("Rows:", nrow(m32), "\n")
cat("Columns:", ncol(m32), "\n")

cat("\nM3.2 columns:\n")
print(names(m32))


# ======================================================================
# 7. IDENTIFY IMPORTANT COLUMNS
# ======================================================================

clinical_id_col <- find_column(
  clinical,
  c(
    "sample_barcode",
    "barcode",
    "sample_id",
    "patient"
  )
)

m32_id_col <- find_column(
  m32,
  c(
    "sample_id",
    "sample_barcode",
    "barcode",
    "patient"
  )
)

state_col <- find_column(
  m32,
  c(
    "ecosystem_state",
    "state",
    "cluster",
    "cluster_state"
  )
)

if (is.na(clinical_id_col)) {

  stop(
    "Could not identify the clinical patient/sample ID column."
  )
}

if (is.na(m32_id_col)) {

  stop(
    "Could not identify the M3.2 patient/sample ID column."
  )
}

if (is.na(state_col)) {

  stop(
    "Could not identify the M3.2 ecosystem state column."
  )
}

cat("\nClinical ID column:", clinical_id_col, "\n")
cat("M3.2 ID column:", m32_id_col, "\n")
cat("M3.2 state column:", state_col, "\n")


# ======================================================================
# 8. CREATE PATIENT-LEVEL NORMALIZED IDS
# ======================================================================

clinical$M3_patient_id <- normalize_tcga_patient_id(
  clinical[[clinical_id_col]]
)

m32$M3_patient_id <- normalize_tcga_patient_id(
  m32[[m32_id_col]]
)

m32$ecosystem_state <- as.character(
  m32[[state_col]]
)

m32$ecosystem_state <- trimws(
  m32$ecosystem_state
)

# Remove accidental empty state values
m32$ecosystem_state[
  m32$ecosystem_state %in%
    c("", "NA", "N/A", "NULL")
] <- NA_character_


# ======================================================================
# 9. CHECK M3.2 PATIENT UNIQUENESS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("CHECKING M3.2 PATIENT UNIQUENESS\n")
cat("=====================================================================\n")

duplicate_m32 <- m32$M3_patient_id[
  duplicated(m32$M3_patient_id) &
    !is.na(m32$M3_patient_id)
]

if (length(duplicate_m32) > 0) {

  cat(
    "WARNING:",
    length(unique(duplicate_m32)),
    "duplicate patient IDs detected.\n"
  )

  cat(
    "Keeping the first occurrence of each patient.\n"
  )

  m32 <- m32 %>%
    filter(!is.na(M3_patient_id)) %>%
    distinct(
      M3_patient_id,
      .keep_all = TRUE
    )

} else {

  cat("No duplicate M3.2 patient keys detected.\n")
}

cat(
  "Unique frozen M3.2 patients:",
  n_distinct(m32$M3_patient_id, na.rm = TRUE),
  "\n"
)


# ======================================================================
# 10. CHECK CLINICAL DUPLICATES
# ======================================================================

duplicate_clinical <- clinical$M3_patient_id[
  duplicated(clinical$M3_patient_id) &
    !is.na(clinical$M3_patient_id)
]

if (length(duplicate_clinical) > 0) {

  cat(
    "WARNING:",
    length(unique(duplicate_clinical)),
    "duplicate clinical patient IDs detected.\n"
  )

  cat(
    "Keeping the first clinical record for each patient.\n"
  )

  clinical <- clinical %>%
    filter(!is.na(M3_patient_id)) %>%
    distinct(
      M3_patient_id,
      .keep_all = TRUE
    )
}


# ======================================================================
# 11. PATIENT-LEVEL MATCHING
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("MERGING FROZEN STATES WITH CLINICAL DATA\n")
cat("=====================================================================\n")

frozen_patient_ids <- unique(
  m32$M3_patient_id[
    !is.na(m32$M3_patient_id)
  ]
)

clinical_patient_ids <- unique(
  clinical$M3_patient_id[
    !is.na(clinical$M3_patient_id)
  ]
)

matched_ids <- intersect(
  frozen_patient_ids,
  clinical_patient_ids
)

cat(
  "Frozen M3.2 patients:",
  length(frozen_patient_ids),
  "\n"
)

cat(
  "Clinical patients:",
  length(clinical_patient_ids),
  "\n"
)

cat(
  "Matched patients:",
  length(matched_ids),
  "\n"
)

match_percentage <- ifelse(
  length(frozen_patient_ids) > 0,
  100 * length(matched_ids) /
    length(frozen_patient_ids),
  0
)

cat(
  "Match percentage:",
  round(match_percentage, 2),
  "%\n"
)


# ----------------------------------------------------------------------
# Stop if there are no matched patients
# ----------------------------------------------------------------------

if (length(matched_ids) == 0) {

  stop(
    paste0(
      "\nNO PATIENTS MATCHED.\n",
      "The normalized patient-level TCGA IDs did not overlap.\n",
      "Check M1_clinical.csv and the M3.2 state file."
    )
  )
}


# ======================================================================
# 12. MERGE CLINICAL + FROZEN STATE DATA
# ======================================================================

frozen_for_merge <- m32 %>%
  filter(
    !is.na(M3_patient_id),
    !is.na(ecosystem_state)
  ) %>%
  select(
    M3_patient_id,
    ecosystem_state,
    everything()
  )

clinical_for_merge <- clinical %>%
  filter(
    !is.na(M3_patient_id)
  )

final_data <- clinical_for_merge %>%
  inner_join(
    frozen_for_merge,
    by = "M3_patient_id",
    suffix = c(
      "_clinical",
      "_M32"
    )
  )


# ----------------------------------------------------------------------
# Ensure ecosystem_state exists after merge
# ----------------------------------------------------------------------

if (!("ecosystem_state" %in% names(final_data))) {

  state_candidates <- grep(
    "^ecosystem_state",
    names(final_data),
    value = TRUE
  )

  if (length(state_candidates) > 0) {

    final_data$ecosystem_state <-
      final_data[[state_candidates[1]]]

  } else {

    stop(
      "ecosystem_state was lost during merge."
    )
  }
}

final_data$ecosystem_state <- as.character(
  final_data$ecosystem_state
)

final_data <- final_data %>%
  filter(
    !is.na(ecosystem_state),
    ecosystem_state != ""
  )

cat(
  "\nFinal matched dataset:",
  nrow(final_data),
  "patients\n"
)


# ======================================================================
# 13. SAVE PATIENT-ID DIAGNOSTIC
# ======================================================================

diagnostic <- data.frame(
  normalized_patient_id = matched_ids,
  stringsAsFactors = FALSE
)

diagnostic$in_frozen_M3_2 <- diagnostic$normalized_patient_id %in%
  frozen_patient_ids

diagnostic$in_clinical <- diagnostic$normalized_patient_id %in%
  clinical_patient_ids

diagnostic$matched <- diagnostic$in_frozen_M3_2 &
  diagnostic$in_clinical

write_csv(
  diagnostic,
  file.path(
    output_dir,
    "M3.5_patient_ID_matching_diagnostic.csv"
  )
)


# ======================================================================
# 14. SAVE FINAL MATCHED DATASET
# ======================================================================

write_csv(
  final_data,
  file.path(
    output_dir,
    "M3.5_FINAL_frozen_ecosystem_clinical_dataset.csv"
  )
)


# ======================================================================
# 15. FROZEN STATE DISTRIBUTION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("FROZEN ECOSYSTEM STATE DISTRIBUTION\n")
cat("=====================================================================\n")

state_distribution <- final_data %>%
  count(
    ecosystem_state,
    name = "n"
  ) %>%
  mutate(
    percentage = 100 * n / sum(n)
  )

print(state_distribution)

write_csv(
  state_distribution,
  file.path(
    output_dir,
    "M3.5_ecosystem_state_distribution.csv"
  )
)


# ======================================================================
# ======================================================================
# LAYER A — BASELINE CLINICAL CHARACTERIZATION
# ======================================================================
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER A — BASELINE CLINICAL CHARACTERIZATION\n")
cat("=====================================================================\n")

cat(
  "Purpose: test whether frozen ecosystem states are associated with\n"
)

cat(
  "baseline clinical characteristics.\n"
)

cat(
  "Survival/outcome variables are explicitly excluded.\n"
)


# ======================================================================
# 16. DEFINE EXCLUDED VARIABLES FOR LAYER A
# ======================================================================

outcome_patterns <- c(
  "survival",
  "os_",
  "^os$",
  "overall_survival",
  "days_to_death",
  "vital_status",
  "vital_status_raw",
  "death",
  "event",
  "follow_up",
  "followup",
  "last_follow"
)


# ----------------------------------------------------------------------
# Variables that should not be tested as baseline phenotype
# ----------------------------------------------------------------------

technical_patterns <- c(
  "^indexed__",
  "^indexed_",
  "^M3_patient_id$",
  "^sample_barcode$",
  "^barcode$",
  "^sample_id$",
  "^patient$",
  "^ecosystem_state$",
  "^cluster_numeric$",
  "^cluster$",
  "updated_datetime",
  "submitter_id",
  "file_name",
  "uuid"
)


# ----------------------------------------------------------------------
# Remove obvious outcome variables and technical duplicate variables
# ----------------------------------------------------------------------

clinical_candidate_names <- names(final_data)

clinical_candidate_names <- clinical_candidate_names[
  !grepl(
    paste(outcome_patterns, collapse = "|"),
    clinical_candidate_names,
    ignore.case = TRUE
  )
]

clinical_candidate_names <- clinical_candidate_names[
  !grepl(
    paste(technical_patterns, collapse = "|"),
    clinical_candidate_names,
    ignore.case = TRUE
  )
]

clinical_candidate_names <- setdiff(
  clinical_candidate_names,
  c(
    "ecosystem_state",
    "M3_patient_id"
  )
)

cat(
  "\nCandidate baseline clinical variables:",
  length(clinical_candidate_names),
  "\n"
)


# ======================================================================
# 17. CLASSIFY BASELINE VARIABLES
# ======================================================================

baseline_continuous <- character(0)

baseline_categorical <- character(0)

for (v in clinical_candidate_names) {

  x <- final_data[[v]]

  numeric_x <- safe_numeric(x)

  n_numeric <- sum(
    is.finite(numeric_x),
    na.rm = TRUE
  )

  n_nonmissing <- sum(
    !is.na(x) &
      trimws(as.character(x)) != ""
  )

  unique_nonmissing <- length(
    unique(
      as.character(
        x[
          !is.na(x) &
            trimws(as.character(x)) != ""
        ]
      )
    )
  )

  if (
    n_numeric >= 5 &&
      unique_nonmissing >= 3
  ) {

    baseline_continuous <- c(
      baseline_continuous,
      v
    )

  } else if (
    n_nonmissing >= 5 &&
      unique_nonmissing >= 2
  ) {

    baseline_categorical <- c(
      baseline_categorical,
      v
    )
  }
}

cat(
  "Baseline continuous variables:",
  length(baseline_continuous),
  "\n"
)

cat(
  "Baseline categorical variables:",
  length(baseline_categorical),
  "\n"
)


# ======================================================================
# 18. LAYER A — CATEGORICAL VARIABLE ASSOCIATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER A — CATEGORICAL VARIABLE ASSOCIATION\n")
cat("=====================================================================\n")

categorical_results <- list()

counter <- 0

for (v in baseline_categorical) {

  x <- as.character(
    final_data[[v]]
  )

  state <- as.character(
    final_data$ecosystem_state
  )

  keep <- !is.na(x) &
    trimws(x) != "" &
    !is.na(state) &
    trimws(state) != ""

  if (sum(keep) < 5) {
    next
  }

  x <- x[keep]

  state <- state[keep]

  x <- factor(x)

  state <- factor(state)

  if (nlevels(x) < 2) {
    next
  }

  if (nlevels(state) < 2) {
    next
  }

  tab <- table(
    state,
    x
  )

  p <- safe_fisher_test(
    tab
  )

  if (is.na(p)) {
    next
  }

  counter <- counter + 1

  categorical_results[[counter]] <- data.frame(
    variable = v,
    n = sum(tab),
    n_levels = nlevels(x),
    test = ifelse(
      nrow(tab) == 2 &&
        ncol(tab) == 2,
      "Fisher_exact",
      "Fisher_exact_simulated"
    ),
    p_value = p,
    stringsAsFactors = FALSE
  )
}


if (counter > 0) {

  categorical_results <- categorical_results[
    seq_len(counter)
  ]

  categorical_results_df <- bind_rows(
    categorical_results
  )

  categorical_results_df$FDR <- p.adjust(
    categorical_results_df$p_value,
    method = "BH"
  )

  categorical_results_df <- categorical_results_df %>%
    arrange(
      FDR,
      p_value
    )

} else {

  categorical_results_df <- data.frame(
    variable = character(0),
    n = numeric(0),
    n_levels = numeric(0),
    test = character(0),
    p_value = numeric(0),
    FDR = numeric(0),
    stringsAsFactors = FALSE
  )
}

cat(
  "Categorical tests completed:",
  nrow(categorical_results_df),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(
    categorical_results_df$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

if (nrow(categorical_results_df) > 0) {

  print(
    head(
      categorical_results_df,
      20
    )
  )
}

write_csv(
  categorical_results_df,
  file.path(
    output_dir,
    "M3.5_LayerA_categorical_associations.csv"
  )
)


# ======================================================================
# 19. LAYER A — CONTINUOUS VARIABLE ASSOCIATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER A — CONTINUOUS VARIABLE ASSOCIATION\n")
cat("=====================================================================\n")

continuous_results <- list()

counter <- 0

for (v in baseline_continuous) {

  x <- safe_numeric(
    final_data[[v]]
  )

  state <- as.character(
    final_data$ecosystem_state
  )

  result <- safe_kruskal(
    x,
    state
  )

  if (is.na(result$p_value)) {
    next
  }

  counter <- counter + 1

  continuous_results[[counter]] <- data.frame(
    variable = v,
    n = result$n,
    n_states = result$n_states,
    statistic = result$statistic,
    p_value = result$p_value,
    stringsAsFactors = FALSE
  )
}


if (counter > 0) {

  continuous_results <- continuous_results[
    seq_len(counter)
  ]

  continuous_results_df <- bind_rows(
    continuous_results
  )

  continuous_results_df$FDR <- p.adjust(
    continuous_results_df$p_value,
    method = "BH"
  )

  continuous_results_df <- continuous_results_df %>%
    arrange(
      FDR,
      p_value
    )

} else {

  continuous_results_df <- data.frame(
    variable = character(0),
    n = numeric(0),
    n_states = numeric(0),
    statistic = numeric(0),
    p_value = numeric(0),
    FDR = numeric(0),
    stringsAsFactors = FALSE
  )
}

cat(
  "Continuous tests completed:",
  nrow(continuous_results_df),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(
    continuous_results_df$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

if (nrow(continuous_results_df) > 0) {

  print(
    head(
      continuous_results_df,
      20
    )
  )
}

write_csv(
  continuous_results_df,
  file.path(
    output_dir,
    "M3.5_LayerA_continuous_associations.csv"
  )
)


# ======================================================================
# 20. LAYER A — EXPLORATORY PAIRWISE COMPARISONS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER A — EXPLORATORY PAIRWISE CONTINUOUS COMPARISONS\n")
cat("=====================================================================\n")


# ----------------------------------------------------------------------
# Pairwise testing is performed ONLY for variables with global
# FDR < 0.05.
# ----------------------------------------------------------------------

significant_continuous <- character(0)

if (nrow(continuous_results_df) > 0) {

  significant_continuous <- continuous_results_df$variable[
    !is.na(continuous_results_df$FDR) &
      continuous_results_df$FDR < 0.05
  ]
}

cat(
  "Exploratory candidate variables:",
  length(significant_continuous),
  "\n"
)

pairwise_list <- list()

counter <- 0

if (length(significant_continuous) > 0) {

  for (v in significant_continuous) {

    pw <- safe_pairwise_wilcox(
      final_data,
      v
    )

    if (!is.null(pw)) {

      counter <- counter + 1

      pairwise_list[[counter]] <- pw
    }
  }
}


if (counter > 0) {

  pairwise_df <- bind_rows(
    pairwise_list
  )

  pairwise_df$FDR <- p.adjust(
    pairwise_df$p_value,
    method = "BH"
  )

} else {

  pairwise_df <- data.frame(
    variable = character(0),
    state_1 = character(0),
    state_2 = character(0),
    n_state_1 = numeric(0),
    n_state_2 = numeric(0),
    median_state_1 = numeric(0),
    median_state_2 = numeric(0),
    p_value = numeric(0),
    FDR = numeric(0),
    stringsAsFactors = FALSE
  )
}

cat(
  "Pairwise comparisons completed:",
  nrow(pairwise_df),
  "\n"
)

write_csv(
  pairwise_df,
  file.path(
    output_dir,
    "M3.5_LayerA_pairwise_continuous_comparisons.csv"
  )
)


# ======================================================================
# 21. LAYER A — STATE-WISE CONTINUOUS PROFILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER A — STATE-WISE CONTINUOUS CLINICAL PROFILES\n")
cat("=====================================================================\n")

profile_list <- list()

counter <- 0

for (v in baseline_continuous) {

  x <- safe_numeric(
    final_data[[v]]
  )

  state <- as.character(
    final_data$ecosystem_state
  )

  temp <- data.frame(
    variable = v,
    value = x,
    ecosystem_state = state,
    stringsAsFactors = FALSE
  )

  temp <- temp[
    is.finite(temp$value) &
      !is.na(temp$ecosystem_state),
    ,
    drop = FALSE
  ]

  if (nrow(temp) == 0) {
    next
  }

  temp$ecosystem_state <- factor(
    temp$ecosystem_state
  )

  temp_profile <- temp %>%
    group_by(
      variable,
      ecosystem_state
    ) %>%
    summarise(
      n = n(),
      mean = mean(value),
      median = median(value),
      sd = ifelse(
        n() > 1,
        sd(value),
        NA_real_
      ),
      min = safe_min(value),
      max = safe_max(value),
      .groups = "drop"
    )

  counter <- counter + 1

  profile_list[[counter]] <- temp_profile
}


if (counter > 0) {

  continuous_profiles <- bind_rows(
    profile_list
  )

} else {

  continuous_profiles <- data.frame(
    variable = character(0),
    ecosystem_state = character(0),
    n = numeric(0),
    mean = numeric(0),
    median = numeric(0),
    sd = numeric(0),
    min = numeric(0),
    max = numeric(0),
    stringsAsFactors = FALSE
  )
}

cat(
  "Continuous state profiles generated:",
  nrow(continuous_profiles),
  "\n"
)

write_csv(
  continuous_profiles,
  file.path(
    output_dir,
    "M3.5_LayerA_state_wise_continuous_profiles.csv"
  )
)


# ======================================================================
# 22. LAYER A — STATE-WISE CATEGORICAL PROFILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER A — STATE-WISE CATEGORICAL PROFILES\n")
cat("=====================================================================\n")

categorical_profile_list <- list()

counter <- 0

for (v in baseline_categorical) {

  x <- as.character(
    final_data[[v]]
  )

  state <- as.character(
    final_data$ecosystem_state
  )

  temp <- data.frame(
    variable = v,
    category = x,
    ecosystem_state = state,
    stringsAsFactors = FALSE
  )

  temp <- temp[
    !is.na(temp$category) &
      trimws(temp$category) != "" &
      !is.na(temp$ecosystem_state),
    ,
    drop = FALSE
  ]

  if (nrow(temp) == 0) {
    next
  }

  profile <- temp %>%
    count(
      variable,
      ecosystem_state,
      category,
      name = "n"
    ) %>%
    group_by(
      variable,
      ecosystem_state
    ) %>%
    mutate(
      percentage = 100 * n / sum(n)
    ) %>%
    ungroup()

  counter <- counter + 1

  categorical_profile_list[[counter]] <- profile
}


if (counter > 0) {

  categorical_profiles <- bind_rows(
    categorical_profile_list
  )

} else {

  categorical_profiles <- data.frame(
    variable = character(0),
    ecosystem_state = character(0),
    category = character(0),
    n = numeric(0),
    percentage = numeric(0),
    stringsAsFactors = FALSE
  )
}

cat(
  "Categorical state profiles generated:",
  nrow(categorical_profiles),
  "\n"
)

write_csv(
  categorical_profiles,
  file.path(
    output_dir,
    "M3.5_LayerA_state_wise_categorical_profiles.csv"
  )
)


# ======================================================================
# 23. LAYER A — SIGNIFICANT PLOTS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER A — SIGNIFICANT / EXPLORATORY PLOTS\n")
cat("=====================================================================\n")


# ----------------------------------------------------------------------
# Continuous boxplots
# ----------------------------------------------------------------------

if (length(significant_continuous) > 0) {

  for (v in significant_continuous) {

    plot_data <- data.frame(
      value = safe_numeric(
        final_data[[v]]
      ),
      ecosystem_state =
        final_data$ecosystem_state
    )

    plot_data <- plot_data[
      is.finite(plot_data$value) &
        !is.na(plot_data$ecosystem_state),
      ,
      drop = FALSE
    ]

    if (nrow(plot_data) < 3) {
      next
    }

    p <- ggplot(
      plot_data,
      aes(
        x = ecosystem_state,
        y = value
      )
    ) +
      geom_boxplot(
        outlier.shape = NA
      ) +
      geom_jitter(
        width = 0.15,
        alpha = 0.6
      ) +
      theme_bw() +
      labs(
        title = paste(
          "Baseline clinical variable:",
          v
        ),
        x = "Frozen ecosystem state",
        y = v
      )

    safe_name <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      v
    )

    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          "M3.5_LayerA_boxplot_",
          safe_name,
          ".png"
        )
      ),
      plot = p,
      width = 7,
      height = 5,
      dpi = 300
    )
  }

} else {

  cat(
    "No FDR-significant continuous clinical variables.\n"
  )
}


# ----------------------------------------------------------------------
# Significant categorical plots
# ----------------------------------------------------------------------

significant_categorical <- character(0)

if (nrow(categorical_results_df) > 0) {

  significant_categorical <-
    categorical_results_df$variable[
      !is.na(categorical_results_df$FDR) &
        categorical_results_df$FDR < 0.05
    ]
}

if (length(significant_categorical) > 0) {

  for (v in significant_categorical) {

    plot_data <- data.frame(
      category = as.character(
        final_data[[v]]
      ),
      ecosystem_state =
        final_data$ecosystem_state
    )

    plot_data <- plot_data[
      !is.na(plot_data$category) &
        trimws(plot_data$category) != "" &
        !is.na(plot_data$ecosystem_state),
      ,
      drop = FALSE
    ]

    if (nrow(plot_data) < 3) {
      next
    }

    plot_data <- plot_data %>%
      count(
        ecosystem_state,
        category
      ) %>%
      group_by(
        ecosystem_state
      ) %>%
      mutate(
        percentage = 100 * n / sum(n)
      ) %>%
      ungroup()

    p <- ggplot(
      plot_data,
      aes(
        x = ecosystem_state,
        y = percentage,
        fill = category
      )
    ) +
      geom_col(
        position = "stack"
      ) +
      theme_bw() +
      labs(
        title = paste(
          "Baseline categorical variable:",
          v
        ),
        x = "Frozen ecosystem state",
        y = "Percentage"
      )

    safe_name <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      v
    )

    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          "M3.5_LayerA_categorical_",
          safe_name,
          ".png"
        )
      ),
      plot = p,
      width = 8,
      height = 5,
      dpi = 300
    )
  }

} else {

  cat(
    "No FDR-significant categorical clinical variables.\n"
  )
}


# ======================================================================
# ======================================================================
# LAYER B — ECOSYSTEM-STATE SURVIVAL VALIDATION
# ======================================================================
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER B — ECOSYSTEM-STATE SURVIVAL VALIDATION\n")
cat("=====================================================================\n")


# ======================================================================
# 24. IDENTIFY OS TIME AND EVENT COLUMNS
# ======================================================================

os_time_col <- find_column(
  final_data,
  c(
    "OS_time_days",
    "OS_time_years",
    "OS_time",
    "overall_survival_time",
    "days_to_death"
  )
)

os_event_col <- find_column(
  final_data,
  c(
    "OS_event",
    "OS_event_raw",
    "OS_status",
    "event"
  )
)

if (is.na(os_time_col)) {

  cat(
    "No suitable OS time column found.\n"
  )

  survival_available <- FALSE

} else if (is.na(os_event_col)) {

  cat(
    "No suitable OS event column found.\n"
  )

  survival_available <- FALSE

} else {

  survival_available <- TRUE
}


# ======================================================================
# 25. SURVIVAL ANALYSIS
# ======================================================================

survival_p <- NA_real_

survival_n <- 0

survival_events <- 0

if (survival_available) {

  cat(
    "Selected OS time column:",
    os_time_col,
    "\n"
  )

  cat(
    "Selected OS event column:",
    os_event_col,
    "\n"
  )


  # --------------------------------------------------------------------
  # Convert OS time
  # --------------------------------------------------------------------

  os_time <- safe_numeric(
    final_data[[os_time_col]]
  )


  # --------------------------------------------------------------------
  # Convert event variable robustly
  # --------------------------------------------------------------------

  raw_event <- final_data[[os_event_col]]

  if (is.numeric(raw_event)) {

    os_event <- as.numeric(
      raw_event
    )

  } else {

    event_text <- tolower(
      trimws(
        as.character(raw_event)
      )
    )

    os_event <- rep(
      NA_real_,
      length(event_text)
    )

    os_event[
      event_text %in%
        c(
          "1",
          "dead",
          "deceased",
          "death",
          "yes",
          "true",
          "event"
        )
    ] <- 1

    os_event[
      event_text %in%
        c(
          "0",
          "alive",
          "living",
          "no",
          "false",
          "censored"
        )
    ] <- 0

    numeric_event <- suppressWarnings(
      as.numeric(event_text)
    )

    os_event[
      is.na(os_event) &
        numeric_event %in% c(0, 1)
    ] <- numeric_event[
      is.na(os_event) &
        numeric_event %in% c(0, 1)
    ]
  }


  # --------------------------------------------------------------------
  # Build survival dataset
  # --------------------------------------------------------------------

  survival_data <- data.frame(
    OS_time = os_time,
    OS_event = os_event,
    ecosystem_state =
      final_data$ecosystem_state,
    M3_patient_id =
      final_data$M3_patient_id,
    stringsAsFactors = FALSE
  )

  survival_data <- survival_data[
    is.finite(survival_data$OS_time) &
      survival_data$OS_time > 0 &
      survival_data$OS_event %in% c(0, 1) &
      !is.na(survival_data$ecosystem_state),
    ,
    drop = FALSE
  ]


  survival_n <- nrow(
    survival_data
  )

  survival_events <- sum(
    survival_data$OS_event == 1,
    na.rm = TRUE
  )


  cat(
    "\nSurvival patients:",
    survival_n,
    "\n"
  )

  cat(
    "OS events:",
    survival_events,
    "\n"
  )

  cat(
    "Ecosystem states represented:",
    length(
      unique(
        survival_data$ecosystem_state
      )
    ),
    "\n"
  )


  # --------------------------------------------------------------------
  # Require at least two states and events
  # --------------------------------------------------------------------

  if (
    survival_n >= 10 &&
      survival_events >= 2 &&
      length(
        unique(
          survival_data$ecosystem_state
        )
      ) >= 2
  ) {

    survival_data$ecosystem_state <- factor(
      survival_data$ecosystem_state
    )


    # ------------------------------------------------------------------
    # Kaplan-Meier
    # ------------------------------------------------------------------

    km_fit <- survfit(
      Surv(
        OS_time,
        OS_event
      ) ~ ecosystem_state,
      data = survival_data
    )


    # ------------------------------------------------------------------
    # Global log-rank
    # ------------------------------------------------------------------

    logrank <- survdiff(
      Surv(
        OS_time,
        OS_event
      ) ~ ecosystem_state,
      data = survival_data
    )

    survival_p <- 1 - pchisq(
      logrank$chisq,
      df = length(logrank$n) - 1
    )


    cat(
      "\nGlobal log-rank P-value:",
      format.pval(
        survival_p,
        digits = 6
      ),
      "\n"
    )


    # ------------------------------------------------------------------
    # Save KM plot
    # ------------------------------------------------------------------

    km_plot <- ggsurvplot(
      km_fit,
      data = survival_data,
      risk.table = TRUE,
      pval = TRUE,
      conf.int = TRUE,
      xlab = "Overall survival time",
      ylab = "Survival probability",
      title = "Overall survival by frozen ecosystem state",
      legend.title = "Ecosystem state",
      ggtheme = theme_bw()
    )

    ggsave(
      filename = file.path(
        output_dir,
        "M3.5_LayerB_overall_survival_KM.png"
      ),
      plot = km_plot$plot,
      width = 8,
      height = 6,
      dpi = 300
    )


    # ------------------------------------------------------------------
    # Save survival dataset
    # ------------------------------------------------------------------

    write_csv(
      survival_data,
      file.path(
        output_dir,
        "M3.5_LayerB_survival_dataset.csv"
      )
    )


    # ------------------------------------------------------------------
    # Save log-rank result
    # ------------------------------------------------------------------

    survival_result <- data.frame(
      analysis = "Overall survival by frozen ecosystem state",
      time_column = os_time_col,
      event_column = os_event_col,
      n_patients = survival_n,
      n_events = survival_events,
      n_states = length(
        unique(
          survival_data$ecosystem_state
        )
      ),
      logrank_chisq = unname(
        logrank$chisq
      ),
      logrank_p_value = survival_p,
      stringsAsFactors = FALSE
    )

  } else {

    cat(
      "\nInsufficient survival data for KM analysis.\n"
    )

    survival_result <- data.frame(
      analysis = "Overall survival by frozen ecosystem state",
      time_column = os_time_col,
      event_column = os_event_col,
      n_patients = survival_n,
      n_events = survival_events,
      n_states = length(
        unique(
          survival_data$ecosystem_state
        )
      ),
      logrank_chisq = NA_real_,
      logrank_p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }

} else {

  survival_result <- data.frame(
    analysis = "Overall survival by frozen ecosystem state",
    time_column = os_time_col,
    event_column = os_event_col,
    n_patients = 0,
    n_events = 0,
    n_states = 0,
    logrank_chisq = NA_real_,
    logrank_p_value = NA_real_,
    stringsAsFactors = FALSE
  )
}


write_csv(
  survival_result,
  file.path(
    output_dir,
    "M3.5_LayerB_survival_result.csv"
  )
)


# ======================================================================
# ======================================================================
# LAYER C — BIOLOGICAL / ECOSYSTEM CHARACTERIZATION
# ======================================================================
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C — BIOLOGICAL / ECOSYSTEM CHARACTERIZATION\n")
cat("=====================================================================\n")

cat(
  "This layer describes molecular/ecosystem variables already present\n"
)

cat(
  "in the frozen M3.2 data. It does not redefine ecosystem states.\n"
)


# ======================================================================
# 26. IDENTIFY TRUE M3.2 BIOLOGICAL VARIABLES
# ======================================================================
#
# IMPORTANT:
#
# The current M3.2 patient state file contains:
#
# sample_id
# ecosystem_state
# cluster_numeric
#
# Therefore cluster_numeric is NOT considered a biological molecular
# axis. It is simply the numeric encoding of the already frozen state.
#
# We explicitly exclude:
#
#   sample_id
#   ecosystem_state
#   cluster_numeric
#   M3_patient_id
#
# This prevents the previous Layer C error.
# ======================================================================

m32_base_columns <- c(
  "sample_id",
  "sample_barcode",
  "barcode",
  "patient",
  "M3_patient_id",
  "ecosystem_state",
  "cluster_numeric",
  "cluster",
  "state"
)

m32_extra_columns <- setdiff(
  names(m32),
  m32_base_columns
)

# Remove any accidental duplicated ID column
m32_extra_columns <- m32_extra_columns[
  !grepl(
    "patient_id|sample_id|barcode|cluster_numeric|ecosystem_state",
    m32_extra_columns,
    ignore.case = TRUE
  )
]

cat(
  "\nPotential M3.2 molecular/ecosystem variables:",
  length(m32_extra_columns),
  "\n"
)


# ======================================================================
# 27. IF NO MOLECULAR VARIABLES EXIST, REPORT AND SKIP SAFELY
# ======================================================================

if (length(m32_extra_columns) == 0) {

  cat(
    "\nM3.2 frozen patient-state file contains no additional molecular\n"
  )

  cat(
    "or ecosystem-axis variables beyond the frozen state itself.\n"
  )

  cat(
    "Therefore Layer C molecular association testing is skipped safely.\n"
  )

  biological_continuous_df <- data.frame(
    variable = character(0),
    n = numeric(0),
    n_states = numeric(0),
    statistic = numeric(0),
    p_value = numeric(0),
    FDR = numeric(0),
    stringsAsFactors = FALSE
  )

  biological_categorical_df <- data.frame(
    variable = character(0),
    n = numeric(0),
    n_levels = numeric(0),
    test = character(0),
    p_value = numeric(0),
    FDR = numeric(0),
    stringsAsFactors = FALSE
  )

  biological_profiles_df <- data.frame(
    variable = character(0),
    ecosystem_state = character(0),
    n = numeric(0),
    mean = numeric(0),
    median = numeric(0),
    sd = numeric(0),
    min = numeric(0),
    max = numeric(0),
    stringsAsFactors = FALSE
  )

  cat(
    "\nLayer C molecular variables tested: 0\n"
  )

} else {


  # ====================================================================
  # 28. MERGE TRUE M3.2 BIOLOGICAL VARIABLES
  # ====================================================================

  bio_source <- m32 %>%
    select(
      M3_patient_id,
      ecosystem_state,
      all_of(m32_extra_columns)
    )

  bio_data <- final_data %>%
    select(
      M3_patient_id,
      ecosystem_state
    ) %>%
    distinct(
      M3_patient_id,
      .keep_all = TRUE
    ) %>%
    inner_join(
      bio_source,
      by = "M3_patient_id",
      suffix = c(
        "_clinical",
        "_M32"
      )
    )


  # --------------------------------------------------------------------
  # Guarantee state vector exists
  # --------------------------------------------------------------------

  if (!("ecosystem_state" %in% names(bio_data))) {

    state_candidates <- grep(
      "^ecosystem_state",
      names(bio_data),
      value = TRUE
    )

    if (length(state_candidates) == 0) {

      stop(
        "Layer C could not recover ecosystem_state."
      )
    }

    bio_data$ecosystem_state <-
      bio_data[[state_candidates[1]]]
  }


  bio_data$ecosystem_state <- as.character(
    bio_data$ecosystem_state
  )


  # ====================================================================
  # 29. CLASSIFY LAYER C VARIABLES
  # ====================================================================

  bio_continuous <- character(0)

  bio_categorical <- character(0)

  for (v in m32_extra_columns) {

    if (!(v %in% names(bio_data))) {
      next
    }

    x <- bio_data[[v]]

    numeric_x <- safe_numeric(
      x
    )

    n_numeric <- sum(
      is.finite(numeric_x),
      na.rm = TRUE
    )

    n_nonmissing <- sum(
      !is.na(x)
    )

    unique_nonmissing <- length(
      unique(
        x[
          !is.na(x)
        ]
      )
    )

    if (
      n_numeric >= 5 &&
        unique_nonmissing >= 3
    ) {

      bio_continuous <- c(
        bio_continuous,
        v
      )

    } else if (
      n_nonmissing >= 5 &&
        unique_nonmissing >= 2
    ) {

      bio_categorical <- c(
        bio_categorical,
        v
      )
    }
  }


  cat(
    "M3.2 continuous variables:",
    length(bio_continuous),
    "\n"
  )

  cat(
    "M3.2 categorical variables:",
    length(bio_categorical),
    "\n"
  )


  # ====================================================================
  # 30. LAYER C CONTINUOUS VARIABLES
  # ====================================================================

  biological_results <- list()

  counter <- 0

  for (v in bio_continuous) {

    value <- safe_numeric(
      bio_data[[v]]
    )

    state <- as.character(
      bio_data$ecosystem_state
    )

    # ---------------------------------------------------------------
    # CRITICAL SAFETY CHECK
    #
    # Never construct a data.frame if lengths differ.
    # ---------------------------------------------------------------

    if (
      length(value) != length(state)
    ) {

      cat(
        "Skipping variable due to length mismatch:",
        v,
        "\n"
      )

      next
    }

    temp <- data.frame(
      value = value,
      ecosystem_state = state,
      stringsAsFactors = FALSE
    )

    temp <- temp[
      is.finite(temp$value) &
        !is.na(temp$ecosystem_state),
      ,
      drop = FALSE
    ]

    if (nrow(temp) < 3) {
      next
    }

    result <- safe_kruskal(
      temp$value,
      temp$ecosystem_state
    )

    if (is.na(result$p_value)) {
      next
    }

    counter <- counter + 1

    biological_results[[counter]] <- data.frame(
      variable = v,
      n = result$n,
      n_states = result$n_states,
      statistic = result$statistic,
      p_value = result$p_value,
      stringsAsFactors = FALSE
    )
  }


  if (counter > 0) {

    biological_continuous_df <- bind_rows(
      biological_results
    )

    biological_continuous_df$FDR <- p.adjust(
      biological_continuous_df$p_value,
      method = "BH"
    )

    biological_continuous_df <- biological_continuous_df %>%
      arrange(
        FDR,
        p_value
      )

  } else {

    biological_continuous_df <- data.frame(
      variable = character(0),
      n = numeric(0),
      n_states = numeric(0),
      statistic = numeric(0),
      p_value = numeric(0),
      FDR = numeric(0),
      stringsAsFactors = FALSE
    )
  }


  # ====================================================================
  # 31. LAYER C CATEGORICAL VARIABLES
  # ====================================================================

  biological_categorical_results <- list()

  counter <- 0

  for (v in bio_categorical) {

    x <- as.character(
      bio_data[[v]]
    )

    state <- as.character(
      bio_data$ecosystem_state
    )

    keep <- !is.na(x) &
      trimws(x) != "" &
      !is.na(state)

    if (sum(keep) < 5) {
      next
    }

    x <- factor(
      x[keep]
    )

    state <- factor(
      state[keep]
    )

    if (
      nlevels(x) < 2 ||
        nlevels(state) < 2
    ) {
      next
    }

    tab <- table(
      state,
      x
    )

    p <- safe_fisher_test(
      tab
    )

    if (is.na(p)) {
      next
    }

    counter <- counter + 1

    biological_categorical_results[[counter]] <-
      data.frame(
        variable = v,
        n = sum(tab),
        n_levels = nlevels(x),
        test = ifelse(
          nrow(tab) == 2 &&
            ncol(tab) == 2,
          "Fisher_exact",
          "Fisher_exact_simulated"
        ),
        p_value = p,
        stringsAsFactors = FALSE
      )
  }


  if (counter > 0) {

    biological_categorical_df <- bind_rows(
      biological_categorical_results
    )

    biological_categorical_df$FDR <- p.adjust(
      biological_categorical_df$p_value,
      method = "BH"
    )

    biological_categorical_df <- biological_categorical_df %>%
      arrange(
        FDR,
        p_value
      )

  } else {

    biological_categorical_df <- data.frame(
      variable = character(0),
      n = numeric(0),
      n_levels = numeric(0),
      test = character(0),
      p_value = numeric(0),
      FDR = numeric(0),
      stringsAsFactors = FALSE
    )
  }


  # ====================================================================
  # 32. LAYER C STATE-WISE PROFILES
  # ====================================================================

  bio_profile_list <- list()

  counter <- 0

  for (v in bio_continuous) {

    value <- safe_numeric(
      bio_data[[v]]
    )

    state <- as.character(
      bio_data$ecosystem_state
    )

    if (
      length(value) != length(state)
    ) {
      next
    }

    temp <- data.frame(
      value = value,
      ecosystem_state = state,
      stringsAsFactors = FALSE
    )

    temp <- temp[
      is.finite(temp$value) &
        !is.na(temp$ecosystem_state),
      ,
      drop = FALSE
    ]

    if (nrow(temp) == 0) {
      next
    }

    profile <- temp %>%
      group_by(
        ecosystem_state
      ) %>%
      summarise(
        n = n(),
        mean = mean(value),
        median = median(value),
        sd = ifelse(
          n() > 1,
          sd(value),
          NA_real_
        ),
        min = safe_min(value),
        max = safe_max(value),
        .groups = "drop"
      ) %>%
      mutate(
        variable = v,
        .before = 1
      )

    counter <- counter + 1

    bio_profile_list[[counter]] <- profile
  }


  if (counter > 0) {

    biological_profiles_df <- bind_rows(
      bio_profile_list
    )

  } else {

    biological_profiles_df <- data.frame(
      variable = character(0),
      ecosystem_state = character(0),
      n = numeric(0),
      mean = numeric(0),
      median = numeric(0),
      sd = numeric(0),
      min = numeric(0),
      max = numeric(0),
      stringsAsFactors = FALSE
    )
  }
}


# ======================================================================
# 33. SAVE LAYER C RESULTS
# ======================================================================

write_csv(
  biological_continuous_df,
  file.path(
    output_dir,
    "M3.5_LayerC_biological_continuous_associations.csv"
  )
)

write_csv(
  biological_categorical_df,
  file.path(
    output_dir,
    "M3.5_LayerC_biological_categorical_associations.csv"
  )
)

write_csv(
  biological_profiles_df,
  file.path(
    output_dir,
    "M3.5_LayerC_biological_state_profiles.csv"
  )
)


cat(
  "\nLayer C continuous tests:",
  nrow(biological_continuous_df),
  "\n"
)

cat(
  "Layer C categorical tests:",
  nrow(biological_categorical_df),
  "\n"
)


# ======================================================================
# ======================================================================
# 34. ECOSYSTEM STATE CHARACTERIZATION
# ======================================================================
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("ECOSYSTEM STATE CHARACTERIZATION\n")
cat("=====================================================================\n")

ecosystem_characterization <- final_data %>%
  count(
    ecosystem_state,
    name = "n"
  ) %>%
  mutate(
    percentage = 100 * n / sum(n)
  )

print(
  ecosystem_characterization
)

write_csv(
  ecosystem_characterization,
  file.path(
    output_dir,
    "M3.5_ecosystem_state_characterization.csv"
  )
)


# ======================================================================
# 35. OVERALL SIGNIFICANT ASSOCIATIONS
# ======================================================================

significant_continuous_count <- sum(
  continuous_results_df$FDR < 0.05,
  na.rm = TRUE
)

significant_categorical_count <- sum(
  categorical_results_df$FDR < 0.05,
  na.rm = TRUE
)

significant_bio_continuous_count <- sum(
  biological_continuous_df$FDR < 0.05,
  na.rm = TRUE
)

significant_bio_categorical_count <- sum(
  biological_categorical_df$FDR < 0.05,
  na.rm = TRUE
)


significant_associations <- data.frame(
  layer = c(
    "Layer_A_Baseline_Continuous",
    "Layer_A_Baseline_Categorical",
    "Layer_C_Biological_Continuous",
    "Layer_C_Biological_Categorical"
  ),
  significant_FDR_lt_0_05 = c(
    significant_continuous_count,
    significant_categorical_count,
    significant_bio_continuous_count,
    significant_bio_categorical_count
  ),
  stringsAsFactors = FALSE
)

write_csv(
  significant_associations,
  file.path(
    output_dir,
    "M3.5_significant_clinical_associations.csv"
  )
)

cat("\n")
cat("Significant associations by layer:\n")
print(
  significant_associations
)


# ======================================================================
# 36. OVERALL VALIDATION SUMMARY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.5 OVERALL VALIDATION SUMMARY\n")
cat("=====================================================================\n")

summary_df <- data.frame(
  metric = c(

    "clinical_file",

    "M3.2_frozen_state_file",

    "clinical_ID_column",

    "M3.2_ID_column",

    "M3.2_state_column",

    "n_frozen_patients",

    "n_clinical_patients",

    "n_matched_patients",

    "match_percentage",

    "n_final_patients",

    "n_ecosystem_states",

    "LayerA_baseline_continuous_tested",

    "LayerA_baseline_categorical_tested",

    "LayerA_continuous_FDR_lt_0.05",

    "LayerA_categorical_FDR_lt_0.05",

    "LayerA_pairwise_tests",

    "LayerB_survival_available",

    "LayerB_survival_patients",

    "LayerB_OS_events",

    "LayerB_logrank_p_value",

    "LayerC_biological_continuous_tested",

    "LayerC_biological_categorical_tested",

    "LayerC_continuous_FDR_lt_0.05",

    "LayerC_categorical_FDR_lt_0.05"

  ),

  value = c(

    clinical_file,

    m32_state_file,

    clinical_id_col,

    m32_id_col,

    state_col,

    length(frozen_patient_ids),

    length(clinical_patient_ids),

    length(matched_ids),

    round(match_percentage, 4),

    nrow(final_data),

    length(
      unique(
        final_data$ecosystem_state
      )
    ),

    length(baseline_continuous),

    length(baseline_categorical),

    significant_continuous_count,

    significant_categorical_count,

    nrow(pairwise_df),

    survival_available,

    survival_n,

    survival_events,

    ifelse(
      is.na(survival_p),
      NA_character_,
      format(
        survival_p,
        scientific = TRUE
      )
    ),

    nrow(biological_continuous_df),

    nrow(biological_categorical_df),

    significant_bio_continuous_count,

    significant_bio_categorical_count

  ),

  stringsAsFactors = FALSE
)

print(
  summary_df
)

write_csv(
  summary_df,
  file.path(
    output_dir,
    "M3.5_overall_validation_summary.csv"
  )
)


# ======================================================================
# 37. WRITE HUMAN-READABLE DIAGNOSTIC FILE
# ======================================================================

diagnostic_file <- file.path(
  output_dir,
  "M3.5_run_diagnostic.txt"
)

sink(
  diagnostic_file
)

cat(
  "M3.5 FROZEN ECOSYSTEM STATE CLINICAL / PHENOTYPIC VALIDATION\n"
)

cat(
  "=====================================================================\n\n"
)

cat(
  "Clinical file:",
  clinical_file,
  "\n"
)

cat(
  "M3.2 frozen state file:",
  m32_state_file,
  "\n\n"
)

cat(
  "Clinical ID column:",
  clinical_id_col,
  "\n"
)

cat(
  "M3.2 ID column:",
  m32_id_col,
  "\n"
)

cat(
  "M3.2 state column:",
  state_col,
  "\n\n"
)

cat(
  "Frozen M3.2 patients:",
  length(frozen_patient_ids),
  "\n"
)

cat(
  "Clinical patients:",
  length(clinical_patient_ids),
  "\n"
)

cat(
  "Matched patients:",
  length(matched_ids),
  "\n"
)

cat(
  "Match percentage:",
  round(match_percentage, 2),
  "%\n\n"
)

cat(
  "Frozen ecosystem states:\n"
)

print(
  ecosystem_characterization
)

cat(
  "\nLAYER A\n"
)

cat(
  "Baseline continuous variables tested:",
  length(baseline_continuous),
  "\n"
)

cat(
  "Baseline categorical variables tested:",
  length(baseline_categorical),
  "\n"
)

cat(
  "Significant continuous variables:",
  significant_continuous_count,
  "\n"
)

cat(
  "Significant categorical variables:",
  significant_categorical_count,
  "\n\n"
)

cat(
  "LAYER B\n"
)

cat(
  "Survival available:",
  survival_available,
  "\n"
)

cat(
  "Survival patients:",
  survival_n,
  "\n"
)

cat(
  "OS events:",
  survival_events,
  "\n"
)

cat(
  "Log-rank p-value:",
  survival_p,
  "\n\n"
)

cat(
  "LAYER C\n"
)

cat(
  "Biological continuous variables tested:",
  nrow(
    biological_continuous_df
  ),
  "\n"
)

cat(
  "Biological categorical variables tested:",
  nrow(
    biological_categorical_df
  ),
  "\n"
)

cat(
  "Biological continuous FDR < 0.05:",
  significant_bio_continuous_count,
  "\n"
)

cat(
  "Biological categorical FDR < 0.05:",
  significant_bio_categorical_count,
  "\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "Ecosystem states were treated as frozen labels from M3.2.\n"
)

cat(
  "No clustering or state redefinition was performed in M3.5.\n"
)

cat(
  "Survival variables were excluded from Layer A baseline testing.\n"
)

cat(
  "Layer C does not treat cluster_numeric as a biological variable.\n"
)

cat(
  "Layer C safely skips variables that are absent or all-missing.\n"
)

sink()


# ======================================================================
# 38. FINAL OUTPUT FILE LIST
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.5 COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n")

cat(
  "Frozen patients:",
  length(frozen_patient_ids),
  "\n"
)

cat(
  "Clinical patients:",
  length(clinical_patient_ids),
  "\n"
)

cat(
  "Matched clinical patients:",
  length(matched_ids),
  "\n"
)

cat(
  "Match percentage:",
  round(match_percentage, 2),
  "%\n"
)

cat(
  "Ecosystem states:",
  length(
    unique(
      final_data$ecosystem_state
    )
  ),
  "\n"
)

cat(
  "Layer A categorical significant variables:",
  significant_categorical_count,
  "\n"
)

cat(
  "Layer A continuous significant variables:",
  significant_continuous_count,
  "\n"
)

cat(
  "Layer A pairwise tests:",
  nrow(pairwise_df),
  "\n"
)

cat(
  "Layer B survival analysis:",
  survival_available,
  "\n"
)

if (survival_available) {

  cat(
    "Layer B OS log-rank p-value:",
    ifelse(
      is.na(survival_p),
      "NA",
      format(
        survival_p,
        digits = 6
      )
    ),
    "\n"
  )
}

cat(
  "Layer C biological continuous variables:",
  nrow(
    biological_continuous_df
  ),
  "\n"
)

cat(
  "Layer C biological categorical variables:",
  nrow(
    biological_categorical_df
  ),
  "\n"
)

cat("\n")
cat("Output directory:\n")
cat(output_dir, "\n")


cat("\n")
cat("Important primary files:\n\n")

cat(
  "1. Final matched dataset:\n",
  file.path(
    output_dir,
    "M3.5_FINAL_frozen_ecosystem_clinical_dataset.csv"
  ),
  "\n\n"
)

cat(
  "2. Patient-ID diagnostic:\n",
  file.path(
    output_dir,
    "M3.5_patient_ID_matching_diagnostic.csv"
  ),
  "\n\n"
)

cat(
  "3. Layer A categorical associations:\n",
  file.path(
    output_dir,
    "M3.5_LayerA_categorical_associations.csv"
  ),
  "\n\n"
)

cat(
  "4. Layer A continuous associations:\n",
  file.path(
    output_dir,
    "M3.5_LayerA_continuous_associations.csv"
  ),
  "\n\n"
)

cat(
  "5. Layer A pairwise comparisons:\n",
  file.path(
    output_dir,
    "M3.5_LayerA_pairwise_continuous_comparisons.csv"
  ),
  "\n\n"
)

cat(
  "6. Layer A continuous state profiles:\n",
  file.path(
    output_dir,
    "M3.5_LayerA_state_wise_continuous_profiles.csv"
  ),
  "\n\n"
)

cat(
  "7. Layer A categorical state profiles:\n",
  file.path(
    output_dir,
    "M3.5_LayerA_state_wise_categorical_profiles.csv"
  ),
  "\n\n"
)

cat(
  "8. Layer B survival result:\n",
  file.path(
    output_dir,
    "M3.5_LayerB_survival_result.csv"
  ),
  "\n\n"
)

cat(
  "9. Layer B Kaplan-Meier plot:\n",
  file.path(
    output_dir,
    "M3.5_LayerB_overall_survival_KM.png"
  ),
  "\n\n"
)

cat(
  "10. Layer C biological continuous associations:\n",
  file.path(
    output_dir,
    "M3.5_LayerC_biological_continuous_associations.csv"
  ),
  "\n\n"
)

cat(
  "11. Layer C biological categorical associations:\n",
  file.path(
    output_dir,
    "M3.5_LayerC_biological_categorical_associations.csv"
  ),
  "\n\n"
)

cat(
  "12. Layer C biological state profiles:\n",
  file.path(
    output_dir,
    "M3.5_LayerC_biological_state_profiles.csv"
  ),
  "\n\n"
)

cat(
  "13. Ecosystem state characterization:\n",
  file.path(
    output_dir,
    "M3.5_ecosystem_state_characterization.csv"
  ),
  "\n\n"
)

cat(
  "14. Overall validation summary:\n",
  file.path(
    output_dir,
    "M3.5_overall_validation_summary.csv"
  ),
  "\n\n"
)

cat(
  "15. Run diagnostic:\n",
  diagnostic_file,
  "\n\n"
)


# ======================================================================
# 39. ELAPSED TIME
# ======================================================================

end_time <- Sys.time()

elapsed_minutes <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "mins"
  )
)

cat(
  "Elapsed time:",
  round(
    elapsed_minutes,
    2
  ),
  "minutes\n"
)

cat("\n")
cat("=====================================================================\n")
cat("M3.5 CLINICAL / PHENOTYPIC VALIDATION COMPLETE\n")
cat("=====================================================================\n")