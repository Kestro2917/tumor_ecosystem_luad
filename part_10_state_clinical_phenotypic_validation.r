# ======================================================================
# M3.5.1-C2 — INDEPENDENT BIOLOGICAL VALIDATION
#
# ANALYTIC FRAMEWORK
# ------------------
# Layer C2:
#
#   Frozen ecosystem state
#              |
#              v
#   independent biological features
#              |
#              +--> Kruskal-Wallis
#              |
#              +--> FDR correction
#              |
#              +--> effect sizes
#              |
#              +--> pairwise exploratory tests
#              |
#              +--> state-wise biological profiles
#              |
#              +--> biological heatmap
#              |
#              +--> PCA visualization
#
# ======================================================================


# ======================================================================
# 0. START
# ======================================================================

options(stringsAsFactors = FALSE)

cat("\n")
cat("=====================================================================\n")
cat("M3.5.1-C2 — INDEPENDENT BIOLOGICAL VALIDATION\n")
cat("=====================================================================\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("\n")


# ======================================================================
# 1. REQUIRED PACKAGES
# ======================================================================

required_packages <- c(
  "tidyverse",
  "readr",
  "dplyr",
  "ggplot2",
  "stringr",
  "purrr",
  "tidyr",
  "tibble",
  "rstatix",
  "scales"
)

cat("=====================================================================\n")
cat("CHECKING R PACKAGES\n")
cat("=====================================================================\n")

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    cat("Installing:", pkg, "\n")

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org",
      quiet = TRUE
    )
  }

  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}


# ======================================================================
# 2. INPUT / OUTPUT PATHS
# ======================================================================

M3_ROOT <- "/content/M3_LIONESS_entropy"

M3_FROZEN_FILE <- file.path(
  M3_ROOT,
  "M3.2_module_axes",
  "M3.2_patient_ecosystem_states.csv"
)

OUTPUT_DIR <- file.path(
  M3_ROOT,
  "M3.5_clinical_validation",
  "M3.5.1_C2_independent_biological_validation"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n")
cat("=====================================================================\n")
cat("INPUT / OUTPUT PATHS\n")
cat("=====================================================================\n")

cat("Frozen M3.2 file:\n")
cat(M3_FROZEN_FILE, "\n\n")

cat("Output directory:\n")
cat(OUTPUT_DIR, "\n\n")


# ======================================================================
# 3. CHECK FROZEN M3.2 FILE
# ======================================================================

cat("=====================================================================\n")
cat("CHECKING FROZEN M3.2 FILE\n")
cat("=====================================================================\n")

if (!file.exists(M3_FROZEN_FILE)) {
  stop(
    paste0(
      "ERROR: Frozen M3.2 file was not found:\n",
      M3_FROZEN_FILE
    )
  )
}

cat("Frozen M3.2 file found.\n\n")


# ======================================================================
# 4. READ FROZEN M3.2 STATES
# ======================================================================

cat("=====================================================================\n")
cat("READING FROZEN M3.2 STATES\n")
cat("=====================================================================\n")

frozen <- readr::read_csv(
  M3_FROZEN_FILE,
  show_col_types = FALSE
)

cat("Rows:", nrow(frozen), "\n")
cat("Columns:", ncol(frozen), "\n\n")

cat("Columns:\n")
print(colnames(frozen))

required_frozen_columns <- c(
  "sample_id",
  "ecosystem_state"
)

missing_frozen <- setdiff(
  required_frozen_columns,
  colnames(frozen)
)

if (length(missing_frozen) > 0) {

  stop(
    paste0(
      "ERROR: Required frozen columns missing: ",
      paste(missing_frozen, collapse = ", ")
    )
  )
}


# ======================================================================
# 5. SAFE PATIENT-ID NORMALIZATION
#
# IMPORTANT
# ----------
# We deliberately avoid the fragile regular expression that caused:
#
# Missing closing bracket on a bracket expression
#
# TCGA barcodes can contain:
#
# TCGA-XX-XXXX-01A-...
#
# The safest approach here is to normalize case/whitespace and then
# retain the first three barcode components.
# ======================================================================

standardize_patient_id <- function(x) {

  x <- as.character(x)

  x <- trimws(x)

  x <- toupper(x)

  # Remove whitespace
  x <- gsub("[[:space:]]+", "", x)

  # Remove quotation marks if present
  x <- gsub('"', "", x)
  x <- gsub("'", "", x)

  # Convert underscores to hyphens
  x <- gsub("_", "-", x)

  # TCGA patient identifier = first 3 barcode components
  #
  # Example:
  # TCGA-55-7574-01A
  #
  # becomes:
  #
  # TCGA-55-7574

  parts <- strsplit(x, "-", fixed = TRUE)

  result <- vapply(
    parts,
    function(z) {

      z <- z[z != ""]

      if (length(z) >= 3) {

        paste(z[1:3], collapse = "-")

      } else {

        paste(z, collapse = "-")
      }
    },
    character(1)
  )

  result
}


# ======================================================================
# 6. NORMALIZE FROZEN PATIENT IDS
# ======================================================================

frozen <- frozen %>%
  mutate(
    sample_id_original = as.character(sample_id),
    sample_id_clean = standardize_patient_id(sample_id),
    ecosystem_state = as.character(ecosystem_state)
  )

frozen <- frozen %>%
  filter(
    !is.na(sample_id_clean),
    sample_id_clean != "",
    !is.na(ecosystem_state),
    ecosystem_state != ""
  )

cat("\n")
cat("Valid frozen patients:", nrow(frozen), "\n")

duplicate_frozen <- frozen %>%
  count(sample_id_clean) %>%
  filter(n > 1)

if (nrow(duplicate_frozen) > 0) {

  cat(
    "WARNING: Duplicate normalized frozen patient IDs detected:",
    nrow(duplicate_frozen),
    "\n"
  )

  write_csv(
    duplicate_frozen,
    file.path(
      OUTPUT_DIR,
      "M3.5.1_C2_duplicate_frozen_patient_IDs.csv"
    )
  )

  stop(
    "Frozen M3.2 file contains duplicate normalized patient IDs."
  )

} else {

  cat("No duplicate frozen patient IDs detected.\n")
}


# ======================================================================
# 7. FROZEN STATE DISTRIBUTION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("FROZEN ECOSYSTEM STATE DISTRIBUTION\n")
cat("=====================================================================\n")

state_distribution <- frozen %>%
  count(ecosystem_state, name = "n") %>%
  mutate(
    percentage = 100 * n / sum(n)
  )

print(state_distribution)

write_csv(
  state_distribution,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_frozen_state_distribution.csv"
  )
)


# ======================================================================
# 8. SEARCH M3 DIRECTORY FOR CSV FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("SEARCHING M3 DIRECTORY FOR BIOLOGICAL DATA\n")
cat("=====================================================================\n")

all_csv_files <- list.files(
  M3_ROOT,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("CSV files found:", length(all_csv_files), "\n")


# ======================================================================
# 9. EXCLUSION RULES
#
# We exclude:
#
#   - M3.5 clinical files
#   - M3.5.1 files
#   - frozen state files
#   - state characterization files
#   - files whose primary purpose is to reproduce the ecosystem state
#
# This prevents obvious circularity.
# ======================================================================

exclude_patterns <- c(

  "M3\\.5",

  "M3\\.5\\.1",

  "patient_ecosystem_states",

  "state_sizes",

  "state_module_profiles",

  "state_standardized_module_profiles",

  "state_module_entropy_profiles",

  "state_directional_biology",

  "top_modules_by_state",

  "top_modules_by_ecosystem_state",

  "cluster_quality",

  "cluster_numeric",

  "bootstrap",

  "confusion_matrix",

  "classification_accuracy",

  "leave_one_out",

  "silhouette",

  "PERMANOVA",

  "pairwise_patient_distances",

  "state_centroid_distances",

  "state_module_effects",

  "module_contribution_to_state_separation",

  "validation_dataset",

  "merged_validation",

  "clinical"
)


exclude_regex <- paste(
  exclude_patterns,
  collapse = "|"
)

candidate_files <- all_csv_files[
  !grepl(
    exclude_regex,
    all_csv_files,
    ignore.case = TRUE
  )
]

cat(
  "Files remaining after exclusion:",
  length(candidate_files),
  "\n"
)


# ======================================================================
# 10. IDENTIFY PATIENT-LEVEL FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("IDENTIFYING PATIENT-LEVEL BIOLOGICAL FILES\n")
cat("=====================================================================\n")


inspect_file <- function(f) {

  result <- tryCatch({

    dat <- readr::read_csv(
      f,
      n_max = 5,
      show_col_types = FALSE
    )

    nm <- colnames(dat)

    has_sample_id <- any(
      tolower(nm) %in% c(
        "sample_id",
        "sampleid",
        "sample",
        "barcode",
        "patient_id",
        "patient"
      )
    )

    tibble(
      file = f,
      n_rows = NA_integer_,
      n_columns = ncol(dat),
      has_sample_id = has_sample_id,
      columns = paste(
        nm,
        collapse = ";"
      )
    )

  }, error = function(e) {

    tibble(
      file = f,
      n_rows = NA_integer_,
      n_columns = NA_integer_,
      has_sample_id = FALSE,
      columns = NA_character_
    )
  })

  result
}


file_inventory <- purrr::map_dfr(
  candidate_files,
  inspect_file
)

patient_candidates <- file_inventory %>%
  filter(
    has_sample_id == TRUE
  )

cat(
  "Patient-level candidate files:",
  nrow(patient_candidates),
  "\n\n"
)

print(patient_candidates)


# ======================================================================
# 11. FUNCTION TO FIND PATIENT ID COLUMN
# ======================================================================

find_patient_id_column <- function(dat) {

  nm <- colnames(dat)

  preferred <- c(
    "sample_id",
    "sampleid",
    "sample",
    "barcode",
    "patient_id",
    "patient"
  )

  idx <- match(
    tolower(preferred),
    tolower(nm)
  )

  idx <- idx[!is.na(idx)]

  if (length(idx) > 0) {
    return(nm[idx[1]])
  }

  return(NULL)
}


# ======================================================================
# 12. READ PATIENT-LEVEL BIOLOGICAL FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING PATIENT-LEVEL BIOLOGICAL FEATURES\n")
cat("=====================================================================\n")


biological_tables <- list()


for (i in seq_len(nrow(patient_candidates))) {

  f <- patient_candidates$file[i]

  cat("\nReading:\n")
  cat(f, "\n")

  dat <- tryCatch({

    readr::read_csv(
      f,
      show_col_types = FALSE,
      progress = FALSE
    )

  }, error = function(e) {

    cat(
      "FAILED:",
      conditionMessage(e),
      "\n"
    )

    NULL
  })

  if (is.null(dat)) {
    next
  }

  id_col <- find_patient_id_column(dat)

  if (is.null(id_col)) {
    next
  }

  dat <- dat %>%
    mutate(
      sample_id_raw_for_C2 = as.character(.data[[id_col]]),
      sample_id_clean = standardize_patient_id(
        .data[[id_col]]
      )
    )

  # Remove duplicate columns caused by repeated joins
  dat <- dat %>%
    select(
      -any_of(
        c(
          "sample_id_clean.1",
          "sample_id_clean.2"
        )
      )
    )

  # Keep one row per patient
  dat <- dat %>%
    distinct(
      sample_id_clean,
      .keep_all = TRUE
    )

  # Prefix columns by file name to avoid accidental collisions
  base_name <- basename(f)

  prefix <- tools::file_path_sans_ext(
    base_name
  )

  prefix <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    prefix
  )

  names(dat) <- make.unique(
    names(dat),
    sep = "_"
  )

  biological_tables[[length(biological_tables) + 1]] <- list(
    file = f,
    prefix = prefix,
    data = dat
  )
}


cat("\n")
cat(
  "Successfully read patient-level biological tables:",
  length(biological_tables),
  "\n"
)


# ======================================================================
# 13. MERGE BIOLOGICAL TABLES
# ======================================================================

if (length(biological_tables) == 0) {

  stop(
    "No patient-level biological files could be read."
  )
}


cat("\n")
cat("=====================================================================\n")
cat("MERGING PATIENT-LEVEL BIOLOGICAL FEATURES\n")
cat("=====================================================================\n")


# Start with frozen patients
biological_data <- frozen %>%
  select(
    sample_id_clean,
    ecosystem_state
  )


for (obj in biological_tables) {

  dat <- obj$data

  # Remove redundant state column if present
  remove_cols <- intersect(
    colnames(dat),
    c(
      "ecosystem_state",
      "cluster_numeric",
      "actual_state",
      "predicted_state"
    )
  )

  dat <- dat %>%
    select(
      -any_of(remove_cols)
    )

  # Keep patient identifier + features
  dat <- dat %>%
    select(
      sample_id_clean,
      everything()
    )

  # Remove duplicate patient IDs
  dat <- dat %>%
    distinct(
      sample_id_clean,
      .keep_all = TRUE
    )

  biological_data <- biological_data %>%
    left_join(
      dat,
      by = "sample_id_clean"
    )
}


cat(
  "Merged biological dataset dimensions:\n"
)

cat(
  "Rows:",
  nrow(biological_data),
  "\n"
)

cat(
  "Columns:",
  ncol(biological_data),
  "\n"
)


# ======================================================================
# 14. CHECK PATIENT RETENTION
# ======================================================================

retention <- tibble(
  metric = c(
    "Frozen patients",
    "Patients after biological merge",
    "Patients lost"
  ),
  value = c(
    nrow(frozen),
    nrow(biological_data),
    nrow(frozen) - nrow(biological_data)
  )
)

print(retention)

write_csv(
  retention,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_patient_retention.csv"
  )
)


# ======================================================================
# 15. IDENTIFY BIOLOGICAL VARIABLES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("IDENTIFYING BIOLOGICAL VARIABLES\n")
cat("=====================================================================\n")


non_biological_names <- c(

  "sample_id_clean",

  "sample_id",
  "sampleid",
  "sample",
  "barcode",
  "patient_id",
  "patient",

  "sample_id_original",
  "sample_id_raw_for_C2",

  "ecosystem_state",

  "cluster_numeric",
  "cluster",
  "cluster_id",

  "actual_state",
  "predicted_state",

  "state",
  "state_id",

  "k2",
  "k3",
  "k4",
  "k5",

  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5"
)


candidate_biological_names <- setdiff(
  colnames(biological_data),
  non_biological_names
)


# ======================================================================
# 16. REMOVE VARIABLES THAT ARE OBVIOUSLY STATE-DEFINITION VARIABLES
# ======================================================================

state_definition_patterns <- c(

  "ecosystem_state",

  "cluster",

  "predicted_state",

  "actual_state",

  "state_label",

  "state_id",

  "state_numeric",

  "^k[2-9]$",

  "^PC[0-9]+$",

  "silhouette",

  "classification",

  "confusion",

  "permanova",

  "centroid",

  "distance_to_state",

  "state_effect",

  "state_score",

  "state_profile",

  "state_rank"
)

state_definition_regex <- paste(
  state_definition_patterns,
  collapse = "|"
)

candidate_biological_names <- candidate_biological_names[
  !grepl(
    state_definition_regex,
    candidate_biological_names,
    ignore.case = TRUE
  )
]


cat(
  "Candidate biological variables:",
  length(candidate_biological_names),
  "\n")


# ======================================================================
# 17. IDENTIFY NUMERIC AND CATEGORICAL VARIABLES
# ======================================================================

is_numeric_variable <- function(x) {

  is.numeric(x) ||
    is.integer(x)
}


continuous_variables <- candidate_biological_names[
  vapply(
    biological_data[candidate_biological_names],
    is_numeric_variable,
    logical(1)
  )
]


categorical_variables <- candidate_biological_names[
  !candidate_biological_names %in%
    continuous_variables
]


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
# 18. INFORMATION FILTERING
#
# Continuous:
#   - at least 10 non-missing observations
#   - at least 3 unique values
#
# Categorical:
#   - at least 10 observations
#   - at least 2 levels
# ======================================================================

continuous_variables_filtered <- continuous_variables[
  vapply(
    biological_data[continuous_variables],
    function(x) {

      x <- suppressWarnings(
        as.numeric(x)
      )

      sum(is.finite(x)) >= 10 &&
        length(unique(x[is.finite(x)])) >= 3

    },
    logical(1)
  )
]


categorical_variables_filtered <- categorical_variables[
  vapply(
    biological_data[categorical_variables],
    function(x) {

      x <- as.character(x)

      x <- x[
        !is.na(x) &
          trimws(x) != ""
      ]

      length(x) >= 10 &&
        length(unique(x)) >= 2

    },
    logical(1)
  )
]


cat("\n")
cat("After information filtering:\n")

cat(
  "Continuous variables:",
  length(continuous_variables_filtered),
  "\n"
)

cat(
  "Categorical variables:",
  length(categorical_variables_filtered),
  "\n"
)


# ======================================================================
# 19. LAYER C2 — CONTINUOUS BIOLOGICAL ASSOCIATIONS
#
# Kruskal-Wallis test
#
# This tests whether the biological feature differs across the
# FROZEN ecosystem states.
#
# NO survival.
# NO Cox.
# NO hazard ratio.
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C2 — BIOLOGICAL CONTINUOUS ASSOCIATIONS\n")
cat("=====================================================================\n")


continuous_results <- purrr::map_dfr(
  continuous_variables_filtered,
  function(v) {

    dat <- biological_data %>%
      select(
        ecosystem_state,
        value = all_of(v)
      ) %>%
      mutate(
        value = suppressWarnings(
          as.numeric(value)
        )
      ) %>%
      filter(
        !is.na(value),
        is.finite(value),
        !is.na(ecosystem_state)
      )

    if (
      nrow(dat) < 10 ||
      n_distinct(dat$ecosystem_state) < 2
    ) {
      return(NULL)
    }

    kw <- tryCatch({

      kruskal.test(
        value ~ ecosystem_state,
        data = dat
      )

    }, error = function(e) {

      NULL
    })

    if (is.null(kw)) {
      return(NULL)
    }

    tibble(
      variable = v,
      n = nrow(dat),
      n_states = n_distinct(
        dat$ecosystem_state
      ),
      statistic = unname(
        kw$statistic
      ),
      p_value = kw$p.value
    )
  }
)


if (nrow(continuous_results) > 0) {

  continuous_results <- continuous_results %>%
    mutate(
      FDR = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    arrange(
      FDR,
      p_value
    )

} else {

  continuous_results <- tibble(
    variable = character(),
    n = integer(),
    n_states = integer(),
    statistic = numeric(),
    p_value = numeric(),
    FDR = numeric()
  )
}


cat(
  "Continuous tests completed:",
  nrow(continuous_results),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(
    continuous_results$FDR < 0.05,
    na.rm = TRUE
  ),
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
    OUTPUT_DIR,
    "M3.5.1_C2_biological_continuous_associations.csv"
  )
)


# ======================================================================
# 20. BIOLOGICAL EFFECT SIZE
#
# Epsilon-squared approximation for Kruskal-Wallis
# ======================================================================

continuous_effects <- continuous_results %>%
  mutate(
    epsilon_squared = case_when(

      n > 1 &
        n_states > 1 ~
        (statistic - n_states + 1) /
        (n - n_states),

      TRUE ~ NA_real_
    )
  ) %>%
  arrange(
    desc(epsilon_squared)
  )


write_csv(
  continuous_effects,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_biological_continuous_effect_sizes.csv"
  )
)


# ======================================================================
# 21. TOP BIOLOGICAL FEATURES
# ======================================================================

top_biological_features <- continuous_effects %>%
  filter(
    is.finite(epsilon_squared)
  ) %>%
  arrange(
    desc(epsilon_squared)
  ) %>%
  slice_head(
    n = 30
  )


write_csv(
  top_biological_features,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_top_biological_features.csv"
  )
)


# ======================================================================
# 22. PAIRWISE BIOLOGICAL COMPARISONS
#
# Only performed for FDR-significant continuous biological variables.
#
# These are exploratory and are NOT used to redefine the ecosystem
# states.
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C2 — PAIRWISE BIOLOGICAL COMPARISONS\n")
cat("=====================================================================\n")


significant_continuous_variables <- continuous_results %>%
  filter(
    FDR < 0.05
  ) %>%
  pull(
    variable
  )


pairwise_results <- purrr::map_dfr(
  significant_continuous_variables,
  function(v) {

    dat <- biological_data %>%
      select(
        ecosystem_state,
        value = all_of(v)
      ) %>%
      mutate(
        value = suppressWarnings(
          as.numeric(value)
        )
      ) %>%
      filter(
        is.finite(value),
        !is.na(ecosystem_state)
      )

    if (
      nrow(dat) < 10 ||
      n_distinct(dat$ecosystem_state) < 2
    ) {
      return(NULL)
    }

    out <- tryCatch({

      rstatix::wilcox_test(
        dat,
        value ~ ecosystem_state,
        p.adjust.method = "BH"
      )

    }, error = function(e) {

      NULL
    })

    if (is.null(out)) {
      return(NULL)
    }

    out %>%
      mutate(
        variable = v
      ) %>%
      select(
        variable,
        everything()
      )
  }
)


if (nrow(pairwise_results) == 0) {

  pairwise_results <- tibble(
    variable = character()
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
    OUTPUT_DIR,
    "M3.5.1_C2_pairwise_biological_comparisons.csv"
  )
)


# ======================================================================
# 23. LAYER C2 — CATEGORICAL BIOLOGICAL ASSOCIATIONS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C2 — BIOLOGICAL CATEGORICAL ASSOCIATIONS\n")
cat("=====================================================================\n")


categorical_results <- purrr::map_dfr(
  categorical_variables_filtered,
  function(v) {

    dat <- biological_data %>%
      select(
        ecosystem_state,
        value = all_of(v)
      ) %>%
      mutate(
        value = as.character(value),
        ecosystem_state = as.character(
          ecosystem_state
        )
      ) %>%
      filter(
        !is.na(value),
        trimws(value) != "",
        !is.na(ecosystem_state)
      )

    if (
      nrow(dat) < 10 ||
      n_distinct(dat$ecosystem_state) < 2 ||
      n_distinct(dat$value) < 2
    ) {
      return(NULL)
    }

    tab <- table(
      dat$ecosystem_state,
      dat$value
    )

    pval <- tryCatch({

      fisher.test(
        tab,
        simulate.p.value = TRUE,
        B = 10000
      )$p.value

    }, error = function(e) {

      NA_real_
    })

    tibble(
      variable = v,
      n = nrow(dat),
      n_levels = n_distinct(
        dat$value
      ),
      p_value = pval
    )
  }
)


if (nrow(categorical_results) > 0) {

  categorical_results <- categorical_results %>%
    mutate(
      FDR = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    arrange(
      FDR,
      p_value
    )

} else {

  categorical_results <- tibble(
    variable = character(),
    n = integer(),
    n_levels = integer(),
    p_value = numeric(),
    FDR = numeric()
  )
}


cat(
  "Categorical tests completed:",
  nrow(categorical_results),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(
    categorical_results$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


write_csv(
  categorical_results,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_biological_categorical_associations.csv"
  )
)


# ======================================================================
# 24. STATE-WISE CONTINUOUS BIOLOGICAL PROFILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("STATE-WISE BIOLOGICAL CONTINUOUS PROFILES\n")
cat("=====================================================================\n")


statewise_continuous <- purrr::map_dfr(
  continuous_variables_filtered,
  function(v) {

    biological_data %>%
      select(
        ecosystem_state,
        value = all_of(v)
      ) %>%
      mutate(
        value = suppressWarnings(
          as.numeric(value)
        )
      ) %>%
      filter(
        is.finite(value),
        !is.na(ecosystem_state)
      ) %>%
      group_by(
        ecosystem_state
      ) %>%
      summarise(

        variable = v,

        n = n(),

        mean = mean(
          value,
          na.rm = TRUE
        ),

        median = median(
          value,
          na.rm = TRUE
        ),

        sd = sd(
          value,
          na.rm = TRUE
        ),

        q25 = quantile(
          value,
          0.25,
          na.rm = TRUE
        ),

        q75 = quantile(
          value,
          0.75,
          na.rm = TRUE
        ),

        .groups = "drop"
      ) %>%
      select(
        variable,
        everything()
      )
  }
)


cat(
  "State-wise continuous profiles generated:",
  nrow(statewise_continuous),
  "\n"
)


write_csv(
  statewise_continuous,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_state_wise_biological_continuous_profiles.csv"
  )
)


# ======================================================================
# 25. STATE-WISE CATEGORICAL PROFILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("STATE-WISE BIOLOGICAL CATEGORICAL PROFILES\n")
cat("=====================================================================\n")


statewise_categorical <- purrr::map_dfr(
  categorical_variables_filtered,
  function(v) {

    biological_data %>%
      select(
        ecosystem_state,
        value = all_of(v)
      ) %>%
      mutate(
        value = as.character(value)
      ) %>%
      filter(
        !is.na(value),
        trimws(value) != ""
      ) %>%
      count(
        ecosystem_state,
        value,
        name = "n"
      ) %>%
      group_by(
        ecosystem_state
      ) %>%
      mutate(
        percentage = 100 * n / sum(n),
        variable = v
      ) %>%
      ungroup() %>%
      select(
        variable,
        everything()
      )
  }
)


cat(
  "State-wise categorical profiles generated:",
  nrow(statewise_categorical),
  "\n"
)


write_csv(
  statewise_categorical,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_state_wise_biological_categorical_profiles.csv"
  )
)


# ======================================================================
# 26. BIOLOGICAL HEATMAP
#
# Use the strongest FDR-significant continuous biological features.
#
# Maximum 25 features to keep figure readable.
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("BIOLOGICAL HEATMAP\n")
cat("=====================================================================\n")


heatmap_features <- continuous_effects %>%
  filter(
    FDR < 0.05,
    is.finite(epsilon_squared)
  ) %>%
  arrange(
    desc(epsilon_squared)
  ) %>%
  slice_head(
    n = 25
  ) %>%
  pull(
    variable
  )


if (length(heatmap_features) >= 2) {

  heatmap_data <- biological_data %>%
    select(
      ecosystem_state,
      all_of(heatmap_features)
    )

  # Convert to numeric matrix
  mat <- as.matrix(
    heatmap_data[
      heatmap_features
    ]
  )

  mat <- apply(
    mat,
    2,
    function(x) {

      x <- suppressWarnings(
        as.numeric(x)
      )

      if (all(is.na(x))) {
        return(rep(0, length(x)))
      }

      med <- median(
        x,
        na.rm = TRUE
      )

      x[is.na(x)] <- med

      as.numeric(
        scale(x)
      )
    }
  )

  rownames(mat) <- biological_data$sample_id_clean

  # Order patients by ecosystem state
  ord <- order(
    biological_data$ecosystem_state
  )

  mat <- mat[ord, , drop = FALSE]

  png(
    filename = file.path(
      OUTPUT_DIR,
      "M3.5.1_C2_biological_heatmap.png"
    ),
    width = 1800,
    height = 2200,
    res = 220
  )

  heatmap(
    mat,
    Rowv = NA,
    Colv = NA,
    scale = "none",
    margins = c(
      12,
      8
    ),
    labRow = NA,
    main = "Frozen Ecosystem States — Biological Features"
  )

  dev.off()

  cat(
    "Heatmap generated.\n"
  )

} else {

  cat(
    "Fewer than two significant biological features.\n",
    "Heatmap not generated.\n"
  )
}


# ======================================================================
# 27. PCA OF BIOLOGICAL FEATURE SPACE
#
# IMPORTANT:
# PCA is descriptive only.
# It does NOT create or modify ecosystem states.
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("BIOLOGICAL PCA\n")
cat("=====================================================================\n")


pca_features <- continuous_effects %>%
  filter(
    FDR < 0.05
  ) %>%
  arrange(
    FDR
  ) %>%
  slice_head(
    n = 50
  ) %>%
  pull(
    variable
  )


if (length(pca_features) >= 3) {

  pca_matrix <- biological_data %>%
    select(
      all_of(pca_features)
    ) %>%
    mutate(
      across(
        everything(),
        ~ suppressWarnings(
          as.numeric(.x)
        )
      )
    )

  # Remove zero-variance features
  keep_features <- names(
    pca_matrix
  )[
    vapply(
      pca_matrix,
      function(x) {

        x <- x[
          is.finite(x)
        ]

        length(x) > 1 &&
          sd(x, na.rm = TRUE) > 0
      },
      logical(1)
    )
  ]

  pca_matrix <- pca_matrix[
    keep_features
  ]

  # Median imputation
  for (j in seq_len(ncol(pca_matrix))) {

    x <- pca_matrix[[j]]

    med <- median(
      x[
        is.finite(x)
      ],
      na.rm = TRUE
    )

    x[
      !is.finite(x)
    ] <- med

    pca_matrix[[j]] <- x
  }

  pca <- prcomp(
    pca_matrix,
    center = TRUE,
    scale. = TRUE
  )

  pca_scores <- as.data.frame(
    pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE]
  )

  pca_scores$sample_id_clean <-
    biological_data$sample_id_clean

  pca_scores$ecosystem_state <-
    biological_data$ecosystem_state

  write_csv(
    pca_scores,
    file.path(
      OUTPUT_DIR,
      "M3.5.1_C2_biological_PCA_scores.csv"
    )
  )

  variance <- tibble(
    PC = paste0(
      "PC",
      seq_along(
        pca$sdev
      )
    ),
    variance_explained =
      (pca$sdev^2) /
      sum(
        pca$sdev^2
      ),
    cumulative_variance =
      cumsum(
        (pca$sdev^2) /
        sum(
          pca$sdev^2
        )
      )
  )

  write_csv(
    variance,
    file.path(
      OUTPUT_DIR,
      "M3.5.1_C2_biological_PCA_variance.csv"
    )
  )

  p <- ggplot(
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
    theme_bw() +
    labs(
      title =
        "Biological Feature Space by Frozen Ecosystem State",
      subtitle =
        "PCA is descriptive; ecosystem states were not recalculated",
      x = paste0(
        "PC1 (",
        round(
          100 *
            variance$variance_explained[1],
          1
        ),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(
          100 *
            variance$variance_explained[2],
          1
        ),
        "%)"
      ),
      shape = "Ecosystem state"
    )

  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      "M3.5.1_C2_biological_PCA.png"
    ),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )

  cat(
    "Biological PCA generated.\n"
  )

} else {

  cat(
    "Insufficient significant biological variables for PCA.\n"
  )
}


# ======================================================================
# 28. SAVE MERGED BIOLOGICAL DATASET
# ======================================================================

write_csv(
  biological_data,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_frozen_states_with_biological_features.csv"
  )
)


# ======================================================================
# 29. FINAL SUMMARY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.5.1-C2 FINAL SUMMARY\n")
cat("=====================================================================\n")


n_continuous_sig <- sum(
  continuous_results$FDR < 0.05,
  na.rm = TRUE
)

n_categorical_sig <- sum(
  categorical_results$FDR < 0.05,
  na.rm = TRUE
)


summary_table <- tibble(

  metric = c(

    "Frozen patients",

    "Frozen ecosystem states",

    "Patient-level biological files",

    "Candidate biological variables",

    "Continuous variables tested",

    "Categorical variables tested",

    "Continuous FDR < 0.05",

    "Categorical FDR < 0.05",

    "Pairwise biological tests",

    "Survival analysis",

    "Univariate hazard ratio",

    "Multivariate hazard ratio",

    "Cox regression",

    "State re-clustering"

  ),

  value = c(

    nrow(frozen),

    n_distinct(
      frozen$ecosystem_state
    ),

    length(
      biological_tables
    ),

    length(
      candidate_biological_names
    ),

    length(
      continuous_variables_filtered
    ),

    length(
      categorical_variables_filtered
    ),

    n_continuous_sig,

    n_categorical_sig,

    nrow(
      pairwise_results
    ),

    "NOT PERFORMED",

    "NOT PERFORMED",

    "NOT PERFORMED",

    "NOT PERFORMED",

    "NOT PERFORMED"

  )
)


print(
  summary_table
)


write_csv(
  summary_table,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_overall_summary.csv"
  )
)


# ======================================================================
# 30. WRITE DIAGNOSTIC LOG
# ======================================================================

diagnostic_file <- file.path(
  OUTPUT_DIR,
  "M3.5.1_C2_run_diagnostic.txt"
)


sink(
  diagnostic_file
)


cat(
  "M3.5.1-C2 INDEPENDENT BIOLOGICAL VALIDATION\n"
)

cat(
  "Run time:",
  as.character(Sys.time()),
  "\n\n"
)

cat(
  "Frozen M3.2 file:\n",
  M3_FROZEN_FILE,
  "\n\n"
)

cat(
  "Frozen patients:",
  nrow(frozen),
  "\n"
)

cat(
  "Ecosystem states:",
  n_distinct(
    frozen$ecosystem_state
  ),
  "\n\n"
)

cat(
  "Continuous variables tested:",
  length(
    continuous_variables_filtered
  ),
  "\n"
)

cat(
  "Continuous FDR < 0.05:",
  n_continuous_sig,
  "\n\n"
)

cat(
  "Categorical variables tested:",
  length(
    categorical_variables_filtered
  ),
  "\n"
)

cat(
  "Categorical FDR < 0.05:",
  n_categorical_sig,
  "\n\n"
)

cat(
  "Pairwise tests:",
  nrow(
    pairwise_results
  ),
  "\n\n"
)

cat(
  "SURVIVAL ANALYSIS: NOT PERFORMED\n"
)

cat(
  "UNIVARIATE HAZARD RATIO: NOT PERFORMED\n"
)

cat(
  "MULTIVARIATE HAZARD RATIO: NOT PERFORMED\n"
)

cat(
  "COX REGRESSION: NOT PERFORMED\n"
)

cat(
  "STATE RECLUSTERING: NOT PERFORMED\n"
)

cat(
  "STATE LABELS MODIFIED: NO\n"
)

cat(
  "FROZEN M3.2 STATES MODIFIED: NO\n"
)

sink()


# ======================================================================
# 31. FINAL MESSAGE
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.5.1-C2 COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n")

cat(
  "Frozen patients:",
  nrow(frozen),
  "\n"
)

cat(
  "Ecosystem states:",
  n_distinct(
    frozen$ecosystem_state
  ),
  "\n"
)

cat(
  "Continuous biological variables tested:",
  length(
    continuous_variables_filtered
  ),
  "\n"
)

cat(
  "Continuous FDR < 0.05:",
  n_continuous_sig,
  "\n"
)

cat(
  "Categorical biological variables tested:",
  length(
    categorical_variables_filtered
  ),
  "\n"
)

cat(
  "Categorical FDR < 0.05:",
  n_categorical_sig,
  "\n"
)

cat("\n")
cat(
  "Survival analysis: NOT PERFORMED\n"
)

cat(
  "Univariate hazard ratios: NOT PERFORMED\n"
)

cat(
  "Multivariate hazard ratios: NOT PERFORMED\n"
)

cat(
  "Cox regression: NOT PERFORMED\n"
)

cat(
  "Ecosystem states were NOT redefined.\n"
)

cat(
  "M3.2 frozen states were NOT modified.\n"
)

cat("\n")
cat("Output directory:\n")
cat(
  OUTPUT_DIR,
  "\n"
)

cat("\n")
cat("Important output files:\n\n")

cat(
  "1. Biological continuous associations:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_biological_continuous_associations.csv"
  ),
  "\n\n"
)

cat(
  "2. Biological continuous effect sizes:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_biological_continuous_effect_sizes.csv"
  ),
  "\n\n"
)

cat(
  "3. Top biological features:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_top_biological_features.csv"
  ),
  "\n\n"
)

cat(
  "4. Pairwise biological comparisons:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_pairwise_biological_comparisons.csv"
  ),
  "\n\n"
)

cat(
  "5. Biological categorical associations:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_biological_categorical_associations.csv"
  ),
  "\n\n"
)

cat(
  "6. State-wise continuous profiles:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_state_wise_biological_continuous_profiles.csv"
  ),
  "\n\n"
)

cat(
  "7. State-wise categorical profiles:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_state_wise_biological_categorical_profiles.csv"
  ),
  "\n\n"
)

cat(
  "8. Biological PCA:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_biological_PCA.png"
  ),
  "\n\n"
)

cat(
  "9. Biological heatmap:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_biological_heatmap.png"
  ),
  "\n\n"
)

cat(
  "10. Overall summary:\n",
  file.path(
    OUTPUT_DIR,
    "M3.5.1_C2_overall_summary.csv"
  ),
  "\n\n"
)

cat(
  "11. Diagnostic log:\n",
  diagnostic_file,
  "\n"
)

cat("\n")
cat("=====================================================================\n")
cat("M3.5.1-C2 FINISHED\n")
cat("=====================================================================\n")