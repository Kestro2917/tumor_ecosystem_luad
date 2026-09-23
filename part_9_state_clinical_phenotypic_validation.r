# ======================================================================
# M3.5.1 — FROZEN ECOSYSTEM STATE BIOLOGICAL CHARACTERIZATION
#
# PURPOSE
# -------
# Characterize the already-frozen M3.2 ecosystem states using
# molecular / ecosystem-axis variables generated upstream.
#
# IMPORTANT:
#   1. Ecosystem states are NOT re-clustered.
#   2. Ecosystem states are NOT redefined.
#   3. No molecular feature is used to create the states.
#   4. This is a characterization / validation analysis only.
#
# ACTUAL M3.2 FROZEN STATE FILE:
#   /content/M3_LIONESS_entropy/M3.2_module_axes/
#       M3.2_patient_ecosystem_states.csv
#
# EXPECTED FROZEN STATE COLUMNS:
#   sample_id
#   ecosystem_state
#   cluster_numeric
#
# OUTPUT:
#   /content/M3_LIONESS_entropy/M3.5_clinical_validation/M3.5.1_LayerC
# ======================================================================


rm(list = ls())

options(stringsAsFactors = FALSE)
options(warn = 1)

start_time <- Sys.time()

cat("\n")
cat("=====================================================================\n")
cat("M3.5.1 — FROZEN ECOSYSTEM STATE BIOLOGICAL CHARACTERIZATION\n")
cat("=====================================================================\n")
cat("Started:", as.character(start_time), "\n")


# ======================================================================
# 1. PACKAGES
# ======================================================================

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "stringr"
)

cat("\n")
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
  "M3.5.1_LayerC"
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

cat("M3.2 frozen state file:\n")
cat(M3_FROZEN_FILE, "\n\n")

cat("M3.5.1 output directory:\n")
cat(OUTPUT_DIR, "\n")


# ======================================================================
# 3. CHECK FROZEN STATE FILE
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("CHECKING FROZEN M3.2 FILE\n")
cat("=====================================================================\n")

if (!file.exists(M3_FROZEN_FILE)) {

  stop(
    paste0(
      "\nERROR: M3.2 frozen state file was not found:\n",
      M3_FROZEN_FILE
    )
  )
}

cat("Frozen M3.2 file found.\n")


# ======================================================================
# 4. READ FROZEN STATES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING FROZEN M3.2 STATES\n")
cat("=====================================================================\n")

frozen <- data.table::fread(
  M3_FROZEN_FILE,
  data.table = FALSE,
  check.names = FALSE
)

cat("Rows:", nrow(frozen), "\n")
cat("Columns:", ncol(frozen), "\n")

cat("\nM3.2 columns:\n")
print(colnames(frozen))


# ======================================================================
# 5. IDENTIFY REQUIRED FROZEN COLUMNS
# ======================================================================

required_frozen <- c(
  "sample_id",
  "ecosystem_state"
)

missing_frozen <- setdiff(
  required_frozen,
  colnames(frozen)
)

if (length(missing_frozen) > 0) {

  stop(
    paste0(
      "\nERROR: Required frozen columns are missing:\n",
      paste(missing_frozen, collapse = ", ")
    )
  )
}


# ======================================================================
# 6. CLEAN FROZEN STATE TABLE
# ======================================================================

frozen <- frozen %>%
  dplyr::select(
    sample_id,
    ecosystem_state,
    dplyr::everything()
  ) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    ecosystem_state = as.character(ecosystem_state)
  )

frozen <- frozen %>%
  dplyr::filter(
    !is.na(sample_id),
    sample_id != "",
    !is.na(ecosystem_state),
    ecosystem_state != ""
  )

cat("\n")
cat("Valid frozen patients:", nrow(frozen), "\n")


# ======================================================================
# 7. CHECK DUPLICATES
# ======================================================================

duplicate_ids <- frozen %>%
  dplyr::count(sample_id) %>%
  dplyr::filter(n > 1)

if (nrow(duplicate_ids) > 0) {

  cat("\nWARNING: Duplicate frozen patient IDs detected.\n")
  print(duplicate_ids)

  stop(
    "M3.5.1 cannot continue because frozen sample_id values are duplicated."
  )

} else {

  cat("No duplicate frozen patient IDs detected.\n")
}


# ======================================================================
# 8. FROZEN STATE DISTRIBUTION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("FROZEN ECOSYSTEM STATE DISTRIBUTION\n")
cat("=====================================================================\n")

state_distribution <- frozen %>%
  dplyr::count(ecosystem_state, name = "n") %>%
  dplyr::mutate(
    percentage = 100 * n / sum(n)
  )

print(state_distribution)

write.csv(
  state_distribution,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_frozen_state_distribution.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 9. SEARCH M3.2 OUTPUT DIRECTORY FOR MOLECULAR / ECOSYSTEM FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("SEARCHING M3.2 OUTPUT DIRECTORY FOR BIOLOGICAL FEATURES\n")
cat("=====================================================================\n")

m3_files <- list.files(
  M3_ROOT,
  recursive = TRUE,
  full.names = TRUE
)

m3_csv_files <- m3_files[
  grepl(
    "\\.csv$",
    m3_files,
    ignore.case = TRUE
  )
]

cat("CSV files found under M3 directory:", length(m3_csv_files), "\n")


# ======================================================================
# 10. EXCLUDE KNOWN NON-BIOLOGICAL FILES
# ======================================================================

excluded_patterns <- c(
  "patient_ecosystem_states",
  "clinical_validation",
  "M3.5",
  "diagnostic",
  "summary",
  "state_distribution",
  "survival",
  "matching"
)

is_excluded_file <- function(x) {

  any(
    stringr::str_detect(
      basename(x),
      stringr::regex(
        paste(excluded_patterns, collapse = "|"),
        ignore_case = TRUE
      )
    )
  )
}

candidate_files <- m3_csv_files[
  !sapply(m3_csv_files, is_excluded_file)
]

cat(
  "Candidate biological CSV files:",
  length(candidate_files),
  "\n"
)


# ======================================================================
# 11. IDENTIFY PATIENT-LEVEL BIOLOGICAL FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("INSPECTING CANDIDATE BIOLOGICAL FILES\n")
cat("=====================================================================\n")

file_inventory <- data.frame(
  file = character(),
  n_rows = integer(),
  n_columns = integer(),
  has_sample_id = logical(),
  has_common_patient_id = logical(),
  stringsAsFactors = FALSE
)

for (f in candidate_files) {

  tmp <- tryCatch(

    data.table::fread(
      f,
      data.table = FALSE,
      check.names = FALSE,
      nrows = 5
    ),

    error = function(e) NULL
  )

  if (is.null(tmp)) {
    next
  }

  cn <- colnames(tmp)

  possible_id <- any(
    cn %in% c(
      "sample_id",
      "sample",
      "patient_id",
      "patient",
      "barcode",
      "sample_barcode"
    )
  )

  file_inventory <- rbind(
    file_inventory,
    data.frame(
      file = f,
      n_rows = NA_integer_,
      n_columns = ncol(tmp),
      has_sample_id = "sample_id" %in% cn,
      has_common_patient_id = possible_id,
      stringsAsFactors = FALSE
    )
  )
}

print(file_inventory)

write.csv(
  file_inventory,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_biological_file_inventory.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 12. PRIORITIZE FILES THAT CONTAIN sample_id
# ======================================================================

biological_files <- file_inventory %>%
  dplyr::filter(
    has_sample_id
  ) %>%
  dplyr::pull(file)

cat(
  "\nPatient-level files containing sample_id:",
  length(biological_files),
  "\n"
)

if (length(biological_files) > 0) {

  print(biological_files)

}


# ======================================================================
# 13. IF NO BIOLOGICAL FILE HAS sample_id
# ======================================================================

if (length(biological_files) == 0) {

  cat("\n")
  cat("=====================================================================\n")
  cat("NO PATIENT-LEVEL BIOLOGICAL FILE FOUND\n")
  cat("=====================================================================\n")

  cat(
    "\nThe frozen M3.2 patient-state file contains only:\n"
  )

  print(colnames(frozen))

  cat(
    "\nTherefore M3.5.1 cannot perform molecular association testing\n",
    "without an upstream M3.2 feature/axis file.\n",
    sep = ""
  )

  no_feature_summary <- data.frame(
    status = "NO_BIOLOGICAL_FEATURE_FILE_FOUND",
    frozen_patients = nrow(frozen),
    ecosystem_states = length(unique(frozen$ecosystem_state)),
    molecular_files_found = 0,
    message = paste(
      "M3.2 frozen state file contains no molecular/ecosystem-axis",
      "variables. M3.5.1 biological association testing was safely skipped."
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    no_feature_summary,
    file.path(
      OUTPUT_DIR,
      "M3.5.1_no_biological_features_summary.csv"
    ),
    row.names = FALSE
  )

  writeLines(
    c(
      "M3.5.1 STATUS",
      "==============================",
      "No patient-level molecular/ecosystem feature file was found.",
      "",
      paste(
        "Frozen patients:",
        nrow(frozen)
      ),
      paste(
        "Frozen ecosystem states:",
        length(unique(frozen$ecosystem_state))
      ),
      "",
      "Biological association testing was skipped safely.",
      "",
      "The M3.2 frozen states were NOT modified."
    ),
    file.path(
      OUTPUT_DIR,
      "M3.5.1_diagnostic.txt"
    )
  )

  cat("\n")
  cat("M3.5.1 completed safely with no molecular analysis.\n")
  cat("No frozen states were changed.\n")

  quit(
    save = "no",
    status = 0
  )
}


# ======================================================================
# 14. READ AND MERGE BIOLOGICAL FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING PATIENT-LEVEL BIOLOGICAL FEATURES\n")
cat("=====================================================================\n")

biological_tables <- list()

for (f in biological_files) {

  cat("\nReading:\n", f, "\n")

  tmp <- tryCatch(

    data.table::fread(
      f,
      data.table = FALSE,
      check.names = FALSE
    ),

    error = function(e) {

      cat(
        "Could not read:",
        basename(f),
        "\n"
      )

      NULL
    }
  )

  if (is.null(tmp)) {
    next
  }

  if (!"sample_id" %in% colnames(tmp)) {
    next
  }

  tmp$sample_id <- as.character(
    tmp$sample_id
  )

  tmp <- tmp %>%
    dplyr::filter(
      !is.na(sample_id),
      sample_id != ""
    )

  if (nrow(tmp) == 0) {
    next
  }

  biological_tables[[basename(f)]] <- tmp
}


# ======================================================================
# 15. MERGE BIOLOGICAL FEATURES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("MERGING BIOLOGICAL FEATURES WITH FROZEN STATES\n")
cat("=====================================================================\n")

bio_merged <- frozen

if (length(biological_tables) > 0) {

  for (nm in names(biological_tables)) {

    tmp <- biological_tables[[nm]]

    # Keep only unique patient-level rows
    tmp <- tmp %>%
      dplyr::group_by(sample_id) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    # Prevent duplicate column names
    common_cols <- intersect(
      colnames(bio_merged),
      colnames(tmp)
    )

    common_cols <- setdiff(
      common_cols,
      "sample_id"
    )

    if (length(common_cols) > 0) {

      tmp <- tmp %>%
        dplyr::select(
          -dplyr::all_of(common_cols)
        )
    }

    bio_merged <- bio_merged %>%
      dplyr::left_join(
        tmp,
        by = "sample_id"
      )
  }
}

cat(
  "Merged biological dataset dimensions:\n"
)

cat(
  "Rows:",
  nrow(bio_merged),
  "\n"
)

cat(
  "Columns:",
  ncol(bio_merged),
  "\n"
)


# ======================================================================
# 16. CHECK MATCHING
# ======================================================================

biological_match <- bio_merged %>%
  dplyr::filter(
    !is.na(sample_id)
  )

cat(
  "\nFrozen patients:",
  nrow(frozen),
  "\n"
)

cat(
  "Patients retained after biological merge:",
  nrow(biological_match),
  "\n"
)


# ======================================================================
# 17. IDENTIFY MOLECULAR / ECOSYSTEM VARIABLES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("IDENTIFYING BIOLOGICAL VARIABLES\n")
cat("=====================================================================\n")

protected_columns <- c(
  "sample_id",
  "ecosystem_state",
  "cluster_numeric"
)

candidate_bio_columns <- setdiff(
  colnames(bio_merged),
  protected_columns
)

cat(
  "Candidate biological variables:",
  length(candidate_bio_columns),
  "\n"
)


# ======================================================================
# 18. CONVERT POSSIBLE NUMERIC VARIABLES
# ======================================================================

bio_numeric <- character()

for (v in candidate_bio_columns) {

  x <- bio_merged[[v]]

  if (is.numeric(x)) {

    bio_numeric <- c(
      bio_numeric,
      v
    )

  } else {

    x_num <- suppressWarnings(
      as.numeric(as.character(x))
    )

    non_missing_original <- sum(
      !is.na(x)
    )

    non_missing_numeric <- sum(
      !is.na(x_num)
    )

    if (
      non_missing_original > 0 &&
      non_missing_numeric >=
      0.8 * non_missing_original
    ) {

      bio_merged[[v]] <- x_num

      bio_numeric <- c(
        bio_numeric,
        v
      )
    }
  }
}

bio_numeric <- unique(
  bio_numeric
)

bio_categorical <- setdiff(
  candidate_bio_columns,
  bio_numeric
)


cat(
  "Biological continuous variables:",
  length(bio_numeric),
  "\n"
)

cat(
  "Biological categorical variables:",
  length(bio_categorical),
  "\n"
)


# ======================================================================
# 19. REMOVE VARIABLES WITH VERY LOW INFORMATION
# ======================================================================

valid_numeric <- character()

for (v in bio_numeric) {

  x <- bio_merged[[v]]

  x <- x[
    is.finite(x)
  ]

  if (length(x) < 5) {
    next
  }

  if (length(unique(x)) < 2) {
    next
  }

  valid_numeric <- c(
    valid_numeric,
    v
  )
}

bio_numeric <- valid_numeric


valid_categorical <- character()

for (v in bio_categorical) {

  x <- as.character(
    bio_merged[[v]]
  )

  x <- x[
    !is.na(x) &
      x != ""
  ]

  if (length(x) < 5) {
    next
  }

  if (length(unique(x)) < 2) {
    next
  }

  valid_categorical <- c(
    valid_categorical,
    v
  )
}

bio_categorical <- valid_categorical


cat(
  "\nAfter information filtering:\n"
)

cat(
  "Continuous variables:",
  length(bio_numeric),
  "\n"
)

cat(
  "Categorical variables:",
  length(bio_categorical),
  "\n"
)


# ======================================================================
# 20. BIOLOGICAL CONTINUOUS ASSOCIATION TESTS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C — BIOLOGICAL CONTINUOUS ASSOCIATIONS\n")
cat("=====================================================================\n")

continuous_results <- data.frame(
  variable = character(),
  n = integer(),
  n_states = integer(),
  statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (v in bio_numeric) {

  tmp <- data.frame(
    value = bio_merged[[v]],
    ecosystem_state =
      bio_merged$ecosystem_state
  )

  tmp <- tmp[
    is.finite(tmp$value) &
      !is.na(tmp$ecosystem_state),
    ,
    drop = FALSE
  ]

  if (nrow(tmp) < 10) {
    next
  }

  n_states <- length(
    unique(
      tmp$ecosystem_state
    )
  )

  if (n_states < 2) {
    next
  }

  test_result <- tryCatch(

    kruskal.test(
      value ~ ecosystem_state,
      data = tmp
    ),

    error = function(e) NULL
  )

  if (is.null(test_result)) {
    next
  }

  continuous_results <- rbind(
    continuous_results,
    data.frame(
      variable = v,
      n = nrow(tmp),
      n_states = n_states,
      statistic = unname(
        test_result$statistic
      ),
      p_value = test_result$p.value,
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
      continuous_results$p_value
    ),
    ,
    drop = FALSE
  ]
}

cat(
  "Continuous tests completed:",
  nrow(continuous_results),
  "\n"
)

if (nrow(continuous_results) > 0) {

  cat(
    "FDR < 0.05:",
    sum(
      continuous_results$FDR < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )

  print(
    head(
      continuous_results,
      20
    )
  )
}


write.csv(
  continuous_results,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_LayerC_biological_continuous_associations.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 21. BIOLOGICAL CATEGORICAL ASSOCIATION TESTS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C — BIOLOGICAL CATEGORICAL ASSOCIATIONS\n")
cat("=====================================================================\n")

categorical_results <- data.frame(
  variable = character(),
  n = integer(),
  n_levels = integer(),
  test = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)


for (v in bio_categorical) {

  tmp <- data.frame(
    value = as.character(
      bio_merged[[v]]
    ),
    ecosystem_state =
      bio_merged$ecosystem_state
  )

  tmp <- tmp[
    !is.na(tmp$value) &
      tmp$value != "" &
      !is.na(tmp$ecosystem_state),
    ,
    drop = FALSE
  ]

  if (nrow(tmp) < 10) {
    next
  }

  tmp$value <- droplevels(
    factor(tmp$value)
  )

  tmp$ecosystem_state <- droplevels(
    factor(tmp$ecosystem_state)
  )

  if (nlevels(tmp$value) < 2) {
    next
  }

  if (nlevels(tmp$ecosystem_state) < 2) {
    next
  }

  tab <- table(
    tmp$ecosystem_state,
    tmp$value
  )

  test_p <- NA_real_
  test_name <- NA_character_

  test_result <- tryCatch(

    fisher.test(
      tab,
      simulate.p.value = TRUE,
      B = 10000
    ),

    error = function(e) NULL
  )

  if (!is.null(test_result)) {

    test_p <- test_result$p.value
    test_name <- "Fisher_exact_simulated"

  } else {

    chi_result <- tryCatch(

      suppressWarnings(
        chisq.test(
          tab
        )
      ),

      error = function(e) NULL
    )

    if (!is.null(chi_result)) {

      test_p <- chi_result$p.value
      test_name <- "Chi_square"
    }
  }

  if (!is.na(test_p)) {

    categorical_results <- rbind(
      categorical_results,
      data.frame(
        variable = v,
        n = nrow(tmp),
        n_levels = nlevels(tmp$value),
        test = test_name,
        p_value = test_p,
        stringsAsFactors = FALSE
      )
    )
  }
}


if (nrow(categorical_results) > 0) {

  categorical_results$FDR <- p.adjust(
    categorical_results$p_value,
    method = "BH"
  )

  categorical_results <- categorical_results[
    order(
      categorical_results$FDR,
      categorical_results$p_value
    ),
    ,
    drop = FALSE
  ]
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

  print(
    head(
      categorical_results,
      20
    )
  )
}


write.csv(
  categorical_results,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_LayerC_biological_categorical_associations.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 22. STATE-WISE CONTINUOUS PROFILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C — STATE-WISE BIOLOGICAL CONTINUOUS PROFILES\n")
cat("=====================================================================\n")

continuous_profiles <- data.frame()

if (length(bio_numeric) > 0) {

  profile_list <- list()

  for (v in bio_numeric) {

    tmp <- data.frame(
      value = bio_merged[[v]],
      ecosystem_state =
        bio_merged$ecosystem_state
    )

    tmp <- tmp[
      is.finite(tmp$value) &
        !is.na(tmp$ecosystem_state),
      ,
      drop = FALSE
    ]

    if (nrow(tmp) == 0) {
      next
    }

    profile <- tmp %>%
      dplyr::group_by(
        ecosystem_state
      ) %>%
      dplyr::summarise(

        variable = v,

        n = sum(
          is.finite(value)
        ),

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
          na.rm = TRUE,
          names = FALSE
        ),

        q75 = quantile(
          value,
          0.75,
          na.rm = TRUE,
          names = FALSE
        ),

        min = ifelse(
          any(is.finite(value)),
          min(
            value,
            na.rm = TRUE
          ),
          NA_real_
        ),

        max = ifelse(
          any(is.finite(value)),
          max(
            value,
            na.rm = TRUE
          ),
          NA_real_
        ),

        .groups = "drop"
      )

    profile_list[[length(profile_list) + 1]] <- profile
  }

  if (length(profile_list) > 0) {

    continuous_profiles <- dplyr::bind_rows(
      profile_list
    )
  }
}


cat(
  "Continuous state profiles generated:",
  nrow(continuous_profiles),
  "\n"
)

write.csv(
  continuous_profiles,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_LayerC_biological_continuous_state_profiles.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 23. STATE-WISE CATEGORICAL PROFILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C — STATE-WISE BIOLOGICAL CATEGORICAL PROFILES\n")
cat("=====================================================================\n")

categorical_profiles <- data.frame()

if (length(bio_categorical) > 0) {

  cat_profile_list <- list()

  for (v in bio_categorical) {

    tmp <- data.frame(
      value = as.character(
        bio_merged[[v]]
      ),
      ecosystem_state =
        bio_merged$ecosystem_state
    )

    tmp <- tmp[
      !is.na(tmp$value) &
        tmp$value != "" &
        !is.na(tmp$ecosystem_state),
      ,
      drop = FALSE
    ]

    if (nrow(tmp) == 0) {
      next
    }

    profile <- tmp %>%
      dplyr::count(
        ecosystem_state,
        value,
        name = "n"
      ) %>%
      dplyr::group_by(
        ecosystem_state
      ) %>%
      dplyr::mutate(
        percentage =
          100 * n / sum(n),
        variable = v
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        variable,
        ecosystem_state,
        value,
        n,
        percentage
      )

    cat_profile_list[[length(cat_profile_list) + 1]] <-
      profile
  }

  if (length(cat_profile_list) > 0) {

    categorical_profiles <- dplyr::bind_rows(
      cat_profile_list
    )
  }
}


cat(
  "Categorical state profiles generated:",
  nrow(categorical_profiles),
  "\n"
)

write.csv(
  categorical_profiles,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_LayerC_biological_categorical_state_profiles.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 24. SIGNIFICANT CONTINUOUS PLOTS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C — SIGNIFICANT CONTINUOUS VARIABLE PLOTS\n")
cat("=====================================================================\n")

significant_continuous <- character()

if (nrow(continuous_results) > 0) {

  significant_continuous <-
    continuous_results %>%
    dplyr::filter(
      FDR < 0.05
    ) %>%
    dplyr::pull(
      variable
    )
}

if (length(significant_continuous) == 0) {

  cat(
    "No FDR-significant biological continuous variables.\n"
  )

} else {

  for (v in significant_continuous) {

    plot_data <- data.frame(
      value = bio_merged[[v]],
      ecosystem_state =
        bio_merged$ecosystem_state
    )

    plot_data <- plot_data[
      is.finite(plot_data$value) &
        !is.na(plot_data$ecosystem_state),
      ,
      drop = FALSE
    ]

    if (nrow(plot_data) < 5) {
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
      labs(
        title = paste(
          "M3.5.1 Layer C:",
          v
        ),
        x = "Frozen ecosystem state",
        y = v
      ) +
      theme_bw()

    safe_name <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      v
    )

    ggsave(
      filename = file.path(
        OUTPUT_DIR,
        paste0(
          "M3.5.1_LayerC_",
          safe_name,
          "_boxplot.png"
        )
      ),
      plot = p,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
}


# ======================================================================
# 25. BIOLOGICAL FEATURE HEATMAP DATA
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("LAYER C — BIOLOGICAL FEATURE STATE SUMMARY\n")
cat("=====================================================================\n")

if (length(bio_numeric) > 0) {

  heatmap_summary <- bio_merged %>%
    dplyr::group_by(
      ecosystem_state
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(bio_numeric),
        ~ ifelse(
          all(is.na(.x)),
          NA_real_,
          mean(
            .x,
            na.rm = TRUE
          )
        )
      ),
      .groups = "drop"
    )

  write.csv(
    heatmap_summary,
    file.path(
      OUTPUT_DIR,
      "M3.5.1_LayerC_biological_state_mean_matrix.csv"
    ),
    row.names = FALSE
  )

} else {

  cat(
    "No continuous biological variables available for feature matrix.\n"
  )
}


# ======================================================================
# 26. STANDARDIZED STATE-WISE BIOLOGICAL MATRIX
# ======================================================================

if (length(bio_numeric) > 0) {

  standardized_matrix <- bio_merged %>%
    dplyr::select(
      ecosystem_state,
      dplyr::all_of(bio_numeric)
    ) %>%
    dplyr::group_by(
      ecosystem_state
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(bio_numeric),
        ~ ifelse(
          all(is.na(.x)),
          NA_real_,
          mean(
            .x,
            na.rm = TRUE
          )
        )
      ),
      .groups = "drop"
    )

  numeric_matrix <- as.matrix(
    standardized_matrix[
      ,
      bio_numeric,
      drop = FALSE
    ]
  )

  for (j in seq_len(ncol(numeric_matrix))) {

    column_values <- numeric_matrix[, j]

    if (
      all(
        is.na(column_values)
      )
    ) {
      next
    }

    med <- median(
      column_values,
      na.rm = TRUE
    )

    if (!is.finite(med)) {
      next
    }

    column_values[
      is.na(column_values)
    ] <- med

    numeric_matrix[, j] <- column_values
  }

  if (
    nrow(numeric_matrix) >= 2 &&
    ncol(numeric_matrix) >= 2
  ) {

    scaled_matrix <- scale(
      numeric_matrix
    )

    scaled_matrix[
      !is.finite(scaled_matrix)
    ] <- 0

    standardized_output <- data.frame(
      ecosystem_state =
        standardized_matrix$ecosystem_state,
      as.data.frame(
        scaled_matrix,
        check.names = FALSE
      ),
      check.names = FALSE
    )

    write.csv(
      standardized_output,
      file.path(
        OUTPUT_DIR,
        "M3.5.1_LayerC_standardized_state_biological_matrix.csv"
      ),
      row.names = FALSE
    )
  }
}


# ======================================================================
# 27. SAVE FINAL BIOLOGICAL DATASET
# ======================================================================

write.csv(
  bio_merged,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_LayerC_frozen_states_with_biological_features.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 28. SIGNIFICANCE SUMMARY
# ======================================================================

n_cont_sig <- 0

n_cat_sig <- 0

if (nrow(continuous_results) > 0) {

  n_cont_sig <- sum(
    continuous_results$FDR < 0.05,
    na.rm = TRUE
  )
}

if (nrow(categorical_results) > 0) {

  n_cat_sig <- sum(
    categorical_results$FDR < 0.05,
    na.rm = TRUE
  )
}


# ======================================================================
# 29. OVERALL M3.5.1 SUMMARY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3.5.1 OVERALL SUMMARY\n")
cat("=====================================================================\n")

summary_table <- data.frame(

  metric = c(
    "M3.2_frozen_state_file",
    "n_frozen_patients",
    "n_ecosystem_states",
    "biological_files_found",
    "biological_continuous_variables",
    "biological_categorical_variables",
    "LayerC_continuous_tests",
    "LayerC_categorical_tests",
    "LayerC_continuous_FDR_lt_0.05",
    "LayerC_categorical_FDR_lt_0.05"
  ),

  value = c(
    M3_FROZEN_FILE,
    nrow(frozen),
    length(
      unique(
        frozen$ecosystem_state
      )
    ),
    length(biological_files),
    length(bio_numeric),
    length(bio_categorical),
    nrow(continuous_results),
    nrow(categorical_results),
    n_cont_sig,
    n_cat_sig
  ),

  stringsAsFactors = FALSE
)

print(summary_table)

write.csv(
  summary_table,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_LayerC_overall_summary.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 30. DIAGNOSTIC FILE
# ======================================================================

diagnostic_lines <- c(

  "M3.5.1 — LAYER C BIOLOGICAL CHARACTERIZATION",

  "==============================================================",

  paste(
    "Started:",
    as.character(start_time)
  ),

  paste(
    "Finished:",
    as.character(Sys.time())
  ),

  "",

  paste(
    "Frozen M3.2 file:",
    M3_FROZEN_FILE
  ),

  paste(
    "Frozen patients:",
    nrow(frozen)
  ),

  paste(
    "Frozen ecosystem states:",
    length(
      unique(
        frozen$ecosystem_state
      )
    )
  ),

  "",

  paste(
    "Biological files found:",
    length(biological_files)
  ),

  paste(
    "Biological continuous variables:",
    length(bio_numeric)
  ),

  paste(
    "Biological categorical variables:",
    length(bio_categorical)
  ),

  "",

  paste(
    "Continuous tests:",
    nrow(continuous_results)
  ),

  paste(
    "Categorical tests:",
    nrow(categorical_results)
  ),

  paste(
    "Continuous FDR < 0.05:",
    n_cont_sig
  ),

  paste(
    "Categorical FDR < 0.05:",
    n_cat_sig
  ),

  "",

  "IMPORTANT:",

  "Frozen ecosystem states were NOT redefined.",

  "No clustering was performed.",

  "No molecular variable was used to construct the frozen states.",

  "Layer C is characterization/validation only."
)

writeLines(
  diagnostic_lines,
  file.path(
    OUTPUT_DIR,
    "M3.5.1_LayerC_diagnostic.txt"
  )
)


# ======================================================================
# 31. FINAL MESSAGE
# ======================================================================

elapsed <- difftime(
  Sys.time(),
  start_time,
  units = "mins"
)

cat("\n")
cat("=====================================================================\n")
cat("M3.5.1 COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n")

cat(
  "Frozen patients:",
  nrow(frozen),
  "\n"
)

cat(
  "Ecosystem states:",
  length(
    unique(
      frozen$ecosystem_state
    )
  ),
  "\n"
)

cat(
  "Biological continuous variables:",
  length(bio_numeric),
  "\n"
)

cat(
  "Biological categorical variables:",
  length(bio_categorical),
  "\n"
)

cat(
  "Continuous FDR < 0.05:",
  n_cont_sig,
  "\n"
)

cat(
  "Categorical FDR < 0.05:",
  n_cat_sig,
  "\n"
)

cat(
  "\nOutput directory:\n",
  OUTPUT_DIR,
  "\n"
)

cat(
  "\nElapsed time:",
  round(
    as.numeric(elapsed),
    2
  ),
  "minutes\n"
)

cat("\n")
cat("IMPORTANT:\n")
cat(
  "M3.2 ecosystem states remain completely frozen.\n"
)

cat(
  "M3.5.1 only characterizes biological/ecosystem features.\n"
)

cat("=====================================================================\n")