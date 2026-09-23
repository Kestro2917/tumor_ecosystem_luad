# ======================================================================
# M3.7 (STREAMLINED) — PRESPECIFIED-ONLY COX PROPORTIONAL HAZARDS
# PURPOSE
# -------
# Auditable prognostic analysis of the FROZEN M3.2 ecosystem states
# using the ORIGINAL patient-level M1 clinical dataset, restricted
# to EXACTLY six prespecified clinical variables:
#
#   1. ecosystem_state   (frozen M3.2 assignment; State_1 = reference)
#   2. age_numeric
#   3. sex_factor        (from sex_at_birth, cleaned to Male/Female)
#   4. stage_factor
#   5. race_factor
#   6. smoking_factor
#
# ======================================================================
#
# ======================================================================
# ======================================================================


# ======================================================================
# SECTION 1 — PACKAGE SETUP
# ======================================================================

required_packages <- c(
  "survival",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "purrr",
  "tibble",
  "ggplot2",
  "broom",
  "forcats"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {

  if (!pkg %in% installed) {

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({

  library(survival)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(broom)
  library(forcats)

})


# ======================================================================
# SECTION 2 — FILE PATHS
# ======================================================================

M32_FILE <- "M3.2_patient_ecosystem_states.csv"
M1_FILE  <- "M1_clinical.csv"

OUT_DIR <- "M3.7_hazard_ratio_analysis"


# ======================================================================
# SECTION 3 — FILE RESOLUTION
# ======================================================================

resolve_file <- function(filename) {

  candidates <- unique(
    c(
      filename,
      file.path("/content", filename)
    )
  )

  hits <- candidates[file.exists(candidates)]

  if (length(hits) == 0) {

    stop(
      paste0(
        "\n============================================================\n",
        "FILE NOT FOUND\n",
        "============================================================\n\n",
        "Could not find:\n",
        filename,
        "\n\nChecked:\n",
        paste(candidates, collapse = "\n"),
        "\n\nUpload the file to Google Colab and rerun."
      )
    )
  }

  hits[1]
}


M32_FILE <- resolve_file(M32_FILE)
M1_FILE  <- resolve_file(M1_FILE)


if (!dir.exists(OUT_DIR)) {

  dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


cat("\n")
cat("============================================================\n")
cat("M3.7 PRESPECIFIED-ONLY COX ANALYSIS\n")
cat("============================================================\n")
cat("M3.2 input :", M32_FILE, "\n")
cat("M1 input   :", M1_FILE, "\n")
cat("Output dir :", OUT_DIR, "\n")
cat("============================================================\n")


# ======================================================================
# SECTION 4 — HELPER FUNCTIONS
# ======================================================================

safe_chr <- function(x) {

  x <- as.character(x)
  x <- trimws(x)

  x[x %in% c(
    "",
    "NA",
    "N/A",
    "NULL",
    "null",
    "NaN"
  )] <- NA_character_

  x
}


safe_numeric <- function(x) {

  if (is.numeric(x)) {

    return(as.numeric(x))
  }

  x <- as.character(x)
  x <- trimws(x)

  x[x %in% c(
    "",
    "NA",
    "N/A",
    "NULL",
    "null",
    "NaN",
    "Inf",
    "-Inf"
  )] <- NA_character_

  x <- gsub(",", "", x)

  suppressWarnings(
    as.numeric(x)
  )
}


find_first_column <- function(data, candidates) {

  hits <- candidates[
    candidates %in% names(data)
  ]

  if (length(hits) == 0) {

    return(NA_character_)
  }

  hits[1]
}


clean_factor <- function(x) {

  x <- safe_chr(x)

  factor(x)
}


# ------------------------------------------------------------------
# Dedicated sex/gender cleaner (see prior discussion). Collapses
# "male"/"MALE"/"Male " etc. onto a clean two-level factor with
# "Male" as the reference level. Anything not recognizably
# male/female becomes NA rather than a spurious third level.
# ------------------------------------------------------------------

clean_sex_factor <- function(x) {

  x <- safe_chr(x)

  x <- tolower(trimws(x))

  x <- dplyr::case_when(

    x %in% c("male", "m")   ~ "Male",

    x %in% c("female", "f") ~ "Female",

    TRUE                    ~ NA_character_
  )

  factor(
    x,
    levels = c("Male", "Female")
  )
}


make_empty_cox_result <- function(
    variable_label,
    analysis_label,
    status_value = "FAILED",
    warning_message = NA_character_,
    n_value = NA_integer_,
    events_value = NA_integer_) {

  tibble(

    variable = variable_label,

    term_original = NA_character_,

    analysis = analysis_label,

    HR = NA_real_,

    lower_95_CI = NA_real_,

    upper_95_CI = NA_real_,

    coefficient = NA_real_,

    standard_error = NA_real_,

    z = NA_real_,

    p_value = NA_real_,

    n = n_value,

    events = events_value,

    status = status_value,

    warning_message = warning_message
  )
}


# ======================================================================
# SECTION 5 — EXACT M3.2 → M1 TCGA PATIENT ID NORMALIZATION
# ======================================================================

normalize_m32_to_patient <- function(x) {

  x <- safe_chr(x)

  x <- toupper(x)

  x <- sub(
    "_LIONESS$",
    "",
    x,
    ignore.case = TRUE
  )

  x <- gsub(
    "[._:/]",
    "-",
    x
  )

  x <- gsub(
    "-+",
    "-",
    x
  )

  x <- trimws(x)

  pieces <- strsplit(
    x,
    "-",
    fixed = TRUE
  )

  patient_id <- vapply(

    pieces,

    function(z) {

      if (
        length(z) >= 3 &&
        toupper(z[1]) == "TCGA"
      ) {

        paste(
          z[1:3],
          collapse = "-"
        )

      } else {

        NA_character_
      }
    },

    character(1)
  )

  toupper(patient_id)
}


normalize_m1_patient <- function(x) {

  x <- safe_chr(x)

  x <- toupper(x)

  x <- gsub(
    "[._:/]",
    "-",
    x
  )

  x <- gsub(
    "-+",
    "-",
    x
  )

  x <- trimws(x)

  pieces <- strsplit(
    x,
    "-",
    fixed = TRUE
  )

  patient_id <- vapply(

    pieces,

    function(z) {

      if (
        length(z) >= 3 &&
        toupper(z[1]) == "TCGA"
      ) {

        paste(
          z[1:3],
          collapse = "-"
        )

      } else if (
        length(z) == 1 &&
        grepl(
          "^TCGA[0-9A-Z]+$",
          z[1]
        )
      ) {

        z[1]

      } else {

        NA_character_
      }
    },

    character(1)
  )

  toupper(patient_id)
}


# ======================================================================
# SECTION 6 — ROBUST COX MODEL FUNCTION
# ======================================================================

run_cox <- function(
    data,
    formula,
    variable_label = NA_character_,
    analysis_label = NA_character_) {

  warning_messages <- character(0)

  fit_error <- NA_character_

  fit <- tryCatch(

    withCallingHandlers(

      survival::coxph(
        formula = formula,
        data = data,
        ties = "efron",
        model = TRUE,
        x = TRUE,
        y = TRUE
      ),

      warning = function(w) {

        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )

        invokeRestart("muffleWarning")
      }
    ),

    error = function(e) {

      fit_error <<- conditionMessage(e)

      NULL
    }
  )


  if (is.null(fit)) {

    return(
      make_empty_cox_result(
        variable_label = variable_label,
        analysis_label = analysis_label,
        status_value = "FAILED",
        warning_message = fit_error
      )
    )
  }


  coefficients <- stats::coef(fit)


  vcov_matrix <- tryCatch(

    stats::vcov(fit),

    error = function(e) {

      NULL
    }
  )


  if (is.null(vcov_matrix)) {

    return(
      make_empty_cox_result(
        variable_label = variable_label,
        analysis_label = analysis_label,
        status_value = "WARNING",
        warning_message = paste(
          c(
            warning_messages,
            "Variance-covariance matrix unavailable"
          ),
          collapse = " | "
        )
      )
    )
  }


  standard_errors <- sqrt(
    diag(vcov_matrix)
  )

  names(standard_errors) <- names(coefficients)


  lower_coef <- coefficients -
    1.96 * standard_errors

  upper_coef <- coefficients +
    1.96 * standard_errors


  model_n <- tryCatch(

    nrow(
      stats::model.frame(fit)
    ),

    error = function(e) {

      NA_integer_
    }
  )


  model_events <- tryCatch({

    y <- fit$y

    if (is.null(y)) {

      NA_integer_

    } else {

      y_matrix <- as.matrix(y)

      if (ncol(y_matrix) >= 2) {

        sum(
          y_matrix[, 2] == 1,
          na.rm = TRUE
        )

      } else {

        NA_integer_
      }
    }

  }, error = function(e) {

    NA_integer_
  })


  tidy_fit <- tryCatch(

    broom::tidy(
      fit,
      exponentiate = FALSE,
      conf.int = FALSE
    ),

    error = function(e) {

      NULL
    }
  )


  if (
    is.null(tidy_fit) ||
    nrow(tidy_fit) == 0
  ) {

    return(
      make_empty_cox_result(
        variable_label = variable_label,
        analysis_label = analysis_label,
        status_value = "FAILED",
        warning_message = paste(
          c(
            warning_messages,
            "No coefficient table returned"
          ),
          collapse = " | "
        ),
        n_value = model_n,
        events_value = model_events
      )
    )
  }


  matched_index <- match(
    tidy_fit$term,
    names(coefficients)
  )


  HR <- exp(
    tidy_fit$estimate
  )


  lower_HR <- exp(
    lower_coef[matched_index]
  )


  upper_HR <- exp(
    upper_coef[matched_index]
  )


  finite_ok <- all(

    is.finite(HR),

    is.finite(lower_HR),

    is.finite(upper_HR),

    is.finite(tidy_fit$estimate),

    is.finite(tidy_fit$std.error),

    is.finite(tidy_fit$statistic),

    is.finite(tidy_fit$p.value)
  )


  if (
    length(warning_messages) == 0 &&
    finite_ok
  ) {

    status_value <- "SUCCESS"

  } else {

    status_value <- "WARNING"
  }


  warning_text <- if (
    length(warning_messages) > 0
  ) {

    paste(
      unique(warning_messages),
      collapse = " | "
    )

  } else if (!finite_ok) {

    "Non-finite coefficient, SE, CI, statistic, or p-value"

  } else {

    NA_character_
  }


  tibble(

    variable = variable_label,

    term_original = tidy_fit$term,

    analysis = analysis_label,

    HR = HR,

    lower_95_CI = lower_HR,

    upper_95_CI = upper_HR,

    coefficient = tidy_fit$estimate,

    standard_error = tidy_fit$std.error,

    z = tidy_fit$statistic,

    p_value = tidy_fit$p.value,

    n = model_n,

    events = model_events,

    status = status_value,

    warning_message = warning_text
  )
}


# ======================================================================
# SECTION 7 — PREPARE COX PREDICTOR
# ======================================================================

prepare_cox_variable <- function(
    data,
    variable) {

  x <- data[[variable]]

  if (
    is.numeric(x) ||
    is.integer(x)
  ) {

    return(
      safe_numeric(x)
    )
  }

  x <- safe_chr(x)

  factor(x)
}


# ======================================================================
# SECTION 8 — UNIVARIATE COX FUNCTION
# ======================================================================

run_univariate_variable <- function(
    data,
    variable,
    analysis_label = "Prespecified univariate Cox") {


  if (!variable %in% names(data)) {

    return(
      make_empty_cox_result(
        variable_label = variable,
        analysis_label = analysis_label,
        status_value = "SKIPPED",
        warning_message = "Variable not found"
      )
    )
  }


  x <- prepare_cox_variable(
    data,
    variable
  )


  valid <- (

    is.finite(data$OS_time_days) &

    !is.na(data$OS_event) &

    !is.na(x)
  )


  n_valid <- sum(valid)


  events_valid <- sum(
    data$OS_event[valid] == 1,
    na.rm = TRUE
  )


  if (n_valid < 5) {

    return(
      make_empty_cox_result(
        variable_label = variable,
        analysis_label = analysis_label,
        status_value = "SKIPPED",
        warning_message = "Fewer than 5 complete observations",
        n_value = n_valid,
        events_value = events_valid
      )
    )
  }


  if (
    length(
      unique(
        data$OS_event[valid]
      )
    ) < 2
  ) {

    return(
      make_empty_cox_result(
        variable_label = variable,
        analysis_label = analysis_label,
        status_value = "SKIPPED",
        warning_message = "Only one event-status category",
        n_value = n_valid,
        events_value = events_valid
      )
    )
  }


  if (is.factor(x)) {

    x_valid <- droplevels(
      x[valid]
    )

    if (
      nlevels(x_valid) < 2
    ) {

      return(
        make_empty_cox_result(
          variable_label = variable,
          analysis_label = analysis_label,
          status_value = "SKIPPED",
          warning_message = "Fewer than 2 predictor levels",
          n_value = n_valid,
          events_value = events_valid
        )
      )
    }

  } else {

    x_valid <- x[valid]

    if (
      length(
        unique(
          x_valid
        )
      ) < 2
    ) {

      return(
        make_empty_cox_result(
          variable_label = variable,
          analysis_label = analysis_label,
          status_value = "SKIPPED",
          warning_message = "Predictor has no variation",
          n_value = n_valid,
          events_value = events_valid
        )
      )
    }
  }


  analysis_data <- tibble(

    OS_time_days =
      data$OS_time_days[valid],

    OS_event =
      data$OS_event[valid],

    predictor =
      x_valid
  )


  if (
    is.factor(
      analysis_data$predictor
    )
  ) {

    analysis_data$predictor <-
      droplevels(
        analysis_data$predictor
      )
  }


  formula <- survival::Surv(
    OS_time_days,
    OS_event
  ) ~ predictor


  run_cox(

    data = analysis_data,

    formula = formula,

    variable_label = variable,

    analysis_label = analysis_label
  )
}


# ======================================================================
# SECTION 9 — READ M3.2
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("SECTION 9 — READING M3.2\n")
cat("============================================================\n")


m32_raw <- readr::read_csv(
  M32_FILE,
  show_col_types = FALSE,
  name_repair = "unique"
)


cat(
  "Raw M3.2 rows    :",
  nrow(m32_raw),
  "\n"
)

cat(
  "Raw M3.2 columns :",
  ncol(m32_raw),
  "\n"
)


required_m32 <- c(
  "sample_id",
  "ecosystem_state"
)


missing_m32 <- setdiff(
  required_m32,
  names(m32_raw)
)


if (
  length(missing_m32) > 0
) {

  stop(
    paste0(
      "\nM3.2 is missing required columns:\n",
      paste(
        missing_m32,
        collapse = ", "
      )
    )
  )
}


# ======================================================================
# SECTION 10 — PREPARE M3.2
# ======================================================================

m32 <- m32_raw %>%

  select(
    any_of(
      c(
        "sample_id",
        "ecosystem_state",
        "cluster_numeric"
      )
    )
  ) %>%

  mutate(

    sample_id_original =
      safe_chr(sample_id),

    patient_id =
      normalize_m32_to_patient(sample_id),

    ecosystem_state =
      safe_chr(ecosystem_state)
  ) %>%

  filter(

    !is.na(patient_id),

    patient_id != "",

    !is.na(ecosystem_state),

    ecosystem_state != ""
  )


cat(
  "\nM3.2 usable rows:",
  nrow(m32),
  "\n"
)

cat(
  "M3.2 unique patient IDs:",
  n_distinct(m32$patient_id),
  "\n"
)


m32_duplicate_ids <- m32 %>%

  count(
    patient_id,
    name = "n"
  ) %>%

  filter(
    n > 1
  )


if (
  nrow(m32_duplicate_ids) > 0
) {

  write_csv(
    m32_duplicate_ids,
    file.path(
      OUT_DIR,
      "M3.7_M3.2_duplicate_patient_ids.csv"
    )
  )

  stop(
    paste0(
      "\nCRITICAL ERROR:\n",
      "M3.2 contains multiple ecosystem-state rows for the same ",
      "patient after normalization.\n\n",
      "The script will NOT arbitrarily select one state.\n\n",
      "Inspect:\n",
      OUT_DIR,
      "/M3.7_M3.2_duplicate_patient_ids.csv"
    )
  )
}


# ======================================================================
# SECTION 11 — FREEZE ECOSYSTEM STATES
# ======================================================================

m32$ecosystem_state <- factor(
  m32$ecosystem_state
)


if (
  "State_1" %in%
  levels(m32$ecosystem_state)
) {

  m32$ecosystem_state <-
    forcats::fct_relevel(
      m32$ecosystem_state,
      "State_1"
    )
}


cat("\nFrozen ecosystem-state distribution:\n")

print(
  table(
    m32$ecosystem_state,
    useNA = "ifany"
  )
)


write_csv(
  m32 %>%
    select(
      sample_id_original,
      patient_id,
      ecosystem_state,
      any_of("cluster_numeric")
    ),
  file.path(
    OUT_DIR,
    "M3.7_frozen_M3.2_ecosystem_states_audit.csv"
  )
)


# ======================================================================
# SECTION 12 — READ M1 CLINICAL DATA
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("SECTION 12 — READING M1 CLINICAL DATA\n")
cat("============================================================\n")


m1_raw <- readr::read_csv(
  M1_FILE,
  show_col_types = FALSE,
  name_repair = "unique"
)


cat(
  "M1 rows    :",
  nrow(m1_raw),
  "\n"
)

cat(
  "M1 columns :",
  ncol(m1_raw),
  "\n"
)


if (
  !"patient_id" %in%
  names(m1_raw)
) {

  stop(
    paste0(
      "\nCRITICAL ERROR:\n",
      "M1_clinical.csv does not contain patient_id."
    )
  )
}


m1 <- m1_raw %>%

  mutate(

    patient_id_original =
      safe_chr(patient_id),

    patient_id =
      normalize_m1_patient(patient_id)
  ) %>%

  filter(

    !is.na(patient_id),

    patient_id != ""
  )


cat(
  "\nM1 normalized unique patient IDs:",
  n_distinct(m1$patient_id),
  "\n"
)


m1_duplicate_ids <- m1 %>%

  count(
    patient_id,
    name = "n"
  ) %>%

  filter(
    n > 1
  )


if (
  nrow(m1_duplicate_ids) > 0
) {

  write_csv(
    m1_duplicate_ids,
    file.path(
      OUT_DIR,
      "M3.7_M1_duplicate_patient_ids.csv"
    )
  )

  stop(
    paste0(
      "\nCRITICAL ERROR:\n",
      "M1 contains multiple clinical records for the same patient.\n\n",
      "The script will not arbitrarily select a clinical record.\n\n",
      "Inspect:\n",
      OUT_DIR,
      "/M3.7_M1_duplicate_patient_ids.csv"
    )
  )
}


# ======================================================================
# SECTION 13 — EXACT PATIENT MATCHING
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("SECTION 13 — EXACT M3.2 → M1 PATIENT MATCHING\n")
cat("============================================================\n")


m32_ids <- unique(m32$patient_id)
m1_ids  <- unique(m1$patient_id)

matched_ids   <- intersect(m32_ids, m1_ids)
m32_only_ids  <- setdiff(m32_ids, m1_ids)
m1_only_ids   <- setdiff(m1_ids, m32_ids)


match_percentage <- if (
  length(m32_ids) > 0
) {

  100 * length(matched_ids) / length(m32_ids)

} else {

  0
}


cat("M3.2 unique patients:", length(m32_ids), "\n")
cat("M1 unique patients:", length(m1_ids), "\n")
cat("Matched patients:", length(matched_ids), "\n")
cat("Match percentage:", round(match_percentage, 2), "%\n")


write_csv(
  tibble(
    metric = c(
      "M3.2_unique_patients",
      "M1_unique_patients",
      "matched_patients",
      "M3.2_only_patients",
      "M1_only_patients",
      "match_percentage"
    ),
    value = c(
      length(m32_ids),
      length(m1_ids),
      length(matched_ids),
      length(m32_only_ids),
      length(m1_only_ids),
      match_percentage
    )
  ),
  file.path(
    OUT_DIR,
    "M3.7_patient_matching_summary.csv"
  )
)


if (
  length(matched_ids) == 0
) {

  stop(
    "\nCRITICAL ERROR — ZERO M3.2/M1 MATCHES.\n"
  )
}


if (
  match_percentage < 95
) {

  warning(
    paste0(
      "\nWARNING: M3.2/M1 match rate is only ",
      round(match_percentage, 2),
      "%. Unmatched patients should be investigated."
    )
  )
}


# ======================================================================
# SECTION 14 — INNER JOIN
# ======================================================================

matched <- m32 %>%

  inner_join(
    m1,
    by = "patient_id",
    suffix = c("_m32", "_m1")
  )


if (nrow(matched) == 0) {

  stop("\nCRITICAL ERROR: inner_join returned zero rows.")
}


cat("\nMATCHING SUCCESSFUL — matched patients:", nrow(matched), "\n")


write_csv(
  matched,
  file.path(
    OUT_DIR,
    "M3.7_matched_M3.2_M1_dataset.csv"
  )
)


# ======================================================================
# SECTION 15 — IDENTIFY CLINICAL SOURCE COLUMNS
# ======================================================================

vital_candidates <- c(
  "vital_status_raw",
  "vital_status",
  "indexed__vital_status"
)

death_candidates <- c(
  "days_to_death",
  "days_to_death_numeric",
  "indexed__days_to_death"
)

followup_candidates <- c(
  "days_to_last_follow_up",
  "days_to_last_followup",
  "days_to_last_follow_up_numeric",
  "indexed__days_to_last_follow_up",
  "days_to_last_known_disease_status",
  "indexed__days_to_last_known_disease_status"
)

age_candidates <- c(
  "age_at_diagnosis_years",
  "age_at_diagnosis",
  "age_at_index",
  "age",
  "age_numeric"
)

# sex_at_birth is tried FIRST: GDC-harmonized, populated for nearly
# every case. paper_Sex (from the marker-paper supplementary table)
# is a last-resort fallback only, since it is far sparser.
sex_candidates <- c(
  "sex_at_birth",
  "sex",
  "sex_raw",
  "gender",
  "paper_Sex"
)

stage_candidates <- c(
  "pathologic_stage",
  "pathologic_stage_raw",
  "pathologic_stage_numeric",
  "ajcc_pathologic_stage"
)

race_candidates <- c(
  "race",
  "race_raw",
  "race_category",
  "indexed__race"
)

smoking_candidates <- c(
  "smoking_group",
  "smoking_status",
  "smoking_status_raw",
  "smoking_history",
  "tobacco_smoking_status",
  "tobacco_smoking_history",
  "indexed__tobacco_smoking_status"
)


vital_col    <- find_first_column(matched, vital_candidates)
death_col    <- find_first_column(matched, death_candidates)
followup_col <- find_first_column(matched, followup_candidates)
age_col      <- find_first_column(matched, age_candidates)
sex_col      <- find_first_column(matched, sex_candidates)
stage_col    <- find_first_column(matched, stage_candidates)
race_col     <- find_first_column(matched, race_candidates)
smoking_col  <- find_first_column(matched, smoking_candidates)


cat("\nClinical source columns selected:\n")
cat("---------------------------------\n")
cat("Vital status :", vital_col, "\n")
cat("Death time   :", death_col, "\n")
cat("Follow-up    :", followup_col, "\n")
cat("Age          :", age_col, "\n")
cat("Sex          :", sex_col, "\n")
cat("Stage        :", stage_col, "\n")
cat("Race         :", race_col, "\n")
cat("Smoking      :", smoking_col, "\n")


# ======================================================================
# SECTION 16 — PRESERVE ORIGINAL M1 SURVIVAL VARIABLES
# ======================================================================

original_M1_OS_event <- if ("OS_event" %in% names(matched)) {
  safe_numeric(matched$OS_event)
} else {
  rep(NA_real_, nrow(matched))
}

original_M1_OS_time_days <- if ("OS_time_days" %in% names(matched)) {
  safe_numeric(matched$OS_time_days)
} else {
  rep(NA_real_, nrow(matched))
}

matched$M1_original_OS_event     <- original_M1_OS_event
matched$M1_original_OS_time_days <- original_M1_OS_time_days


# ======================================================================
# SECTION 17 — OUTCOME RECONSTRUCTION
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("SECTION 17 — RECONSTRUCTING OVERALL SURVIVAL\n")
cat("============================================================\n")


matched$OS_event_reconstructed      <- NA_integer_
matched$OS_time_days_reconstructed  <- NA_real_
matched$OS_source                   <- NA_character_


vital <- rep(NA_character_, nrow(matched))

if (!is.na(vital_col)) {
  vital <- tolower(safe_chr(matched[[vital_col]]))
}

dead_values  <- c("dead", "deceased", "death", "dead/deceased", "1")
alive_values <- c("alive", "living", "0")

matched$OS_event_reconstructed[vital %in% dead_values]  <- 1L
matched$OS_event_reconstructed[vital %in% alive_values] <- 0L


death_time <- rep(NA_real_, nrow(matched))

if (!is.na(death_col)) {
  death_time <- safe_numeric(matched[[death_col]])
}


followup_time <- rep(NA_real_, nrow(matched))

if (!is.na(followup_col)) {
  followup_time <- safe_numeric(matched[[followup_col]])
}


# Explicit death cases
death_valid <-
  matched$OS_event_reconstructed == 1 &
  is.finite(death_time) &
  death_time > 0

matched$OS_time_days_reconstructed[death_valid] <- death_time[death_valid]
matched$OS_source[death_valid] <- "death_time"


# Explicit alive cases
alive_valid <-
  matched$OS_event_reconstructed == 0 &
  is.finite(followup_time) &
  followup_time > 0

matched$OS_time_days_reconstructed[alive_valid] <- followup_time[alive_valid]
matched$OS_source[alive_valid] <- "followup_time"


# Missing vital status but valid death time
fallback_death <-
  is.na(matched$OS_event_reconstructed) &
  is.finite(death_time) &
  death_time > 0

matched$OS_event_reconstructed[fallback_death]     <- 1L
matched$OS_time_days_reconstructed[fallback_death] <- death_time[fallback_death]
matched$OS_source[fallback_death] <- "death_time_vital_missing"


# Missing vital status/death time: valid follow-up implies censoring
fallback_alive <-
  is.na(matched$OS_event_reconstructed) &
  is.finite(followup_time) &
  followup_time > 0

matched$OS_event_reconstructed[fallback_alive]     <- 0L
matched$OS_time_days_reconstructed[fallback_alive] <- followup_time[fallback_alive]
matched$OS_source[fallback_alive] <- "followup_time_vital_missing"


# Final fallback: existing M1 OS_time_days (only when event status known)
existing_time_fallback <-
  is.na(matched$OS_time_days_reconstructed) &
  is.finite(original_M1_OS_time_days) &
  original_M1_OS_time_days > 0 &
  !is.na(matched$OS_event_reconstructed)

matched$OS_time_days_reconstructed[existing_time_fallback] <-
  original_M1_OS_time_days[existing_time_fallback]

matched$OS_source[existing_time_fallback] <- "existing_M1_OS_time_days_fallback"


# Invalid values
invalid_time <-
  !is.finite(matched$OS_time_days_reconstructed) |
  matched$OS_time_days_reconstructed <= 0

matched$OS_time_days_reconstructed[invalid_time] <- NA_real_


# Official M3.7 endpoint
matched$OS_event     <- matched$OS_event_reconstructed
matched$OS_time_days <- matched$OS_time_days_reconstructed
matched$OS_time_years <- matched$OS_time_days / 365.25


cat("\nOS reconstruction sources:\n")
print(table(matched$OS_source, useNA = "ifany"))


write_csv(
  matched %>%
    transmute(
      patient_id,
      M1_original_OS_event,
      M3.7_reconstructed_OS_event = OS_event,
      M1_original_OS_time_days,
      M3.7_reconstructed_OS_time_days = OS_time_days,
      OS_source
    ),
  file.path(
    OUT_DIR,
    "M3.7_survival_outcome_reconstruction_audit.csv"
  )
)


# ======================================================================
# SECTION 18 — SURVIVAL DATASET
# ======================================================================

survival_data <- matched %>%

  filter(
    is.finite(OS_time_days),
    OS_time_days > 0,
    !is.na(OS_event),
    OS_event %in% c(0, 1)
  )


cat("\n============================================================\n")
cat("OVERALL SURVIVAL\n")
cat("============================================================\n")
cat("Matched patients:", nrow(matched), "\n")
cat("Usable OS:", nrow(survival_data), "\n")
cat("Deaths:", sum(survival_data$OS_event == 1, na.rm = TRUE), "\n")
cat("Censored:", sum(survival_data$OS_event == 0, na.rm = TRUE), "\n")


if (nrow(survival_data) == 0) {

  stop(
    paste0(
      "\nCRITICAL ERROR: No patients have usable overall survival.\n",
      "Inspect:\n",
      OUT_DIR,
      "/M3.7_survival_outcome_reconstruction_audit.csv"
    )
  )
}


# ======================================================================
# SECTION 19 — BUILD THE SIX PRESPECIFIED VARIABLES ONLY
# ======================================================================

survival_data$ecosystem_state <-
  factor(safe_chr(survival_data$ecosystem_state))

if ("State_1" %in% levels(survival_data$ecosystem_state)) {

  survival_data$ecosystem_state <-
    forcats::fct_relevel(survival_data$ecosystem_state, "State_1")
}


survival_data$age_numeric <- if (!is.na(age_col)) {
  safe_numeric(survival_data[[age_col]])
} else {
  NA_real_
}


# Sex — clean_sex_factor(), not the generic clean_factor()
survival_data$sex_factor <- if (!is.na(sex_col)) {
  clean_sex_factor(survival_data[[sex_col]])
} else {
  factor(rep(NA_character_, nrow(survival_data)), levels = c("Male", "Female"))
}

cat(
  "\nSex distribution (", sex_col, ") after cleaning:\n",
  sep = ""
)
print(table(survival_data$sex_factor, useNA = "ifany"))


survival_data$stage_factor <- if (!is.na(stage_col)) {
  clean_factor(survival_data[[stage_col]])
} else {
  factor(rep(NA_character_, nrow(survival_data)))
}


survival_data$race_factor <- if (!is.na(race_col)) {
  clean_factor(survival_data[[race_col]])
} else {
  factor(rep(NA_character_, nrow(survival_data)))
}


survival_data$smoking_factor <- if (!is.na(smoking_col)) {
  clean_factor(survival_data[[smoking_col]])
} else {
  factor(rep(NA_character_, nrow(survival_data)))
}


write_csv(
  survival_data,
  file.path(
    OUT_DIR,
    "M3.7_matched_clinical_ecosystem_dataset.csv"
  )
)


prespecified_variables <- c(
  "ecosystem_state",
  "age_numeric",
  "sex_factor",
  "stage_factor",
  "race_factor",
  "smoking_factor"
)


availability_table <- map_dfr(

  prespecified_variables,

  function(v) {

    x <- survival_data[[v]]

    tibble(
      variable = v,
      n_total = length(x),
      n_nonmissing = sum(!is.na(x)),
      n_missing = sum(is.na(x)),
      n_unique = dplyr::n_distinct(x, na.rm = TRUE)
    )
  }
)

print(availability_table)

write_csv(
  availability_table,
  file.path(
    OUT_DIR,
    "M3.7_prespecified_covariate_availability.csv"
  )
)


# ======================================================================
# SECTION 20 — UNIVARIATE COX (ALL SIX PRESPECIFIED VARIABLES ONLY)
# ======================================================================

cat("\n============================================================\n")
cat("UNIVARIATE COX — PRESPECIFIED VARIABLES ONLY\n")
cat("============================================================\n")


prespecified_results <- map_dfr(

  prespecified_variables,

  function(v) {

    cat("Testing:", v, "\n")

    run_univariate_variable(
      survival_data,
      v,
      analysis_label = "Prespecified univariate Cox"
    )
  }
)


M3.7_UNIVARIATE_Cox_hazard_ratios <-

  prespecified_results %>%

  mutate(model_class = "Prespecified univariate") %>%

  select(
    model_class,
    variable,
    term_original,
    analysis,
    HR,
    lower_95_CI,
    upper_95_CI,
    coefficient,
    standard_error,
    z,
    p_value,
    n,
    events,
    status,
    warning_message
  ) %>%

  arrange(is.na(p_value), p_value)


write_csv(
  M3.7_UNIVARIATE_Cox_hazard_ratios,
  file.path(
    OUT_DIR,
    "M3.7_UNIVARIATE_Cox_hazard_ratios.csv"
  )
)


cat(
  "\nUnivariate HR file written. Rows:",
  nrow(M3.7_UNIVARIATE_Cox_hazard_ratios),
  "\n"
)


# ======================================================================
# SECTION 21 — PRIMARY MULTIVARIABLE COX MODEL (ALL SIX VARIABLES)
# ======================================================================

cat("\n============================================================\n")
cat("PRIMARY MULTIVARIABLE COX MODEL\n")
cat("============================================================\n")


primary_formula <-

  survival::Surv(OS_time_days, OS_event) ~

  ecosystem_state +
  age_numeric +
  sex_factor +
  stage_factor +
  race_factor +
  smoking_factor


primary_vars <- c(
  "OS_time_days",
  "OS_event",
  prespecified_variables
)


primary_complete <-

  survival_data %>%

  select(all_of(primary_vars)) %>%

  filter(complete.cases(.)) %>%

  droplevels()


primary_N      <- nrow(primary_complete)
primary_events <- sum(primary_complete$OS_event == 1, na.rm = TRUE)

cat("Complete-case N:", primary_N, "\n")
cat("Events:", primary_events, "\n")


if (primary_N >= 20 && primary_events >= 5) {

  primary_model_result <-

    run_cox(
      data = primary_complete,
      formula = primary_formula,
      variable_label = "Primary multivariable model",
      analysis_label = "Primary multivariable Cox"
    )

} else {

  primary_model_result <-

    make_empty_cox_result(
      variable_label = "Primary multivariable model",
      analysis_label = "Primary multivariable Cox",
      status_value = "SKIPPED",
      warning_message = "Insufficient complete cases/events",
      n_value = primary_N,
      events_value = primary_events
    )
}


write_csv(
  primary_model_result,
  file.path(
    OUT_DIR,
    "M3.7_PRIMARY_multivariate_Cox_model.csv"
  )
)


# ======================================================================
# SECTION 22 — PROPORTIONAL HAZARDS DIAGNOSTIC
# ======================================================================

cat("\n============================================================\n")
cat("PROPORTIONAL HAZARDS DIAGNOSTIC\n")
cat("============================================================\n")


ph_results <- tibble(
  term = character(),
  rho = numeric(),
  chisq = numeric(),
  p_value = numeric(),
  status = character(),
  warning_message = character()
)


if (primary_N >= 20 && primary_events >= 5) {

  ph_warning_messages <- character(0)

  primary_fit <- tryCatch(

    withCallingHandlers(

      survival::coxph(
        formula = primary_formula,
        data = primary_complete,
        ties = "efron"
      ),

      warning = function(w) {
        ph_warning_messages <<- c(ph_warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),

    error = function(e) NULL
  )


  if (!is.null(primary_fit)) {

    zph <- tryCatch(
      survival::cox.zph(primary_fit),
      error = function(e) NULL
    )

    if (!is.null(zph)) {

      zph_table <- as.data.frame(zph$table)
      zph_table$term <- rownames(zph_table)
      rownames(zph_table) <- NULL
      names(zph_table) <- tolower(names(zph_table))

      if (all(c("rho", "chisq", "p", "term") %in% names(zph_table))) {

        ph_results <- tibble(
          term = zph_table$term,
          rho = zph_table$rho,
          chisq = zph_table$chisq,
          p_value = zph_table$p,
          status = ifelse(
            zph_table$p < 0.05,
            "PH_violation_signal",
            "No_PH_violation_signal"
          ),
          warning_message = if (length(ph_warning_messages) > 0) {
            paste(unique(ph_warning_messages), collapse = " | ")
          } else {
            NA_character_
          }
        )
      }
    }
  }
}


write_csv(
  ph_results,
  file.path(
    OUT_DIR,
    "M3.7_primary_cox_PH_test.csv"
  )
)


# ======================================================================
# SECTION 23 — STAGE I-II SENSITIVITY ANALYSIS
# ======================================================================

cat("\n============================================================\n")
cat("STAGE I-II SENSITIVITY ANALYSIS\n")
cat("============================================================\n")


stage_sensitivity_result <-

  make_empty_cox_result(
    variable_label = "Stage I-II sensitivity",
    analysis_label = "Stage I-II sensitivity",
    status_value = "SKIPPED",
    warning_message = "Stage variable unavailable"
  )


if (!is.na(stage_col)) {

  stage_text <- toupper(safe_chr(survival_data[[stage_col]]))
  stage_text <- gsub("[_-]+", " ", stage_text)
  stage_text <- gsub("\\s+", " ", stage_text)
  stage_text <- trimws(stage_text)

  stage_I_II_indicator <-
    str_detect(stage_text, "^(STAGE[ ]*)?(I|IA|IB|II|IIA|IIB)$") |
    str_detect(stage_text, "^(1|1A|1B|2|2A|2B)$")

  stage_I_II <-
    survival_data[
      !is.na(stage_I_II_indicator) & stage_I_II_indicator,
      ,
      drop = FALSE
    ] %>%
    filter(
      is.finite(OS_time_days),
      !is.na(OS_event)
    ) %>%
    droplevels()

  stage_N <- nrow(stage_I_II)
  stage_events <- sum(stage_I_II$OS_event == 1, na.rm = TRUE)

  cat("Stage I-II N:", stage_N, "\n")
  cat("Stage I-II events:", stage_events, "\n")

  if (stage_N >= 20 && stage_events >= 5) {

    sensitivity_formula <-

      survival::Surv(OS_time_days, OS_event) ~

      ecosystem_state +
      age_numeric +
      sex_factor

    stage_sensitivity_result <-

      run_cox(
        data = stage_I_II,
        formula = sensitivity_formula,
        variable_label = "Stage I-II sensitivity",
        analysis_label = "Stage I-II sensitivity"
      )
  }
}


write_csv(
  stage_sensitivity_result,
  file.path(
    OUT_DIR,
    "M3.7_STAGE_I_II_ecosystem_state_Cox_model.csv"
  )
)


# ======================================================================
# SECTION 24 — FOREST PLOT: PRESPECIFIED UNIVARIATE RESULTS
# ======================================================================

forest_univariate <-

  prespecified_results %>%

  filter(
    status %in% c("SUCCESS", "WARNING"),
    is.finite(HR),
    is.finite(lower_95_CI),
    is.finite(upper_95_CI),
    HR > 0,
    lower_95_CI > 0,
    upper_95_CI > 0
  ) %>%

  mutate(
    display_term = paste(variable, term_original, sep = " : ")
  )


if (nrow(forest_univariate) > 0) {

  forest_univariate <-

    forest_univariate %>%

    mutate(
      display_term = forcats::fct_reorder(display_term, HR)
    )

  p_univariate <-

    ggplot(
      forest_univariate,
      aes(x = HR, y = display_term)
    ) +

    geom_point(size = 2.8) +

    geom_errorbar(
      aes(xmin = lower_95_CI, xmax = upper_95_CI),
      width = 0.20
    ) +

    geom_vline(xintercept = 1, linetype = "dashed") +

    scale_x_log10() +

    labs(
      title = "Prespecified Univariate Cox Associations",
      x = "Hazard ratio (log scale)",
      y = NULL
    ) +

    coord_flip() +

    theme_bw(base_size = 12)

  ggsave(
    filename = file.path(
      OUT_DIR,
      "M3.7_prespecified_univariate_forest.png"
    ),
    plot = p_univariate,
    width = 9,
    height = 7,
    dpi = 600
  )
}


# ======================================================================
# SECTION 25 — FOREST PLOT: PRIMARY MULTIVARIABLE MODEL
# ======================================================================

forest_primary <-

  primary_model_result %>%

  filter(
    status %in% c("SUCCESS", "WARNING"),
    is.finite(HR),
    is.finite(lower_95_CI),
    is.finite(upper_95_CI),
    HR > 0,
    lower_95_CI > 0,
    upper_95_CI > 0
  ) %>%

  mutate(display_term = term_original)


if (nrow(forest_primary) > 0) {

  forest_primary <-

    forest_primary %>%

    mutate(
      display_term = forcats::fct_reorder(display_term, HR)
    )

  p_primary <-

    ggplot(
      forest_primary,
      aes(x = HR, y = display_term)
    ) +

    geom_point(size = 3) +

    geom_errorbar(
      aes(xmin = lower_95_CI, xmax = upper_95_CI),
      width = 0.20
    ) +

    geom_vline(xintercept = 1, linetype = "dashed") +

    scale_x_log10() +

    labs(
      title = "Primary Multivariable Cox Model",
      x = "Hazard ratio (log scale)",
      y = NULL
    ) +

    coord_flip() +

    theme_bw(base_size = 12)

  ggsave(
    filename = file.path(
      OUT_DIR,
      "M3.7_primary_multivariable_forest.png"
    ),
    plot = p_primary,
    width = 9,
    height = 7,
    dpi = 600
  )
}


# ======================================================================
# SECTION 26 — MASTER HAZARD-RATIO TABLE
# ======================================================================

master_hazard_ratio_table <-

  bind_rows(

    prespecified_results %>%
      mutate(model_class = "Prespecified univariate"),

    primary_model_result %>%
      mutate(model_class = "Primary multivariable"),

    stage_sensitivity_result %>%
      mutate(model_class = "Stage I-II sensitivity")
  ) %>%

  select(
    model_class,
    variable,
    term_original,
    analysis,
    HR,
    lower_95_CI,
    upper_95_CI,
    coefficient,
    standard_error,
    z,
    p_value,
    n,
    events,
    status,
    warning_message
  )


write_csv(
  master_hazard_ratio_table,
  file.path(
    OUT_DIR,
    "M3.7_MASTER_hazard_ratio_table.csv"
  )
)


# ======================================================================
# SECTION 27 — MODEL STATUS AUDIT
# ======================================================================

model_status_audit <-

  bind_rows(

    prespecified_results %>%
      count(status, name = "n") %>%
      mutate(analysis = "Prespecified univariate"),

    primary_model_result %>%
      count(status, name = "n") %>%
      mutate(analysis = "Primary multivariable"),

    stage_sensitivity_result %>%
      count(status, name = "n") %>%
      mutate(analysis = "Stage I-II sensitivity")
  ) %>%

  select(analysis, status, n)


write_csv(
  model_status_audit,
  file.path(
    OUT_DIR,
    "M3.7_model_status_audit.csv"
  )
)


# ======================================================================
# SECTION 28 — SESSION INFORMATION
# ======================================================================

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "M3.7_sessionInfo.txt")
)


# ======================================================================
# SECTION 29 — FINAL CONSOLE SUMMARY
# ======================================================================

cat("\n============================================================\n")
cat("M3.7 PRESPECIFIED-ONLY ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("Matched patients:", nrow(matched), "\n")
cat("Usable OS:", nrow(survival_data), "\n")
cat("Deaths:", sum(survival_data$OS_event == 1, na.rm = TRUE), "\n")
cat("Censored:", sum(survival_data$OS_event == 0, na.rm = TRUE), "\n")
cat("Primary model status:", primary_model_result$status[1], "\n")
cat("Primary complete-case N:", primary_N, "\n")
cat("Primary events:", primary_events, "\n")

cat("\nKEY OUTPUT FILES\n")
cat("----------------\n")
cat("1. M3.7_UNIVARIATE_Cox_hazard_ratios.csv\n")
cat("2. M3.7_PRIMARY_multivariate_Cox_model.csv\n")
cat("3. M3.7_STAGE_I_II_ecosystem_state_Cox_model.csv\n")
cat("4. M3.7_primary_cox_PH_test.csv\n")
cat("5. M3.7_prespecified_univariate_forest.png\n")
cat("6. M3.7_primary_multivariable_forest.png\n")
cat("7. M3.7_MASTER_hazard_ratio_table.csv\n")
cat("8. M3.7_model_status_audit.csv\n")

cat("\nOutput directory:\n")
cat(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE), "\n")

cat("\nEnd of M3.7 (prespecified-only).\n")
cat("============================================================\n")