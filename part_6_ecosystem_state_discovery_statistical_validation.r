# ======================================================================
# M3.3 — BIOLOGICAL VALIDATION OF LIONESS ECOSYSTEM STATES
# PURPOSE
# -------
# Biologically validate the FROZEN LIONESS ecosystem states from M3.2
# using patient-level standardized module entropy.
#
# INPUTS
# -------
# M3.2_patient_ecosystem_states.csv
# M3.2_standardized_module_entropy.csv
#
# OUTPUT
# -------
# /content/M3_LIONESS_entropy/M3.3_biological_validation/
#
# ANALYSES
# --------
# 1. Frozen ecosystem-state structure
# 2. Module entropy profiles by state
# 3. Module-wise standardized state profiles
# 4. Kruskal-Wallis tests
# 5. Pairwise Wilcoxon tests
# 6. Effect-size summaries
# 7. Top biological modules per state
# 8. Heatmap
# 9. Boxplots
# 10. PCA visualization
# 11. Biological interpretation tables
#
# IMPORTANT
# ---------
# Ecosystem states are NOT reclustered.
# M3.2 ecosystem states are treated as FROZEN.
# ======================================================================


options(stringsAsFactors = FALSE)

cat("\n")
cat("=====================================================================\n")
cat("M3.3 — BIOLOGICAL VALIDATION OF LIONESS ECOSYSTEM STATES\n")
cat("=====================================================================\n\n")

cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# 1. INSTALL / LOAD PACKAGES
# ======================================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "purrr",
  "tibble",
  "stringr",
  "pheatmap"
)

cat("Checking required R packages...\n")

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    cat("Installing:", pkg, "\n")

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org",
      dependencies = TRUE
    )
  }
}

suppressPackageStartupMessages({

  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(stringr)
  library(pheatmap)

})

cat("Packages ready.\n\n")


# ======================================================================
# 2. DIRECTORIES
# ======================================================================

base_dir <- "/content/M3_LIONESS_entropy"

input_dir <- file.path(
  base_dir,
  "M3.2_module_axes"
)

output_dir <- file.path(
  base_dir,
  "M3.3_biological_validation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Output directory:\n")
cat(output_dir, "\n\n")


# ======================================================================
# 3. INPUT FILES
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
    "ERROR: Missing ecosystem-state file:\n",
    state_file
  )

}

if (!file.exists(entropy_file)) {

  stop(
    "ERROR: Missing standardized module entropy file:\n",
    entropy_file
  )

}

cat("FOUND:", state_file, "\n")
cat("FOUND:", entropy_file, "\n\n")

cat("All required files found.\n\n")


# ======================================================================
# 4. READ FROZEN ECOSYSTEM STATES
# ======================================================================

cat("=====================================================================\n")
cat("READING FROZEN ECOSYSTEM STATES\n")
cat("=====================================================================\n")

state_df_raw <- read_csv(
  state_file,
  show_col_types = FALSE
)

cat(
  "Rows:",
  nrow(state_df_raw),
  "\n"
)

cat(
  "Columns:",
  ncol(state_df_raw),
  "\n\n"
)

cat("Columns detected:\n")
print(names(state_df_raw))
cat("\n")


# ======================================================================
# 5. IDENTIFY STATE SAMPLE-ID COLUMN
# ======================================================================

sample_candidates <- c(
  "sample_id",
  "Sample",
  "sample",
  "patient_id",
  "Patient",
  "patient",
  "barcode",
  "Barcode"
)

state_id_candidates_found <- intersect(
  sample_candidates,
  names(state_df_raw)
)

if (length(state_id_candidates_found) == 0) {

  char_cols <- names(state_df_raw)[
    vapply(
      state_df_raw,
      function(x) {
        is.character(x) || is.factor(x)
      },
      logical(1)
    )
  ]

  if (length(char_cols) == 0) {

    stop(
      "Could not identify sample ID in ecosystem-state file."
    )

  }

  state_id_col <- char_cols[1]

} else {

  state_id_col <- state_id_candidates_found[1]

}

cat(
  "Sample ID column:",
  state_id_col,
  "\n"
)


# ======================================================================
# 6. IDENTIFY ECOSYSTEM STATE COLUMN
# ======================================================================

state_candidates <- c(
  "ecosystem_state",
  "ecosystem",
  "state",
  "cluster",
  "cluster_label",
  "Cluster"
)

state_candidates_found <- intersect(
  state_candidates,
  names(state_df_raw)
)

if (length(state_candidates_found) == 0) {

  stop(
    "Could not identify ecosystem-state column."
  )

}

state_col <- state_candidates_found[1]

cat(
  "Ecosystem-state column:",
  state_col,
  "\n\n"
)


# ======================================================================
# 7. CREATE CLEAN FROZEN STATE TABLE
# ======================================================================

state_df <- state_df_raw %>%
  transmute(
    sample_id = as.character(.data[[state_id_col]]),
    ecosystem_state = as.character(.data[[state_col]])
  ) %>%
  mutate(
    sample_id = str_trim(sample_id),
    ecosystem_state = str_trim(ecosystem_state)
  ) %>%
  filter(
    !is.na(sample_id),
    sample_id != "",
    !is.na(ecosystem_state),
    ecosystem_state != ""
  ) %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  )


# ======================================================================
# 8. FREEZE STATE LABELS
# ======================================================================

state_levels <- sort(
  unique(
    state_df$ecosystem_state
  )
)

state_df$ecosystem_state <- factor(
  state_df$ecosystem_state,
  levels = state_levels
)

cat("Frozen ecosystem-state structure:\n")

state_summary <- state_df %>%
  count(
    ecosystem_state,
    name = "n"
  ) %>%
  mutate(
    fraction = n / sum(n)
  )

print(state_summary)

cat("\n")


# ======================================================================
# 9. READ STANDARDIZED MODULE ENTROPY
# ======================================================================

cat("=====================================================================\n")
cat("READING STANDARDIZED MODULE ENTROPY\n")
cat("=====================================================================\n")

entropy_df_raw <- read_csv(
  entropy_file,
  show_col_types = FALSE,
  name_repair = "unique"
)

cat(
  "Rows:",
  nrow(entropy_df_raw),
  "\n"
)

cat(
  "Columns:",
  ncol(entropy_df_raw),
  "\n\n"
)

cat("Detected columns:\n")
print(names(entropy_df_raw))
cat("\n")


# ======================================================================
# 10. IDENTIFY ENTROPY SAMPLE-ID COLUMN
# ======================================================================

entropy_id_candidates_found <- intersect(
  sample_candidates,
  names(entropy_df_raw)
)

if (length(entropy_id_candidates_found) == 0) {

  char_cols <- names(entropy_df_raw)[
    vapply(
      entropy_df_raw,
      function(x) {
        is.character(x) || is.factor(x)
      },
      logical(1)
    )
  ]

  if (length(char_cols) == 0) {

    stop(
      "Could not identify sample ID in standardized module entropy file."
    )

  }

  entropy_id_col <- char_cols[1]

} else {

  entropy_id_col <- entropy_id_candidates_found[1]

}

cat(
  "Sample ID column in module entropy file:",
  entropy_id_col,
  "\n\n"
)


# ======================================================================
# 11. IDENTIFY MODULE COLUMNS
# ======================================================================

module_cols <- names(entropy_df_raw)

module_cols <- module_cols[
  grepl(
    "^S_module_",
    module_cols
  )
]


# ----------------------------------------------------------------------
# Fallback if S_module_ names are not detected
# ----------------------------------------------------------------------

if (length(module_cols) < 2) {

  numeric_cols <- names(entropy_df_raw)[
    vapply(
      entropy_df_raw,
      is.numeric,
      logical(1)
    )
  ]

  module_cols <- setdiff(
    numeric_cols,
    entropy_id_col
  )

}

if (length(module_cols) == 0) {

  stop(
    "No module entropy columns detected."
  )

}

cat(
  "Module entropy features detected:",
  length(module_cols),
  "\n\n"
)


# ======================================================================
# 12. CLEAN MODULE ENTROPY DATA
# ======================================================================

entropy_analysis <- entropy_df_raw %>%
  transmute(
    sample_id = as.character(
      .data[[entropy_id_col]]
    ),
    across(
      all_of(module_cols),
      ~ suppressWarnings(
        as.numeric(.x)
      )
    )
  ) %>%
  mutate(
    sample_id = str_trim(sample_id)
  ) %>%
  filter(
    !is.na(sample_id),
    sample_id != ""
  ) %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  )


# ======================================================================
# 13. MODULE QC
# ======================================================================

cat("=====================================================================\n")
cat("MODULE QC\n")
cat("=====================================================================\n")

module_qc <- tibble(
  module = module_cols
) %>%
  mutate(

    n = nrow(entropy_analysis),

    missing = map_int(
      module,
      function(m) {

        sum(
          is.na(
            entropy_analysis[[m]]
          )
        )

      }
    ),

    sd = map_dbl(
      module,
      function(m) {

        sd(
          entropy_analysis[[m]],
          na.rm = TRUE
        )

      }
    )

  )

print(module_qc)

cat("\n")

cat(
  "Modules before QC:",
  length(module_cols),
  "\n"
)


# ======================================================================
# 14. RETAIN USABLE MODULES
# ======================================================================

valid_modules <- module_qc %>%
  filter(
    missing < n,
    !is.na(sd),
    sd > 0
  ) %>%
  pull(module)

if (length(valid_modules) < 2) {

  stop(
    "Fewer than 2 usable module entropy features remain after QC."
  )

}

cat(
  "Modules after QC:",
  length(valid_modules),
  "\n\n"
)


entropy_analysis <- entropy_analysis %>%
  select(
    sample_id,
    all_of(valid_modules)
  )


# ======================================================================
# 15. MERGE FROZEN STATES WITH MODULE ENTROPY
# ======================================================================

cat("=====================================================================\n")
cat("MERGING STATES WITH MODULE ENTROPY\n")
cat("=====================================================================\n")

cat(
  "Patients in frozen state table:",
  nrow(state_df),
  "\n"
)

cat(
  "Patients in module entropy:",
  nrow(entropy_analysis),
  "\n\n"
)


analysis_df <- state_df %>%
  inner_join(
    entropy_analysis,
    by = "sample_id"
  )


cat(
  "Patients after merge:",
  nrow(analysis_df),
  "\n\n"
)


if (nrow(analysis_df) == 0) {

  stop(
    "No patients matched between state and entropy files."
  )

}


# ======================================================================
# 16. VERIFY MERGED STATE STRUCTURE
# ======================================================================

merged_state_summary <- analysis_df %>%
  count(
    ecosystem_state,
    name = "n"
  ) %>%
  mutate(
    fraction = n / sum(n)
  )

cat("Frozen ecosystem-state structure:\n")
print(merged_state_summary)

cat("\n")


# ======================================================================
# 17. MODULE ENTROPY PROFILES BY STATE
# ======================================================================

cat("=====================================================================\n")
cat("MODULE ENTROPY PROFILES BY ECOSYSTEM STATE\n")
cat("=====================================================================\n")

state_profiles <- analysis_df %>%
  select(
    ecosystem_state,
    all_of(valid_modules)
  ) %>%
  pivot_longer(
    cols = all_of(valid_modules),
    names_to = "module",
    values_to = "entropy"
  ) %>%
  group_by(
    ecosystem_state,
    module
  ) %>%
  summarise(

    n = sum(
      !is.na(entropy)
    ),

    mean = ifelse(
      all(is.na(entropy)),
      NA_real_,
      mean(
        entropy,
        na.rm = TRUE
      )
    ),

    median = ifelse(
      all(is.na(entropy)),
      NA_real_,
      median(
        entropy,
        na.rm = TRUE
      )
    ),

    sd = ifelse(
      sum(!is.na(entropy)) <= 1,
      NA_real_,
      sd(
        entropy,
        na.rm = TRUE
      )
    ),

    .groups = "drop"

  )


cat(
  "State-module profiles calculated successfully.\n\n"
)


# ======================================================================
# 18. MODULE-WISE GLOBAL REFERENCE
# ======================================================================
#
# For each module:
#
# global_mean = mean entropy across all 75 patients
# global_sd   = SD entropy across all 75 patients
#
# State standardized score:
#
#     Z = (state_mean - global_mean) / global_sd
#
# This is calculated independently for every module.
# ======================================================================

cat("=====================================================================\n")
cat("CALCULATING MODULE-WISE REFERENCE STATISTICS\n")
cat("=====================================================================\n")

module_reference <- analysis_df %>%
  select(
    all_of(valid_modules)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "module",
    values_to = "entropy"
  ) %>%
  group_by(
    module
  ) %>%
  summarise(

    global_n = sum(
      !is.na(entropy)
    ),

    global_mean = ifelse(
      all(is.na(entropy)),
      NA_real_,
      mean(
        entropy,
        na.rm = TRUE
      )
    ),

    global_sd = ifelse(
      sum(!is.na(entropy)) <= 1,
      NA_real_,
      sd(
        entropy,
        na.rm = TRUE
      )
    ),

    .groups = "drop"

  )


print(
  head(
    module_reference,
    10
  )
)

cat("\n")


# ======================================================================
# 19. STANDARDIZE STATE MODULE PROFILES
# ======================================================================
#
# IMPORTANT:
# ----------
# This section deliberately avoids:
#
#     entropy_analysis[[module]]
#
# inside mutate().
#
# The state profile is already summarized.
# The module-specific global mean and SD are joined by module name.
# Therefore every row has exactly one reference mean and SD.
# ======================================================================

cat("=====================================================================\n")
cat("STANDARDIZED STATE MODULE PROFILES\n")
cat("=====================================================================\n")

state_profiles_standardized <- state_profiles %>%

  left_join(
    module_reference,
    by = "module"
  ) %>%

  mutate(

    standardized_entropy = case_when(

      is.na(mean) ~ NA_real_,

      is.na(global_mean) ~ NA_real_,

      is.na(global_sd) ~ 0,

      global_sd == 0 ~ 0,

      TRUE ~ (
        mean - global_mean
      ) / global_sd

    )

  )


# ----------------------------------------------------------------------
# Verify standardization
# ----------------------------------------------------------------------

if (
  nrow(state_profiles_standardized) == 0
) {

  stop(
    "Standardized state profiles contain zero rows."
  )

}

cat(
  "Standardized profiles calculated successfully.\n"
)

cat(
  "Rows:",
  nrow(state_profiles_standardized),
  "\n"
)

cat(
  "Modules:",
  length(unique(state_profiles_standardized$module)),
  "\n"
)

cat(
  "States:",
  length(unique(state_profiles_standardized$ecosystem_state)),
  "\n\n"
)


# ======================================================================
# 20. SAVE PROFILE TABLES
# ======================================================================

write_csv(
  state_profiles,
  file.path(
    output_dir,
    "M3.3_state_module_entropy_profiles.csv"
  )
)

write_csv(
  state_profiles_standardized,
  file.path(
    output_dir,
    "M3.3_standardized_state_module_profiles.csv"
  )
)

write_csv(
  module_reference,
  file.path(
    output_dir,
    "M3.3_module_reference_statistics.csv"
  )
)


# ======================================================================
# 21. KRUSKAL-WALLIS TESTS
# ======================================================================

cat("=====================================================================\n")
cat("MODULE DIFFERENCES BETWEEN ECOSYSTEM STATES\n")
cat("=====================================================================\n")

kw_results <- map_dfr(
  valid_modules,
  function(mod) {

    tmp <- analysis_df %>%
      select(
        ecosystem_state,
        value = all_of(mod)
      ) %>%
      filter(
        !is.na(value),
        !is.na(ecosystem_state)
      )

    n_states <- length(
      unique(
        tmp$ecosystem_state
      )
    )

    if (
      nrow(tmp) < 3 ||
      n_states < 2
    ) {

      return(
        tibble(
          module = mod,
          statistic = NA_real_,
          p_value = NA_real_
        )
      )

    }

    kt <- tryCatch(

      kruskal.test(
        value ~ ecosystem_state,
        data = tmp
      ),

      error = function(e) {
        NULL
      }

    )

    if (is.null(kt)) {

      return(
        tibble(
          module = mod,
          statistic = NA_real_,
          p_value = NA_real_
        )
      )

    }

    tibble(
      module = mod,
      statistic = unname(
        kt$statistic
      ),
      p_value = kt$p.value
    )

  }
) %>%

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


cat("\nTop modules by FDR:\n\n")

print(
  head(
    kw_results,
    20
  )
)

cat("\n")


write_csv(
  kw_results,
  file.path(
    output_dir,
    "M3.3_kruskal_wallis_module_tests.csv"
  )
)


# ======================================================================
# 22. PAIRWISE WILCOXON TESTS
# ======================================================================

cat("=====================================================================\n")
cat("PAIRWISE STATE COMPARISONS\n")
cat("=====================================================================\n")

available_states <- levels(
  droplevels(
    analysis_df$ecosystem_state
  )
)

if (length(available_states) >= 2) {

  state_pairs <- combn(
    available_states,
    2,
    simplify = FALSE
  )

} else {

  state_pairs <- list()

}


pairwise_results <- map_dfr(
  valid_modules,
  function(mod) {

    if (length(state_pairs) == 0) {

      return(
        tibble(
          module = mod,
          state_1 = character(),
          state_2 = character(),
          statistic = numeric(),
          p_value = numeric(),
          median_state_1 = numeric(),
          median_state_2 = numeric(),
          median_difference = numeric()
        )
      )

    }

    map_dfr(
      state_pairs,
      function(pair) {

        tmp <- analysis_df %>%
          filter(
            ecosystem_state %in% pair
          ) %>%
          select(
            ecosystem_state,
            value = all_of(mod)
          ) %>%
          filter(
            !is.na(value)
          )

        if (
          nrow(tmp) < 3 ||
          length(
            unique(
              tmp$ecosystem_state
            )
          ) < 2
        ) {

          return(
            tibble(
              module = mod,
              state_1 = pair[1],
              state_2 = pair[2],
              statistic = NA_real_,
              p_value = NA_real_,
              median_state_1 = NA_real_,
              median_state_2 = NA_real_,
              median_difference = NA_real_
            )
          )

        }

        wt <- tryCatch(

          suppressWarnings(
            wilcox.test(
              value ~ ecosystem_state,
              data = tmp,
              exact = FALSE
            )
          ),

          error = function(e) {
            NULL
          }

        )

        med1 <- median(
          tmp$value[
            tmp$ecosystem_state == pair[1]
          ],
          na.rm = TRUE
        )

        med2 <- median(
          tmp$value[
            tmp$ecosystem_state == pair[2]
          ],
          na.rm = TRUE
        )

        if (is.null(wt)) {

          return(
            tibble(
              module = mod,
              state_1 = pair[1],
              state_2 = pair[2],
              statistic = NA_real_,
              p_value = NA_real_,
              median_state_1 = med1,
              median_state_2 = med2,
              median_difference = med1 - med2
            )
          )

        }

        tibble(
          module = mod,
          state_1 = pair[1],
          state_2 = pair[2],
          statistic = unname(
            wt$statistic
          ),
          p_value = wt$p.value,
          median_state_1 = med1,
          median_state_2 = med2,
          median_difference = med1 - med2
        )

      }
    )

  }
) %>%

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


cat(
  "Pairwise testing completed.\n\n"
)


write_csv(
  pairwise_results,
  file.path(
    output_dir,
    "M3.3_pairwise_module_tests.csv"
  )
)


# ======================================================================
# 23. EFFECT-SIZE SUMMARY
# ======================================================================
#
# Kruskal-Wallis epsilon-squared:
#
#     epsilon² = (H - k + 1) / (N - k)
#
# where:
#
# H = Kruskal-Wallis statistic
# k = number of states
# N = number of observations
#
# Negative values caused by small H are truncated to zero.
# ======================================================================

cat("=====================================================================\n")
cat("EFFECT-SIZE SUMMARY\n")
cat("=====================================================================\n")

effect_results <- map_dfr(
  valid_modules,
  function(mod) {

    tmp <- analysis_df %>%
      select(
        ecosystem_state,
        value = all_of(mod)
      ) %>%
      filter(
        !is.na(value),
        !is.na(ecosystem_state)
      )

    N <- nrow(tmp)

    k <- length(
      unique(
        tmp$ecosystem_state
      )
    )

    if (
      N <= k ||
      k < 2
    ) {

      return(
        tibble(
          module = mod,
          N = N,
          states = k,
          H = NA_real_,
          epsilon_squared = NA_real_
        )
      )

    }

    kt <- tryCatch(

      kruskal.test(
        value ~ ecosystem_state,
        data = tmp
      ),

      error = function(e) {
        NULL
      }

    )

    if (is.null(kt)) {

      return(
        tibble(
          module = mod,
          N = N,
          states = k,
          H = NA_real_,
          epsilon_squared = NA_real_
        )
      )

    }

    H <- as.numeric(
      kt$statistic
    )

    epsilon_squared <- max(
      0,
      (H - k + 1) / (N - k)
    )

    tibble(
      module = mod,
      N = N,
      states = k,
      H = H,
      epsilon_squared = epsilon_squared
    )

  }
) %>%

  arrange(
    desc(epsilon_squared)
  )


write_csv(
  effect_results,
  file.path(
    output_dir,
    "M3.3_module_effect_sizes.csv"
  )
)


cat("Effect-size analysis completed.\n\n")


# ======================================================================
# 24. TOP MODULES CHARACTERIZING EACH STATE
# ======================================================================

cat("=====================================================================\n")
cat("TOP MODULES CHARACTERIZING EACH ECOSYSTEM STATE\n")
cat("=====================================================================\n")


state_top_modules <- state_profiles_standardized %>%

  filter(
    !is.na(standardized_entropy)
  ) %>%

  group_by(
    ecosystem_state
  ) %>%

  arrange(
    desc(
      abs(
        standardized_entropy
      )
    ),
    .by_group = TRUE
  ) %>%

  slice_head(
    n = 10
  ) %>%

  ungroup()


for (
  st in available_states
) {

  cat("\n")
  cat(st, ":\n\n")

  tmp <- state_top_modules %>%
    filter(
      ecosystem_state == st
    ) %>%
    select(
      module,
      mean,
      median,
      global_mean,
      global_sd,
      standardized_entropy
    )

  print(tmp)

}


cat("\n")


write_csv(
  state_top_modules,
  file.path(
    output_dir,
    "M3.3_top_modules_by_ecosystem_state.csv"
  )
)


# ======================================================================
# 25. STATE-SPECIFIC DIRECTIONAL BIOLOGY
# ======================================================================
#
# +1 SD or greater:
#     Strongly elevated
#
# +0.5 to <1 SD:
#     Elevated
#
# -0.5 to >-1 SD:
#     Reduced
#
# <= -1 SD:
#     Strongly reduced
# ======================================================================

state_directional <- state_profiles_standardized %>%

  mutate(

    direction = case_when(

      standardized_entropy >= 1 ~
        "Strongly elevated",

      standardized_entropy >= 0.5 ~
        "Elevated",

      standardized_entropy <= -1 ~
        "Strongly reduced",

      standardized_entropy <= -0.5 ~
        "Reduced",

      TRUE ~
        "Near cohort average"

    )

  )


write_csv(
  state_directional,
  file.path(
    output_dir,
    "M3.3_state_directional_biology.csv"
  )
)


# ======================================================================
# 26. HEATMAP
# ======================================================================

cat("=====================================================================\n")
cat("GENERATING STATE MODULE HEATMAP\n")
cat("=====================================================================\n")


heat_df <- state_profiles_standardized %>%

  select(
    ecosystem_state,
    module,
    standardized_entropy
  ) %>%

  pivot_wider(
    names_from = ecosystem_state,
    values_from = standardized_entropy
  )


heat_matrix <- heat_df %>%

  column_to_rownames(
    "module"
  ) %>%

  as.matrix()


heat_matrix <- heat_matrix[
  apply(
    heat_matrix,
    1,
    function(x) {

      all(
        is.finite(x)
      )

    }
  ),
  ,
  drop = FALSE
]


if (
  nrow(heat_matrix) >= 2 &&
  ncol(heat_matrix) >= 2
) {

  pdf(
    file.path(
      output_dir,
      "M3.3_state_module_heatmap.pdf"
    ),
    width = 9,
    height = 12
  )

  pheatmap(
    heat_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "none",
    border_color = NA,
    main = "Standardized Module Entropy by Ecosystem State"
  )

  dev.off()


  png(
    file.path(
      output_dir,
      "M3.3_state_module_heatmap.png"
    ),
    width = 1800,
    height = 2400,
    res = 220
  )

  pheatmap(
    heat_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "none",
    border_color = NA,
    main = "Standardized Module Entropy by Ecosystem State"
  )

  dev.off()


  cat(
    "Heatmap generated successfully.\n\n"
  )

} else {

  cat(
    "Heatmap skipped because matrix dimensions are insufficient.\n\n"
  )

}


# ======================================================================
# 27. TOP MODULE BOXPLOTS
# ======================================================================

cat("=====================================================================\n")
cat("GENERATING TOP MODULE BOXPLOTS\n")
cat("=====================================================================\n")


# ----------------------------------------------------------------------
# Select modules ranked by Kruskal-Wallis FDR
# ----------------------------------------------------------------------

kw_results_filtered <- kw_results %>%
  filter(
    !is.na(FDR),
    !is.na(p_value)
  ) %>%
  arrange(
    FDR,
    p_value
  )


# ----------------------------------------------------------------------
# Determine number of modules to plot
#
# IMPORTANT:
# Do NOT use:
#
#     slice_head(n = min(9, n()))
#
# because newer dplyr versions require slice_head(n=...)
# to receive a constant integer.
# ----------------------------------------------------------------------

top_n <- min(
  9L,
  nrow(kw_results_filtered)
)


# ----------------------------------------------------------------------
# Select top modules
# ----------------------------------------------------------------------

if (top_n > 0L) {

  top_plot_modules <- kw_results_filtered %>%
    slice_head(
      n = top_n
    ) %>%
    pull(
      module
    )

} else {

  top_plot_modules <- character(0)

}


cat(
  "Number of modules selected for boxplots:",
  length(top_plot_modules),
  "\n\n"
)


# ----------------------------------------------------------------------
# Generate boxplots
# ----------------------------------------------------------------------

if (length(top_plot_modules) > 0L) {

  plot_long <- analysis_df %>%
    select(
      sample_id,
      ecosystem_state,
      all_of(top_plot_modules)
    ) %>%
    pivot_longer(
      cols = all_of(top_plot_modules),
      names_to = "module",
      values_to = "entropy"
    ) %>%
    filter(
      !is.na(entropy),
      !is.na(ecosystem_state)
    )


  # --------------------------------------------------------------------
  # Boxplot
  # --------------------------------------------------------------------

  p_box <- ggplot(
    plot_long,
    aes(
      x = ecosystem_state,
      y = entropy
    )
  ) +
    geom_boxplot(
      outlier.alpha = 0.35
    ) +
    facet_wrap(
      ~ module,
      scales = "free_y",
      ncol = 3
    ) +
    theme_bw() +
    labs(
      title =
        "Top Differential Module Entropy by Ecosystem State",

      x =
        "Frozen ecosystem state",

      y =
        "Standardized module entropy"
    )


  # --------------------------------------------------------------------
  # Save PDF
  # --------------------------------------------------------------------

  ggsave(
    file.path(
      output_dir,
      "M3.3_top_module_boxplots.pdf"
    ),
    p_box,
    width = 12,
    height = 10
  )


  # --------------------------------------------------------------------
  # Save PNG
  # --------------------------------------------------------------------

  ggsave(
    file.path(
      output_dir,
      "M3.3_top_module_boxplots.png"
    ),
    p_box,
    width = 12,
    height = 10,
    dpi = 300
  )


  cat(
    "Boxplots generated successfully.\n\n"
  )


} else {

  cat(
    "No modules were available for boxplot generation.\n\n"
  )

}


# ======================================================================
# END OF SECTION 27
# ======================================================================

# ======================================================================
# 28. PCA OF PATIENT MODULE ENTROPY
# ======================================================================

cat("=====================================================================\n")
cat("PCA OF PATIENT MODULE ENTROPY\n")
cat("=====================================================================\n")


pca_complete <- analysis_df %>%

  select(
    sample_id,
    ecosystem_state,
    all_of(valid_modules)
  ) %>%

  filter(
    if_all(
      all_of(valid_modules),
      ~ is.finite(.x)
    )
  )


if (
  nrow(pca_complete) >= 3 &&
  length(valid_modules) >= 2
) {

  pca_matrix <- pca_complete %>%

    select(
      all_of(valid_modules)
    ) %>%

    as.matrix()


  pca_fit <- prcomp(
    pca_matrix,
    center = TRUE,
    scale. = TRUE
  )


  n_pcs <- min(
    5,
    ncol(
      pca_fit$x
    )
  )


  pca_scores <- as.data.frame(
    pca_fit$x[
      ,
      1:n_pcs,
      drop = FALSE
    ]
  )


  pca_scores$sample_id <-
    pca_complete$sample_id

  pca_scores$ecosystem_state <-
    pca_complete$ecosystem_state


  variance_explained <- (
    pca_fit$sdev^2
  ) / sum(
    pca_fit$sdev^2
  )


  if (
    ncol(pca_scores) >= 2
  ) {

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

      theme_bw() +

      labs(

        title =
          "PCA of Patient-Level Module Entropy",

        x = paste0(
          "PC1 (",
          round(
            variance_explained[1] * 100,
            1
          ),
          "%)"
        ),

        y = paste0(
          "PC2 (",
          round(
            variance_explained[2] * 100,
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
        "M3.3_module_entropy_PCA.pdf"
      ),
      pca_plot,
      width = 8,
      height = 6
    )


    ggsave(
      file.path(
        output_dir,
        "M3.3_module_entropy_PCA.png"
      ),
      pca_plot,
      width = 8,
      height = 6,
      dpi = 300
    )

  }


  write_csv(
    pca_scores,
    file.path(
      output_dir,
      "M3.3_module_entropy_PCA_scores.csv"
    )
  )


  pca_variance <- tibble(

    PC = paste0(
      "PC",
      seq_along(
        variance_explained
      )
    ),

    variance_explained =
      variance_explained

  )


  write_csv(
    pca_variance,
    file.path(
      output_dir,
      "M3.3_PCA_variance_explained.csv"
    )
  )


  cat(
    "PCA completed successfully.\n\n"
  )

} else {

  cat(
    "PCA skipped because insufficient complete observations/features.\n\n"
  )

}


# ======================================================================
# 29. STATE SEPARATION SUMMARY
# ======================================================================

cat("=====================================================================\n")
cat("ECOSYSTEM-STATE BIOLOGICAL VALIDATION SUMMARY\n")
cat("=====================================================================\n")


significant_modules <- kw_results %>%

  filter(
    !is.na(FDR),
    FDR < 0.05
  )


cat(
  "Total modules tested:",
  nrow(kw_results),
  "\n"
)


cat(
  "Modules with FDR < 0.05:",
  nrow(significant_modules),
  "\n"
)


cat(
  "Modules with FDR < 0.01:",
  sum(
    kw_results$FDR < 0.01,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Modules with FDR < 0.001:",
  sum(
    kw_results$FDR < 0.001,
    na.rm = TRUE
  ),
  "\n\n"
)


# ======================================================================
# 30. SAVE VALIDATION SUMMARY
# ======================================================================

validation_summary <- tibble(

  metric = c(

    "Total patients in frozen states",

    "Patients included in biological validation",

    "Number of ecosystem states",

    "Number of module features",

    "Modules FDR < 0.05",

    "Modules FDR < 0.01",

    "Modules FDR < 0.001"

  ),

  value = c(

    nrow(state_df),

    nrow(analysis_df),

    nlevels(
      droplevels(
        analysis_df$ecosystem_state
      )
    ),

    length(valid_modules),

    sum(
      kw_results$FDR < 0.05,
      na.rm = TRUE
    ),

    sum(
      kw_results$FDR < 0.01,
      na.rm = TRUE
    ),

    sum(
      kw_results$FDR < 0.001,
      na.rm = TRUE
    )

  )

)


write_csv(
  validation_summary,
  file.path(
    output_dir,
    "M3.3_validation_summary.csv"
  )
)


# ======================================================================
# 31. SAVE MERGED VALIDATION DATASET
# ======================================================================

write_csv(
  analysis_df,
  file.path(
    output_dir,
    "M3.3_merged_validation_dataset.csv"
  )
)


# ======================================================================
# 32. SAVE STATE SUMMARY
# ======================================================================

write_csv(
  state_summary,
  file.path(
    output_dir,
    "M3.3_frozen_state_summary.csv"
  )
)


# ======================================================================
# 33. FINAL OUTPUT INVENTORY
# ======================================================================

cat("=====================================================================\n")
cat("M3.3 COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Output files:\n\n")

output_files <- list.files(
  output_dir,
  full.names = FALSE
)

for (
  f in output_files
) {

  cat(
    " -",
    f,
    "\n"
  )

}


cat("\n")

cat(
  "Finished:",
  as.character(Sys.time()),
  "\n\n"
)

cat("=====================================================================\n")
cat("BIOLOGICAL VALIDATION COMPLETE\n")
cat("=====================================================================\n")