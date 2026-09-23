# ======================================================================
# V6b — GSE72094 SURVIVAL VALIDATION UNDER THE LOCKED STATE ASSIGNMENT
#
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

output_dir <- "/content/GSE72094_validation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

states_file <- file.path(output_dir, "GSE72094_ecosystem_states.csv")
if (!file.exists(states_file)) {
  stop("ERROR: ", states_file, " not found. Run V5 first.")
}
if (!file.exists(GSE72094_CLINICAL_FILE)) {
  stop("ERROR: GSE72094_CLINICAL_FILE not found:\n", GSE72094_CLINICAL_FILE,
       "\nEdit the path at the top of this cell.")
}

states   <- read_csv(states_file, show_col_types = FALSE)
clinical <- read_csv(GSE72094_CLINICAL_FILE, show_col_types = FALSE)

cat("Clinical file columns (use these to fix auto-detection below if needed):\n")
print(colnames(clinical))
cat("\n")

# ---- identify sample id column ----
sample_candidates <- c("sample_id", "Sample_ID", "sample", "patient",
                        "Patient", "geo_accession", "GSM")
existing_candidates <- sample_candidates[sample_candidates %in% colnames(clinical)]
clinical_id_col <- if (length(existing_candidates) > 0) existing_candidates[1] else colnames(clinical)[1]

clinical$sample_id <- trimws(as.character(clinical[[clinical_id_col]]))
cat("Using clinical sample-ID column:", clinical_id_col, "\n")

# ---- identify survival-time column ----
time_candidates <- c("OS_time_days", "os_time_days", "OS_MONTHS", "os_months",
                      "overall_survival_months", "time", "OS_time",
                      "survival_time", "follow_up_months", "follow_up_days",
                      "days before death/censor", "months before death/censor",
                      "survival time", "survival (months)", "survival (days)",
                      "survival_time_in_days")
# ^ "survival_time_in_days" is GSE72094's exact GEO field name, confirmed
#   from your own tested External_Validation_GSE72094.ipynb - already in
#   days, unlike GSE31210 (days, different field name) or GSE50081 (years).
time_col <- intersect(time_candidates, colnames(clinical))

if (length(time_col) == 0) {
  stop("Could not auto-detect a survival-time column.\n",
       "Available columns are printed above — add your column name to\n",
       "'time_candidates' above (or hardcode time_col <- \"your_column\"),\n",
       "then re-run this cell.")
}
time_col <- time_col[1]
cat("Using survival-time column:", time_col, "\n")

# ---- identify vital-status column ----
event_candidates <- c("OS_event", "os_event", "vital_status", "OS_STATUS",
                       "status", "death", "event", "death/censor")
# ^ "vital_status" (already listed) is GSE72094's exact GEO field name.
event_col <- intersect(event_candidates, colnames(clinical))

if (length(event_col) == 0) {
  stop("Could not auto-detect a vital-status column.\n",
       "Available columns are printed above — add your column name to\n",
       "'event_candidates' above (or hardcode event_col <- \"your_column\"),\n",
       "then re-run this cell.")
}
event_col <- event_col[1]
cat("Using vital-status column:", event_col, "\n\n")

clinical$OS_time_raw <- suppressWarnings(as.numeric(clinical[[time_col]]))

# ---- normalize vital status to 1 = dead / 0 = censored ----
raw_status <- clinical[[event_col]]

if (is.numeric(raw_status)) {
  clinical$OS_event <- as.integer(raw_status)
} else {
  raw_status_lower <- tolower(trimws(as.character(raw_status)))
  clinical$OS_event <- case_when(
    raw_status_lower %in% c("dead", "death", "deceased", "1", "yes") ~ 1L,
    raw_status_lower %in% c("alive", "living", "0", "no")            ~ 0L,
    TRUE ~ NA_integer_
  )
}

cat("Vital-status coding check (verify this looks correct before continuing):\n")
print(table(raw_original = clinical[[event_col]], OS_event = clinical$OS_event, useNA = "ifany"))

# ---- unit detection / conversion ----
if (GSE72094_TIME_UNIT == "auto") {
  time_col_lower <- tolower(time_col)
  if (grepl("day", time_col_lower)) {
    resolved_time_unit <- "days"
  } else if (grepl("month", time_col_lower)) {
    resolved_time_unit <- "months"
  } else {
    stop("Could not auto-detect whether '", time_col, "' is in days or ",
         "months from its name. Set GSE72094_TIME_UNIT to \"days\" or ",
         "\"months\" manually above, then re-run this cell.")
  }
  cat("Auto-detected time unit for column '", time_col, "': ", resolved_time_unit, "\n", sep = "")
} else {
  resolved_time_unit <- GSE72094_TIME_UNIT
  cat("Using manually set time unit:", resolved_time_unit, "\n")
}

if (resolved_time_unit == "months") {
  clinical$OS_time_days <- clinical$OS_time_raw * 30.44
} else {
  clinical$OS_time_days <- clinical$OS_time_raw
}

cat("\nSurvival-time range (days) after unit conversion:\n")
print(summary(clinical$OS_time_days))

# ---- merge locked states with clinical survival data ----
survival_data <- states %>%
  inner_join(clinical %>% select(sample_id, OS_time_days, OS_event), by = "sample_id") %>%
  filter(is.finite(OS_time_days), OS_time_days > 0, !is.na(OS_event))

cat("\nPatients with both a frozen-state assignment and usable survival data:",
    nrow(survival_data), "\n")

if (nrow(survival_data) < 10) {
  stop("Fewer than 10 patients have usable survival data after merging — ",
       "check that sample IDs match between the expression file and the ",
       "clinical file (print both and compare a few values).")
}

survival_data$ecosystem_state <- factor(survival_data$ecosystem_state)
cat("\nState sizes in the merged survival dataset:\n")
print(table(survival_data$ecosystem_state))

write_csv(survival_data, file.path(output_dir, "V6_GSE72094_survival_data_used.csv"))

# ======================================================================
# KAPLAN-MEIER
# ======================================================================

km_fit <- survfit(Surv(OS_time_days, OS_event) ~ ecosystem_state, data = survival_data)

km_plot <- ggsurvplot(
  km_fit,
  data = survival_data,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "GSE72094 — Overall Survival by Frozen TCGA Ecosystem State",
  xlab = "Days"
)

pdf(file.path(output_dir, "V6_GSE72094_KM_plot.pdf"), width = 8, height = 8)
print(km_plot)
dev.off()

# ======================================================================
# LOG-RANK TEST
# ======================================================================

logrank <- survdiff(Surv(OS_time_days, OS_event) ~ ecosystem_state, data = survival_data)
logrank_p <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)

cat("\nLog-rank p-value:", signif(logrank_p, 4), "\n")

# ======================================================================
# UNIVARIATE COX MODEL
# ======================================================================

cox_fit <- coxph(Surv(OS_time_days, OS_event) ~ ecosystem_state, data = survival_data)
cox_summary <- summary(cox_fit)

cox_result <- data.frame(
  analysis    = "GSE72094 validation (frozen TCGA classifier)",
  n           = cox_summary$n,
  events      = cox_summary$nevent,
  HR          = cox_summary$coefficients[, "exp(coef)"],
  lower_95_CI = cox_summary$conf.int[, "lower .95"],
  upper_95_CI = cox_summary$conf.int[, "upper .95"],
  p_value     = cox_summary$coefficients[, "Pr(>|z|)"],
  logrank_p   = logrank_p,
  stringsAsFactors = FALSE
)

cat("\nCox proportional-hazards result:\n")
print(cox_result)

write_csv(cox_result, file.path(output_dir, "V6_GSE72094_Cox_result.csv"))

# ======================================================================
# PROPORTIONAL-HAZARDS ASSUMPTION CHECK
# ======================================================================

ph_test <- cox.zph(cox_fit)
cat("\nProportional-hazards assumption test:\n")
print(ph_test)

write_csv(as.data.frame(ph_test$table), file.path(output_dir, "V6_GSE72094_PH_test.csv"))

cat("\n=====================================================================\n")
cat("GSE72094 VALIDATION COMPLETE\n")
cat("=====================================================================\n")
cat("If the log-rank p-value is < 0.05 AND the Cox HR direction matches\n")
cat("TCGA's State_1-vs-State_2 direction, this is genuine external\n")
cat("replication of the frozen classifier — not a re-discovered result.\n")
cat("If the PH test above is significant (p < 0.05) for ecosystem_state,\n")
cat("report the HR with that caveat, as the TCGA pipeline does for its\n")
cat("own PH violations.\n")
