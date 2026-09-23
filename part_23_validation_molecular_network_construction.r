# ======================================================================
# M3 — LIONESS PATIENT-SPECIFIC MOLECULAR NETWORK CONSTRUCTION
#
#
# METHOD
# ------
# 1. Read expression matrix
# 2. Map genes using annotation
# 3. Remove duplicated genes
# 4. Remove low-expression / low-variance genes
# 5. Select top variable genes
# 6. Construct cohort-wide Pearson correlation network
# 7. Calculate leave-one-patient-out networks
# 8. Apply LIONESS:
#
#       LIONESS_ij,k =
#           N * r_ij,all
#           -
#           (N-1) * r_ij,-k
#
# 9. Convert each patient's LIONESS network into an edge list
# 10. Apply edge threshold
# 11. Save one CSV network per patient
#
# INPUTS
# -------
# /content/M1_expression_matrix.csv
# /content/02_ProteinCoding_GeneAnnotation.csv
#
# OUTPUTS
# --------
# /content/M3_LIONESS/
#
#   M3_selected_genes.csv
#   M3_cohort_correlation_network.csv
#   M3_network_summary.csv
#   patient_networks/
#       SAMPLE1.csv
#       SAMPLE2.csv
#       ...
#
# ======================================================================


# ======================================================================
# 0. FILE SETTINGS
# ======================================================================

EXPRESSION_FILE <- "/content/GSE72094/M1_expression_matrix.csv"

ANNOTATION_FILE <- "/content/GSE72094/02_ProteinCoding_GeneAnnotation.csv"

OUTPUT_DIR <- "/content/GSE72094_M3_LIONESS"

PATIENT_NETWORK_DIR <- file.path(
  OUTPUT_DIR,
  "patient_networks"
)


# ======================================================================
# 1. PARAMETERS
# ======================================================================

# Number of highly variable genes used for LIONESS.
#
# Recommended starting point:
# 500
#
# You can later test:
# 300
# 500
# 750
# 1000

N_TOP_GENES <- 500


# Minimum number of non-zero expression values required.
MIN_NONZERO_FRACTION <- 0.10


# Minimum standard deviation.
MIN_SD <- 0.05


# Correlation threshold for the cohort network.
#
# This is only used for the cohort reference network.
# LIONESS patient networks use the PATIENT_EDGE_THRESHOLD below.

COHORT_CORRELATION_THRESHOLD <- 0.30


# Patient-specific LIONESS edge threshold.
#
# Edges with:
#
# abs(LIONESS edge) >= 0.30
#
# will be retained.

PATIENT_EDGE_THRESHOLD <- 0.30


# Minimum number of patient-specific edges required
# for a network to be considered usable.

MIN_PATIENT_EDGES <- 10


# Whether to keep positive and negative edges.

KEEP_POSITIVE_EDGES <- TRUE
KEEP_NEGATIVE_EDGES <- TRUE


# ======================================================================
# 2. START
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3 - LIONESS PATIENT-SPECIFIC NETWORK CONSTRUCTION\n")
cat("=====================================================================\n")
cat("\n")

cat("Started:", as.character(Sys.time()), "\n")
cat("\n")


# ======================================================================
# 3. PACKAGES
# ======================================================================

cat("Checking required R packages...\n")

if (!requireNamespace("data.table", quietly = TRUE)) {

  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )

}

suppressPackageStartupMessages({
  library(data.table)
})

cat("Packages ready.\n")


# ======================================================================
# 4. CREATE OUTPUT DIRECTORIES
# ======================================================================

if (!dir.exists(OUTPUT_DIR)) {

  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )

}

if (!dir.exists(PATIENT_NETWORK_DIR)) {

  dir.create(
    PATIENT_NETWORK_DIR,
    recursive = TRUE
  )

}

cat("\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Patient network directory:", PATIENT_NETWORK_DIR, "\n")


# ======================================================================
# 5. CHECK INPUT FILES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("CHECKING INPUT FILES\n")
cat("=====================================================================\n")

required_files <- c(
  EXPRESSION_FILE,
  ANNOTATION_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("\nMissing files:\n")

  print(missing_files)

  stop(
    "One or more required input files are missing."
  )

}

cat("All required files found.\n")


# ======================================================================
# 6. READ EXPRESSION MATRIX
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING EXPRESSION MATRIX\n")
cat("=====================================================================\n")

expr_df <- data.table::fread(
  EXPRESSION_FILE,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "Rows:",
  nrow(expr_df),
  "\n"
)

cat(
  "Columns:",
  ncol(expr_df),
  "\n"
)

if (!("matrix_gene_id" %in% colnames(expr_df))) {

  stop(
    "The expression file must contain a column named 'matrix_gene_id'."
  )

}


# ======================================================================
# 7. IDENTIFY PATIENT COLUMNS
# ======================================================================

expression_samples <- setdiff(
  colnames(expr_df),
  "matrix_gene_id"
)

if (length(expression_samples) < 10) {

  stop(
    "Fewer than 10 expression samples were detected."
  )

}

cat(
  "Number of samples:",
  length(expression_samples),
  "\n"
)


# ======================================================================
# 8. CONVERT EXPRESSION TO NUMERIC
# ======================================================================

cat("\n")
cat("Converting expression values to numeric...\n")

for (sample_name in expression_samples) {

  expr_df[[sample_name]] <- suppressWarnings(
    as.numeric(
      expr_df[[sample_name]]
    )
  )

}


# ======================================================================
# 9. READ ANNOTATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("READING GENE ANNOTATION\n")
cat("=====================================================================\n")

annotation <- data.table::fread(
  ANNOTATION_FILE,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "Annotation rows:",
  nrow(annotation),
  "\n"
)

cat(
  "Annotation columns:",
  ncol(annotation),
  "\n"
)

if (!("matrix_gene_id" %in% colnames(annotation))) {

  stop(
    "Annotation file does not contain 'matrix_gene_id'."
  )

}

if (!("gene_symbol" %in% colnames(annotation))) {

  stop(
    "Annotation file does not contain 'gene_symbol'."
  )

}


# ======================================================================
# 10. CLEAN ANNOTATION
# ======================================================================

annotation$matrix_gene_id <- trimws(
  as.character(
    annotation$matrix_gene_id
  )
)

annotation$gene_symbol <- trimws(
  as.character(
    annotation$gene_symbol
  )
)

annotation$gene_symbol <- toupper(
  annotation$gene_symbol
)


# ======================================================================
# 11. MAP GENES
# ======================================================================

cat("\n")
cat("Mapping expression genes to gene symbols...\n")

expression_gene_ids <- trimws(
  as.character(
    expr_df$matrix_gene_id
  )
)

annotation_index <- match(
  expression_gene_ids,
  annotation$matrix_gene_id
)

mapped_symbols <- annotation$gene_symbol[
  annotation_index
]

matched_count <- sum(
  !is.na(mapped_symbols) &
    mapped_symbols != ""
)

cat(
  "Expression genes:",
  length(expression_gene_ids),
  "\n"
)

cat(
  "Mapped genes:",
  matched_count,
  "\n"
)

expr_df$gene_symbol <- mapped_symbols


# ======================================================================
# 12. REMOVE UNMAPPED GENES
# ======================================================================

expr_df <- expr_df[
  !is.na(expr_df$gene_symbol) &
    expr_df$gene_symbol != "",
  ,
  drop = FALSE
]


# ======================================================================
# 13. REMOVE DUPLICATE GENE SYMBOLS
# ======================================================================

cat("\n")
cat("Removing duplicated gene symbols...\n")

expr_df$gene_symbol <- toupper(
  trimws(
    expr_df$gene_symbol
  )
)

duplicate_symbols <- duplicated(
  expr_df$gene_symbol
)

cat(
  "Duplicate genes:",
  sum(duplicate_symbols),
  "\n"
)

expr_df <- expr_df[
  !duplicate_symbols,
  ,
  drop = FALSE
]

cat(
  "Unique genes:",
  nrow(expr_df),
  "\n"
)


# ======================================================================
# 14. CREATE EXPRESSION MATRIX
# ======================================================================

cat("\n")
cat("Creating expression matrix...\n")

expr_matrix <- as.matrix(
  expr_df[
    ,
    expression_samples,
    drop = FALSE
  ]
)

rownames(expr_matrix) <- expr_df$gene_symbol

storage.mode(expr_matrix) <- "numeric"


# ======================================================================
# 15. REMOVE INVALID VALUES
# ======================================================================

expr_matrix[
  !is.finite(expr_matrix)
] <- 0

expr_matrix[
  expr_matrix < 0
] <- 0


# ======================================================================
# 16. LOG TRANSFORMATION
# ======================================================================

cat("\n")
cat("Applying log2(x + 1) transformation...\n")

expr_matrix <- log2(
  expr_matrix + 1
)


# ======================================================================
# 17. FILTER LOW-EXPRESSION GENES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("FILTERING LOW-EXPRESSION GENES\n")
cat("=====================================================================\n")

nonzero_fraction <- rowMeans(
  expr_matrix > 0
)

keep_expression <- (
  nonzero_fraction >=
    MIN_NONZERO_FRACTION
)

cat(
  "Genes before expression filtering:",
  nrow(expr_matrix),
  "\n"
)

cat(
  "Genes retained:",
  sum(keep_expression),
  "\n"
)

expr_matrix <- expr_matrix[
  keep_expression,
  ,
  drop = FALSE
]


# ======================================================================
# 18. VARIANCE FILTERING
# ======================================================================

cat("\n")
cat("Calculating gene variability...\n")

gene_sd <- apply(
  expr_matrix,
  1,
  sd,
  na.rm = TRUE
)

keep_sd <- (
  is.finite(gene_sd) &
    gene_sd >= MIN_SD
)

cat(
  "Genes after SD filtering:",
  sum(keep_sd),
  "\n"
)

expr_matrix <- expr_matrix[
  keep_sd,
  ,
  drop = FALSE
]

gene_sd <- gene_sd[
  keep_sd
]


# ======================================================================
# 19. SELECT TOP VARIABLE GENES
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("SELECTING HIGHLY VARIABLE GENES\n")
cat("=====================================================================\n")

gene_order <- order(
  gene_sd,
  decreasing = TRUE
)

N_TOP_GENES_ACTUAL <- min(
  N_TOP_GENES,
  length(gene_order)
)

selected_indices <- gene_order[
  seq_len(
    N_TOP_GENES_ACTUAL
  )
]

selected_genes <- rownames(
  expr_matrix
)[
  selected_indices
]

selected_sd <- gene_sd[
  selected_indices
]

expr_selected <- expr_matrix[
  selected_genes,
  ,
  drop = FALSE
]

cat(
  "Requested genes:",
  N_TOP_GENES,
  "\n"
)

cat(
  "Genes selected:",
  nrow(expr_selected),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_selected),
  "\n"
)


# ======================================================================
# 20. SAVE SELECTED GENES
# ======================================================================

selected_gene_table <- data.frame(
  gene = selected_genes,
  standard_deviation = selected_sd,
  rank = seq_along(selected_genes),
  stringsAsFactors = FALSE
)

write.csv(
  selected_gene_table,
  file.path(
    OUTPUT_DIR,
    "M3_selected_genes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Saved selected gene list.\n"
)


# ======================================================================
# 21. TRANSPOSE EXPRESSION
# ======================================================================

# Rows    = patients
# Columns = genes

X <- t(
  expr_selected
)

N <- nrow(X)

G <- ncol(X)

cat("\n")
cat(
  "LIONESS matrix:",
  N,
  "patients x",
  G,
  "genes\n"
)


# ======================================================================
# 22. CHECK SAMPLE NUMBER
# ======================================================================

if (N < 20) {

  stop(
    "Too few patients for reliable cohort-level co-expression analysis."
  )

}


# ======================================================================
# 23. REMOVE ZERO-VARIANCE COLUMNS
# ======================================================================

column_sd <- apply(
  X,
  2,
  sd,
  na.rm = TRUE
)

valid_columns <- (
  is.finite(column_sd) &
    column_sd > 0
)

if (!all(valid_columns)) {

  cat(
    "Removing",
    sum(!valid_columns),
    "zero-variance genes.\n"
  )

  X <- X[
    ,
    valid_columns,
    drop = FALSE
  ]

  selected_genes <- selected_genes[
    valid_columns
  ]

}

G <- ncol(X)


# ======================================================================
# 24. COHORT CORRELATION MATRIX
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("CALCULATING COHORT-LEVEL CORRELATION NETWORK\n")
cat("=====================================================================\n")

cat(
  "Genes:",
  G,
  "\n"
)

cat(
  "Patients:",
  N,
  "\n"
)

cat(
  "Calculating Pearson correlation matrix...\n"
)

cohort_cor <- suppressWarnings(
  cor(
    X,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
)

cohort_cor[
  !is.finite(cohort_cor)
] <- 0

diag(
  cohort_cor
) <- 0


# ======================================================================
# 25. SAVE COHORT NETWORK
# ======================================================================

cat("\n")
cat("Creating cohort correlation edge list...\n")

cohort_edge_indices <- which(
  upper.tri(cohort_cor) &
    abs(cohort_cor) >=
      COHORT_CORRELATION_THRESHOLD,
  arr.ind = TRUE
)

if (nrow(cohort_edge_indices) > 0) {

  cohort_edges <- data.frame(

    gene1 = selected_genes[
      cohort_edge_indices[, 1]
    ],

    gene2 = selected_genes[
      cohort_edge_indices[, 2]
    ],

    correlation = cohort_cor[
      cohort_edge_indices
    ],

    stringsAsFactors = FALSE

  )

  cohort_edges <- cohort_edges[
    order(
      -abs(
        cohort_edges$correlation
      )
    ),
    ,
    drop = FALSE
  ]

} else {

  cohort_edges <- data.frame(
    gene1 = character(0),
    gene2 = character(0),
    correlation = numeric(0),
    stringsAsFactors = FALSE
  )

}

write.csv(
  cohort_edges,
  file.path(
    OUTPUT_DIR,
    "M3_cohort_correlation_network.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Cohort edges retained:",
  nrow(cohort_edges),
  "\n"
)


# ======================================================================
# 26. PREPARE LIONESS CALCULATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("PREPARING LIONESS CALCULATION\n")
cat("=====================================================================\n")

cat(
  "Formula:\n"
)

cat(
  "LIONESS = N * r_all - (N - 1) * r_minus_patient\n"
)

cat(
  "Patients:",
  N,
  "\n"
)

cat(
  "Genes:",
  G,
  "\n"
)

cat(
  "Possible edges:",
  G * (G - 1) / 2,
  "\n"
)


# ======================================================================
# 27. CREATE SAMPLE NAMES
# ======================================================================

sample_ids <- rownames(
  X
)

if (is.null(sample_ids)) {

  sample_ids <- expression_samples

}

if (length(sample_ids) != N) {

  sample_ids <- expression_samples[
    seq_len(N)
  ]

}


# ======================================================================
# 28. PROGRESS SETTINGS
# ======================================================================

progress_step <- max(
  1,
  floor(
    N / 10
  )
)


# ======================================================================
# 29. NETWORK SUMMARY STORAGE
# ======================================================================

network_summary <- data.frame(

  sample_id = sample_ids,

  n_genes = integer(N),

  n_edges = integer(N),

  network_density = numeric(N),

  mean_abs_edge = numeric(N),

  mean_signed_edge = numeric(N),

  positive_edges = integer(N),

  negative_edges = integer(N),

  stringsAsFactors = FALSE

)


# ======================================================================
# 30. LIONESS PATIENT NETWORKS
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("CONSTRUCTING PATIENT-SPECIFIC LIONESS NETWORKS\n")
cat("=====================================================================\n")

for (k in seq_len(N)) {

  patient_id <- sample_ids[k]

  cat(
    "\nProcessing patient",
    k,
    "/",
    N,
    ":",
    patient_id,
    "\n"
  )


  # --------------------------------------------------------------------
  # Remove patient k
  # --------------------------------------------------------------------

  X_minus <- X[
    -k,
    ,
    drop = FALSE
  ]


  # --------------------------------------------------------------------
  # Calculate leave-one-out correlation
  # --------------------------------------------------------------------

  correlation_minus <- suppressWarnings(
    cor(
      X_minus,
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  )

  correlation_minus[
    !is.finite(
      correlation_minus
    )
  ] <- 0

  diag(
    correlation_minus
  ) <- 0


  # --------------------------------------------------------------------
  # LIONESS equation
  # --------------------------------------------------------------------

  lioness_matrix <- (
    N * cohort_cor
  ) -
    (
      (N - 1) *
        correlation_minus
    )


  # --------------------------------------------------------------------
  # Numerical protection
  # --------------------------------------------------------------------

  lioness_matrix[
    !is.finite(
      lioness_matrix
    )
  ] <- 0

  diag(
    lioness_matrix
  ) <- 0


  # --------------------------------------------------------------------
  # Extract upper triangle
  # --------------------------------------------------------------------

  edge_indices <- which(
    upper.tri(
      lioness_matrix
    ),
    arr.ind = TRUE
  )


  edge_values <- lioness_matrix[
    edge_indices
  ]


  # --------------------------------------------------------------------
  # Create edge table
  # --------------------------------------------------------------------

  patient_edges <- data.frame(

    gene1 = selected_genes[
      edge_indices[, 1]
    ],

    gene2 = selected_genes[
      edge_indices[, 2]
    ],

    lioness_weight = edge_values,

    stringsAsFactors = FALSE

  )


  # --------------------------------------------------------------------
  # Remove weak edges
  # --------------------------------------------------------------------

  patient_edges <- patient_edges[
    is.finite(
      patient_edges$lioness_weight
    ),
    ,
    drop = FALSE
  ]

  patient_edges <- patient_edges[
    abs(
      patient_edges$lioness_weight
    ) >= PATIENT_EDGE_THRESHOLD,
    ,
    drop = FALSE
  ]


  # --------------------------------------------------------------------
  # Keep positive / negative edges according to settings
  # --------------------------------------------------------------------

  if (!KEEP_POSITIVE_EDGES) {

    patient_edges <- patient_edges[
      patient_edges$lioness_weight < 0,
      ,
      drop = FALSE
    ]

  }

  if (!KEEP_NEGATIVE_EDGES) {

    patient_edges <- patient_edges[
      patient_edges$lioness_weight > 0,
      ,
      drop = FALSE
    ]

  }


  # --------------------------------------------------------------------
  # Sort edges
  # --------------------------------------------------------------------

  if (nrow(patient_edges) > 0) {

    patient_edges <- patient_edges[
      order(
        -abs(
          patient_edges$lioness_weight
        )
      ),
      ,
      drop = FALSE
    ]

  }


  # --------------------------------------------------------------------
  # Network statistics
  # --------------------------------------------------------------------

  n_edges_patient <- nrow(
    patient_edges
  )

  possible_edges <- (
    G * (G - 1)
  ) / 2

  if (n_edges_patient > 0) {

    density_patient <- (
      n_edges_patient /
        possible_edges
    )

    mean_abs_weight <- mean(
      abs(
        patient_edges$lioness_weight
      ),
      na.rm = TRUE
    )

    mean_signed_weight <- mean(
      patient_edges$lioness_weight,
      na.rm = TRUE
    )

    positive_count <- sum(
      patient_edges$lioness_weight > 0
    )

    negative_count <- sum(
      patient_edges$lioness_weight < 0
    )

  } else {

    density_patient <- 0

    mean_abs_weight <- 0

    mean_signed_weight <- 0

    positive_count <- 0

    negative_count <- 0

  }


  # --------------------------------------------------------------------
  # Save patient network
  # --------------------------------------------------------------------

  safe_patient_id <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    patient_id
  )

  patient_file <- file.path(
    PATIENT_NETWORK_DIR,
    paste0(
      safe_patient_id,
      "_LIONESS_network.csv"
    )
  )


  # --------------------------------------------------------------------
  # Add patient ID
  # --------------------------------------------------------------------

  if (nrow(patient_edges) > 0) {

    patient_edges$sample_id <- patient_id

    patient_edges <- patient_edges[
      ,
      c(
        "sample_id",
        "gene1",
        "gene2",
        "lioness_weight"
      ),
      drop = FALSE
    ]

  } else {

    patient_edges <- data.frame(

      sample_id = character(0),

      gene1 = character(0),

      gene2 = character(0),

      lioness_weight = numeric(0),

      stringsAsFactors = FALSE

    )

  }


  write.csv(
    patient_edges,
    patient_file,
    row.names = FALSE,
    quote = FALSE
  )


  # --------------------------------------------------------------------
  # Store summary
  # --------------------------------------------------------------------

  network_summary$n_genes[k] <- G

  network_summary$n_edges[k] <- n_edges_patient

  network_summary$network_density[k] <- density_patient

  network_summary$mean_abs_edge[k] <- mean_abs_weight

  network_summary$mean_signed_edge[k] <- mean_signed_weight

  network_summary$positive_edges[k] <- positive_count

  network_summary$negative_edges[k] <- negative_count


  # --------------------------------------------------------------------
  # Progress
  # --------------------------------------------------------------------

  cat(
    "  Edges retained:",
    n_edges_patient,
    "\n"
  )

  cat(
    "  Density:",
    round(
      density_patient,
      6
    ),
    "\n"
  )

  if (
    k == 1 ||
      k == N ||
      k %% progress_step == 0
  ) {

    cat(
      "Completed:",
      k,
      "/",
      N,
      "patients\n"
    )

  }

}


# ======================================================================
# 31. SAVE NETWORK SUMMARY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("SAVING NETWORK SUMMARY\n")
cat("=====================================================================\n")

write.csv(
  network_summary,
  file.path(
    OUTPUT_DIR,
    "M3_network_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Saved:",
  file.path(
    OUTPUT_DIR,
    "M3_network_summary.csv"
  ),
  "\n"
)


# ======================================================================
# 32. QUALITY CONTROL
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("NETWORK QUALITY CONTROL\n")
cat("=====================================================================\n")

cat(
  "Number of patients:",
  N,
  "\n"
)

cat(
  "Genes per network:",
  G,
  "\n"
)

cat(
  "LIONESS edge threshold:",
  PATIENT_EDGE_THRESHOLD,
  "\n"
)

cat(
  "Minimum patient edges:",
  MIN_PATIENT_EDGES,
  "\n"
)

cat(
  "Mean number of edges:",
  mean(
    network_summary$n_edges,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median number of edges:",
  median(
    network_summary$n_edges,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Minimum number of edges:",
  min(
    network_summary$n_edges,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Maximum number of edges:",
  max(
    network_summary$n_edges,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Mean network density:",
  mean(
    network_summary$network_density,
    na.rm = TRUE
  ),
  "\n"
)


# ======================================================================
# 33. IDENTIFY LOW-EDGE NETWORKS
# ======================================================================

low_edge_patients <- network_summary[
  network_summary$n_edges <
    MIN_PATIENT_EDGES,
  ,
  drop = FALSE
]

cat("\n")

cat(
  "Patients with fewer than",
  MIN_PATIENT_EDGES,
  "edges:",
  nrow(low_edge_patients),
  "\n"
)

if (nrow(low_edge_patients) > 0) {

  print(
    low_edge_patients
  )

  write.csv(
    low_edge_patients,
    file.path(
      OUTPUT_DIR,
      "M3_low_edge_patients.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )

}


# ======================================================================
# 34. EDGE DISTRIBUTION
# ======================================================================

cat("\n")
cat("Positive/negative edge statistics:\n")

cat(
  "Mean positive edges:",
  mean(
    network_summary$positive_edges,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Mean negative edges:",
  mean(
    network_summary$negative_edges,
    na.rm = TRUE
  ),
  "\n"
)


# ======================================================================
# 35. FINAL OUTPUT INFORMATION
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("OUTPUT FILES\n")
cat("=====================================================================\n")

cat(
  "1. Selected genes:\n   ",
  file.path(
    OUTPUT_DIR,
    "M3_selected_genes.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "2. Cohort correlation network:\n   ",
  file.path(
    OUTPUT_DIR,
    "M3_cohort_correlation_network.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "3. Network summary:\n   ",
  file.path(
    OUTPUT_DIR,
    "M3_network_summary.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "4. Patient-specific networks:\n   ",
  PATIENT_NETWORK_DIR,
  "\n",
  sep = ""
)


# ======================================================================
# 36. SHOW EXAMPLE NETWORK
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("EXAMPLE PATIENT NETWORK\n")
cat("=====================================================================\n")

example_patient <- network_summary$sample_id[1]

example_safe_id <- gsub(
  "[^A-Za-z0-9_.-]",
  "_",
  example_patient
)

example_file <- file.path(
  PATIENT_NETWORK_DIR,
  paste0(
    example_safe_id,
    "_LIONESS_network.csv"
  )
)

if (file.exists(example_file)) {

  example_network <- read.csv(
    example_file,
    stringsAsFactors = FALSE
  )

  cat(
    "Example patient:",
    example_patient,
    "\n"
  )

  cat(
    "Number of edges:",
    nrow(example_network),
    "\n"
  )

  if (nrow(example_network) > 0) {

    print(
      head(
        example_network,
        10
      )
    )

  }

}


# ======================================================================
# 37. FINISHED
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("M3 LIONESS COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n")

cat("\n")

cat(
  "The resulting networks are PATIENT-SPECIFIC.\n"
)

cat(
  "No STRING/PPI interaction database was used.\n"
)

cat(
  "Each patient has an independently estimated LIONESS edge profile.\n"
)

cat(
  "These networks can now be used for patient-specific entropy analysis.\n"
)

cat("\n")

cat(
  "Finished:",
  as.character(Sys.time()),
  "\n"
)

cat("\n")
cat("=====================================================================\n")
cat("END OF M3\n")
cat("=====================================================================\n")