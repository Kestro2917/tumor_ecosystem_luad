# ======================================================================
# M3 — LIONESS PATIENT-SPECIFIC MULTISCALE NETWORK ENTROPY
#
# PURPOSE
# -------
# Calculate patient-specific molecular network entropy using the
# ALREADY-GENERATED LIONESS patient-specific networks.
#
# THREE SCALES
# ------------
# 1. LOCAL
#    Gene/node-level entropy within each patient's network.
#
# 2. MESO
#    MSigDB module/category-level entropy.
#    Each module is defined by genes belonging to the corresponding
#    MSigDB category from the M2b libraries.
#
# 3. GLOBAL
#    Whole-network entropy for each individual patient.
#
#
# INPUTS
# -------
#
# 1. Already-generated LIONESS patient-specific networks.
#
# 2. M2b MSigDB libraries:
#
#    /content/M2_MSigDB_LUAD/M2_CellType_MSigDB.csv
#    /content/M2_MSigDB_LUAD/M2_CellState_MSigDB.csv
#    /content/M2_MSigDB_LUAD/M2_TumorStroma_MSigDB.csv
#
# 3. Gene annotation:
#
#    /content/02_ProteinCoding_GeneAnnotation.csv
#
#
# EXPECTED LIONESS FORMAT
# -----------------------
#
# The script attempts to support common LIONESS formats:
#
# A) One CSV per patient:
#
#    /content/LIONESS_networks/
#       patient1.csv
#       patient2.csv
#       patient3.csv
#
# Each file should contain:
#
#    gene1
#    gene2
#    weight
#
#
# B) A single directory containing patient-specific edge files.
#
# The patient ID is inferred from the filename.
#
#
# OUTPUTS
# -------
#
# M3_LIONESS_entropy.csv
#
#    One row per patient containing:
#
#    sample_id
#    S_global
#    S_local_mean
#    S_local_sd
#    S_local_median
#    S_mesoscale_mean
#    S_mesoscale_sd
#    S_mesoscale_median
#    S_module_<category>
#
#
# M3_LIONESS_local_entropy.csv
#
#    Patient x gene local entropy matrix.
#
#
# M3_LIONESS_module_entropy.csv
#
#    Patient x MSigDB module/category entropy matrix.
#
#
# M3_LIONESS_network_summary.csv
#
#    Patient-specific network statistics.
#
#
# M3_LIONESS_module_definitions.csv
#
#    MSigDB category -> gene mapping.
#
#
# M3_LIONESS_run_summary.csv
#
#    Run-level quality-control summary.
#
# ======================================================================


# ======================================================================
# 0. FILE SETTINGS
# ======================================================================

# ----------------------------------------------------------------------
# Directory containing ALREADY-GENERATED LIONESS networks
# ----------------------------------------------------------------------

LIONESS_DIR <- "/content/GSE72094_M3_LIONESS/patient_networks"


# ----------------------------------------------------------------------
# M2b MSigDB directory
# ----------------------------------------------------------------------

MSIGDB_DIR <- "/content/M2_MSigDB_LUAD"


# ----------------------------------------------------------------------
# Gene annotation
# ----------------------------------------------------------------------

ANNOTATION_FILE <- "/content/GSE72094/02_ProteinCoding_GeneAnnotation.csv"


# ----------------------------------------------------------------------
# Output directory
# ----------------------------------------------------------------------

OUTPUT_DIR <- "/content/GSE72094_M3_LIONESS_entropy"


# ----------------------------------------------------------------------
# M2b MSigDB libraries
# ----------------------------------------------------------------------

LIBRARY_FILES <- c(

  file.path(
    MSIGDB_DIR,
    "M2_CellType_MSigDB.csv"
  ),

  file.path(
    MSIGDB_DIR,
    "M2_CellState_MSigDB.csv"
  ),

  file.path(
    MSIGDB_DIR,
    "M2_TumorStroma_MSigDB.csv"
  )

)


# ----------------------------------------------------------------------
# Minimum genes required for a module
# ----------------------------------------------------------------------

MIN_MODULE_GENES <- 3


# ----------------------------------------------------------------------
# Minimum genes required in a patient network
# ----------------------------------------------------------------------

MIN_NETWORK_GENES <- 20


# ----------------------------------------------------------------------
# Whether to normalize local entropy by log(degree)
# ----------------------------------------------------------------------

NORMALIZE_LOCAL_ENTROPY <- TRUE


# ======================================================================
# 1. OUTPUT FILES
# ======================================================================

ENTROPY_FILE <- file.path(
  OUTPUT_DIR,
  "M3_LIONESS_entropy.csv"
)

LOCAL_FILE <- file.path(
  OUTPUT_DIR,
  "M3_LIONESS_local_entropy.csv"
)

MODULE_FILE <- file.path(
  OUTPUT_DIR,
  "M3_LIONESS_module_entropy.csv"
)

NETWORK_SUMMARY_FILE <- file.path(
  OUTPUT_DIR,
  "M3_LIONESS_network_summary.csv"
)

MODULE_DEFINITION_FILE <- file.path(
  OUTPUT_DIR,
  "M3_LIONESS_module_definitions.csv"
)

RUN_SUMMARY_FILE <- file.path(
  OUTPUT_DIR,
  "M3_LIONESS_run_summary.csv"
)


# ======================================================================
# 2. START
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("M3 — LIONESS PATIENT-SPECIFIC NETWORK ENTROPY\n")
cat("============================================================\n")
cat("\n")

cat(
  "Started:",
  as.character(Sys.time()),
  "\n\n"
)


# ======================================================================
# 3. CREATE OUTPUT DIRECTORY
# ======================================================================

if (!dir.exists(OUTPUT_DIR)) {

  dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )

}


# ======================================================================
# 4. INSTALL / LOAD PACKAGES
# ======================================================================

cat("Checking required R packages...\n")


if (!requireNamespace(
  "data.table",
  quietly = TRUE
)) {

  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )

}


if (!requireNamespace(
  "igraph",
  quietly = TRUE
)) {

  install.packages(
    "igraph",
    repos = "https://cloud.r-project.org"
  )

}


suppressPackageStartupMessages({

  library(data.table)

  library(igraph)

})


cat("Packages ready.\n")


# ======================================================================
# 5. CHECK INPUT DIRECTORIES / FILES
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("CHECKING INPUTS\n")
cat("============================================================\n")


if (!dir.exists(LIONESS_DIR)) {

  stop(
    paste0(
      "LIONESS directory does not exist:\n",
      LIONESS_DIR,
      "\n\n",
      "Please change LIONESS_DIR to the directory containing ",
      "your already-generated patient-specific networks."
    )
  )

}


missing_libraries <- LIBRARY_FILES[
  !file.exists(LIBRARY_FILES)
]


if (length(missing_libraries) > 0) {

  cat("\nMissing MSigDB library files:\n")

  print(missing_libraries)

  stop(
    "One or more M2b MSigDB libraries are missing."
  )

}


if (!file.exists(ANNOTATION_FILE)) {

  stop(
    paste0(
      "Gene annotation file not found:\n",
      ANNOTATION_FILE
    )
  )

}


cat(
  "LIONESS directory found:",
  LIONESS_DIR,
  "\n"
)

cat(
  "MSigDB libraries found:",
  length(LIBRARY_FILES),
  "\n"
)

cat(
  "Annotation file found.\n"
)


# ======================================================================
# 6. FIND LIONESS NETWORK FILES
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("DISCOVERING LIONESS NETWORK FILES\n")
cat("============================================================\n")


lioness_files <- list.files(

  LIONESS_DIR,

  pattern = "\\.(csv|CSV|tsv|TSV)$",

  full.names = TRUE,

  recursive = TRUE

)


if (length(lioness_files) == 0) {

  stop(
    paste0(
      "No CSV/TSV LIONESS network files were found in:\n",
      LIONESS_DIR
    )
  )

}


cat(
  "LIONESS network files found:",
  length(lioness_files),
  "\n"
)


print(
  head(
    basename(lioness_files),
    10
  )
)


if (length(lioness_files) > 10) {

  cat(
    "... and",
    length(lioness_files) - 10,
    "more files.\n"
  )

}


# ======================================================================
# 7. READ GENE ANNOTATION
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("READING GENE ANNOTATION\n")
cat("============================================================\n")


annotation <- fread(

  ANNOTATION_FILE,

  data.table = FALSE,

  check.names = FALSE

)


if (!(
  "matrix_gene_id" %in%
    colnames(annotation)
)) {

  stop(
    "matrix_gene_id column missing from annotation."
  )

}


if (!(
  "gene_symbol" %in%
    colnames(annotation)
)) {

  stop(
    "gene_symbol column missing from annotation."
  )

}


annotation$matrix_gene_id <- trimws(
  as.character(
    annotation$matrix_gene_id
  )
)


annotation$gene_symbol <- toupper(
  trimws(
    as.character(
      annotation$gene_symbol
    )
  )
)


annotation <- annotation[
  !is.na(annotation$gene_symbol) &
    annotation$gene_symbol != "",
  ,
  drop = FALSE
]


cat(
  "Annotation genes:",
  nrow(annotation),
  "\n"
)


# ======================================================================
# 8. READ MSigDB LIBRARIES
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("READING MSigDB MODULE LIBRARIES\n")
cat("============================================================\n")


all_modules <- data.frame(

  library = character(),

  signature = character(),

  category = character(),

  gene = character(),

  stringsAsFactors = FALSE

)


for (lib_file in LIBRARY_FILES) {

  cat(
    "\nReading:",
    basename(lib_file),
    "\n"
  )


  lib <- fread(

    lib_file,

    data.table = FALSE,

    check.names = FALSE

  )


  required_columns <- c(

    "signature",

    "category",

    "gene"

  )


  missing_columns <- setdiff(

    required_columns,

    colnames(lib)

  )


  if (length(missing_columns) > 0) {

    stop(
      paste0(
        "Library ",
        basename(lib_file),
        " is missing: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )

  }


  lib$signature <- trimws(
    as.character(
      lib$signature
    )
  )


  lib$category <- trimws(
    as.character(
      lib$category
    )
  )


  lib$gene <- toupper(
    trimws(
      as.character(
        lib$gene
      )
    )
  )


  lib$library <- tools::file_path_sans_ext(
    basename(
      lib_file
    )
  )


  lib <- lib[
    !is.na(lib$gene) &
      lib$gene != "" &
      !is.na(lib$category) &
      lib$category != "",
    ,
    drop = FALSE
  ]


  all_modules <- rbind(

    all_modules,

    lib[
      ,
      c(
        "library",
        "signature",
        "category",
        "gene"
      ),
      drop = FALSE
    ]

  )

}


all_modules <- unique(
  all_modules
)


cat(
  "\nTotal MSigDB records:",
  nrow(all_modules),
  "\n"
)


cat(
  "Libraries:",
  length(
    unique(
      all_modules$library
    )
  ),
  "\n"
)


cat(
  "Signatures:",
  length(
    unique(
      all_modules$signature
    )
  ),
  "\n"
)


cat(
  "Categories/modules:",
  length(
    unique(
      all_modules$category
    )
  ),
  "\n"
)


# ======================================================================
# 9. DETERMINE MODULE GENES
#
# IMPORTANT
# ---------
# A "module" here means an MSigDB CATEGORY.
#
# Example:
#
# category = Macrophage
#
# all genes belonging to the MSigDB signatures assigned to the
# Macrophage category become the Macrophage module.
#
# This gives the MESO SCALE.
# ======================================================================


cat("\n")
cat("============================================================\n")
cat("BUILDING MSigDB MODULES\n")
cat("============================================================\n")


module_names <- sort(
  unique(
    all_modules$category
  )
)


module_gene_list <- list()


module_definition_table <- data.frame(

  module = character(),

  genes_total = integer(),

  genes_in_annotation = integer(),

  stringsAsFactors = FALSE

)


for (module_name in module_names) {

  genes <- unique(

    all_modules$gene[
      all_modules$category ==
        module_name
    ]

  )


  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]


  genes_annotation <- intersect(

    genes,

    annotation$gene_symbol

  )


  if (
    length(genes_annotation) >=
      MIN_MODULE_GENES
  ) {

    module_gene_list[[module_name]] <-
      genes_annotation


    module_definition_table <-
      rbind(

        module_definition_table,

        data.frame(

          module = module_name,

          genes_total = length(genes),

          genes_in_annotation =
            length(
              genes_annotation
            ),

          stringsAsFactors =
            FALSE

        )

      )

  }

}


cat(
  "Usable MSigDB modules:",
  length(module_gene_list),
  "\n"
)


if (
  length(module_gene_list) == 0
) {

  stop(
    "No MSigDB modules passed the minimum gene threshold."
  )

}


# ======================================================================
# 10. SAVE MODULE DEFINITIONS
# ======================================================================

module_definition_long <- do.call(

  rbind,

  lapply(

    names(module_gene_list),

    function(m) {

      data.frame(

        module = m,

        gene = module_gene_list[[m]],

        stringsAsFactors = FALSE

      )

    }

  )

)


module_definition_long <- merge(

  module_definition_long,

  module_definition_table,

  by = "module",

  all.x = TRUE

)


module_definition_long <-
  module_definition_long[
    order(
      module_definition_long$module,
      module_definition_long$gene
    ),
    ,
    drop = FALSE
  ]


write.csv(

  module_definition_long,

  MODULE_DEFINITION_FILE,

  row.names = FALSE,

  quote = FALSE

)


cat(
  "Saved:",
  MODULE_DEFINITION_FILE,
  "\n"
)


# ======================================================================
# 11. HELPER FUNCTION
#     IDENTIFY PATIENT ID FROM LIONESS FILE
# ======================================================================

get_patient_id <- function(file_path) {

  patient_id <- tools::file_path_sans_ext(
    basename(file_path)
  )


  # Remove common LIONESS naming prefixes

  patient_id <- sub(
    "^LIONESS_",
    "",
    patient_id,
    ignore.case = TRUE
  )


  patient_id <- sub(
    "_network$",
    "",
    patient_id,
    ignore.case = TRUE
  )


  patient_id <- sub(
    "_network_edges$",
    "",
    patient_id,
    ignore.case = TRUE
  )


  patient_id <- trimws(
    patient_id
  )


  patient_id

}


# ======================================================================
# 12. HELPER FUNCTION
#     READ ONE LIONESS NETWORK
# ======================================================================

read_lioness_network <- function(file_path) {

  cat(
    "\nReading:",
    basename(file_path),
    "\n"
  )


  network_df <- fread(

    file_path,

    data.table = FALSE,

    check.names = FALSE

  )


  cat(
    "Rows:",
    nrow(network_df),
    "\n"
  )


  cat(
    "Columns:",
    ncol(network_df),
    "\n"
  )


  # ---------------------------------------------------------------
  # Normalize column names
  # ---------------------------------------------------------------

  original_names <- colnames(
    network_df
  )


  normalized_names <- tolower(
    trimws(
      original_names
    )
  )


  normalized_names <- gsub(
    "[ .-]+",
    "_",
    normalized_names
  )


  colnames(network_df) <-
    normalized_names


  # ---------------------------------------------------------------
  # Detect gene 1 column
  # ---------------------------------------------------------------

  gene1_candidates <- c(

    "gene1",

    "gene_1",

    "from",

    "source",

    "node1",

    "node_1"

  )


  gene1_column <- intersect(

    gene1_candidates,

    colnames(network_df)

  )


  # ---------------------------------------------------------------
  # Detect gene 2 column
  # ---------------------------------------------------------------

  gene2_candidates <- c(

    "gene2",

    "gene_2",

    "to",

    "target",

    "node2",

    "node_2"

  )


  gene2_column <- intersect(

    gene2_candidates,

    colnames(network_df)

  )


  # ---------------------------------------------------------------
  # Detect weight column
  # ---------------------------------------------------------------

  weight_candidates <- c(

    "weight",

    "weights",

    "edge_weight",

    "edgeweight",

    "lioness_weight",

    "value",

    "correlation",

    "corr",

    "score"

  )


  weight_column <- intersect(

    weight_candidates,

    colnames(network_df)

  )


  if (
    length(gene1_column) == 0 ||
      length(gene2_column) == 0
  ) {

    stop(
      paste0(
        "\nCould not identify gene1/gene2 columns in:\n",
        file_path,
        "\n\nAvailable columns:\n",
        paste(
          colnames(network_df),
          collapse = ", "
        )
      )
    )

  }


  if (
    length(weight_column) == 0
  ) {

    stop(
      paste0(
        "\nCould not identify edge weight column in:\n",
        file_path,
        "\n\nAvailable columns:\n",
        paste(
          colnames(network_df),
          collapse = ", "
        )
      )
    )

  }


  gene1_column <- gene1_column[1]

  gene2_column <- gene2_column[1]

  weight_column <- weight_column[1]


  network <- data.frame(

    gene1 = toupper(
      trimws(
        as.character(
          network_df[[gene1_column]]
        )
      )
    ),

    gene2 = toupper(
      trimws(
        as.character(
          network_df[[gene2_column]]
        )
      )
    ),

    weight = suppressWarnings(
      as.numeric(
        network_df[[weight_column]]
      )
    ),

    stringsAsFactors = FALSE

  )


  # ---------------------------------------------------------------
  # Remove invalid edges
  # ---------------------------------------------------------------

  network <- network[

    !is.na(network$gene1) &

      !is.na(network$gene2) &

      !is.na(network$weight) &

      is.finite(network$weight) &

      network$gene1 != "" &

      network$gene2 != "" &

      network$gene1 != network$gene2,

    ,

    drop = FALSE

  ]


  # ---------------------------------------------------------------
  # Convert negative weights
  #
  # Entropy requires non-negative transition weights.
  #
  # LIONESS networks may contain positive and negative associations.
  #
  # We therefore use absolute edge strength:
  #
  #       |weight|
  #
  # This preserves interaction strength but removes direction sign.
  # ---------------------------------------------------------------

  network$weight <- abs(
    network$weight
  )


  network <- network[
    network$weight > 0,
    ,
    drop = FALSE
  ]


  # ---------------------------------------------------------------
  # Remove duplicate undirected edges
  # ---------------------------------------------------------------

  network$edge_a <- pmin(

    network$gene1,

    network$gene2

  )


  network$edge_b <- pmax(

    network$gene1,

    network$gene2

  )


  network <- network[
    !duplicated(
      paste(
        network$edge_a,
        network$edge_b,
        sep = "__"
      )
    ),
    ,
    drop = FALSE
  ]


  network <- network[
    ,
    c(
      "edge_a",
      "edge_b",
      "weight"
    ),
    drop = FALSE
  ]


  colnames(network) <- c(

    "gene1",

    "gene2",

    "weight"

  )


  network

}


# ======================================================================
# 13. ENTROPY FUNCTION
# ======================================================================

calculate_lioness_entropy <- function(

  network,

  module_gene_list,

  normalize_local = TRUE

) {


  # ---------------------------------------------------------------
  # Network genes
  # ---------------------------------------------------------------

  network_genes <- sort(

    unique(
      c(
        network$gene1,
        network$gene2
      )
    )

  )


  n_genes <- length(
    network_genes
  )


  if (
    n_genes <
      MIN_NETWORK_GENES
  ) {

    return(NULL)

  }


  # ---------------------------------------------------------------
  # Gene index
  # ---------------------------------------------------------------

  gene_index <- setNames(

    seq_along(network_genes),

    network_genes

  )


  edge_from <- unname(

    gene_index[
      network$gene1
    ]

  )


  edge_to <- unname(

    gene_index[
      network$gene2
    ]

  )


  edge_weights <- as.numeric(
    network$weight
  )


  # ---------------------------------------------------------------
  # Safety
  # ---------------------------------------------------------------

  valid <- (

    is.finite(edge_weights) &

      edge_weights > 0 &

      !is.na(edge_from) &

      !is.na(edge_to)

  )


  edge_from <- edge_from[
    valid
  ]

  edge_to <- edge_to[
    valid
  ]

  edge_weights <- edge_weights[
    valid
  ]


  if (
    length(edge_weights) == 0
  ) {

    return(NULL)

  }


  # ---------------------------------------------------------------
  # Weighted degree
  # ---------------------------------------------------------------

  weighted_degree <- numeric(
    n_genes
  )


  for (i in seq_along(edge_weights)) {

    a <- edge_from[i]

    b <- edge_to[i]

    w <- edge_weights[i]


    weighted_degree[a] <-
      weighted_degree[a] + w


    weighted_degree[b] <-
      weighted_degree[b] + w

  }


  # ---------------------------------------------------------------
  # Unweighted degree
  # ---------------------------------------------------------------

  degree <- numeric(
    n_genes
  )


  for (i in seq_along(edge_weights)) {

    a <- edge_from[i]

    b <- edge_to[i]


    degree[a] <-
      degree[a] + 1


    degree[b] <-
      degree[b] + 1

  }


  names(degree) <-
    network_genes


  # ---------------------------------------------------------------
  # LOCAL ENTROPY
  #
  # For every node:
  #
  # p_ij = w_ij / sum_j(w_ij)
  #
  # H_i = -sum_j p_ij log(p_ij)
  #
  # Optional normalization:
  #
  # H_i / log(k_i)
  #
  # resulting in approximately [0,1].
  # ---------------------------------------------------------------

  local_entropy <- numeric(
    n_genes
  )


  for (i in seq_len(n_genes)) {

    incident_edges <- which(

      edge_from == i |
        edge_to == i

    )


    if (
      length(incident_edges) == 0
    ) {

      local_entropy[i] <- 0

      next

    }


    incident_weights <-
      edge_weights[
        incident_edges
      ]


    total_weight <- sum(
      incident_weights
    )


    if (
      total_weight <= 0 ||
        !is.finite(total_weight)
    ) {

      local_entropy[i] <- 0

      next

    }


    probabilities <-
      incident_weights /
      total_weight


    probabilities <- probabilities[
      is.finite(probabilities) &
        probabilities > 0
    ]


    H <- -sum(

      probabilities *
        log(
          probabilities
        )

    )


    if (
      normalize_local &
        degree[i] > 1
    ) {

      H <- H /
        log(
          degree[i]
        )

    }


    if (
      !is.finite(H)
    ) {

      H <- 0

    }


    local_entropy[i] <- H

  }


  local_entropy[
    !is.finite(local_entropy)
  ] <- 0


  if (normalize_local) {

    local_entropy <- pmin(

      pmax(
        local_entropy,
        0
      ),

      1

    )

  }


  names(local_entropy) <-
    network_genes


  # ---------------------------------------------------------------
  # GLOBAL ENTROPY
  #
  # Degree-weighted mean of local node entropy.
  # ---------------------------------------------------------------

  if (
    sum(weighted_degree) > 0
  ) {

    degree_weights <-
      weighted_degree /
      sum(weighted_degree)


    global_entropy <- sum(

      degree_weights *
        local_entropy

    )

  } else {

    global_entropy <- 0

  }


  # ---------------------------------------------------------------
  # MESO / MODULE ENTROPY
  #
  # For every MSigDB module:
  #
  # module entropy =
  # mean(local entropy of genes in that module)
  #
  # Only genes actually present in the patient's network are used.
  # ---------------------------------------------------------------

  module_entropy <- numeric(

    length(
      module_gene_list
    )

  )


  names(module_entropy) <-
    names(module_gene_list)


  module_gene_counts <- integer(

    length(
      module_gene_list
    )

  )


  names(module_gene_counts) <-
    names(module_gene_list)


  for (
    m in seq_along(module_gene_list)
  ) {

    module_genes <-
      module_gene_list[[m]]


    genes_present <- intersect(

      module_genes,

      network_genes

    )


    module_gene_counts[m] <-
      length(
        genes_present
      )


    if (
      length(genes_present) >=
        MIN_MODULE_GENES
    ) {

      module_entropy[m] <-
        mean(
          local_entropy[
            genes_present
          ],
          na.rm = TRUE
        )

    } else {

      module_entropy[m] <-
        NA_real_

    }

  }


  # ---------------------------------------------------------------
  # LOCAL SUMMARY
  # ---------------------------------------------------------------

  local_mean <- mean(

    local_entropy,

    na.rm = TRUE

  )


  local_sd <- sd(

    local_entropy,

    na.rm = TRUE

  )


  local_median <- median(

    local_entropy,

    na.rm = TRUE

  )


  # ---------------------------------------------------------------
  # MESO SUMMARY
  # ---------------------------------------------------------------

  valid_modules <-
    is.finite(
      module_entropy
    )


  if (
    any(valid_modules)
  ) {

    meso_mean <- mean(

      module_entropy[
        valid_modules
      ],

      na.rm = TRUE

    )


    meso_sd <- sd(

      module_entropy[
        valid_modules
      ],

      na.rm = TRUE

    )


    meso_median <- median(

      module_entropy[
        valid_modules
      ],

      na.rm = TRUE

    )

  } else {

    meso_mean <- NA_real_

    meso_sd <- NA_real_

    meso_median <- NA_real_

  }


  # ---------------------------------------------------------------
  # RETURN
  # ---------------------------------------------------------------

  list(

    local_entropy =
      local_entropy,

    global_entropy =
      global_entropy,

    local_mean =
      local_mean,

    local_sd =
      local_sd,

    local_median =
      local_median,

    module_entropy =
      module_entropy,

    module_gene_counts =
      module_gene_counts,

    meso_mean =
      meso_mean,

    meso_sd =
      meso_sd,

    meso_median =
      meso_median,

    n_genes =
      n_genes,

    n_edges =
      length(
        edge_weights
      ),

    total_edge_weight =
      sum(
        edge_weights
      )

  )

}


# ======================================================================
# 14. INITIALIZE OUTPUT CONTAINERS
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("INITIALIZING PATIENT-SPECIFIC ENTROPY ANALYSIS\n")
cat("============================================================\n")


n_samples <- length(
  lioness_files
)


n_modules <- length(
  module_gene_list
)


all_sample_ids <- character(
  n_samples
)


global_entropy_values <- rep(

  NA_real_,

  n_samples

)


local_mean_values <- rep(

  NA_real_,

  n_samples

)


local_sd_values <- rep(

  NA_real_,

  n_samples

)


local_median_values <- rep(

  NA_real_,

  n_samples

)


meso_mean_values <- rep(

  NA_real_,

  n_samples

)


meso_sd_values <- rep(

  NA_real_,

  n_samples

)


meso_median_values <- rep(

  NA_real_,

  n_samples

)


network_gene_counts <- rep(

  NA_integer_,

  n_samples

)


network_edge_counts <- rep(

  NA_integer_,

  n_samples

)


network_total_weights <- rep(

  NA_real_,

  n_samples

)


module_entropy_matrix <- matrix(

  NA_real_,

  nrow = n_samples,

  ncol = n_modules

)


colnames(
  module_entropy_matrix
) <- paste0(

  "S_module_",

  make.names(
    names(
      module_gene_list
    ),
    unique = TRUE
  )

)


rownames(
  module_entropy_matrix
) <-
  all_sample_ids


# ======================================================================
# 15. LOCAL ENTROPY STORAGE
# ======================================================================

local_entropy_list <- vector(

  "list",

  n_samples

)


# ======================================================================
# 16. NETWORK SUMMARY STORAGE
# ======================================================================

network_summary_list <- vector(

  "list",

  n_samples

)


# ======================================================================
# 17. PROCESS EACH PATIENT NETWORK
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("CALCULATING PATIENT-SPECIFIC ENTROPY\n")
cat("============================================================\n")


progress_step <- max(

  1,

  floor(
    n_samples / 10
  )

)


successful_samples <- 0


failed_samples <- 0


for (
  s in seq_len(n_samples)
) {


  file_path <-
    lioness_files[s]


  sample_id <-
    get_patient_id(
      file_path
    )


  all_sample_ids[s] <-
    sample_id


  cat("\n")
  cat(
    "------------------------------------------------------------\n"
  )


  cat(
    "Patient",
    s,
    "/",
    n_samples,
    ":",
    sample_id,
    "\n"
  )


  # ---------------------------------------------------------------
  # Read network
  # ---------------------------------------------------------------

  network_result <- tryCatch(

    {

      read_lioness_network(
        file_path
      )

    },

    error = function(e) {

      cat(
        "ERROR reading network:",
        conditionMessage(e),
        "\n"
      )

      NULL

    }

  )


  if (
    is.null(network_result)
  ) {

    failed_samples <-

      failed_samples + 1

    next

  }


  # ---------------------------------------------------------------
  # Calculate entropy
  # ---------------------------------------------------------------

  result <- tryCatch(

    {

      calculate_lioness_entropy(

        network =
          network_result,

        module_gene_list =
          module_gene_list,

        normalize_local =
          NORMALIZE_LOCAL_ENTROPY

      )

    },

    error = function(e) {

      cat(
        "ERROR calculating entropy:",
        conditionMessage(e),
        "\n"
      )

      NULL

    }

  )


  if (
    is.null(result)
  ) {

    failed_samples <-

      failed_samples + 1

    next

  }


  # ---------------------------------------------------------------
  # Store global entropy
  # ---------------------------------------------------------------

  global_entropy_values[s] <-
    result$global_entropy


  # ---------------------------------------------------------------
  # Store local summaries
  # ---------------------------------------------------------------

  local_mean_values[s] <-
    result$local_mean


  local_sd_values[s] <-
    result$local_sd


  local_median_values[s] <-
    result$local_median


  # ---------------------------------------------------------------
  # Store meso summaries
  # ---------------------------------------------------------------

  meso_mean_values[s] <-
    result$meso_mean


  meso_sd_values[s] <-
    result$meso_sd


  meso_median_values[s] <-
    result$meso_median


  # ---------------------------------------------------------------
  # Store network statistics
  # ---------------------------------------------------------------

  network_gene_counts[s] <-
    result$n_genes


  network_edge_counts[s] <-
    result$n_edges


  network_total_weights[s] <-
    result$total_edge_weight


  # ---------------------------------------------------------------
  # Store module entropy
  # ---------------------------------------------------------------

  module_entropy_matrix[s, ] <-

    as.numeric(
      result$module_entropy
    )


  # ---------------------------------------------------------------
  # Store local entropy
  # ---------------------------------------------------------------

  local_entropy_list[[s]] <-

    result$local_entropy


  # ---------------------------------------------------------------
  # Network summary
  # ---------------------------------------------------------------

  network_summary_list[[s]] <-

    data.frame(

      sample_id =
        sample_id,

      network_genes =
        result$n_genes,

      network_edges =
        result$n_edges,

      total_edge_weight =
        result$total_edge_weight,

      S_global =
        result$global_entropy,

      S_local_mean =
        result$local_mean,

      S_local_sd =
        result$local_sd,

      S_local_median =
        result$local_median,

      S_meso_mean =
        result$meso_mean,

      S_meso_sd =
        result$meso_sd,

      S_meso_median =
        result$meso_median,

      stringsAsFactors =
        FALSE

    )


  successful_samples <-

    successful_samples + 1


  # ---------------------------------------------------------------
  # Progress
  # ---------------------------------------------------------------

  if (

    s == 1 ||

      s == n_samples ||

      s %% progress_step == 0

  ) {

    cat(
      "Completed:",
      s,
      "/",
      n_samples,
      "\n"
    )

    cat(
      "Global entropy:",
      round(
        result$global_entropy,
        6
      ),
      "\n"
    )

  }

}


# ======================================================================
# 18. UPDATE MODULE MATRIX ROW NAMES
# ======================================================================

rownames(
  module_entropy_matrix
) <-
  all_sample_ids


# ======================================================================
# 19. CREATE LOCAL ENTROPY MATRIX
# ======================================================================

cat("\n")
cat("Creating local entropy matrix...\n")


all_network_genes <- sort(

  unique(

    unlist(

      lapply(

        local_entropy_list,

        names

      )

    )

  )

)


local_entropy_matrix <- matrix(

  NA_real_,

  nrow = n_samples,

  ncol = length(
    all_network_genes
  )

)


colnames(
  local_entropy_matrix
) <-
  all_network_genes


rownames(
  local_entropy_matrix
) <-
  all_sample_ids


for (
  s in seq_len(n_samples)
) {

  if (
    !is.null(
      local_entropy_list[[s]]
    )
  ) {

    genes <- names(
      local_entropy_list[[s]]
    )


    local_entropy_matrix[
      s,
      genes
    ] <-
      local_entropy_list[[s]]

  }

}


# ======================================================================
# 20. SAVE LOCAL ENTROPY
# ======================================================================

local_entropy_df <-

  as.data.frame(

    local_entropy_matrix,

    check.names = FALSE

  )


local_entropy_df$sample_id <-
  rownames(
    local_entropy_df
  )


local_entropy_df <-

  local_entropy_df[

    ,

    c(

      "sample_id",

      setdiff(

        colnames(
          local_entropy_df
        ),

        "sample_id"

      )

    ),

    drop = FALSE

  ]


write.csv(

  local_entropy_df,

  LOCAL_FILE,

  row.names = FALSE,

  quote = FALSE

)


cat(
  "Saved:",
  LOCAL_FILE,
  "\n"
)


# ======================================================================
# 21. SAVE MODULE ENTROPY
# ======================================================================

module_entropy_df <-

  as.data.frame(

    module_entropy_matrix,

    check.names = FALSE

  )


module_entropy_df$sample_id <-
  rownames(
    module_entropy_df
  )


module_entropy_df <-

  module_entropy_df[

    ,

    c(

      "sample_id",

      setdiff(

        colnames(
          module_entropy_df
        ),

        "sample_id"

      )

    ),

    drop = FALSE

  ]


write.csv(

  module_entropy_df,

  MODULE_FILE,

  row.names = FALSE,

  quote = FALSE

)


cat(
  "Saved:",
  MODULE_FILE,
  "\n"
)


# ======================================================================
# 22. CREATE FINAL MULTISCALE ENTROPY TABLE
# ======================================================================

cat("\n")
cat("Creating final entropy table...\n")


entropy_table <- data.frame(

  sample_id =
    all_sample_ids,

  S_global =
    global_entropy_values,

  S_local_mean =
    local_mean_values,

  S_local_sd =
    local_sd_values,

  S_local_median =
    local_median_values,

  S_meso_mean =
    meso_mean_values,

  S_meso_sd =
    meso_sd_values,

  S_meso_median =
    meso_median_values,

  network_genes =
    network_gene_counts,

  network_edges =
    network_edge_counts,

  total_edge_weight =
    network_total_weights,

  stringsAsFactors =
    FALSE

)


# ----------------------------------------------------------------------
# Add module-level entropy
# ----------------------------------------------------------------------

entropy_table <-

  cbind(

    entropy_table,

    module_entropy_matrix

  )


# ======================================================================
# 23. SAVE FINAL ENTROPY TABLE
# ======================================================================

write.csv(

  entropy_table,

  ENTROPY_FILE,

  row.names = FALSE,

  quote = FALSE

)


cat(
  "Saved:",
  ENTROPY_FILE,
  "\n"
)


# ======================================================================
# 24. SAVE NETWORK SUMMARY
# ======================================================================

network_summary_df <-

  do.call(

    rbind,

    network_summary_list[
      !sapply(
        network_summary_list,
        is.null
      )
    ]

  )


write.csv(

  network_summary_df,

  NETWORK_SUMMARY_FILE,

  row.names = FALSE,

  quote = FALSE

)


cat(
  "Saved:",
  NETWORK_SUMMARY_FILE,
  "\n"
)


# ======================================================================
# 25. RUN SUMMARY
# ======================================================================

run_summary <- data.frame(

  parameter = c(

    "LIONESS_directory",

    "LIONESS_files_found",

    "Successful_networks",

    "Failed_networks",

    "MSigDB_libraries",

    "MSigDB_modules",

    "Minimum_module_genes",

    "Minimum_network_genes",

    "Local_entropy_normalization",

    "Total_output_samples"

  ),

  value = c(

    LIONESS_DIR,

    n_samples,

    successful_samples,

    failed_samples,

    length(
      LIBRARY_FILES
    ),

    length(
      module_gene_list
    ),

    MIN_MODULE_GENES,

    MIN_NETWORK_GENES,

    NORMALIZE_LOCAL_ENTROPY,

    n_samples

  ),

  stringsAsFactors =
    FALSE

)


write.csv(

  run_summary,

  RUN_SUMMARY_FILE,

  row.names = FALSE,

  quote = FALSE

)


# ======================================================================
# 26. FINAL QUALITY CONTROL
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("FINAL QUALITY CONTROL\n")
cat("============================================================\n")


cat(
  "LIONESS networks discovered:",
  n_samples,
  "\n"
)


cat(
  "Successfully processed:",
  successful_samples,
  "\n"
)


cat(
  "Failed networks:",
  failed_samples,
  "\n"
)


cat(
  "MSigDB modules:",
  length(
    module_gene_list
  ),
  "\n"
)


cat(
  "Genes represented in local matrix:",
  length(
    all_network_genes
  ),
  "\n"
)


# ----------------------------------------------------------------------
# Global entropy statistics
# ----------------------------------------------------------------------

valid_global <- is.finite(

  entropy_table$S_global

)


if (
  any(valid_global)
) {

  cat(
    "\nGlobal entropy statistics:\n"
  )


  cat(
    "Mean:",
    mean(
      entropy_table$S_global[
        valid_global
      ]
    ),
    "\n"
  )


  cat(
    "SD:",
    sd(
      entropy_table$S_global[
        valid_global
      ]
    ),
    "\n"
  )


  cat(
    "Median:",
    median(
      entropy_table$S_global[
        valid_global
      ]
    ),
    "\n"
  )


  cat(
    "Minimum:",
    min(
      entropy_table$S_global[
        valid_global
      ]
    ),
    "\n"
  )


  cat(
    "Maximum:",
    max(
      entropy_table$S_global[
        valid_global
      ]
    ),
    "\n"
  )

}


# ----------------------------------------------------------------------
# Local entropy statistics
# ----------------------------------------------------------------------

valid_local <- is.finite(

  entropy_table$S_local_mean

)


if (
  any(valid_local)
) {

  cat(
    "\nLocal entropy mean:\n"
  )


  cat(
    "Mean:",
    mean(
      entropy_table$S_local_mean[
        valid_local
      ]
    ),
    "\n"
  )


  cat(
    "SD:",
    sd(
      entropy_table$S_local_mean[
        valid_local
      ]
    ),
    "\n"
  )

}


# ----------------------------------------------------------------------
# Meso entropy statistics
# ----------------------------------------------------------------------

valid_meso <- is.finite(

  entropy_table$S_meso_mean

)


if (
  any(valid_meso)
) {

  cat(
    "\nMeso entropy mean:\n"
  )


  cat(
    "Mean:",
    mean(
      entropy_table$S_meso_mean[
        valid_meso
      ]
    ),
    "\n"
  )


  cat(
    "SD:",
    sd(
      entropy_table$S_meso_mean[
        valid_meso
      ]
    ),
    "\n"
  )

}


# ======================================================================
# 27. CHECK GLOBAL ENTROPY RANGE
# ======================================================================

if (
  any(valid_global)
) {

  if (

    all(

      entropy_table$S_global[
        valid_global
      ] >= 0

    ) &&

      all(

        entropy_table$S_global[
          valid_global
        ] <= 1

      )

  ) {

    cat(
      "\nS_global range check: PASSED\n"
    )

  } else {

    cat(
      "\nWARNING: Some S_global values are outside [0,1].\n"
    )

  }

}


# ======================================================================
# 28. CHECK MISSING VALUES
# ======================================================================

cat("\n")
cat("Missing-value checks:\n")


cat(

  "Missing S_global:",

  sum(
    is.na(
      entropy_table$S_global
    )
  ),

  "\n"

)


cat(

  "Missing S_local_mean:",

  sum(
    is.na(
      entropy_table$S_local_mean
    )
  ),

  "\n"

)


cat(

  "Missing S_meso_mean:",

  sum(
    is.na(
      entropy_table$S_meso_mean
    )
  ),

  "\n"

)


# ======================================================================
# 29. SHOW FIRST FIVE PATIENTS
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("FIRST FIVE PATIENTS\n")
cat("============================================================\n")


print(

  head(

    entropy_table[

      ,

      c(

        "sample_id",

        "S_global",

        "S_local_mean",

        "S_local_sd",

        "S_local_median",

        "S_meso_mean",

        "S_meso_sd",

        "S_meso_median",

        "network_genes",

        "network_edges"

      ),

      drop = FALSE

    ],

    5

  )

)


# ======================================================================
# 30. OUTPUT FILES
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("OUTPUT FILES\n")
cat("============================================================\n")


cat(
  "1. ",
  ENTROPY_FILE,
  "\n",
  sep = ""
)


cat(
  "2. ",
  LOCAL_FILE,
  "\n",
  sep = ""
)


cat(
  "3. ",
  MODULE_FILE,
  "\n",
  sep = ""
)


cat(
  "4. ",
  NETWORK_SUMMARY_FILE,
  "\n",
  sep = ""
)


cat(
  "5. ",
  MODULE_DEFINITION_FILE,
  "\n",
  sep = ""
)


cat(
  "6. ",
  RUN_SUMMARY_FILE,
  "\n",
  sep = ""
)


# ======================================================================
# 31. INTERPRETATION SUMMARY
# ======================================================================

cat("\n")
cat("============================================================\n")
cat("WHAT THESE SCORES REPRESENT\n")
cat("============================================================\n")


cat(
  "\nLOCAL:\n"
)


cat(
  "S_local = gene-level entropy.\n"
)


cat(
  "It measures how evenly a patient's network connectivity\n"
)


cat(
  "weight is distributed among the neighbors of each gene.\n"
)


cat(
  "\nMESO:\n"
)


cat(
  "S_module_* = MSigDB module/category entropy.\n"
)


cat(
  "Each module summarizes the local entropy of genes belonging\n"
)


cat(
  "to an MSigDB-defined biological category.\n"
)


cat(
  "\nGLOBAL:\n"
)


cat(
  "S_global = degree-weighted whole-network entropy.\n"
)


cat(
  "It summarizes the overall structural organization of the\n"
)


cat(
  "patient-specific LIONESS network.\n"
)


cat("\n")


# ======================================================================
# 32. FINISHED
# ======================================================================

cat("============================================================\n")
cat("M3 LIONESS ENTROPY COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")


cat(
  "Finished:",
  as.character(Sys.time()),
  "\n"
)


cat("\n")
cat(
  "The analysis now provides patient-specific:\n"
)


cat(
  "  Local  -> gene-level entropy\n"
)


cat(
  "  Meso   -> MSigDB module/category entropy\n"
)


cat(
  "  Global -> whole-network entropy\n"
)


cat("\n")
cat("============================================================\n")
cat("M3 FINISHED\n")
cat("============================================================\n")