# ======================================================================
# M8 — NETWORK TOPOLOGY COMPARISON BETWEEN ECOSYSTEM STATES
# PURPOSE
# -------
# Compare the structural (topological) properties of patient-specific
# LIONESS networks between the FROZEN M3.2 ecosystem states, going
# beyond the entropy summaries already computed in M3: network
# density, weighted degree, community structure (Louvain modularity),
# and hub-gene composition.
#
# IMPORTANT
# ---------
#   - M3.2 ecosystem states are FROZEN and are not recomputed here.
#   - Per-patient LIONESS edge lists (already written by M3) are read
#     as-is; no network is re-estimated in this script.
#   - Sample identifiers are compared via a NORMALIZED key (case,
#     dash/dot unification, and a stripped "_LIONESS" suffix), since
#     the ecosystem-state file's sample_id and the per-patient network
#     filenames use slightly different ID conventions in this
#     pipeline (confirmed while debugging M4/M5/M6/M7).
#
# METHOD
# ------
# 1. Load the FROZEN M3.2 ecosystem-state assignment and the per-
#    patient LIONESS edge lists written by M3.
# 2. For each patient network: build an igraph object (undirected,
#    weighted by absolute edge weight), then compute
#      - network density
#      - mean weighted degree
#      - number of Louvain communities and modularity
#      - top-decile-degree "hub" genes
# 3. Compare density, mean weighted degree, community count, and
#    modularity between states (Wilcoxon rank-sum tests).
# 4. Compare hub-gene composition between states via the Jaccard
#    overlap of each state's most frequently recurring hub genes.
# 5. Build one representative consensus network per state (mean edge
#    weight across that state's patients, thresholded) purely for
#    visualization.
# 6. Figures: state-specific consensus network diagrams, boxplots of
#    density / modularity / mean degree by state, and a heatmap of
#    top hub-gene recurrence frequency by state.
#
# INPUTS
# -------
# /content/M3_LIONESS_entropy/M3.2_module_axes/
#     M3.2_patient_ecosystem_states.csv
# /content/M3_LIONESS/patient_networks/<safe_patient_id>_LIONESS_network.csv
#     (one file per patient, written by M3; columns: sample_id, gene1,
#     gene2, lioness_weight)
#
# OUTPUTS
# --------
# /content/M8_network_topology/
#     M8_patient_topology_metrics.csv
#     M8_topology_state_comparison.csv
#     M8_hub_gene_frequency.csv
#     M8_topology_boxplots.pdf
#     M8_hub_gene_heatmap.pdf
#     M8_state_consensus_networks.pdf
#
# ======================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

cat("\n")
cat("=====================================================================\n")
cat("M8 — NETWORK TOPOLOGY COMPARISON BETWEEN ECOSYSTEM STATES\n")
cat("=====================================================================\n\n")

cat("Started:", as.character(Sys.time()), "\n\n")


# ======================================================================
# 1. PACKAGES
# ======================================================================

required_packages <- c("igraph", "readr", "dplyr", "tidyr", "ggplot2",
                        "pheatmap", "tibble", "purrr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(igraph)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
  library(purrr)
})

cat("Packages ready.\n\n")


# ======================================================================
# 2. FILE SETTINGS
# ======================================================================

STATE_FILE <- "/content/M3_LIONESS_entropy/M3.2_module_axes/M3.2_patient_ecosystem_states.csv"

PATIENT_NETWORK_DIR <- "/content/M3_LIONESS/patient_networks"

OUTPUT_DIR <- "/content/M8_network_topology"

HUB_DEGREE_PERCENTILE <- 0.90   # top-decile weighted degree = "hub"
CONSENSUS_EDGE_THRESHOLD <- 0.30  # matches PATIENT_EDGE_THRESHOLD in M3
TOP_N_HUB_GENES_FOR_HEATMAP <- 25

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

cat("Output directory:", normalizePath(OUTPUT_DIR), "\n\n")

kv <- function(label, value) cat(label, ":", value, "\n")

# Normalize any sample/file-derived identifier to a common key: upper
# case, unify "-"/"." to "_", and strip a trailing "_LIONESS" suffix
# (present on the M3.2 ecosystem-state file's sample_id but not on the
# per-patient network filenames, confirmed while debugging M4-M7).
normalize_id <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("[.\\-]", "_", x)
  x <- sub("_LIONESS$", "", x)
  x
}


# ======================================================================
# 3. LOAD FROZEN ECOSYSTEM STATE AND LOCATE PATIENT NETWORK FILES
# ======================================================================

if (!file.exists(STATE_FILE)) {
  stop(paste0("FROZEN M3.2 ecosystem-state file not found: ", STATE_FILE,
              "\nRun M3.2 before M8."))
}

if (!dir.exists(PATIENT_NETWORK_DIR)) {
  stop(paste0("Patient LIONESS network directory not found: ",
              PATIENT_NETWORK_DIR, "\nRun M3 before M8."))
}

ecosystem_states <- read_csv(STATE_FILE, show_col_types = FALSE)

if (!"sample_id" %in% colnames(ecosystem_states)) {
  colnames(ecosystem_states)[1] <- "sample_id"
}

ecosystem_states <- ecosystem_states %>%
  select(sample_id, ecosystem_state) %>%
  mutate(.join_key = normalize_id(sample_id))

network_files <- list.files(
  PATIENT_NETWORK_DIR, pattern = "_LIONESS_network\\.csv$", full.names = TRUE
)

kv("Frozen-state patients", nrow(ecosystem_states))
kv("Patient network files found", length(network_files))

if (length(network_files) == 0) {
  stop("No per-patient LIONESS network files found.")
}

get_sample_id_from_filename <- function(path) {
  sub("_LIONESS_network\\.csv$", "", basename(path))
}

# Diagnostic: confirm ID formats line up before doing any heavy work.
cat("Sample ID format check:\n")
cat("  Network filenames (normalized, first 3): ",
    paste(head(normalize_id(get_sample_id_from_filename(network_files)), 3), collapse = " | "), "\n")
cat("  Ecosystem-state sample_id (normalized, first 3): ",
    paste(head(ecosystem_states$.join_key, 3), collapse = " | "), "\n\n")


# ======================================================================
# 4. PER-PATIENT TOPOLOGY METRICS
# ======================================================================

cat("Computing per-patient network topology metrics...\n")

hub_gene_records <- list()

topology_metrics <- map_dfr(network_files, function(path) {

  raw_sample_id <- get_sample_id_from_filename(path)
  join_key <- normalize_id(raw_sample_id)

  edges <- tryCatch(
    read_csv(path, show_col_types = FALSE),
    error = function(e) NULL
  )

  if (is.null(edges) || nrow(edges) == 0) {
    return(NULL)
  }

  edge_cols_needed <- c("gene1", "gene2", "lioness_weight")
  if (!all(edge_cols_needed %in% colnames(edges))) {
    return(NULL)
  }
  edges <- edges %>% select(gene1, gene2, lioness_weight) %>%
    rename(gene_1 = gene1, gene_2 = gene2, weight = lioness_weight)

  edges$abs_weight <- abs(edges$weight)

  g <- graph_from_data_frame(
    edges[, c("gene_1", "gene_2", "abs_weight")],
    directed = FALSE
  )
  E(g)$weight <- edges$abs_weight

  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  max_possible_edges <- n_nodes * (n_nodes - 1) / 2

  density_val <- if (max_possible_edges > 0) n_edges / max_possible_edges else NA_real_

  weighted_degree <- strength(g, weights = E(g)$weight)
  mean_degree <- mean(weighted_degree)

  louvain <- tryCatch(
    cluster_louvain(g, weights = E(g)$weight),
    error = function(e) NULL
  )

  n_communities <- if (!is.null(louvain)) length(unique(membership(louvain))) else NA_integer_
  modularity_val <- if (!is.null(louvain)) modularity(louvain) else NA_real_

  hub_cutoff <- quantile(weighted_degree, HUB_DEGREE_PERCENTILE, na.rm = TRUE)
  hub_genes <- names(weighted_degree[weighted_degree >= hub_cutoff])

  hub_gene_records[[join_key]] <<- hub_genes

  tibble(
    sample_id       = raw_sample_id,
    .join_key       = join_key,
    n_nodes         = n_nodes,
    n_edges         = n_edges,
    density         = density_val,
    mean_weighted_degree = mean_degree,
    n_communities   = n_communities,
    modularity      = modularity_val,
    n_hub_genes     = length(hub_genes)
  )
})

kv("Patients with computed topology metrics", nrow(topology_metrics))

write_csv(
  topology_metrics %>% select(-.join_key),
  file.path(OUTPUT_DIR, "M8_patient_topology_metrics.csv")
)


# ======================================================================
# 5. MERGE WITH FROZEN ECOSYSTEM STATE (ID-NORMALIZED JOIN)
# ======================================================================

merged <- ecosystem_states %>%
  inner_join(
    topology_metrics %>% select(-sample_id),
    by = ".join_key"
  )

kv("Patients with both topology metrics and frozen state", nrow(merged))

if (nrow(merged) == 0) {
  cat("\nStill 0 matches after normalization. Example IDs:\n")
  cat("  network files (normalized):", paste(head(topology_metrics$.join_key, 5), collapse = " | "), "\n")
  cat("  states (normalized):       ", paste(head(ecosystem_states$.join_key, 5), collapse = " | "), "\n")
  stop("No overlapping sample IDs between per-patient network files and frozen ecosystem states. ",
       "Compare the printed ID formats above and adjust normalize_id() accordingly.")
}

merged <- merged %>% select(-.join_key)

if (n_distinct(merged$ecosystem_state) != 2) {
  stop("Expected exactly two frozen ecosystem states (State_1, State_2).")
}

merged$ecosystem_state <- factor(merged$ecosystem_state)
state_levels <- levels(merged$ecosystem_state)


# ======================================================================
# 6. STATISTICAL COMPARISON: TOPOLOGY METRICS BY STATE
# ======================================================================

cat("\nTesting topology metrics by frozen ecosystem state...\n")

metrics_to_test <- c("density", "mean_weighted_degree", "n_communities", "modularity")

topology_comparison <- map_dfr(metrics_to_test, function(metric) {

  values <- merged[[metric]]
  group  <- merged$ecosystem_state

  test <- tryCatch(wilcox.test(values ~ group), error = function(e) NULL)
  medians <- tapply(values, group, median, na.rm = TRUE)

  tibble(
    metric        = metric,
    median_state1 = medians[state_levels[1]],
    median_state2 = medians[state_levels[2]],
    p_value       = if (is.null(test)) NA_real_ else test$p.value
  )
})

topology_comparison$FDR <- p.adjust(topology_comparison$p_value, method = "BH")
topology_comparison <- topology_comparison %>% arrange(FDR, p_value)

print(topology_comparison)

write_csv(
  topology_comparison,
  file.path(OUTPUT_DIR, "M8_topology_state_comparison.csv")
)


# ======================================================================
# 7. HUB-GENE FREQUENCY BY STATE (ID-NORMALIZED JOIN)
# ======================================================================

cat("\nSummarizing hub-gene recurrence frequency by state...\n")

# hub_gene_records is keyed by the NORMALIZED join key (set in Section 4).
hub_df <- imap_dfr(hub_gene_records, function(genes, join_key) {
  if (length(genes) == 0) return(NULL)
  tibble(.join_key = join_key, gene = genes)
})

hub_df <- hub_df %>%
  inner_join(
    ecosystem_states %>% select(.join_key, ecosystem_state),
    by = ".join_key"
  )

hub_frequency <- hub_df %>%
  group_by(ecosystem_state, gene) %>%
  summarise(n_patients_hub = n(), .groups = "drop") %>%
  left_join(
    ecosystem_states %>% count(ecosystem_state, name = "n_patients_total"),
    by = "ecosystem_state"
  ) %>%
  mutate(frequency = n_patients_hub / n_patients_total) %>%
  arrange(desc(frequency))

write_csv(
  hub_frequency,
  file.path(OUTPUT_DIR, "M8_hub_gene_frequency.csv")
)

# Jaccard overlap of each state's top hub-gene set.
top_hubs_by_state <- hub_frequency %>%
  group_by(ecosystem_state) %>%
  slice_max(order_by = frequency, n = TOP_N_HUB_GENES_FOR_HEATMAP, with_ties = FALSE) %>%
  ungroup()

hub_sets <- split(top_hubs_by_state$gene, top_hubs_by_state$ecosystem_state)

if (length(hub_sets) == 2) {
  intersection_n <- length(intersect(hub_sets[[1]], hub_sets[[2]]))
  union_n <- length(union(hub_sets[[1]], hub_sets[[2]]))
  jaccard <- if (union_n > 0) intersection_n / union_n else NA_real_
  kv("Jaccard overlap of top hub genes between states", round(jaccard, 3))
}


# ======================================================================
# 8. FIGURE — TOPOLOGY METRIC BOXPLOTS BY STATE
# ======================================================================

cat("\nGenerating topology metric boxplots...\n")

topology_long <- merged %>%
  select(sample_id, ecosystem_state, all_of(metrics_to_test)) %>%
  pivot_longer(cols = all_of(metrics_to_test), names_to = "metric", values_to = "value")

topology_boxplots <- ggplot(topology_long, aes(x = ecosystem_state, y = value,
                                                fill = ecosystem_state)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "LIONESS Network Topology Metrics by Frozen Ecosystem State",
    x = "Ecosystem state", y = "Value"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  file.path(OUTPUT_DIR, "M8_topology_boxplots.pdf"),
  topology_boxplots, width = 10, height = 8
)


# ======================================================================
# 9. FIGURE — HUB-GENE FREQUENCY HEATMAP
# ======================================================================

cat("Generating hub-gene frequency heatmap...\n")

heatmap_genes <- top_hubs_by_state %>%
  distinct(gene) %>%
  pull(gene)

if (length(heatmap_genes) >= 2) {

  heatmap_matrix <- hub_frequency %>%
    filter(gene %in% heatmap_genes) %>%
    select(ecosystem_state, gene, frequency) %>%
    pivot_wider(names_from = ecosystem_state, values_from = frequency, values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  pdf(file.path(OUTPUT_DIR, "M8_hub_gene_heatmap.pdf"), width = 6, height = 10)

  pheatmap(
    heatmap_matrix,
    cluster_cols = FALSE,
    main = "Hub-Gene Recurrence Frequency by Ecosystem State",
    fontsize_row = 6
  )

  dev.off()

} else {
  cat("Fewer than two hub genes identified; heatmap skipped.\n")
}


# ======================================================================
# 10. FIGURE — STATE-SPECIFIC CONSENSUS NETWORK DIAGRAMS
#     (ID-NORMALIZED FILTERING OF STATE MEMBERSHIP)
# ======================================================================

cat("Building and plotting state-specific consensus networks...\n")

build_consensus_network <- function(state_label, edge_threshold) {

  state_keys <- ecosystem_states %>%
    filter(ecosystem_state == state_label) %>%
    pull(.join_key)

  file_keys <- normalize_id(get_sample_id_from_filename(network_files))

  state_files <- network_files[file_keys %in% state_keys]

  if (length(state_files) == 0) return(NULL)

  all_edges <- map_dfr(state_files, function(path) {
    e <- read_csv(path, show_col_types = FALSE)
    if (!all(c("gene1", "gene2", "lioness_weight") %in% colnames(e))) {
      return(NULL)
    }
    e %>%
      select(gene1, gene2, lioness_weight) %>%
      rename(gene_1 = gene1, gene_2 = gene2, weight = lioness_weight)
  })

  consensus <- all_edges %>%
    mutate(pair = paste(pmin(gene_1, gene_2), pmax(gene_1, gene_2), sep = "__")) %>%
    group_by(pair) %>%
    summarise(
      gene_1 = first(gene_1), gene_2 = first(gene_2),
      mean_weight = mean(weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(abs(mean_weight) >= edge_threshold)

  graph_from_data_frame(
    consensus[, c("gene_1", "gene_2", "mean_weight")],
    directed = FALSE
  )
}

pdf(file.path(OUTPUT_DIR, "M8_state_consensus_networks.pdf"), width = 14, height = 7)
par(mfrow = c(1, 2))

for (state_label in state_levels) {

  g_state <- build_consensus_network(state_label, CONSENSUS_EDGE_THRESHOLD)

  if (is.null(g_state) || vcount(g_state) == 0) {
    plot.new()
    title(main = paste0(state_label, " (no consensus edges at threshold)"))
    next
  }

  edge_col <- ifelse(E(g_state)$mean_weight > 0, "firebrick", "steelblue")

  plot(
    g_state,
    vertex.size = 3,
    vertex.label = NA,
    edge.width = abs(E(g_state)$mean_weight) * 3,
    edge.color = edge_col,
    layout = layout_with_fr(g_state),
    main = paste0(state_label, " consensus network (|mean r| \u2265 ",
                  CONSENSUS_EDGE_THRESHOLD, ")")
  )
}

dev.off()

cat("\nM8 network topology comparison complete.\n")
cat("Finished:", as.character(Sys.time()), "\n")