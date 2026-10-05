#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(HDF5Array)
  library(rhdf5)
  library(Seurat)
  library(Matrix)
})

# read_neighbors(): load a CSR neighbors graph from a *_neighbors.h5 file.
source("src/read_neighbors.R")
# flavor_algorithm_id(): map a --flavor name to a Seurat FindClusters algorithm id.
source("src/flavor.R")

# arg parsing
source("src/common/cli.R")
p <- arg_parser("CLUST module")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "CLUST")  # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--flavor", type = "character", help = "Clustering algorithm (louvain_original|louvain_multilevel_refinement|slm|leiden)")
p <- add_argument(p, "--resolution", type = "numeric", help = "Clustering resolution")
p <- add_argument(p, "--random_seed", type = "integer", help = "Random seed")
# One job, many resolutions: graph built once, FindClusters at every point of a
# linear grid. Writes {name}_clusters_sweep.tsv (cell_id + one column per
# resolution) and {name}_sweep.json instead of {name}_clusters.tsv.
# leidenbase num_iter: runs EXACTLY this many iterations (must be >= 1; no early
# stop, no until-convergence mode). Default 10 is Seurat's FindClusters default.
# On tm-facs modularity is 0.96261 at 2, 0.96293 at 10, 0.96292 at 50 -- 2 is
# near-converged and 2.8x faster. Not the same unit as scanpy's n_iterations
# (igraph, negative = until no improvement) or rapids' max_iter (cuGraph cap).
p <- add_argument(p, "--n_iterations", type = "integer", default = 10L,
                  help = "Leiden iterations (leidenbase num_iter, >= 1)")
p <- add_argument(p, "--sweep", type = "character",
                  help = "MIN:MAX:STEP resolution grid (inclusive); replaces --resolution")
args <- parse_args(p)                      # argparser's own parser
if (is.na(args$resolution) == is.na(args$sweep))
  stop("give exactly one of --resolution or --sweep", call. = FALSE)

sweep_grid <- function(spec) {
  v <- as.numeric(strsplit(spec, ":", fixed = TRUE)[[1]])
  if (length(v) != 3 || anyNA(v) || !(v[1] > 0 && v[1] <= v[2] && v[3] > 0))
    stop(sprintf("bad --sweep '%s': need 0 < MIN <= MAX and STEP > 0", spec), call. = FALSE)
  n <- floor((v[2] - v[1]) / v[3] + 1e-9) + 1
  if (n > 500) stop(sprintf("--sweep '%s' gives %d resolutions; cap is 500", spec, n), call. = FALSE)
  round(v[1] + (seq_len(n) - 1) * v[3], 10)
}

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


# Reproducibility
set.seed(args$random_seed)

# Load neighbors graph into Seurat Object. Cluster on the connectivities graph
# (UMAP-style affinities); the distances graph is the flat root layout.
neighbors_mat <- read_neighbors(args$neighbors_h5, group = "connectivities")

neighbors_graph <- as.Graph(neighbors_mat)

cat("Neighbors graph dimensions:\n")
print(dim(neighbors_graph))
cat("\n")

so <- CreateSeuratObject(
  counts = neighbors_graph,
  assay = "RNA"
)

so@graphs$neighbors <- neighbors_graph


# Map the requested --flavor to a Seurat FindClusters algorithm id (errors if
# --flavor is missing or unknown).
algorithm_seurat_id <- flavor_algorithm_id(args$flavor)


if (!is.na(args$sweep)) {
  res <- sweep_grid(args$sweep)
  out <- data.frame(cell_id = colnames(so), check.names = FALSE)
  info <- vector("list", length(res))
  for (i in seq_along(res)) {
    t0 <- Sys.time()
    so <- FindClusters(so, algorithm = algorithm_seurat_id, resolution = res[i], n.iter = args$n_iterations,
                       graph.name = "neighbors", random.seed = args$random_seed, verbose = FALSE)
    lab <- as.character(so$seurat_clusters)
    secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)
    out[[format(res[i], digits = 6)]] <- lab
    info[[i]] <- list(resolution = res[i], n_clusters = length(unique(lab)), seconds = secs)
    cat(sprintf("  resolution %g: %d clusters (%.3f s)\n", res[i], length(unique(lab)), secs))
  }
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(args$output_dir, paste0(args$name, "_clusters_sweep.tsv"))
  write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE)
  jsonlite::write_json(list(sweep = args$sweep, flavor = args$flavor, random_seed = args$random_seed,
                            n_iterations = args$n_iterations,
                            resolutions = info), file.path(args$output_dir, paste0(args$name, "_sweep.json")),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
  cat("wrote:", f, "\n")
  quit(save = "no")
}

# Run clustering
so <- FindClusters(
  so,
  algorithm = algorithm_seurat_id,
  resolution = args$resolution,
  n.iter = args$n_iterations,
  graph.name = "neighbors",
  random.seed = args$random_seed,
  verbose = TRUE
)

cat("Running clustering...\n")
cat("Selected algorithm:", args$flavor, "\n")
cat("Algorithm ID:", algorithm_seurat_id, "\n\n")


# Extract clusters matrix
m_clusters <- data.frame(cell_id = colnames(so),
                         cluster = as.character(so$seurat_clusters))

cat("Cluster matrix dimensions:\n")
print(dim(m_clusters))


# Save cluster matrix as .tsv
output_file <- file.path(
  args$output_dir,
  paste0(args$name, "_clusters.tsv")
)

cat("Writing output to:\n")
cat(output_file, "\n\n")

write.table(
  m_clusters,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(file.info(output_file)[, c("size", "ctime")])
