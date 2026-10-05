#!/usr/bin/env Rscript
# PCA module (Seurat-backed) for omnibenchmark.
#
# Solvers:
#   approximate — Seurat::RunPCA(approx = TRUE)  [IRLBA]
#   exact       — Seurat::RunPCA(approx = FALSE) [full SVD]
#
# Phases (obkit-events.jsonl):
#   load      h5 -> dgCMatrix -> Seurat object (everything that gets the input
#             into the framework's expected in-memory form).
#   densify   ScaleData(do.scale = FALSE) — allocates a dense
#             (n_genes × n_cells) double matrix in the `scale.data` layer
#             because that's what RunPCA reads from. This is the memory
#             elephant for Seurat: ~n_genes * n_cells * 8 bytes.
#   pca       RunPCA (IRLBA or full SVD).
#   write     TSV serialization.
#
# Mean-center only (do.scale = FALSE) to match scanpy/scrapper/rapids,
# which don't divide by per-gene std. Seurat's default ScaleData does both.

suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(Seurat)
  library(data.table)
  library(jsonlite)  # logger emits via toJSON; satisfied transitively by Seurat
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("PCA module")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "PCA")  # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--solver", type = "character", help = "solver type")
p <- add_argument(p, "--n_components", type = "integer", help = "number of PCs")
p <- add_argument(p, "--random_seed", type = "integer", help = "random seed")
args <- parse_args(p)                      # argparser's own parser

if (!(args$solver %in% c("exact", "approximate", "irlba", "random")))
  stop("Invalid --solver: ", args$solver, " (valid: exact, approximate)")

# tracking stuff
source("src/obkit_logger.R")
source("src/phases.R")


run_pca <- function(so, X, args, phase) {
  # so: Seurat object with `counts` + `data` layers populated from X.
  # X:  the same gene-by-cell sparse dgCMatrix (for gene/cell ids + PCA features).
  set.seed(args$random_seed)

  approx <- switch(args$solver,
    approximate = TRUE,
    exact       = FALSE,
    stop("unknown solver: ", args$solver)
  )

  so <- phase("densify", function(attrs) {
    attrs$do_scale        <- FALSE
    attrs$dense_bytes_est <- as.numeric(nrow(X)) * as.numeric(ncol(X)) * 8
    ScaleData(so, features = rownames(X), do.scale = FALSE, verbose = FALSE)
  })

  so <- phase("pca", function(attrs) {
    attrs$solver       <- args$solver
    attrs$approx       <- approx
    attrs$n_components <- args$n_components
    RunPCA(so,
           features = rownames(X),
           npcs = args$n_components,
           approx = approx,
           seed.use = args$random_seed,
           verbose = FALSE)
  })

  dimred    <- so[["pca"]]
  embedding <- Embeddings(dimred)           # (n_cells, n_components)
  loadings  <- Loadings(dimred)             # (n_genes, n_components)
  rownames(embedding) <- colnames(X)
  colnames(embedding) <- paste0("PC", seq_len(ncol(embedding)))
  rownames(loadings)  <- rownames(X)
  colnames(loadings)  <- paste0("PC", seq_len(ncol(loadings)))

  list(embedding = embedding, loadings = loadings)
}

main <- function() {
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "normalized_selected_h5",
              "solver", "n_components", "random_seed")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  logger_init(args$output_dir)

  loaded <- phase("load", function(attrs) {
    # TENxMatrix is lazy/disk-backed; as(., dgCMatrix) materializes it
    # in memory. CreateSeuratObject + SetAssayData wrap the sparse matrix
    # in Seurat's container — sparse-pointer-cheap but allocates metadata
    # (cell IDs as factor levels, etc.). All three steps are "ingest into
    # framework-native form" so they belong in `load`.
    m_lazy <- TENxMatrix(args$normalized_selected_h5, group = "matrix")
    m_mem  <- as(m_lazy, "dgCMatrix")
    so     <- CreateSeuratObject(counts = m_mem)
    so     <- SetAssayData(so, layer = "data", new.data = m_mem)
    attrs$n_genes <- nrow(m_mem)
    attrs$n_cells <- ncol(m_mem)
    attrs$nnz     <- as.integer(length(m_mem@x))
    list(X = m_mem, so = so)
  })
  cat(sprintf("  matrix (genes x cells): %d x %d\n",
              nrow(loaded$X), ncol(loaded$X)))

  res <- run_pca(loaded$so, loaded$X, args, phase)
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
              nrow(res$embedding), ncol(res$embedding),
              nrow(res$loadings),  ncol(res$loadings)))

  phase("write", function(attrs) {
    out_pcas <- file.path(args$output_dir, sprintf("%s_embedding.tsv", args$name))
    fwrite(data.frame(cell_id = rownames(res$embedding), res$embedding), out_pcas,
           sep = "\t", quote = FALSE, row.names = FALSE)
    out_loadings <- file.path(args$output_dir, sprintf("%s_loadings.tsv", args$name))
    fwrite(data.frame(gene = rownames(res$loadings), res$loadings), out_loadings,
           sep = "\t", quote = FALSE, row.names = FALSE)
    attrs$pcas     <- out_pcas
    attrs$loadings <- out_loadings
    cat(sprintf("  wrote: %s\n  wrote: %s\n", out_pcas, out_loadings))
  })
}

if (sys.nframe() == 0L) {
  main()
}
