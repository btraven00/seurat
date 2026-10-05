#!/usr/bin/env Rscript
# NNG module: Seurat FindNeighbors -> SNN graph.
#
# Output HDF5 {output_dir}/{name}_neighbors.h5, same layout as the scanpy knn
# entrypoint, so every CLUST module reads it unchanged:
#   /cell_ids                               barcodes, in (permuted) row order
#   /data /indices /indptr                  CSR kNN distances (self excluded)
#   /connectivities/{data,indices,indptr}   CSR SNN graph (Jaccard, pruned)
#
# Seeds:
#   --random_seed       seeds the Annoy index. Seurat itself never calls
#                       setSeed, and R's set.seed does not reach Annoy's RNG,
#                       so stock FindNeighbors is fixed at Annoy's built-in
#                       seed. Ignored for rann (exact kd-tree, no RNG).
#                       rann is order-invariant (tested) except under exact
#                       ties at the k-th neighbour, i.e. duplicate rows.
#   --permutation_seed  shuffles cell order before the search; 0 = identity.
#                       Leiden walks nodes in index order, so this perturbs
#                       clustering even when the kNN search is exact.
#                       Uses R's sample(), so seed N is NOT the same
#                       permutation as the scanpy knn module's seed N.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(rhdf5)
})

source("src/common/cli.R")
p <- arg_parser("NNG module (Seurat SNN)")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "NNG")              # --embedding_tsv
p <- add_argument(p, "--n_neighbors", type = "integer", help = "k.param (self included, as in Seurat)")
p <- add_argument(p, "--nn_method", type = "character", help = "annoy|rann")
p <- add_argument(p, "--prune_snn", type = "numeric", default = 1/15, help = "drop SNN edges below this Jaccard")
p <- add_argument(p, "--random_seed", type = "integer", help = "Annoy index seed")
p <- add_argument(p, "--permutation_seed", type = "integer", default = 0, help = "shuffle cell order; 0 = identity")
args <- parse_args(p)

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
stopifnot(args$nn_method %in% c("annoy", "rann"))

# skip=1 + header=FALSE parses both TSV layouts (unnamed or named id column);
# column 1 is the cell ids either way.
df <- data.table::fread(args$embedding_tsv, skip = 1, header = FALSE)
X <- as.matrix(df[, -1])
rownames(X) <- df[[1]]
cat(sprintf("  embedding (cells x PCs): %d x %d\n", nrow(X), ncol(X)))

if (args$permutation_seed != 0) {
  set.seed(args$permutation_seed)
  X <- X[sample(nrow(X)), , drop = FALSE]
  cat(sprintf("  permuted %d cells (seed %d)\n", nrow(X), args$permutation_seed))
}

index <- NULL
if (args$nn_method == "annoy") {
  # Same index Seurat:::AnnoyBuildIndex builds (euclidean, 50 trees), plus a seed.
  index <- new(RcppAnnoy::AnnoyEuclidean, ncol(X))
  index$setSeed(args$random_seed)
  for (i in seq_len(nrow(X))) index$addItem(i - 1, X[i, ])
  index$build(50)
}

# One search; the SNN step is exactly what FindNeighbors does after it.
nn <- FindNeighbors(X, k.param = args$n_neighbors, nn.method = args$nn_method,
                    index = index, return.neighbor = TRUE, verbose = FALSE)
idx <- Indices(nn)
snn <- Seurat:::ComputeSNN(nn_ranked = idx, prune = args$prune_snn)
# SNN keeps its unit diagonal (self is in every neighbourhood), as Seurat's own
# FindClusters sees it.

row <- rep(seq_len(nrow(idx)), ncol(idx))
keep <- as.vector(idx) != row               # scanpy's distances exclude self
dist <- sparseMatrix(i = row[keep], j = as.vector(idx)[keep],
                     x = as.vector(Distances(nn))[keep], dims = dim(snn))

# CSR of M == CSC of t(M): dgCMatrix p/i/x are scipy's indptr/indices/data.
write_csr <- function(f, grp, m) {
  m <- as(t(m), "CsparseMatrix")
  h5write(as.numeric(m@x), f, paste0(grp, "data"))
  h5write(m@i, f, paste0(grp, "indices"))
  h5write(m@p, f, paste0(grp, "indptr"))
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(args$output_dir, paste0(args$name, "_neighbors.h5"))
if (file.exists(out)) file.remove(out)
h5createFile(out)
h5write(rownames(X), out, "cell_ids")
write_csr(out, "", dist)
h5createGroup(out, "connectivities")
write_csr(out, "connectivities/", snn)
cat(sprintf("  SNN: %d edges (%.1f per cell)\n", length(snn@x), length(snn@x) / nrow(X)))
cat(sprintf("  wrote: %s\n", out))
