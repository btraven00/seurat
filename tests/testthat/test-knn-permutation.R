suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
})

# knn.R sources src/ relative to the module root, so run it from there.
root <- normalizePath("../..")
rscript <- file.path(R.home("bin"), "Rscript")

# Overlapping blobs: well-separated ones make Annoy look exact.
set.seed(1)
n <- 600
X <- matrix(rnorm(n * 10), n) + 1.5 * rep(1:3, each = n / 3)
emb <- tempfile(fileext = ".tsv")
writeLines(c(paste(c("cell_id", paste0("PC", 1:10)), collapse = "\t"),
             paste(sprintf("c%03d", 1:n), apply(X, 1, paste, collapse = "\t"), sep = "\t")), emb)
ids <- sprintf("c%03d", 1:n)

run <- function(method, perm, seed = 1) {
  out <- tempfile()
  status <- withr::with_dir(root, system2(rscript, c("knn.R", "--output_dir", out, "--name", "t",
    "--embedding_tsv", emb, "--n_neighbors", "20", "--nn_method", method,
    "--random_seed", seed, "--permutation_seed", perm), stdout = FALSE, stderr = FALSE))
  stopifnot(status == 0)
  file.path(out, "t_neighbors.h5")
}

# Graph as a barcode-keyed, sorted triplet table, comparable across permutations.
by_barcode <- function(path, group = "") {
  m <- read_neighbors_csr(path, group)
  t <- summary(m)
  b <- rownames(m)
  d <- data.frame(i = b[t$i], j = b[t$j], x = t$x)
  d <- d[order(d$i, d$j), , drop = FALSE]
  rownames(d) <- NULL
  d
}
read_neighbors_csr <- function(path, group) {
  ids <- as.character(h5read(path, "cell_ids"))
  g <- if (nzchar(group)) paste0(group, "/") else ""
  sparseMatrix(p = h5read(path, paste0(g, "indptr")), j = h5read(path, paste0(g, "indices")),
               x = as.numeric(h5read(path, paste0(g, "data"))), dims = c(length(ids), length(ids)),
               index1 = FALSE, dimnames = list(ids, ids))
}

test_that("same arguments give a byte-identical file", {
  for (m in c("annoy", "rann"))
    expect_identical(unname(tools::md5sum(run(m, 3))), unname(tools::md5sum(run(m, 3))))
})

test_that("permutation_seed 0 keeps the input order", {
  expect_identical(as.character(h5read(run("rann", 0), "cell_ids")), ids)
})

test_that("a permutation reorders rows but keeps the cell set", {
  got <- as.character(h5read(run("rann", 3), "cell_ids"))
  expect_false(identical(got, ids))
  expect_setequal(got, ids)
})

test_that("rann (exact): the graph is unchanged by permutation, weights included", {
  base <- run("rann", 0)
  for (p in c(3, 11)) {
    perm <- run("rann", p)
    for (g in c("", "connectivities"))
      expect_identical(by_barcode(perm, g), by_barcode(base, g), info = paste("perm", p, g))
  }
})

test_that("annoy: permutation changes the graph even at a fixed seed", {
  # Annoy's trees depend on insertion order. If this ever passes as identical,
  # permutation_seed has stopped reaching the search.
  expect_false(identical(by_barcode(run("annoy", 0), "connectivities"),
                         by_barcode(run("annoy", 3), "connectivities")))
})

test_that("annoy: random_seed reaches the index", {
  # Stock Seurat never calls setSeed, so without this the seed axis is inert.
  expect_false(identical(by_barcode(run("annoy", 0, seed = 1), "connectivities"),
                         by_barcode(run("annoy", 0, seed = 2), "connectivities")))
})
