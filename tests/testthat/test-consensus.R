test_that("getModuleCrossGenes works with multi = FALSE (single-omics)", {
  skip_on_cran()

  td <- create_test_data(ngenes = 80, nsamples = 20, nmodules = 3)
  res <- computeWGCNA(
    X = td$X,
    samples = td$samples,
    contrasts = td$contrasts,
    ngenes = 80,
    power = 6,
    minmodsize = 5,
    minKME = 0.1,
    deepsplit = 0,
    mergeCutHeight = 0.25,
    calcMethod = "fast",
    lowrank = 10,
    maxBlockSize = 999,
    sv.tom = 10,
    compute.stats = TRUE,
    merge.dendro = FALSE,
    verbose = 0
  )

  # multi = FALSE used to crash: sapply() collapsed the length-1 layer list
  # into a bare matrix instead of keeping it as a list.
  result <- getModuleCrossGenes(res, multi = FALSE)

  expect_type(result, "list")
  expect_true(length(result) > 0)
  for (df in result) {
    expect_true(is.data.frame(df))
    expect_true(all(c("gene", "rho", "module") %in% colnames(df)))
  }
})

test_that("getModuleCrossGenes works with multi = TRUE and custom (non-gx/px) layer names", {
  skip_on_cran()

  td1 <- create_test_data(ngenes = 80, nsamples = 20, nmodules = 3)
  res1 <- computeWGCNA(
    X = td1$X,
    samples = td1$samples,
    contrasts = td1$contrasts,
    ngenes = 80,
    power = 6,
    minmodsize = 5,
    minKME = 0.1,
    deepsplit = 0,
    mergeCutHeight = 0.25,
    calcMethod = "fast",
    lowrank = 10,
    maxBlockSize = 999,
    sv.tom = 10,
    compute.stats = TRUE,
    merge.dendro = FALSE,
    verbose = 0
  )

  set.seed(103)
  X2 <- td1$X + matrix(rnorm(length(td1$X), sd = 0.3), nrow = nrow(td1$X))
  dimnames(X2) <- dimnames(td1$X)
  res2 <- computeWGCNA(
    X = X2,
    samples = td1$samples,
    contrasts = td1$contrasts,
    ngenes = 80,
    power = 6,
    minmodsize = 5,
    minKME = 0.1,
    deepsplit = 0,
    mergeCutHeight = 0.25,
    calcMethod = "fast",
    lowrank = 10,
    maxBlockSize = 999,
    sv.tom = 10,
    compute.stats = TRUE,
    merge.dendro = FALSE,
    verbose = 0
  )

  # Custom layer names ("act"/"notact" instead of "gx"/"px") with no explicit
  # `ref` used to crash: ref resolution fell through to character(0) and blew
  # up the length-zero `if (!ref %in% names(wgcna))` condition.
  wgcna_list <- list(act = res1, notact = res2)
  result <- getModuleCrossGenes(wgcna_list, multi = TRUE)

  expect_type(result, "list")
  expect_true(length(result) > 0)
  for (df in result) {
    expect_true(is.data.frame(df))
    expect_true(all(c("gene", "rho", "module") %in% colnames(df)))
  }
})
