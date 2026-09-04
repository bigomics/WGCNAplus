test_that("computeWGCNA runs on small data", {
  skip_on_cran()

  td <- create_test_data(ngenes = 80, nsamples = 20, nmodules = 3)

  result <- computeWGCNA(
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

  expect_type(result, "list")
  expect_true("datExpr" %in% names(result))
  expect_true("datTraits" %in% names(result))
  expect_true("net" %in% names(result))
  expect_true("me.genes" %in% names(result))
  expect_true("me.colors" %in% names(result))
  expect_true("W" %in% names(result))
  expect_true("modTraits" %in% names(result))

  # Check net has expected fields
  expect_true("MEs" %in% names(result$net))
  expect_true("colors" %in% names(result$net))

  # Check dimensions
  expect_equal(nrow(result$datExpr), 20)
  expect_true(ncol(result$datExpr) <= 80)
})

test_that("fastTOMsimilarity produces valid similarity matrix", {
  skip_on_cran()

  set.seed(42)
  n <- 30
  X <- matrix(rnorm(n * 15), nrow = 15, ncol = n)
  A <- WGCNA::adjacency(X, power = 6, type = "signed")

  TOM <- WGCNAplus:::fastTOMsimilarity(A, tomtype = "signed", lowrank = 10)

  expect_equal(nrow(TOM), n)
  expect_equal(ncol(TOM), n)
  expect_true(all(TOM >= 0, na.rm = TRUE))
  expect_true(all(TOM <= 1, na.rm = TRUE))
})

test_that("fastTOMsimilarity validates tomtype", {
  skip_on_cran()

  set.seed(42)
  n <- 10
  X <- matrix(rnorm(n * 5), nrow = 5, ncol = n)
  A <- WGCNA::adjacency(X, power = 6, type = "signed")

  expect_error(
    WGCNAplus:::fastTOMsimilarity(A, tomtype = "invalid", lowrank = 5),
    "only works for signed and unsigned tomtype"
  )
})

test_that("computeGeneStats accepts a full wgcna object as net", {
  skip_on_cran()

  td <- create_test_data(ngenes = 30, nsamples = 10, nmodules = 2)
  res <- computeWGCNA(
    X = td$X, samples = td$samples,
    ngenes = 30, power = 6, minmodsize = 3,
    minKME = 0, mergeCutHeight = 0.3,
    calcMethod = "fast", lowrank = 5, sv.tom = 5,
    compute.stats = FALSE, verbose = 0
  )

  # Pass the full wgcna result as 'net' — computeGeneStats should extract internals
  stats <- computeGeneStats(res)
  expect_type(stats, "list")
  expect_true("moduleTraitCor" %in% names(stats))
  expect_true("moduleMembership" %in% names(stats))
  expect_true("traitSignificance" %in% names(stats))
})

test_that("getGeneStats safe_col handles missing traits gracefully", {
  skip_on_cran()

  td <- create_test_data(ngenes = 30, nsamples = 10, nmodules = 2)
  res <- computeWGCNA(
    X = td$X, samples = td$samples,
    ngenes = 30, power = 6, minmodsize = 3,
    minKME = 0, mergeCutHeight = 0.3,
    calcMethod = "fast", lowrank = 5, sv.tom = 5,
    compute.stats = TRUE, verbose = 0
  )

  # Request a non-existent trait — should return NULL (not error)
  df <- getGeneStats(res, trait = "nonexistent_trait", plot = FALSE)
  expect_null(df)
})
