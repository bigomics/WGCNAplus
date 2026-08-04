test_that("runPreservationWGCNA runs on small two-layer data (reference = 1)", {
  skip_on_cran()

  td <- create_test_data(ngenes = 100, nsamples = 30, nmodules = 3)

  # Second layer: same genes/samples as layer 1, lightly perturbed so it's a
  # distinct (but non-degenerate) dataset for module preservation to compare
  # against the reference.
  set.seed(101)
  X2 <- td$X + matrix(rnorm(length(td$X), sd = 0.3), nrow = nrow(td$X))
  dimnames(X2) <- dimnames(td$X)

  exprList <- list(layer1 = td$X, layer2 = X2)

  result <- runPreservationWGCNA(
    exprList,
    phenoData = td$samples,
    contrasts = td$contrasts,
    power = 6,
    reference = 1,
    ngenes = 100,
    minModuleSize = 5,
    deepSplit = 0,
    compute.stats = TRUE,
    compute.enrichment = FALSE
  )

  expect_type(result, "list")
  expect_true("Zsummary" %in% names(result))
  expect_true("medianRank" %in% names(result))
  expect_true("moduleSize" %in% names(result))
  expect_true("modTraits" %in% names(result))
  expect_true("MEs" %in% names(result))

  # Guards the length-1 vector-collapse bug: with a single non-reference
  # comparison, sapply() used to drop dimensions and turn these into plain
  # vectors instead of matrices.
  expect_true(is.matrix(result$Zsummary))
  expect_true(is.matrix(result$medianRank))
})

test_that("runPreservationWGCNA runs with compute.stats = FALSE, compute.enrichment = TRUE (ref-scope guard)", {
  skip_on_cran()

  td <- create_test_data(ngenes = 100, nsamples = 30, nmodules = 3)

  set.seed(102)
  X2 <- td$X + matrix(rnorm(length(td$X), sd = 0.3), nrow = nrow(td$X))
  dimnames(X2) <- dimnames(td$X)

  exprList <- list(layer1 = td$X, layer2 = X2)

  # This is the exact combination (compute.stats = FALSE, compute.enrichment
  # = TRUE) that used to crash with "object 'ref' not found", because `ref`
  # was only assigned inside the compute.stats branch but read from the
  # compute.enrichment branch below it. GMT is left NULL here (a real
  # gene-set matrix isn't available in this fast unit test), which means the
  # enrichment body itself is skipped via `compute.enrichment && !is.null(GMT)`
  # — but the call must still complete without erroring on `ref`.
  result <- runPreservationWGCNA(
    exprList,
    phenoData = td$samples,
    contrasts = td$contrasts,
    power = 6,
    reference = 1,
    ngenes = 100,
    minModuleSize = 5,
    deepSplit = 0,
    compute.stats = FALSE,
    compute.enrichment = TRUE,
    GMT = NULL
  )

  expect_type(result, "list")
  expect_true("Zsummary" %in% names(result))
})
