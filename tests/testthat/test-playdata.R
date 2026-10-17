test_that("getLibxGMT reads the default GMT from lib-x", {
  skip_if_not_installed("omicspgx")
  skip_if(is.null(omicspgx::libx_dir()), "lib-x not available")

  gsets <- omicspgx::libx_read("featuresets/legacy/gsetxgene")
  msets <- omicspgx::libx_read("featuresets/legacy/msetxmetabolite")
  G <- getLibxGMT()

  # one column per set of either source; one row per distinct feature
  expect_equal(ncol(G), nrow(gsets) + nrow(msets))
  expect_setequal(colnames(G), c(rownames(gsets), rownames(msets)))
  genes <- colnames(gsets)
  metabolites <- sub(".*:", "", colnames(msets))
  expect_setequal(rownames(G), union(genes, metabolites))
  expect_false(any(duplicated(rownames(G))))

  # the two sources keep their memberships after the merge
  expect_equal(sum(G[genes, rownames(gsets)]), sum(gsets))
  expect_equal(sum(G[metabolites, rownames(msets)]), sum(msets))
})

test_that("getPlaydataGMT is an alias of getLibxGMT", {
  expect_identical(body(getPlaydataGMT), quote(getLibxGMT()))
})
