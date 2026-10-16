test_that("getPlaydataGMT reads the default GMT from lib-x", {
  skip_if_not_installed("omicspgx")
  skip_if(is.null(omicspgx::libx_dir()), "lib-x not available")

  G <- getPlaydataGMT()
  expect_true(all(dim(G) > 0))

  # gene rownames are kept; metabolite rownames lose their "db:" prefix
  genes <- colnames(omicspgx::libx_read("featuresets/legacy/gsetxgene"))
  metabolites <- colnames(omicspgx::libx_read("featuresets/legacy/msetxmetabolite"))
  expect_true(all(genes %in% rownames(G)))
  expect_true(all(sub(".*:", "", metabolites) %in% rownames(G)))
})
