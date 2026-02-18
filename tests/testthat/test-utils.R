test_that("labels2colors handles numeric labels", {
  colors <- labels2colors(c(1, 2, 3, 0))
  expect_true(all(is.character(colors)))
  expect_equal(colors[4], "grey") # 0 maps to grey
})

test_that("labels2colors handles color labels", {
  colors <- labels2colors(c("blue", "red", "grey"))
  expect_equal(colors, c("blue", "red", "grey"))
})

test_that("labels2colors handles non-standard labels", {
  colors <- labels2colors(c("moduleA", "moduleB", "moduleA"))
  expect_true(all(colors %in% WGCNA::standardColors()))
  expect_equal(colors[1], colors[3]) # same label -> same color
})

test_that("purpleGreyYellow returns correct number of colors", {
  cols <- WGCNAplus:::purpleGreyYellow(5)
  expect_length(cols, 5)
  expect_true(all(grepl("^#", cols)))
})

test_that("rho2bluered handles vector input", {
  cols <- WGCNAplus:::rho2bluered(c(-1, 0, 1), f = 1)
  expect_length(cols, 3)
})

test_that("rho2bluered handles matrix input", {
  R <- matrix(c(-1, 0, 0, 1), nrow = 2)
  cols <- WGCNAplus:::rho2bluered(R)
  expect_equal(dim(cols), dim(R))
})
