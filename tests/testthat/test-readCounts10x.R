test_that("Read h5 with multiple assays", {
  path <- testthat::test_path('test_data', 'cellranger-output.h5')
  data <- gencoreSC::readCounts10x(capID = 'test', filepath = path, format = 'h5')
  # expect_equal(length(names(data@assays)), 2)
  ## expect contains tests both ways, so redundant with expect_equal
  expect_contains(names(data@assays), c('RNA', 'ADT'))
})
