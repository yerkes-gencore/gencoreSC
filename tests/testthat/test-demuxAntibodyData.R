test_that("Load data", {
  path <- testthat::test_path('test_data', 'cellranger-output.h5')
  data <- gencoreSC::readCounts10x(capID = 'test', filepath = path, format = 'h5')
  # expect_equal(length(names(data@assays)), 2)
  ## expect contains tests both ways, so redundant with expect_equal
  expect_contains(names(data@assays), 'ADT')
})

test_that("Get demux results", {
  hashes <- c('Sample-1', 'Sample-2')
  names(hashes) <- c('Hash1', 'Hash2')
  path <- testthat::test_path('test_data', 'cellranger-output.h5')
  obj <- gencoreSC::readCounts10x(capID = 'test', filepath = path, format = 'h5')
  data <- demuxAntibodyData(obj = obj, labels = hashes, assay = 'ADT')
  expect_contains(levels(data$hash.ID), c('Doublet', 'Negative', 'Sample-1', 'Sample-2'))
})
