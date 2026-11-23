library(testthat)
library(MMAPPR2)

BiocParallel::register(BiocParallel::SerialParam())

testthat::test_check("MMAPPR2")

test_that(".avgFiles single-file shortcut returns per-position values", {
  x <- data.table::data.table(CHROM = "chr1", POS = 1:3, FILE_ID = 1L,
    A.FREQ = c(.5, 1, 0), C.FREQ = c(.5, 0, 1), G.FREQ = 0, T.FREQ = 0,
    CVG = c(10L, 20L, 30L))
  for (mode in c("simple", "weighted")) {
    out <- MMAPPR2:::.avgFiles(data.table::copy(x), mode, 1L)
    expect_equal(out$AVE.CVG, c(10, 20, 30))
    expect_equal(out$AVE.A.FREQ, x$A.FREQ)
    expect_equal(names(out), c("CHROM", "POS", "AVE.CVG", "AVE.A.FREQ",
                               "AVE.C.FREQ", "AVE.G.FREQ", "AVE.T.FREQ"))
    expect_identical(data.table::key(out), c("CHROM", "POS"))
  }
})

test_that(".avgFiles with nFiles = 1 still averages duplicated positions", {
  x <- data.table::data.table(CHROM = "chr1", POS = c(1L, 1L), FILE_ID = 1L,
    A.FREQ = c(.5, 0), C.FREQ = c(.5, 1), G.FREQ = 0, T.FREQ = 0,
    CVG = c(10L, 30L))
  out <- MMAPPR2:::.avgFiles(data.table::copy(x), "simple", 1L)
  expect_equal(out$AVE.A.FREQ, .25)
  expect_equal(out$AVE.CVG, 40)
})
