testthat::test_that(
  desc = "Peak deconvolution of test data works",
  code = {
    x <- HiResTEC::DeconvoluteSpectrum(dat = HiResTEC::raw, rt = 1026)
    testthat::expect_true(is.matrix(x))
    testthat::expect_equal(nrow(x), 18L)
    testthat::expect_equal(attr(x, "rt"), 1026.535)
    testthat::expect_equal(max(x[,2]), 3312706)
})
