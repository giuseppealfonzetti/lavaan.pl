if (requireNamespace("tinytest", quietly = TRUE)) {
  library(lavaan)
  tinytest::test_package("lavaan.pl", at_home = identical(Sys.getenv("NOT_CRAN"), "true"))
}
