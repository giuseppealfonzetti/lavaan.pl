if (requireNamespace("tinytest", quietly = TRUE)) {
  tinytest::test_package("lavaan.pl", at_home = identical(Sys.getenv("NOT_CRAN"), "true"))
}
