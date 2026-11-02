library(lavaan.pl)

source("inst/tinytest/fixtures/lavaan-pml-cases.R")

fits <- lapply(lavaan_pml_cases, function(case) suppressWarnings(case(lavaan_pml_fitters)))
reference <- list(
  lavaan_version = as.character(utils::packageVersion("lavaan")),
  cases = Map(function(fit, name) pml_summary(fit, startsWith(name, "grid/")), fits, names(fits)),
  nested = nested_tests(fits, function(...) suppressWarnings(lavaan::lavTestLRT(...)))
)
saveRDS(reference, "inst/tinytest/fixtures/lavaan-pml-reference.rds")
