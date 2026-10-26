.onAttach <- function(libname, pkgname) {

  # Load lavaan
  packageStartupMessage("Loading required package: lavaan")
  attach_lavaan(without = c("cfa", "sem", "lavTestLRT", "anova"))

  # Inform conflicts
  packageStartupMessage("lavaan.pl ", utils::packageVersion("lavaan.pl"), ": lavaan.pl::cfa(), sem(), lavTestLRT() and anova() mask lavaan's")
}

attach_lavaan <- function(without) {
  do.call("library", list(
    package = "lavaan",
    exclude = without,
    character.only = TRUE
  ))
}
