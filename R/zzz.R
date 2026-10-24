.onAttach <- function(libname, pkgname) {

  # Load lavaan
  packageStartupMessage("Loading required package: lavaan")
  attach_lavaan(without = "cfa")

  # Inform conflicts
  packageStartupMessage("lavaan.pl ", utils::packageVersion("lavaan.pl"), ": lavaan.pl::cfa() masks lavaan::cfa()")
}

attach_lavaan <- function(without) {
  do.call("library", list(
    package = "lavaan",
    exclude = without,
    character.only = TRUE
  ))
}
