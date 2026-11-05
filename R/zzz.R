.onLoad <- function(libname, pkgname) {
  ns <- topenv()
  for (table in grep("^[.]__T__", names(ns), value = TRUE)) {
    rm(list = setdiff(names(ns[[table]]), "lavaan.pl"), envir = ns[[table]])
  }
}
