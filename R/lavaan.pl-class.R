#' Fits of lavaan.pl
#'
#' @param object An object of class `lavaan.pl`.
#' @param ... Further arguments for the lavaan method.
#' @return As the lavaan method. `logLik()` returns the pairwise log-likelihood.
#' @seealso [lavaan.pl()], [lavaan::lavaan-class]
#' @importClassesFrom lavaan lavaan
#' @export
setClass("lavaan.pl", contains = "lavaan")
