#' Standardized solution
#'
#' @param object A fitted model of class `lavaan`.
#' @param ... Further arguments for [lavaan::standardizedSolution()].
#' @return The output of [lavaan::standardizedSolution()], with the standard errors of a `lavaan.pl` object.
#' @seealso [lavaan::standardizedSolution()]
#' @export
standardizedSolution <- function(object, ...) {
  lavaan::standardizedSolution(if (inherits(object, "lavaan.pl")) as(object, "lavaan") else object, ...)
}

#' @rdname standardizedSolution
#' @export
standardizedsolution <- standardizedSolution
