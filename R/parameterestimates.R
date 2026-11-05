#' Parameter estimates
#'
#' @param object A fitted model of class `lavaan`.
#' @param ... Further arguments for [lavaan::parameterEstimates()].
#' @return The output of [lavaan::parameterEstimates()], with the standard errors of a `lavaan.pl` object.
#' @seealso [lavaan::parameterEstimates()]
#' @export
parameterEstimates <- function(object, ...) {
  lavaan::parameterEstimates(if (inherits(object, "lavaan.pl")) as(object, "lavaan") else object, ...)
}

#' @rdname parameterEstimates
#' @export
parameterestimates <- parameterEstimates

#' @rdname parameterEstimates
#' @export
lavParameterEstimates <- parameterEstimates
