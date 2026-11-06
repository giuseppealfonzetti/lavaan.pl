#' lavaan.pl fit of a lavaan.pl object
#'
#' @param object An object of class `lavaan.pl`.
#' @return The internal fit, a list of class `latpl_fit`.
#' @seealso [lavaan.pl()]
#' @export
get_latpl_fit <- function(object) {
  if (!inherits(object, "lavaan.pl")) {
    stop("`object` must come from lavaan.pl(), pmlcfa() or pmlsem().", call. = FALSE)
  }
  object@external$lavaan.pl
}
