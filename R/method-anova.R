#' @rdname lavaan.pl-class
#' @importFrom stats anova
#' @export
setMethod("anova", "lavaan.pl", function(object, ...) {
  mc <- as.list(match.call(expand.dots = TRUE))[-1]
  modp <- vapply(list(...), inherits, logical(1), "lavaan")
  lavTestLRT(object, ..., model_names = vapply(mc[c(TRUE, modp)], function(x) deparse1(x), ""))
})
