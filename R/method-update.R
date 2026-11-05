#' @rdname lavaan.pl-class
#' @importFrom stats update
#' @export
setMethod("update", "lavaan.pl", function(object, ...) {
  call <- object@call
  extras <- match.call(expand.dots = FALSE)$...
  for (name in names(extras)) call[[name]] <- extras[[name]]
  eval(call, parent.frame())
})
