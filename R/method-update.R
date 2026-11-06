#' @rdname lavaan.pl-class
#' @importFrom stats update
#' @export
setMethod("update", "lavaan.pl", function(object, ...) {
  call <- object@call
  extras <- match.call(expand.dots = FALSE)$...
  for (name in setdiff(names(extras), "add")) call[[name]] <- extras[[name]]
  if (!is.null(extras$add)) call$model <- call("paste", call$model, extras$add, sep = "\n")
  eval(call, parent.frame())
})
