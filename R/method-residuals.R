#' @rdname lavaan.pl-class
#' @importFrom stats residuals
#' @export
setMethod("residuals", "lavaan.pl", function(object, ...) lavaan::residuals(as(object, "lavaan"), ...))

#' @rdname lavaan.pl-class
#' @importFrom stats resid
#' @export
setMethod("resid", "lavaan.pl", function(object, ...) lavaan::resid(as(object, "lavaan"), ...))
