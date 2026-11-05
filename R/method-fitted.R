#' @rdname lavaan.pl-class
#' @importFrom stats fitted
#' @export
setMethod("fitted", "lavaan.pl", function(object, ...) lavaan::fitted(as(object, "lavaan"), ...))

#' @rdname lavaan.pl-class
#' @importFrom stats fitted.values
#' @export
setMethod("fitted.values", "lavaan.pl", function(object, ...) lavaan::fitted.values(as(object, "lavaan"), ...))
