#' @rdname lavaan.pl-class
#' @importFrom stats vcov
#' @export
setMethod("vcov", "lavaan.pl", function(object, ...) lavaan::vcov(as(object, "lavaan"), ...))
