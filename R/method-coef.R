#' @rdname lavaan.pl-class
#' @importFrom stats coef
#' @export
setMethod("coef", "lavaan.pl", function(object, ...) lavaan::coef(as(object, "lavaan"), ...))
