#' @rdname lavaan.pl-class
#' @importFrom stats predict
#' @export
setMethod("predict", "lavaan.pl", function(object, ...) lavaan::predict(as(object, "lavaan"), ...))
