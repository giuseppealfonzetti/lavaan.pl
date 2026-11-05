#' @rdname lavaan.pl-class
#' @importFrom stats nobs
#' @export
setMethod("nobs", "lavaan.pl", function(object, ...) lavaan::nobs(as(object, "lavaan"), ...))
