#' @rdname lavaan.pl-class
#' @export
setMethod("summary", "lavaan.pl", function(object, ...) {
  cat("lavaan.pl", get_latpl_fit(object)$version, "\n")
  lavaan::summary(as(object, "lavaan"), ...)
})
