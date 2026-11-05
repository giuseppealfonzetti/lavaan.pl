#' @rdname lavaan.pl-class
#' @importFrom methods show
#' @export
setMethod("show", "lavaan.pl", function(object) {
  cat("lavaan.pl", get_latpl_fit(object)$version, "\n")
  show(as(object, "lavaan"))
})
