#' @rdname lavaan.pl-class
#' @importFrom methods show
#' @export
setMethod("show", "lavaan.pl", function(object) {
  cat("lavaan.pl", get_lavaanpl_internal(object)$version, "\n")
  show(as(object, "lavaan"))
})
