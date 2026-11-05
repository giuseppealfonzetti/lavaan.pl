#' @rdname lavaan.pl-class
#' @exportS3Method lavaan::lavInspect
lavInspect.lavaan.pl <- function(object, ...) lavaan::lavInspect(as(object, "lavaan"), ...)

#' @rdname lavaan.pl-class
#' @exportS3Method lavaan::lavTech
lavTech.lavaan.pl <- function(object, ...) lavaan::lavTech(as(object, "lavaan"), ...)
