#' @rdname lavaan.pl-class
#' @importFrom stats logLik
#' @export
setMethod("logLik", "lavaan.pl", function(object, ...) {
  fit <- get_latpl_fit(object)
  vars <- if (is.null(object@external$vars)) compute_var(fit) else object@external$vars
  structure(object@optim$logl, df = sum(diag(solve(vars$H, vars$J))), nobs = fit$dims$n, class = "logLik")
})
