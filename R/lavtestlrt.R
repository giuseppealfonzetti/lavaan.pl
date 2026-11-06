#' Likelihood ratio test for nested models
#'
#' @param object Object of class `lavaan`.
#' @param ... Further objects of class `lavaan`, and arguments for [lavaan::lavTestLRT()].
#' @param model_names Names of the models.
#' @return An object of class `anova`, as [lavaan::lavTestLRT()]. When all models are of class
#'   `lavaan.pl`, the difference tests are adjusted pairwise likelihood ratio tests.
#' @export
lavTestLRT <- function(object, ..., model_names = NULL) {
  dots <- list(...)
  modp <- vapply(dots, inherits, logical(1), "lavaan")
  mods <- c(list(object), dots[modp])
  if (is.null(model_names)) {
    mc <- as.list(match.call(expand.dots = TRUE))[-1]
    model_names <- vapply(mc[c(TRUE, modp)], function(x) deparse1(x), "")
  }
  ours <- length(mods) > 1 &&
    all(vapply(mods, inherits, logical(1), "lavaan.pl")) &&
    (is.null(dots$method) || dots$method %in% c("default", "mean.var.adjusted.PLRT")) &&
    (is.null(dots$type) || tolower(dots$type) == "chisq")
  if (!ours) {
    return(lavaan::lavTestLRT(object, ..., model_names = model_names))
  }

  names(mods) <- make.unique(model_names)
  mods <- mods[order(vapply(mods, function(x) x@test[[1]]$df, numeric(1)))]
  fits <- lapply(mods, get_latpl_fit)
  stat_delta <- df_delta <- rep(NA, length(mods))
  for (m in seq_len(length(mods) - 1)) {
    out <- compute_plrt(fits[[m + 1]], H1 = fits[[m]])$model
    stat_delta[m + 1] <- out$stat.scaled
    df_delta[m + 1] <- out$df.scaled
  }

  val <- data.frame(
    Df = vapply(mods, function(x) x@test[[1]]$df, numeric(1)),
    PL_AIC = vapply(mods, stats::AIC, numeric(1)),
    PL_BIC = vapply(mods, stats::BIC, numeric(1)),
    Chisq = vapply(mods, function(x) x@test[[1]]$stat, numeric(1)),
    `Chisq diff` = stat_delta,
    `Df diff` = df_delta,
    `Pr(>Chisq)` = stats::pchisq(stat_delta, df_delta, lower.tail = FALSE),
    row.names = names(mods),
    check.names = FALSE
  )
  attr(val, "heading") <- paste0("\nMean and Variance Adjusted Chi-Squared Difference Test (method = ",
                                 dQuote("mean.var.adjusted.PLRT"), ", lavaan.pl)\n")
  class(val) <- c("anova", class(val))
  val
}
