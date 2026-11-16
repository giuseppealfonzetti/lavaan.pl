#' Plots of the optimisation
#'
#' Traces of the optimisation across epochs, a pass over all pairs, with a background band for each phase: Adam, then
#' averaging for `SA` or preconditioning and quasi-Newton for `SA-QN`, where each evaluation of the pairwise
#' likelihood counts as one epoch and only accepted steps are shown.
#'
#' @param x A `lavaan.pl` object from [lavaan.pl()].
#' @param type `"trace"` for the estimates by parameter block, `"objective"` for the negative pairwise log-likelihood
#'   per unit (at the running average in the averaged phase, evaluated when plotting), `"noise"` for the optimisation noise relative to the standard
#'   error of each parameter, with stochastic approximation only.
#' @param blocks Parameter blocks to show, such as `"lambda"` or `"tau"`. If `NULL`, all of them.
#' @param ... Unused.
#' @return The plotted data, invisibly. Without the tinyplot package, a warning and the data.
#' @seealso [diagnostics()]
#' @export
plot.lavaan.pl <- function(
  x,
  type = c("trace", "objective", "noise"),
  blocks = NULL,
  ...
) {
  fit <- get_lavaanpl_internal(x)
  type <- match.arg(type)
  sa <- fit$stoFit$path
  qn <- fit$numFit$path
  sa_epoch <- sa$iter / fit$stoFit$epoch_length
  epoch <- c(sa_epoch, max(0, sa_epoch) + qn$epoch - 1)
  phase <- factor(
    c(pmax(sa$phase, 1), rep(4, length(qn$epoch))),
    1:4,
    c("Adam", "Preconditioning", "Averaged", "QN")
  )
  long <- switch(
    type,
    trace = do.call(
      rbind,
      Map(
        function(epoch, phase, vector) {
          data.frame(
            epoch = epoch,
            phase = phase,
            estimate = transform_parameters(fit$model, vector),
            id = seq_along(vector),
            block = fit$blocks
          )
        },
        epoch,
        phase,
        c(sa$parameter_vector, qn$parameter_vector)
      )
    ),
    objective = {
      nll <- c(sa$nll, qn$nll)
      miss <- is.na(nll)
      nll[miss] <- vapply(
        c(sa$parameter_vector, qn$parameter_vector)[miss],
        function(v) {
          -cpp_pl(v, fit$cpp_inputs$model, fit$cpp_inputs$data)$pl / fit$data$n
        },
        numeric(1)
      )
      data.frame(epoch = epoch, phase = phase, objective = nll)
    },
    noise = if (fit$method == "SA") {
      fit_parameters(x)
    } else {
      stop("`x` was not fitted by stochastic approximation.", call. = FALSE)
    }
  )
  if (!is.null(blocks) && type != "objective") {
    long <- long[long$block %in% blocks, ]
  }
  if (!requireNamespace("tinyplot", quietly = TRUE)) {
    warning(
      "{tinyplot} must be installed to plot; the data are returned.",
      call. = FALSE
    )
    return(invisible(long))
  }
  if (type == "noise") {
    do.call(
      tinyplot::tinyplot,
      list(
        ratio ~ seq_along(ratio),
        data = long,
        facet = ~block,
        facet.args = list(free = TRUE),
        ylim = c(0, max(long$ratio)),
        type = "h",
        xlab = "Parameter index",
        ylab = "Optimisation noise / standard error",
        theme = "clean2"
      )
    )
    return(invisible(long))
  }
  right <- tapply(epoch, phase, max)
  shown <- !is.na(right)
  right <- unname(right[shown])
  bands <- data.frame(
    phase = factor(levels(phase)[shown], levels(phase)[shown]),
    left = c(0, right[-length(right)]),
    right = right
  )
  if (type == "trace") {
    bands <- merge(bands, data.frame(block = unique(long$block)))
  }
  y <- if (type == "trace") long$estimate else long$objective
  lim <- range(y) + c(-0.04, 0.04) * diff(range(y))
  tinyplot::tinyplot(
    x = bands$left,
    y = rep(lim[1], nrow(bands)),
    xmin = bands$left,
    xmax = bands$right,
    ymin = lim[1],
    ymax = lim[2],
    by = bands$phase,
    facet = bands$block,
    type = "rect",
    col = NA,
    fill = grDevices::palette.colors(4, "Tableau 10", alpha = 0.2)[shown],
    yaxs = "i",
    legend = list(title = NULL),
    xlab = "Epochs",
    ylab = if (type == "trace") {
      "Estimates"
    } else {
      "Negative pairwise log-likelihood per unit"
    },
    theme = "clean2"
  )
  switch(
    type,
    trace = do.call(
      tinyplot::tinyplot,
      list(
        estimate ~ epoch | factor(id),
        data = long,
        facet = ~block,
        type = "l",
        legend = FALSE,
        col = "black",
        add = TRUE
      )
    ),
    objective = do.call(
      tinyplot::tinyplot,
      list(
        objective ~ epoch,
        data = long,
        type = "b",
        pch = 19,
        col = "black",
        add = TRUE
      )
    )
  )
  invisible(long)
}
