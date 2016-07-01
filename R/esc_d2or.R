#' @title Compute OR from effect size d
#' @name esc_d2or
#'
#' @description Compute effect size \code{OR} from effect size \code{d}.
#'
#' @param d The effect size \code{d}.
#' @param se The standard error of \code{d}. One of \code{se} or \code{v}
#'        must be specified.
#' @param v The variance of \code{d}. One of \code{se} or \code{v} must be
#'        specified.
#' @param info String with information on the transformation. Used for the
#'        print-method. Usually, this argument can be ignored
#' @inheritParams esc_beta
#'
#' @note Effect size is returned as \code{exp(log_values)} (odds ratio), i.e. variance,
#'       standard error and confidence intervals are also exponentiated. To
#'       get the log-values, use \code{\link{esc_d2logit}}.
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower
#'         and upper confidence limits \code{ci.lo} and \code{ci.hi} as well as
#'         the weight factor \code{w}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'
#' @examples
#' esc_d2or(0.7, se = 0.5)
#'
#' @export
esc_d2or <- function(d, se, v, info = NULL) {
  # check if parameter are complete
  if ((missing(se) || is.null(se)) && (missing(v) || is.null(v))) {
    stop("Either `se` or `v` must be specified.", call. = F)
  }

  # do we have se?
  if (!missing(se) && !is.null(se)) v <- se ^ 2

  # do we have a separate info string?
  if (is.null(info)) info <- "effect size d to effect size OR"

  es <- pi / sqrt(3) * d
  v <- (pi  ^ 2) / 3 * v

  # return effect size d
  return(structure(
    class = c("esc", "esc_d2or"),
    list(es = exp(es), se = exp(sqrt(v)), var = exp(v),
         ci.lo = exp(lower_d(es, v)), ci.hi = exp(upper_d(es, v)),
         w = 1 / v, info = info)
  ))
}


#' @title Compute logits from effect size d
#' @name esc_d2logit
#'
#' @description Compute effect size \code{log odds} from effect size \code{d}.
#'
#' @inheritParams esc_beta
#' @inheritParams esc_d2or
#'
#' @note Effect size, variance, standard error and confidence intervals are
#'       returned on the log-scale. To get the odds ratios and exponentiated
#'       values, use \code{\link{esc_d2or}}.
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower
#'         and upper confidence limits \code{ci.lo} and \code{ci.hi} as well as
#'         the weight factor \code{w}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'
#' @examples
#' esc_d2logit(0.7, se = 0.5)
#'
#' @export
esc_d2logit <- function(d, se, v, info = NULL) {
  # check if parameter are complete
  if ((missing(se) || is.null(se)) && (missing(v) || is.null(v))) {
    stop("Either `se` or `v` must be specified.", call. = F)
  }

  # do we have se?
  if (!missing(se) && !is.null(se)) v <- se ^ 2

  # do we have a separate info string?
  if (is.null(info)) info <- "effect size d to effect size logit"

  es <- pi / sqrt(3) * d
  v <- (pi  ^ 2) / 3 * v

  # return effect size d
  return(structure(
    class = c("esc", "esc_d2logit"),
    list(es = es, se = sqrt(v), var = v,
         ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, info = info)
  ))
}


#' @title Compute correlation from effect size d
#' @name esc_d2r
#'
#' @description Compute effect size correlation from effect size \code{d}.
#'
#' @inheritParams esc_beta
#' @inheritParams esc_d2or
#' @inheritParams esc_mean_sd
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower
#'         and upper confidence limits \code{ci.lo} and \code{ci.hi} as well as
#'         the weight factor \code{w}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'
#' @examples
#' esc_d2r(d = 0.7, se = 0.5, grp1n = 70, grp2n = 80)
#'
#' @export
esc_d2r <- function(d, se, v, grp1n, grp2n, info = NULL) {
  # check if parameter are complete
  if ((missing(se) || is.null(se)) && (missing(v) || is.null(v))) {
    stop("Either `se` or `v` must be specified.", call. = F)
  }

  # do we have se?
  if (!missing(se) && !is.null(se)) v <- se ^ 2

  # do we have a separate info string?
  if (is.null(info)) info <- "effect size d to effect size correlation"

  p <- grp1n / (grp1n + grp2n)
  es <- d / sqrt(d ^ 2 + 1 / (p * (1 - p)))
  v <- v / (v + 1 / (p * (1 - p)))

  # return effect size d
  return(structure(
    class = c("esc", "esc_d2r"),
    list(es = es, se = sqrt(v), var = v,
         ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, info = info)
  ))
}
