#' @title Convert effect size d into OR
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
#' @note Effect size is returned as \code{exp(log_values)} (odds ratio),
#'       confidence intervals are also exponentiated. To get the log-values,
#'       use \code{\link{esc_d2logit}}.
#'       \cr \cr
#'       \strong{However}, variance and standard error of this function
#'       are returned on the log-scale!
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
#' esc_or2d(3.56, se = 0.91)
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
    list(es = exp(es), se = sqrt(v), var = v,
         ci.lo = exp(lower_d(es, v)), ci.hi = exp(upper_d(es, v)),
         w = 1 / v, info = info)
  ))
}


#' @title Convert effect size OR from d
#' @name esc_or2d
#'
#' @description Compute effect size \code{d} from effect size \code{OR}.
#'
#' @param or The effect size as odds ratio.
#'
#' @inheritParams esc_beta
#' @inheritParams esc_d2or
#'
#' @note While \code{or} is the exponentiated log odds, the variance or standard
#'       error need to be on the log-scale!
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower
#'         and upper confidence limits \code{ci.lo} and \code{ci.hi} as well as
#'         the weight factor \code{w}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'
#' @examples
#' esc_or2d(3.56, se = 0.91)
#' esc_d2or(0.7, se = 0.5)
#'
#' @export
esc_or2d <- function(or, se, v, info = NULL) {
  # check if parameter are complete
  if ((missing(se) || is.null(se)) && (missing(v) || is.null(v))) {
    stop("Either `se` or `v` must be specified.", call. = F)
  }

  # do we have se?
  if (!missing(se) && !is.null(se)) v <- se ^ 2

  # do we have a separate info string?
  if (is.null(info)) info <- "effect size OR to effect size d"

  es <- log(or) / (pi / sqrt(3))
  v <- v / ((pi  ^ 2) / 3)

  # return effect size d
  return(structure(
    class = c("esc", "esc_or2d"),
    list(es = es, se = sqrt(v), var = v,
         ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, info = info)
  ))
}


#' @title Convert effect size d into log odds
#' @name esc_d2logit
#'
#' @description Compute effect size \code{log odds} from effect size \code{d}.
#'
#' @inheritParams esc_beta
#' @inheritParams esc_d2or
#'
#' @note Effect size, variance, standard error and confidence intervals are
#'       returned on the log-scale. To get the odds ratios and exponentiated
#'       confidence intervals, use \code{\link{esc_d2or}}.
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


#' @title Convert effect size d into correlation
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
#'         the weight factor \code{w}. Furthermore, Fisher's transformation for
#'         the effect size \code{r} and their confidence intervals are returned.
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

  # compute fisher's corrections for r
  es.zr <- esc.zr(es)

  # return effect size d
  return(structure(
    class = c("esc", "esc_d2r"),
    list(es = es, se = sqrt(v), var = v,
         ci.lo = esc.inv.zr(lower_d(es.zr, v)), ci.hi = esc.inv.zr(upper_d(es.zr, v)),
         w = 1 / v, zr = es.zr, ci.lo.zr = lower_d(es.zr, v), ci.hi.zr = upper_d(es.zr, v),
         info = info)
  ))
}
