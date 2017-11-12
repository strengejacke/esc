#' @title Convert effect size d into Eta Squared
#' @name esc_d2etasq
#'
#' @description Compute effect size Eta Squared from effect size \code{d}.
#'
#' @param d The effect size \code{d}.
#' @param se The standard error of \code{d}. One of \code{se} or \code{v}
#'        must be specified.
#' @param v The variance of \code{d}. One of \code{se} or \code{v} must be
#'        specified.
#' @param info String with information on the transformation. Used for the
#'        print-method. Usually, this argument can be ignored
#'
#' @inheritParams esc_d2r
#' @inheritParams esc_d2or
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @references Cohen J. 1988. Statistical Power Analysis for the Behavioral Sciences. 2nd ed. Hillsdale, NJ: Erlbaum
#'
#' @examples
#' # d to eta squared
#' esc_d2etasq(d = 0.7, se = 0.5, grp1n = 70, grp2n = 80)
#'
#' @export
esc_d2etasq <- function(d, se, v, grp1n, grp2n, info = NULL, study = NULL) {

  # check if parameter are complete
  if ((missing(se) || is.null(se) || is.na(se)) && (missing(v) || is.null(v) || is.na(v))) {
    warning("Either `se` or `v` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = "eta", grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # do we have se?
  if (!missing(se) && !is.null(se) && !is.na(se)) v <- se ^ 2

  # do we have a separate info string?
  if (is.null(info)) info <- "effect size d to effect size eta squared"

  p <- grp1n / (grp1n + grp2n)
  es <- eta_squared(d = d)
  v <- v / (v + 1 / (p * (1 - p)))

  # return effect size f
  structure(
    class = c("esc", "esc_d2etasq"),
    list(
      es = es,
      se = sqrt(v),
      var = v,
      ci.lo = exp(lower_d(es, v)),
      ci.hi = exp(upper_d(es, v)),
      w = 1 / v,
      totaln = grp1n + grp2n,
      measure = "etasq",
      info = info,
      study = study
    )
  )
}

