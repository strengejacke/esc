#' @title Convert effect size d into f
#' @name convert_d2f
#'
#' @description Compute effect size \code{f} from effect size \code{d}.
#'
#' @param d The effect size \code{d}.
#' @param se The standard error of \code{d}. One of \code{se} or \code{v}
#'        must be specified.
#' @param v The variance of \code{d}. One of \code{se} or \code{v} must be
#'        specified.
#' @param info String with information on the transformation. Used for the
#'        print-method. Usually, this argument can be ignored
#'
#' @inheritParams convert_d2or
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @references Cohen J. 1988. Statistical Power Analysis for the Behavioral Sciences. 2nd ed. Hillsdale, NJ: Erlbaum
#'
#' @examples
#' # d to f
#' convert_d2f(d = 0.2, se = .1, totaln = 50)
#'
#' @export
convert_d2f <- function(d, se, v, totaln, info = NULL, study = NULL) {
  # compute effect size f
  es <- d / 2

  # check if parameter are complete
  if ((missing(se) || is.null(se) || is.na(se)) && (missing(v) || is.null(v) || is.na(v))) {
    warning("Either `se` or `v` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = "f", grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # do we have se?
  if (!missing(se) && !is.null(se) && !is.na(se)) v <- se^2

  # return effect size f
  structure(
    class = c("esc", "convert_d2f"),
    list(
      es = es,
      se = sqrt(v),
      var = v,
      ci.lo = lower_d(es, v),
      ci.hi = upper_d(es, v),
      w = 1 / v,
      totaln = totaln,
      measure = "f",
      info = info,
      study = study
    )
  )
}
