#' @title Convert effect size d into correlation
#' @name convert_d2r
#'
#' @description Compute effect size correlation from effect size \code{d}.
#'
#' @inheritParams esc_beta
#' @inheritParams convert_d2or
#' @inheritParams esc_mean_sd
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}. Furthermore, Fisher's z and
#'         confidence intervals are returned.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' convert_d2r(d = 0.7, se = 0.5, grp1n = 70, grp2n = 80)
#'
#' @export
convert_d2r <- function(d, se, v, grp1n, grp2n, info = NULL, study = NULL) {
  # check if parameter are complete
  if ((missing(se) || is.null(se) || is.na(se)) && (missing(v) || is.null(v) || is.na(v))) {
    warning("Either `se` or `v` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = NA, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # do we have se?
  if (!missing(se) && !is.null(se) && !is.na(se)) v <- se^2

  # do we have a separate info string?
  if (is.null(info)) info <- "effect size d to effect size correlation"

  p <- grp1n / (grp1n + grp2n)
  es <- d / sqrt(d^2 + 1 / (p * (1 - p)))
  v <- v / (v + 1 / (p * (1 - p)))

  # compute fisher's corrections for r
  es.zr <- convert_r2z(es)

  # return effect size d
  structure(
    class = c("esc", "convert_d2r"),
    list(
      es = es,
      se = sqrt(v),
      var = v,
      ci.lo = convert_z2r(lower_d(es.zr, v)),
      ci.hi = convert_z2r(upper_d(es.zr, v)),
      w = 1 / v,
      zr = es.zr,
      ci.lo.zr = lower_d(es.zr, v),
      ci.hi.zr = upper_d(es.zr, v),
      totaln = (grp1n + grp2n),
      measure = "r",
      info = info,
      study = study
    )
  )
}