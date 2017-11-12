#' @title Compute effect size from Point-Biserial Correlation
#' @name esc_rpb
#'
#' @description Compute effect size from Point-Biserial Correlation.
#'
#' @param r The point-biserial r-value. One of \code{r} or \code{p} must be specified.
#' @param p The p-value of the point-biserial correlation. One of \code{r} or \code{p} must be specified.
#'
#' @inheritParams esc_beta
#' @inheritParams esc_t
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' # unequal sample size
#' esc_rpb(r = .3, grp1n = 100, grp2n = 150)
#'
#' # equal sample size
#' esc_rpb(r = .3, totaln = 200)
#'
#' # unequal sample size, with p-value
#' esc_rpb(p = 0.03, grp1n = 100, grp2n = 150)
#'
#' # equal sample size, with p-value
#' esc_rpb(p = 0.03, totaln = 200)
#'
#' @importFrom stats qt
#' @export
esc_rpb <- function(r, p, totaln, grp1n, grp2n, es.type = c("d", "g", "or", "logit", "f", "eta", "cox.or", "cox.log"), study = NULL) {
  es.type <- match.arg(es.type)

  # check if parameter are complete
  if ((missing(r) || is.null(r) || is.na(r)) && (missing(p) || is.null(p) || is.na(p))) {
    warning("Either `r` or `p` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # check if parameter are complete
  if ((missing(totaln) || is.null(totaln) || is.na(totaln)) &&
      ((missing(grp1n) || is.null(grp1n) || is.na(grp1n)) ||
       (missing(grp2n) || is.null(grp2n) || is.na(grp2n)))) {
    warning("Either `totaln` or both `grp1n` and `grp2n` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # if we have no total sample size, compute it from group sizes
  if (missing(totaln) || is.null(totaln) || is.na(totaln)) {
    totaln <- grp1n + grp2n
    equal.size <- F
  } else {
    equal.size <- T
  }

  # init t, only not NULL, when p is given
  t <- NULL

  # if we have no t-value, compute it from p.
  # divide p by two, because two-tailed.
  if (missing(r) || is.null(r) || is.na(r))
    t <- stats::qt(p = p / 2, df = totaln - 2, lower.tail = F)

  # if t is not NULL, p is given
  if (!is.null(t)) {
    # equal sample size?
    if (equal.size) {
      es <- 2 * t / sqrt(totaln)
      # manually get values for each group
      grp1n <- grp2n <- totaln / 2
    } else {
      # effect size for unequal sample sizes
      es <- t * sqrt((grp1n + grp2n) / (grp1n * grp2n))
    }
  } else {
    # r is given
    # equal sample size?
    if (equal.size) {
      es <- 2 * r / sqrt(1 - r ^ 2)
      # manually get values for each group
      grp1n <- grp2n <- totaln / 2
    } else {
      # effect size for unequal sample sizes
      pr <- grp1n / (grp1n + grp2n)
      es <- r / sqrt((1 - r ^ 2) * (pr * (1 - pr)))
    }
  }

  # get variance
  v <- esc.vd(es, grp1n, grp2n)

  # return effect size
  esc_generic(
    es = es,
    v = v,
    es.type = es.type,
    grp1n = grp1n,
    grp2n = grp2n,
    info = "point-biserial r",
    study = study
  )
}
