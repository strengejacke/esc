#' @title Compute effect size from Student's t-test
#' @name esc_t
#'
#' @description Compute effect size from Student's t-test for \emph{independent samples}.
#'
#' @param t The t-value of the t-test. One of \code{t} or \code{p} must be specified.
#' @param p The p-value of the t-test. One of \code{t} or \code{p} must be specified.
#' @param totaln Total sample size. Either \code{totaln}, or \code{grp1n} and
#'        \code{grp2n} must be specified.
#'
#' @inheritParams esc_beta
#'
#' @note This function only applies to \emph{independent sample} t-tests, either
#'       equal or unequal sample sizes. It can't be used for t-values from
#'       dependent or paired t-tests, or t-values from other statistical procedures
#'       (like regressions).
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower
#'         and upper confidence limits \code{ci.lo} and \code{ci.hi} as well as
#'         the weight factor \code{w}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'
#' @examples
#' esc_t(t = 3.3, grp1n = 100, grp2n = 150)
#'
#' @export
esc_t <- function(t, p, totaln, grp1n, grp2n, es.type = c("d", "OR", "logit", "r")) {
  es.type <- match.arg(es.type)

  # check if parameter are complete
  if ((missing(t) || is.null(t)) && (missing(p) || is.null(p))) {
    stop("Either `t` or `p` must be specified.", call. = F)
  }

  # check if parameter are complete
  if ((missing(totaln) || is.null(totaln)) &&
      ((missing(grp1n) || is.null(grp1n)) ||
       (missing(grp2n) || is.null(grp2n)))) {
    stop("Either `totaln` or both `grp1n` and `grp2n` must be specified.", call. = F)
  }

  # if we have no t-value, compute it from p
  if (missing(t) || is.null(t))
    t <- qt(p = p, df = totaln - 2, lower.tail = F)

  # equal sample size?
  if (!missing(totaln)) {
    es <- 2 * t / sqrt(totaln)
    # manually get values for each group
    grp1n <- grp2n <- totaln / 2
  } else {
    # effect size for unequal sample sizes
    es <- t * sqrt((grp1n + grp2n) / (grp1n * grp2n))
  }
  # get variance
  v <- esc.vd(es, grp1n, grp2n)

  # which es type to be returned?
  if (es.type == "OR") return(esc_d2or(d = es, v = v, info = "t-value to effect size odds ratio"))

  # which es type to be returned?
  if (es.type == "logit") return(esc_d2logit(d = es, v = v, info = "t-value to effect size logit"))

  # which es type to be returned?
  if (es.type == "r") return(esc_d2r(d = es, v = v, grp1n = grp1n, grp2n = grp2n, info = "t-value to effect size correlation"))

  # return effect size d
  return(structure(
    class = c("esc", "esc_t"),
    list(es = es, se = sqrt(v), var = v, ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, info = "t-value to effect size d")
  ))
}