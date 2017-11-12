#' @title Compute effect size from Student's t-test
#' @name esc_t
#'
#' @description Compute effect size from Student's t-test for \emph{independent samples}.
#'
#' @param t The t-value of the t-test. One of \code{t} or \code{p} must be specified.
#' @param p The p-value of the t-test. One of \code{t} or \code{p} must be specified.
#' @param totaln Total sample size. Either \code{totaln}, or \code{grp1n} and
#'        \code{grp2n} must be specified.
#' @param ... Other parameters, passed down to further functions. For internal
#'        use only, can be ignored.
#'
#' @inheritParams esc_beta
#'
#' @note This function only applies to \emph{independent sample} t-tests, either
#'       equal or unequal sample sizes. It can't be used for t-values from
#'       dependent or paired t-tests, or t-values from other statistical procedures
#'       (like regressions).
#'       \cr \cr
#'       If \code{es.type = "r"}, Fisher's transformation for the effect size
#'       \code{r} and their confidence intervals are also returned.
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
#' esc_t(t = 3.3, grp1n = 100, grp2n = 150)
#'
#' # equal sample size
#' esc_t(t = 3.3, totaln = 200)
#'
#' # unequal sample size, with p-value
#' esc_t(p = 0.03, grp1n = 100, grp2n = 150)
#'
#' # equal sample size, with p-value
#' esc_t(p = 0.03, totaln = 200)
#'
#' @importFrom stats qt
#' @export
esc_t <- function(t, p, totaln, grp1n, grp2n, es.type = c("d", "g", "or", "logit", "r", "f", "eta", "cox.or", "cox.log"), study = NULL, ...) {
  es.type <- match.arg(es.type)

  # evaluate ellipses
  ell <- match.call(expand.dots = FALSE)$`...`
  # do we have info string?
  if (!is.null(ell$info))
    info <- ell$info
  else
    # if not, set default
    info <- "t-value"

  # check if parameter are complete
  if ((missing(t) || is.null(t) || is.na(t)) && (missing(p) || is.null(p) || is.na(p))) {
    warning("Either `t` or `p` must be specified.", call. = F)
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

  # if we have no t-value, compute it from p.
  # divide p by two, because two-tailed.
  if (missing(t) || is.null(t) || is.na(t))
    t <- stats::qt(p = p / 2, df = totaln - 2, lower.tail = F)

  # for t-test, directly estimate effect size
  if (es.type == "r") {
    es <- t / sqrt(t ^ 2 + totaln - 2)
    es.zr <- esc_r2z(es)
    v <- 1 / (totaln - 3)
    # return effect size d
    return(structure(
      class = c("esc", "esc_t2r"),
      list(es = es, se = sqrt(v), var = v, totaln = totaln,
           ci.lo = esc_z2r(lower_d(es.zr, v)), ci.hi = esc_z2r(upper_d(es.zr, v)),
           w = 1 / v, zr = es.zr, ci.lo.zr = lower_d(es.zr, v), ci.hi.zr = upper_d(es.zr, v),
           measure = "r", info = paste0(info, " to effect size correlation"), study = study)
    ))
  }

  # equal sample size?
  if (equal.size) {
    es <- 2 * t / sqrt(totaln)
    # manually get values for each group
    grp1n <- grp2n <- totaln / 2
  } else {
    # effect size for unequal sample sizes
    es <- t * sqrt((grp1n + grp2n) / (grp1n * grp2n))
  }
  # get variance
  v <- esc.vd(es, grp1n, grp2n)

  # return effect size
  esc_generic(
    es = es, v = v, es.type = es.type, grp1n = grp1n,
    grp2n = grp2n, info = info, study = study
  )
}
