#' @title Compute effect size from One-way Anova
#' @name esc_f
#'
#' @description Compute effect size from One-way Anova with two independent groups.
#'
#' @param f The F-value of the F-test.
#'
#' @inheritParams esc_beta
#' @inheritParams esc_t
#'
#' @note This function only applies to \emph{one-way Anova} F-tests with
#'       \emph{two independent groups}, either equal or unequal sample sizes.
#'       \cr \cr
#'       If \code{es.type = "r"}, Fisher's transformation for the effect size
#'       \code{r} and their confidence intervals are also returned.
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower
#'         and upper confidence limits \code{ci.lo} and \code{ci.hi} as well as
#'         the weight factor \code{w}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' # unequal sample size
#' esc_f(f = 5.5, grp1n = 100, grp2n = 150)
#'
#' # equal sample size
#' esc_f(f = 5.5, totaln = 200)
#'
#' @export
esc_f <- function(f, totaln, grp1n, grp2n, es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log"), study = NULL) {
  return(esc_t(t = sqrt(f), totaln = totaln, grp1n = grp1n, grp2n = grp2n, es.type = es.type, study = study, info = "F-value (one-way-Anova)"))
}
