#' @title Compute effect size from binary proportions
#' @name esc_bin_prop
#'
#' @description Compute effect size from binary proportions
#'
#' @param prop1event Proportion of successes in treatment group (proportion of outcome = yes).
#' @param prop2event Proportion of successes in control group (proportion of outcome = yes).
#'
#' @inheritParams esc_beta
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @note If \code{es.type = "r"}, Fisher's transformation for the effect size
#'       \code{r} and their confidence intervals are also returned.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' # effect size log odds
#' esc_bin_prop(prop1event = .375, grp1n = 80, prop2event = .47, grp2n = 85)
#'
#' # effect size odds ratio
#' esc_bin_prop(prop1event = .375, grp1n = 80, prop2event = .47, grp2n = 85,
#'              es.type = "or")
#'
#' @export
esc_bin_prop <- function(prop1event, grp1n, prop2event, grp2n,
                         es.type = c("logit", "d", "g", "or", "r", "cox.d"),
                         study = NULL) {
  # match  arguments
  es.type <- match.arg(es.type)

  # compute 2x2 effect size
  return(esc_2x2(grp1yes = round(prop1event * grp1n), grp1no = round((1 - prop1event) * grp1n),
                 grp2yes = round(prop2event * grp2n), grp2no = round((1 - prop2event) * grp2n),
                 es.type = es.type, study = study))
}
