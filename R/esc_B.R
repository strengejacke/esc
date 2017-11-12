#' @title Compute effect size from Unstandardized Regression Coefficient
#' @name esc_B
#'
#' @description Compute effect size from Unstandardized Regression Coefficient.
#'
#' @param b The unstandardized coefficient B.
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
#' esc_B(3.3, 5, 100, 150)
#'
#' @export
esc_B <- function(b, sdy, grp1n, grp2n,
                  es.type = c("d", "g", "or", "logit", "r", "f", "eta", "cox.or", "cox.log"),
                  study = NULL) {
  # match  arguments
  es.type <- match.arg(es.type)

  totaln <- grp1n + grp2n
  sdpooled <- sqrt(abs(((sdy ^ 2 * (totaln - 1)) - (b ^ 2 * ((grp1n * grp2n) / (grp1n + grp2n)))) / (totaln - 2)))
  es <- b / sdpooled
  v <- esc.vd(es, grp1n, grp2n)

  # return effect size
  esc_generic(
    es = es,
    v = v,
    es.type = es.type,
    grp1n = grp1n,
    grp2n = grp2n,
    info = "unstandardized regression coefficient",
    study = study
  )
}
