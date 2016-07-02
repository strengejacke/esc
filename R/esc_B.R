#' @title Compute effect size from Unstandardized Regression Coefficient
#' @name esc_B
#'
#' @description Compute effect size from Unstandardized Regression Coefficient.
#'
#' @param b The unstandardized coefficient B.
#' @inheritParams esc_beta
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower
#'         and upper confidence limits \code{ci.lo} and \code{ci.hi} as well as
#'         the weight factor \code{w}.
#'
#' @note If \code{es.type = "r"}, Fisher's transformation for the effect size
#'       \code{r} and their confidence intervals are also returned.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'
#' @examples
#' esc_B(3.3, 5, 100, 150)
#'
#' @export
esc_B <- function(b, sdy, grp1n, grp2n, es.type = c("d", "OR", "logit", "r")) {
  es.type <- match.arg(es.type)

  totaln <- grp1n + grp2n
  sdpooled <- sqrt(((sdy ^ 2 * (totaln - 1)) - (b ^ 2 * ((grp1n * grp2n) / (grp1n + grp2n)))) / (totaln - 2))
  es <- b / sdpooled
  v <- esc.vd(es, grp1n, grp2n)

  # which es type to be returned?
  if (es.type == "OR") return(esc_d2or(d = es, v = v, info = "unstandardized regression coefficient to effect size odds ratio"))

  # which es type to be returned?
  if (es.type == "logit") return(esc_d2logit(d = es, v = v, info = "unstandardized regression coefficient to effect size logit"))

  # which es type to be returned?
  if (es.type == "r") return(esc_d2r(d = es, v = v, grp1n = grp1n, grp2n = grp2n, info = "unstandardized regression coefficient to effect size correlation"))

  # return effect size d
  return(structure(
    class = c("esc", "esc_B"),
    list(es = es, se = sqrt(v), var = v, ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, info = "unstandardized regression coefficient to effect size d")
  ))
}
