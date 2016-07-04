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
esc_B <- function(b, sdy, grp1n, grp2n, es.type = c("d", "or", "logit", "r", "cox.or", "cox.log")) {
  es.type <- match.arg(es.type)

  totaln <- grp1n + grp2n
  sdpooled <- sqrt(((sdy ^ 2 * (totaln - 1)) - (b ^ 2 * ((grp1n * grp2n) / (grp1n + grp2n)))) / (totaln - 2))
  es <- b / sdpooled
  v <- esc.vd(es, grp1n, grp2n)

  # return effect size
  return(esc_generic(es = es, v = v, es.type = es.type,
                     info = "unstandardized regression coefficient"))
}
