#' @title Compute effect size from Standardized Regression Coefficient
#' @name esc_beta
#'
#' @description Compute effect size from Standardized Regression Coefficient.
#'
#' @param beta The standardized beta coefficient.
#' @param sdy The standard deviation of the dependent variable.
#' @param grp1n Treatment group sample size.
#' @param grp2n Control group sample size.
#' @param es.type Type of effect size that should be returned. By default,
#'        this is \code{"d"}, i.e. effect size \code{d} is returned.
#'        Use \code{es.type = "OR"} to return effect size as odds ratios,
#'        \code{es.type = "logit"} to return effect size as log odds and
#'        \code{es.type = "r"} to return effect size as correlation.
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
#' esc_beta(1.25, 3, 100, 150)
#' esc_beta(1.25, 3, 100, 150, es.type = "logit")
#'
#' @export
esc_beta <- function(beta, sdy, grp1n, grp2n, es.type = c("d", "OR", "logit", "r")) {
  es.type <- match.arg(es.type)

  totaln <- grp1n + grp2n
  sdx <- sqrt((grp1n - (grp1n ^ 2 / totaln)) / (totaln - 1))
  b <- beta * (sdy / sdx)

  sdpooled <- sqrt(((sdy ^ 2 * (totaln - 1)) - (b ^ 2 * ((grp1n * grp2n) / (grp1n + grp2n)))) / (totaln - 2))
  es <- b / sdpooled
  v <- esc.vd(es, grp1n, grp2n)

  # which es type to be returned?
  if (es.type == "OR") return(esc_d2or(d = es, v = v, info = "standardized regression coefficient to effect size odds ratios"))

  # which es type to be returned?
  if (es.type == "logit") return(esc_d2logit(d = es, v = v, info = "standardized regression coefficient to effect size logits"))

  # which es type to be returned?
  if (es.type == "r") return(esc_d2r(d = es, v = v, grp1n = grp1n, grp2n = grp2n, info = "standardized regression coefficient to effect size correlation"))

  return(structure(
    class = c("esc", "esc_beta"),
    list(es = es, se = sqrt(v), var = v, ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, info = "standardized regression coefficient to effect size d")
  ))
}
