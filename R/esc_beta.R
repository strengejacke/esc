#' @title Compute effect size from Standardized Regression Coefficient
#' @name esc_beta
#'
#' @description Compute effect size from Standardized Regression Coefficient.
#'
#' @param beta The standardized beta coefficient.
#' @param sdy The standard deviation of the dependent variable.
#' @param grp1n Treatment group sample size.
#' @param grp2n Control group sample size.
#' @param es.type Type of effect size that should be returned.
#'        \describe{
#'          \item{\code{"d"}}{returns standardized mean difference effect size \code{d}}
#'          \item{\code{"g"}}{returns adjusted standardized mean difference effect size Hedge's \code{g}}
#'          \item{\code{"or"}}{returns effect size as odds ratio}
#'          \item{\code{"cox.or"}}{returns effect size as Cox-odds ratio (see \code{\link{esc_d2or}} for details)}
#'          \item{\code{"logit"}}{returns effect size as log odds}
#'          \item{\code{"cox.log"}}{returns effect size as Cox-log odds (see \code{\link{esc_d2logit}} for details)}
#'          \item{\code{"r"}}{returns correlation effect size \code{r}}
#'        }
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
#'
#' @examples
#' esc_beta(1.25, 3, 100, 150)
#' esc_beta(1.25, 3, 100, 150, es.type = "cox.log")
#'
#' @export
esc_beta <- function(beta, sdy, grp1n, grp2n, es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log")) {
  es.type <- match.arg(es.type)

  totaln <- grp1n + grp2n
  sdx <- sqrt((grp1n - (grp1n ^ 2 / totaln)) / (totaln - 1))
  b <- beta * (sdy / sdx)

  sdpooled <- sqrt(((sdy ^ 2 * (totaln - 1)) - (b ^ 2 * ((grp1n * grp2n) / (grp1n + grp2n)))) / (totaln - 2))
  es <- b / sdpooled
  v <- esc.vd(es, grp1n, grp2n)

  # return effect size
  return(esc_generic(es = es, v = v, es.type = es.type, grp1n = grp1n, grp2n = grp2n,
                     info = "standardized regression coefficient"))
}
