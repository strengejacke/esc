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
#'          \item{\code{"f"}}{returns effect size Cohen's \code{f}}
#'          \item{\code{"g"}}{returns adjusted standardized mean difference effect size Hedges' \code{g}}
#'          \item{\code{"or"}}{returns effect size as odds ratio}
#'          \item{\code{"cox.or"}}{returns effect size as Cox-odds ratio (see \code{\link{esc_d2or}} for details)}
#'          \item{\code{"logit"}}{returns effect size as log odds}
#'          \item{\code{"cox.log"}}{returns effect size as Cox-log odds (see \code{\link{esc_d2logit}} for details)}
#'          \item{\code{"r"}}{returns correlation effect size \code{r}}
#'          \item{\code{"eta"}}{returns effect size eta squared}
#'        }
#' @param study Optional string with the study name. Using \code{\link{combine_esc}} or
#'        \code{as.data.frame} on \code{esc}-objects will add this as column
#'        in the returned data frame.
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
#' esc_beta(.7, 3, 100, 150)
#' esc_beta(.7, 3, 100, 150, es.type = "cox.log")
#'
#' @export
esc_beta <- function(beta, sdy, grp1n, grp2n,
                     es.type = c("d", "g", "or", "logit", "r", "f", "eta", "cox.or", "cox.log"),
                     study = NULL) {
  # match  arguments
  es.type <- match.arg(es.type)

  totaln <- grp1n + grp2n
  sdx <- sqrt(abs((grp1n - (grp1n ^ 2 / totaln)) / (totaln - 1)))
  b <- beta * (sdy / sdx)

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
    info = "standardized regression coefficient",
    study = study
  )
}
