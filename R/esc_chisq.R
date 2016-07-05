#' @title Compute effect size from Chi-Square coefficient
#' @name esc_chisq
#'
#' @description Compute effect size from Chi-Square coefficient
#'
#' @param chisq The chi-squared value. One of \code{chisq} or \code{p} must be reported.
#' @param p The p-value of the chi-squared or phi-value.
#'
#' @inheritParams esc_beta
#' @inheritParams hedges_g
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @note This effect size should only be used for data from 2x2 frequency
#'       tables. Furthermore, use this approximation for the effect size only,
#'       if information about the 2x2 frequencies or proportions are not available.
#'       Else, \code{\link{esc_2x2}} or \code{\link{esc_bin_prop}} provide better
#'       estimates for the effect size.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' # Effect size based on chi-squared value
#' esc_chisq(chisq = 9.9, totaln = 100)
#'
#' # Effect size based on p-value of chi-squared
#' esc_chisq(p = .04, totaln = 100)
#'
#' @importFrom stats qnorm
#' @export
esc_chisq <- function(chisq, p, totaln, es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log"), study = NULL) {
  # match  arguments
  es.type <- match.arg(es.type)

  # check if parameter are complete
  if ((missing(chisq) || is.null(chisq)) && (missing(p) || is.null(p))) {
    stop("Either `chisq` or `p` must be specified.", call. = F)
  }

  # if we have no phi-value, compute it from p.
  # divide p by two, because two-tailed.
  if (missing(chisq) || is.null(chisq)) chisq <- stats::qnorm(p / 2, lower.tail = F) ^ 2

  # compute effect size
  es <- 2 * sqrt(chisq / (totaln - chisq))

  # compute variancr
  v <- es ^ 2 / chisq

  # return effect size
  return(esc_generic(es = es, v = v, es.type = es.type, grp1n = totaln / 2, grp2n = totaln / 2,
                     info = "chi-squared-value", study = study))
}


#' @title Compute effect size from Phi coefficient
#' @name esc_phi
#'
#' @description Compute effect size from phi coefficient
#'
#' @param phi The phi value. One of \code{phi} or \code{p} must be reported.
#' @param p The p-value of the chi-squared or phi-value.
#'
#' @inheritParams esc_beta
#' @inheritParams hedges_g
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @note This effect size should only be used for data from 2x2 frequency
#'       tables. Furthermore, use this approximation for the effect size only,
#'       if information about the 2x2 frequencies or proportions are not available.
#'       Else, \code{\link{esc_2x2}} or \code{\link{esc_bin_prop}} provide better
#'       estimates for the effect size.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' # Effect size based on chi-squared value
#' esc_phi(phi = .67, totaln = 100)
#'
#' # Effect size based on p-value of chi-squared
#' esc_phi(p = .003, totaln = 100)
#'
#' @export
esc_phi <- function(phi, p, totaln, es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log"), study = NULL) {
  es.type <- match.arg(es.type)

  # check if parameter are complete
  if ((missing(phi) || is.null(phi)) && (missing(p) || is.null(p))) {
    stop("Either `phi` or `p` must be specified.", call. = F)
  }

  # if we have no phi-value, compute it from p.
  # divide p by two, because two-tailed.
  if (missing(phi) || is.null(phi)) return(esc_chisq(p = p, totaln = totaln, es.type = es.type, study = study))

  # compute effect size
  es <- (2 * phi ) / sqrt(1 - phi ^ 2)

  # compute chi-squared
  chisq <- phi ^ 2 * totaln

  # compute variance
  v = es ^ 2 / chisq

  # return effect size
  return(esc_generic(es = es, v = v, es.type = es.type, grp1n = totaln / 2, grp2n = totaln / 2,
                     info = "phi-value", study = study))
}
