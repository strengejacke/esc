#' @title Compute effect size from Mean and Standard Deviation
#' @name esc_mean_sd
#'
#' @description Compute effect size from mean and either group-based standard
#'              deviations or full sample standard deviation.
#'
#' @param grp1m The mean of the first group.
#' @param grp1sd The standard deviation of the first group.
#' @param grp1n The sample size of the first group.
#' @param grp2m The mean of the second group.
#' @param grp2sd The standard deviation of the second group.
#' @param grp2n The sample size of the second group.
#' @param totalsd The full sample standard deviation. Either \code{grp1sd} and
#'        \code{grp2sd}, or \code{totalsd} must be specified.
#'
#' @inheritParams esc_beta
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
#' @note If \code{es.type = "r"}, Fisher's transformation for the effect size
#'       \code{r} and their confidence intervals are also returned.
#'
#' @examples
#' # with standard deviations for each group
#' esc_mean_sd(grp1m = 7, grp1sd = 2, grp1n = 50,
#'             grp2m = 9, grp2sd = 3, grp2n = 60, es.type = "logit")
#'
#' # with full sample standard deviations
#' esc_mean_sd(grp1m = 7, grp1n = 50, grp2m = 9, grp2n = 60, totalsd = 4)
#'
#' @export
esc_mean_sd <- function(grp1m, grp1sd, grp1n, grp2m, grp2sd, grp2n, totalsd,
                        es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log"), study = NULL) {
  es.type <- match.arg(es.type)

  # check if parameter are complete
  if ((missing(totalsd) || is.null(totalsd) || is.na(totalsd)) &&
      ((missing(grp1sd) || is.null(grp1sd) || is.na(grp1sd)) ||
       (missing(grp2sd) || is.null(grp2sd) || is.na(grp2sd)))) {
    warning("Either `totalsd` or both `grp1sd` and `grp2sd` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # compute totaln, better overview
  totaln <- grp1n + grp2n

  # compute mean difference
  dm <- grp1m - grp2m

  # compute pooled standard deviation.
  if (!missing(totalsd) && !is.null(totalsd)) {

    # pooled sd from full sample sd, formula from book
    sdp <- ((totalsd ^ 2 * (totaln - 1) - ((dm ^ 2 * grp1n * grp2n) / totaln)) / (totaln - 1))

    # pooled sd from full sample sd, formula from unpublished manuscript. formulas vary,
    # email-correspondence with author suggests that book-formula should be correct
    # however, in some case value might be negative, so sqrt is not possible, use
    # alternative formula then
    if (sdp < 0)
      sdp <- (totalsd ^ 2 * (totaln - 1) - ((grp1m ^ 2 + grp2m ^ 2 - 2 * grp1m * grp2m) / totaln)) / totaln

    sd_pooled <- sqrt(sdp)

  } else {
    # pooled sd from group sd's
    sd_pooled <- sqrt((grp1sd ^ 2 * (grp1n - 1) + grp2sd ^ 2 * (grp2n - 1)) / (grp1n + grp2n - 2))
  }

  # compute effect size
  es <- (grp1m - grp2m) / sd_pooled
  # compute variance
  v <- esc.vd(es, grp1n, grp2n)

  # return effect size
  esc_generic(
    es = es,
    v = v,
    es.type = es.type,
    grp1n = grp1n,
    grp2n = grp2n,
    info = "mean and sd",
    study = study
  )
}


#' @title Compute effect size from Mean and Standard Error
#' @name esc_mean_se
#'
#' @description Compute effect size from Mean and Standard Error.
#'
#' @param grp1se The standard error of the first group.
#' @param grp2se The standard error of the second group.
#' @inheritParams esc_mean_sd
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
#' esc_mean_se(grp1m = 7, grp1se = 1.5, grp1n = 50,
#'             grp2m = 9, grp2se = 1.8, grp2n = 60, es.type = "or")
#'
#' @export
esc_mean_se <- function(grp1m, grp1se, grp1n, grp2m, grp2se, grp2n,
                        es.type = c("d", "g", "or", "logit", "r", "f", "eta", "cox.or", "cox.log"), study = NULL) {
  es.type <- match.arg(es.type)

  grp1sd <- grp1se * sqrt(grp1n - 1)
  grp2sd <- grp2se * sqrt(grp2n - 1)

  sd_pooled <- sqrt((grp1sd ^ 2 * (grp1n - 1) + grp2sd ^ 2 * (grp2n - 1)) / (grp1n + grp2n - 2))
  es <- (grp1m - grp2m) / sd_pooled
  v <- esc.vd(es, grp1n, grp2n)

  # return effect size
  return(esc_generic(es = es, v = v, es.type = es.type, grp1n = grp1n, grp2n = grp2n,
                     info = "mean and se", study = study))
}
