#' @title Compute effect size from Mean and Standard Deviation
#' @name esc_mean_sd
#'
#' @description Compute effect size from Mean and Standard Deviation.
#'
#' @param grp1m The mean of the first group.
#' @param grp1sd The standard deviation of the first group.
#' @param grp1n The sample size of the first group.
#' @param grp2m The mean of the second group.
#' @param grp2sd The standard deviation of the second group.
#' @param grp2n The sample size of the second group.
#' @param es.type Type of effect size that should be returned. By default,
#'        this is \code{"d"}, i.e. effect size \code{d} is returned.
#'        Use \code{es.type = "OR"} to return effect size as odds ratios or
#'        \code{es.type = "logit"} to return effect size as log odds,
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower
#'         and upper confidence limits \code{ci.lo} and \code{ci.hi} as well as
#'         the weight factor \code{w}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'
#' @note If \code{es.type = "r"}, Fisher's transformation for the effect size
#'       \code{r} and their confidence intervals are also returned.
#'
#' @examples
#' esc_mean_sd(grp1m = 7, grp1sd = 2, grp1n = 50,
#'             grp2m = 9, grp2sd = 3, grp2n = 60, es.type = "logit")
#'
#' @export
esc_mean_sd <- function(grp1m, grp1sd, grp1n, grp2m, grp2sd, grp2n, es.type = c("d", "OR", "logit", "r")) {
  es.type <- match.arg(es.type)

  sd_pooled <- sqrt((grp1sd ^ 2 * (grp1n - 1) + grp2sd ^ 2 * (grp2n - 1)) / (grp1n + grp2n - 2))
  es <- (grp1m - grp2m) / sd_pooled
  v <- esc.vd(es, grp1n, grp2n)

  # which es type to be returned?
  if (es.type == "OR") return(esc_d2or(d = es, v = v, info = "mean and sd to effect size odds ratio"))

  # which es type to be returned?
  if (es.type == "logit") return(esc_d2logit(d = es, v = v, info = "mean and sd to effect size logits"))

  # which es type to be returned?
  if (es.type == "r") return(esc_d2r(d = es, v = v, grp1n = grp1n, grp2n = grp2n, info = "mean and sd to effect size correlation"))

  return(structure(
    class = c("esc", "esc_mean_sd"),
    list(es = es, se = sqrt(v), var = v, ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, info = "mean and sd to effect size d")
  ))
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
#' esc_mean_se(grp1m = 7, grp1se = 1.5, grp1n = 50,
#'             grp2m = 9, grp2se = 1.8, grp2n = 60, es.type = "OR")
#'
#' @export
esc_mean_se <- function(grp1m, grp1se, grp1n, grp2m, grp2se, grp2n, es.type = c("d", "OR", "logit", "r")) {
  es.type <- match.arg(es.type)

  grp1sd <- grp1se * sqrt(grp1n - 1)
  grp2sd <- grp2se * sqrt(grp2n - 1)

  sd_pooled <- sqrt((grp1sd ^ 2 * (grp1n - 1) + grp2sd ^ 2 * (grp2n - 1)) / (grp1n + grp2n - 2))
  es <- (grp1m - grp2m) / sd_pooled
  v <- esc.vd(es, grp1n, grp2n)

  # which es type to be returned?
  if (es.type == "OR") return(esc_d2or(d = es, v = v, info = "mean and se to effect size odds ratio"))

  # which es type to be returned?
  if (es.type == "logit") return(esc_d2logit(d = es, v = v, info = "mean and se to effect size logits"))

  # which es type to be returned?
  if (es.type == "r") return(esc_d2r(d = es, v = v, grp1n = grp1n, grp2n = grp2n, info = "mean and se to effect size correlation"))

  return(structure(
    class = c("esc", "esc_mean_se"),
    list(es = es, se = sqrt(v), var = v, ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, info = "mean and se to effect size d")
  ))
}
