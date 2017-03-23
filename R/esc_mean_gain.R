#' @title Compute effect size from Mean Gain Scores and Standard Deviations
#' @name esc_mean_gain
#'
#' @description Compute effect size from Mean Gain Scores and Standard
#'   Deviations for pre-post tests.
#'
#' @param pre1mean The mean of the first group at pre-test.
#' @param pre1sd The standard deviation of the first group at pre-test.
#' @param post1mean The mean of the first group at post-test.
#' @param post1sd The standard deviation of the first group at post-test.
#' @param grp1n The sample size of the first group.
#' @param gain1mean The mean gain between pre and post of the first group.
#' @param gain1sd The standard deviation gain between pre and post of the first
#'   group.
#' @param grp1r The (estimated) correlation of pre-post scores for the first
#'   group.
#' @param pre2mean The mean of the second group at pre-test.
#' @param pre2sd The standard deviation of the second group at pre-test.
#' @param post2mean The mean of the second group at post-test.
#' @param post2sd The standard deviation of the second group at post-test.
#' @param grp2n The sample size of the second group.
#' @param gain2mean The mean gain between pre and post of the second group.
#' @param gain2sd The standard deviation gain between pre and post of the second
#'   group.
#' @param grp2r The (estimated) correlation of pre-post scores for the second
#'   group.
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
#'   \code{r} and their confidence intervals are also returned.
#'
#' @details For this function, either the gain scores of mean and sd
#'   (\code{gain1mean} and \code{gain1sd} for the first group and
#'   \code{gain2mean} and \code{gain2sd} for the second group) must be
#'   specified, or the pre-post values (\code{pre1mean}, \code{post1mean},
#'   \code{pre1sd} and \code{post1sd} and the counterpart arguments for the
#'   second group).
#'   \cr \cr
#'   If the pre-post standard deviations are available, no correlation value
#'   \code{grp1r} resp. \code{grp2r} needs to be specified, because these can
#'   then be computed based on t-value computation. However, if \code{grp1r}
#'   is specified, this value will be used (and no t-test performed).
#'
#' @examples
#' # effect size of mean gain scores, with available pre-post values
#' esc_mean_gain(pre1mean = 13.07, pre1sd = 11.95, post1mean = 6.1,
#'               post1sd = 8.33, grp1n = 78, pre2mean = 10.77, pre2sd = 10.73,
#'               post2mean = 8.83, post2sd = 9.67, grp2n = 83)
#'
#' # same as above, but with assumed correlation of .5
#' # Note that effect size is the same, but variance differs
#' esc_mean_gain(pre1mean = 13.07, pre1sd = 11.95, post1mean = 6.1, grp1r = .5,
#'               post1sd = 8.33, grp1n = 78, pre2mean = 10.77, pre2sd = 10.73,
#'               post2mean = 8.83, post2sd = 9.67, grp2n = 83, grp2r = .5)
#'
#' # effect size based on gain scores for mean and sd. note that the
#' # pre-post correlations must be given
#' esc_mean_gain(gain1mean = 1.5, gain1sd = 1, grp1n = 40, grp1r = .5,
#'               gain2mean = .7, gain2sd = .8, grp2n = 50, grp2r = .5)
#'
#' @export
esc_mean_gain <- function(pre1mean, pre1sd, post1mean, post1sd, grp1n, gain1mean, gain1sd, grp1r,
                          pre2mean, pre2sd, post2mean, post2sd, grp2n, gain2mean, gain2sd, grp2r,
                          es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log"), study = NULL) {
  es.type <- match.arg(es.type)

  # check if arguments are complete
  if ((missing(gain1mean) || is.null(gain1mean) || is.na(gain1mean)) &&
      ((missing(pre1mean) || is.null(pre1mean) || is.na(pre1mean)) ||
       (missing(post1mean) || is.null(post1mean) || is.na(post1mean)))) {
    warning("Either `pre1mean` and `post1mean` or `gain1mean` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # check if arguments are complete
  if ((missing(gain2mean) || is.null(gain2mean) || is.na(gain2mean)) &&
      ((missing(pre2mean) || is.null(pre2mean) || is.na(pre2mean)) ||
       (missing(post2mean) || is.null(post2mean) || is.na(post2mean)))) {
    warning("Either `pre2mean` and `post2mean` or `gain2mean` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # check if arguments are complete
  if ((missing(gain1sd) || is.null(gain1sd) || is.na(gain1sd)) &&
      ((missing(pre1sd) || is.null(pre1sd) || is.na(pre1sd)) ||
       (missing(post1sd) || is.null(post1sd) || is.na(post1sd)))) {
    warning("Either `pre1sd` and `post1sd` or `gain1sd` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # check if arguments are complete
  if ((missing(gain2sd) || is.null(gain2sd) || is.na(gain2sd)) &&
      ((missing(pre2sd) || is.null(pre2sd) || is.na(pre2sd)) ||
       (missing(post2sd) || is.null(post2sd) || is.na(post2sd)))) {
    warning("Either `pre2sd` and `post2sd` or `gain2sd` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # check if arguments are complete
  if (!missing(pre1sd) && !is.null(pre1sd) && !is.na(pre1sd) && !missing(post1sd) && !is.null(post1sd) && !is.na(post1sd) && !missing(grp1r) && !is.null(grp1r) && !is.na(grp1r)) {
    message("Pre-post correlation `grp1r` is specified, although it could be computed from `pre1sd` and `post1sd`. If `grp1r` is missing, correlation of pre-post scores will be computed automatically.", call. = F)
  }

  # check if arguments are complete
  if (!missing(pre2sd) && !is.null(pre2sd) && !is.na(pre2sd) && !missing(post2sd) && !is.null(post2sd) && !is.na(post2sd) && !missing(grp2r) && !is.null(grp2r) && !is.na(grp2r)) {
    message("Pre-post correlation `grp2r` is specified, although it could be computed from `pre2sd` and `post2sd`. If `grp2r` is missing, correlation of pre-post scores will be computed automatically.", call. = F)
  }

  # check if arguments are complete
  if (((missing(pre1sd) || is.null(pre1sd) || is.na(pre1sd)) || (missing(post1sd) || is.null(post1sd) || is.na(post1sd))) && (missing(grp1r) || is.null(grp1r) || is.na(grp1r))) {
    warning("`grp1r` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # check if arguments are complete
  if (((missing(pre2sd) || is.null(pre2sd) || is.na(pre2sd)) || (missing(post2sd) || is.null(post2sd) || is.na(post2sd))) && (missing(grp2r) || is.null(grp2r) || is.na(grp2r))) {
    warning("`grp2r` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # compute r for group 1, based on t-test
  if (missing(grp1r)) {
    # compute t-value
    tgrp1 <- esc_compute_t(pre1mean, post1mean, pre1sd, post1sd, grp1n)
    # compute correlation
    grp1r <- ((pre1sd ^ 2 * tgrp1 ^ 2 + post1sd ^ 2 * tgrp1 ^ 2) - (post1mean - pre1mean) ^ 2 * grp1n) / (2 * pre1sd * post1sd * tgrp1 ^ 2)
  }

  # compute r for group 2, based on t-test
  if (missing(grp2r)) {
    # compute t-value
    tgrp2 <- esc_compute_t(pre2mean, post2mean, pre2sd, post2sd, grp2n)
    # compute correlation
    grp2r <- ((pre2sd ^ 2 * tgrp2 ^ 2 + post2sd ^ 2 * tgrp2 ^ 2) - (post2mean - pre2mean) ^ 2 * grp2n) / (2 * pre2sd * post2sd * tgrp2 ^ 2)
  }

  # compute sd for group 1
  if (missing(pre1sd) || missing(post1sd) || is.null(pre1sd) || is.null(post1sd))
    grp1sd <- gain1sd / sqrt(2 * (1 - grp1r))
  else
    grp1sd <- sqrt((pre1sd ^ 2 + post1sd ^ 2) / 2)

  # compute sd for group 2
  if (missing(pre2sd) || missing(post2sd) || is.null(pre2sd) || is.null(post2sd))
    grp2sd <- gain2sd / sqrt(2 * (1 - grp2r))
  else
    grp2sd <- sqrt((pre2sd ^ 2 + post2sd ^ 2) / 2)


  # compute mean gain scores for groups 1 and 2
  if (missing(gain1mean)) gain1mean <- pre1mean - post1mean
  if (missing(gain2mean)) gain2mean <- pre2mean - post2mean

  # compute pooled standard deviation
  sd_pooled <- sqrt((grp1sd ^ 2 * (grp1n - 1) + grp2sd ^ 2 * (grp2n - 1)) / (grp1n + grp2n - 2))

  # compute effect size d
  es <- (gain1mean - gain2mean) / sd_pooled

  # compute variance
  v <- (2 * (1 - grp1r)) / grp1n + (2 * (1 - grp2r)) / grp2n + es ^ 2 / (2 * (grp1n + grp2n))

  # return effect size
  return(esc_generic(es = es, v = v, grp1n = grp1n, grp2n = grp2n,
                     es.type = es.type, info = "mean gain score", study = study))
}


esc_compute_t <- function(m1, m2, s1, s2, n) {
  se <- sqrt((s1 ^ 2 / (n / 2)) + (s2 ^ 2 / (n / 2)))
  t <- (m1 - m2) / se
  return(t)
}