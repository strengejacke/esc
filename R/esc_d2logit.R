#' @title Convert effect size d into log odds
#' @name esc_d2logit
#'
#' @description Compute effect size \code{log odds} from effect size \code{d}.
#'
#' @inheritParams esc_d2or
#' @inheritParams esc_beta
#' @inheritParams hedges_g
#'
#' @note Effect size, variance, standard error and confidence intervals are
#'       returned on the log-scale. To get the odds ratios and exponentiated
#'       confidence intervals, use \code{\link{esc_d2or}}.
#'
#' @details Conversion from \code{d} to odds ratios can be done with two
#'          methods:
#'          \describe{
#'            \item{\code{es.type = "logit"}}{uses the Hasselblad and Hedges logit method.}
#'            \item{\code{es.type = "cox"}}{uses the modified logit method as proposed by Cox.
#'                  This method performs slightly better for rare or frequent events, i.e.
#'                  if the success rate is close to 0 or 1.}
#'          }
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'             \cr \cr
#'             Cox DR. 1970. Analysis of binary data. New York: Chapman & Hall/CRC
#'             \cr \cr
#'             Hasselblad V, Hedges LV. 1995. Meta-analysis of screening and diagnostic tests. Psychological Bulletin 117(1): 167–178. \doi{10.1037/0033-2909.117.1.167}
#'
#' @examples
#' # to logits
#' esc_d2logit(0.7, se = 0.5)
#'
#' # to Cox-logits
#' esc_d2logit(0.7, v = 0.25, es.type = "cox")
#'
#' @export
esc_d2logit <- function(d, se, v, totaln,
                        es.type = c("logit", "cox"),
                        info = NULL,
                        study = NULL) {
  # match  arguments
  es.type <- match.arg(es.type)

  # check if parameter are complete
  if ((missing(se) || is.null(se) || is.na(se)) && (missing(v) || is.null(v) || is.na(v))) {
    warning("Either `se` or `v` must be specified.", call. = F)
    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
  }

  # do we have se?
  if (!missing(se) && !is.null(se) && !is.na(se)) v <- se ^ 2

  # do we have total n?
  if (missing(totaln) || is.na(totaln)) totaln <- NULL

  # do we have a separate info string?
  if (is.null(info)) {
    info <- "effect size d to effect size logit"
    if (es.type == "cox") info <- paste0(info, "(Cox)")
  }

  if (es.type == "logit") {
    # Hasselblad and Hedges logit
    es <- pi / sqrt(3) * d
    v <- (pi  ^ 2) / 3 * v
    measure <- "logit"
  } else {
    # Cox logit
    es <- d * 1.65
    v <- v / .367
    measure <- "cox-logit"
  }

  # return effect size d
  return(structure(
    class = c("esc", "esc_d2logit"),
    list(es = es, se = sqrt(v), var = v, ci.lo = lower_d(es, v),
         ci.hi = upper_d(es, v), w = 1 / v, totaln = totaln,
         measure = measure, info = info, study = study)
  ))
}
