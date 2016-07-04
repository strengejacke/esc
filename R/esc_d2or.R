#' @title Convert effect size d into OR
#' @name esc_d2or
#'
#' @description Compute effect size \code{OR} from effect size \code{d}.
#'
#' @param d The effect size \code{d}.
#' @param se The standard error of \code{d}. One of \code{se} or \code{v}
#'        must be specified.
#' @param v The variance of \code{d}. One of \code{se} or \code{v} must be
#'        specified.
#' @param es.type Type of effect size odds ratio that should be returned.
#'        May be \code{es.type = "logit"} or \code{es.type = "cox"}
#'        (see 'Details').
#' @param info String with information on the transformation. Used for the
#'        print-method. Usually, this argument can be ignored
#'
#' @inheritParams esc_beta
#' @inheritParams hedges_g
#'
#' @note Effect size is returned as \code{exp(log_values)} (odds ratio),
#'       confidence intervals are also exponentiated. To get the log-values,
#'       use \code{\link{esc_d2logit}}.
#'       \strong{However}, variance and standard error of this function
#'       are returned on the log-scale!
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
#'             Cox DR. 1970. Analysis of binary data. New York: Chapman & Hall/CRC
#'             \cr \cr
#'             Hasselblad V, Hedges LV. 1995. Meta-analysis of screening and diagnostic tests. Psychological Bulletin 117(1): 167–178. \doi{10.1037/0033-2909.117.1.167}
#'
#' @examples
#' # d to odds ratio
#' esc_d2or(0.7, se = 0.5)
#' # odds ratio to d
#' esc_or2d(3.56, se = 0.91)
#'
#' @export
esc_d2or <- function(d, se, v, totaln, es.type = c("logit", "cox"), info = NULL) {
  es.type <- match.arg(es.type)

  # check if parameter are complete
  if ((missing(se) || is.null(se)) && (missing(v) || is.null(v))) {
    stop("Either `se` or `v` must be specified.", call. = F)
  }

  # do we have se?
  if (!missing(se) && !is.null(se)) v <- se ^ 2

  # do we have total n?
  if (missing(totaln)) totaln <- NULL

  # do we have a separate info string?
  if (is.null(info)) {
    info <- "effect size d to effect size OR"
    if (es.type == "cox") info <- paste0(info, "(Cox)")
  }

  if (es.type == "logit") {
    # Hasselblad and Hedges logit
    es <- pi / sqrt(3) * d
    v <- (pi  ^ 2) / 3 * v
    measure <- "or"
  } else {
    # Cox logit
    es <- d * 1.65
    v <- v / .367
    measure <- "cox-or"
  }

  # return effect size d
  return(structure(
    class = c("esc", "esc_d2or"),
    list(es = exp(es), se = sqrt(v), var = v, ci.lo = exp(lower_d(es, v)),
         ci.hi = exp(upper_d(es, v)), w = 1 / v, totaln = totaln,
         measure = measure, info = info)
  ))
}

