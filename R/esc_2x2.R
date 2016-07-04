#' @title Compute effect size from 2 by 2 Contingency Table
#' @name esc_2x2
#'
#' @description Compute effect size from a 2 by 2 frequency table.
#'
#' @param grp1yes Size of treatment group with successes (outcome = yes).
#' @param grp1no Size of treatment group with non-successes (outcome = no).
#' @param grp2yes Size of control group with successes (outcome = yes).
#' @param grp2no Size of control group with non-successes (outcome = no).
#'
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
#' # effect size log odds
#' esc_2x2(grp1yes = 30, grp1no = 50, grp2yes = 40, grp2no = 45)
#'
#' # effect size odds ratio
#' esc_2x2(grp1yes = 30, grp1no = 50, grp2yes = 40, grp2no = 45, es.type = "or")
#'
#' @export
esc_2x2 <- function(grp1yes, grp1no, grp2yes, grp2no, es.type = c("logit", "d", "or", "r", "cox.d")) {
  es.type <- match.arg(es.type)

  es <- (grp1yes * grp2no) / (grp1no * grp2yes)
  v <- 1 / grp1yes + 1 / grp1no + 1 / grp2yes + 1 / grp2no

  # which es type to be returned?
  if (es.type == "or") {
    return(structure(
      class = c("esc", "esc_2x2"),
      list(es = es, se = sqrt(v), var = v, ci.lo = exp(lower_d(log(es), v)), ci.hi = exp(upper_d(log(es), v)),
           w = 1 / v, measure = "or", info = "2x2 table (OR) coefficient to effect size odds ratio")
    ))
  }

  # which es type to be returned?
  if (es.type == "d") return(esc_or2d(or = es, v = v, es.type = "d", info = "2x2 table (OR) to effect size d"))

  # which es type to be returned?
  if (es.type == "cox.d") return(esc_or2d(or = es, v = v, es.type = "cox.d", info = "2x2 table (OR) to effect size Cox d"))

  # convert to plain d
  es <- log(es) * sqrt(3) / pi
  v <- v * 3 / pi ^ 2

  # which es type to be returned?
  if (es.type == "logit") return(esc_d2logit(d = es, v = v, info = "2x2 table (OR) to effect size logits"))

  # which es type to be returned?
  if (es.type == "r") return(esc_d2r(d = es, v = v, grp1n = (grp1yes + grp1no), grp2n = (grp2yes + grp2no), info = "2x2 table (OR) coefficient to effect size correlation"))
}
