# Compute variance of d-type effect size
esc.vd <- function(d, grp1n, grp2n) {
  (grp1n + grp2n) / (grp1n * grp2n) + (d * d) / (2 * (grp1n + grp2n))
}

# 95% confidence interval
lower_d <- function(d, v) d - 1.959964 * sqrt(v)
upper_d <- function(d, v) d + 1.959964 * sqrt(v)

# small sample size bias correction
sssbc <- function(totaln) return(1 - (3 / (4 * totaln - 9)))

# convert correlation coefficient r into fisher's zr
esc.zr <- function(r) return(.5 * log((1 + r) / (1 - r)))

# Inverse Fisher's Zr transformation
esc.inv.zr <- function(zr) return((exp(2 * zr) - 1 ) / (exp(2 * zr) + 1))


# generic conversion function
esc_generic <- function(es, v, grp1n, grp2n, es.type, info) {
  # return effect size as odds ratio
  if (es.type == "or")
    return(esc_d2or(d = es, v = v, es.type = "logit",
                    info = paste0(info, " to effect size odds ratios")))

  # return effect size as cox odds ratio
  if (es.type == "cox.or")
    return(esc_d2or(d = es, v = v, es.type = "cox",
                    info = paste0(info, " to effect size Cox odds ratios")))

  # return effect size as log odds
  if (es.type == "logit")
    return(esc_d2logit(d = es, v = v, es.type = "logit",
                       info = paste0(info, " to effect size logits")))

  # return effect size as cox log odds
  if (es.type == "cox.log")
    return(esc_d2logit(d = es, v = v, es.type = "cox",
                       info = paste0(info, " to effect size Cox logits")))

  # return effect size as correlation
  if (es.type == "r")
    return(esc_d2r(d = es, v = v, grp1n = grp1n, grp2n = grp2n,
                   info = paste0(info, " to effect size correlation")))

  # return effect size as standardized mean difference d
  return(structure(
    class = c("esc", "esc_d"),
    list(es = es, se = sqrt(v), var = v, ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, measure = "d", info = paste0(info, " to effect size d")
  )))
}