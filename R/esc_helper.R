# Compute variance of d-type effect size
esc.vd <- function(d, grp1n, grp2n) {
  (grp1n + grp2n) / (grp1n * grp2n) + (d * d) / (2 * (grp1n + grp2n))
}

# 95% confidence interval
lower_d <- function(d, v) d - 1.959964 * sqrt(v)
upper_d <- function(d, v) d + 1.959964 * sqrt(v)

# small sample size bias correction
sssbc <- function(totaln) return(1 - (3 / (4 * totaln - 9)))
