#' @title Combine one or more 'esc' objects into a data frame
#' @name combine_esc
#'
#' @description This method takes one or more objects of class \code{esc} (which
#'              are returned by each effect size calculation function) and
#'              returns the combined result as a single data frame. This can
#'              then be used for further computation, e.g. with the
#'              \code{\link[metafor]{rma}}-function of the \pkg{metafor}-package.
#'
#' @param ... One or more objects of class \code{esc}
#'
#' @return A data frame with all relevant information from the effect size
#'         calculation.
#'
#' @examples
#' e1 <- esc_2x2(grp1yes = 30, grp1no = 50, grp2yes = 40,
#'               grp2no = 45, study = "Study 1")
#' e2 <- esc_2x2(grp1yes = 30, grp1no = 50, grp2yes = 40, grp2no = 45,
#'               es.type = "or", study = "Study 2")
#' e3 <- esc_t(p = 0.03, grp1n = 100, grp2n = 150, study = "Study 3")
#' e4 <- esc_mean_sd(grp1m = 7, grp1sd = 2, grp1n = 50, grp2m = 9, grp2sd = 3,
#'                   grp2n = 60, es.type = "logit", study = "Study 4")
#'
#' combine_esc(e1, e2, e3, e4)
#'
#' @export
combine_esc <- function(...) {
  # get input
  obj <- list(...)
  # do we have any object with correlation effect size?
  any.res <- unlist(lapply(obj, function(x) any(class(x) == "esc_d2r")))
  # if yes, add this class attribute to each object
  if (any(any.res)) {
    # iterate all esc-objects
    for (i in seq_len(length(obj))) {
      # add class attribute
      if (!any(class(obj[[i]]) == "esc_d2r")) {
        class(obj[[i]]) <- c(class(obj[[i]]), "esc_d2r")
        # fill non existing values
        obj[[i]]$zr <- NA
        obj[[i]]$ci.lo.zr <- NA
        obj[[i]]$ci.hi.zr <- NA
      }
    }
  }
  # create new data frame
  ret.val <- data.frame()
  # iterate objects again
  for (i in seq_len(length(obj))) {
    # bind converted data frame to each row
    ret.val <- rbind(ret.val, as.data.frame(obj[[i]]))
  }
  ret.val
}