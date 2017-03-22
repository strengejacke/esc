#' @title Generate effect size data frame from other data
#' @name effect_sizes
#'
#' @description This method computes any effect size from raw values from a data
#'              frame. Convenient method to compute multiple effect sizes at once.
#'
#' @param ... Named arguments. The name (left-hand side) is the name of one of
#'          \pkg{esc} functions' argument, the argument (right-hand side) is the
#'          name of the column in \code{data} that holds the data values.
#'          See 'Examples'.
#' @param data A data frame with columns that contain the values that are passed
#'          to one of the \pkg{esc}-functions.
#' @param type Name of one of the \pkg{esc}-functions, where arguments in \code{...}
#'          are passed to.
#' @inheritParams esc_beta
#'
#' @note This function rowwise iterates \code{data} and calls the function named
#'         in \code{type} for the values taken from each row of \code{data}. With
#'         this, it is possible to compute multiple effect size results with
#'         just one function call, if the necessary values are stored in a data
#'         frame.
#'
#' @return A data frame with the results of multiple effect size function calls.
#'
#' @examples
#' tmp <- data.frame(
#'   tvalue = c(3.3, 2.9, 2.3),
#'   n = c(250, 200, 210),
#'   studyname = c("Study 1", "Study 2", "Study 3")
#' )
#'
#' effect_sizes(t = tvalue, totaln = n, study = studyname, type = "esc_t", data = tmp)
#'
#'
#' tmp <- data.frame(
#'   tvalue = c(3.3, 2.9, NA, 2.3),
#'   n = c(250, 200, 210, 210),
#'   studyname = c("Study 1", "Study 2", NA, "Study 4")
#' )
#'
#' effect_sizes(t = tvalue, totaln = n, study = studyname, type = "esc_t", data = tmp)
#'
#'
#' tmp <- data.frame(
#'   coefficient = c(0.4, 0.2, 0.6),
#'   se = c(.15, .1, .2),
#'   treat = c(50, 60, 50),
#'   cntrl = c(45, 70, 40),
#'   author = c("Smith 2000", "Smith 2010 2", "Smith 2012")
#' )
#'
#' effect_sizes(beta = coefficient, sdy = se, grp1n = treat, grp2n = cntrl,
#'     study = author, type = "esc_beta", data = tmp, es.type = "or")
#'
#' @importFrom tibble as_tibble
#' @export
effect_sizes <- function(..., data, type, es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log")) {
  es.type <- match.arg(es.type)

  # get arguments
  params <- match.call(expand.dots = FALSE)$`...`

  # get argument names
  argnames <- names(params)

  # and get columns names
  cols <- unname(unlist(lapply(params, function(x) as.character(x))))

  results <- list()

  # iterate data frame rowwise
  for (i in 1:nrow(data)) {
    # get required data row
    datarow <- data[i, cols]
    # make sure we have no missing, else es can't be computed
    if (anyNA(datarow)) {
      warning(sprintf(
        "Row %i has missing values in columns '%s' and was dropped.",
        i,
        paste0(colnames(data)[which(is.na(datarow))], collapse = ", ")
      ),
      call. = F)
    } else {
      # copy values as list, for do.call
      callparams <- as.list(datarow)
      # give argument proper names
      names(callparams) <- argnames
      # call esc-method
      results[[length(results) + 1]] <- do.call(eval(type), args = c(callparams, es.type))
    }
  }

  # returned results as combined tibble
  results <- purrr::map_df(results, function(x) tibble::as_tibble(x))
  purrr::map_df(results, function(x) if (is.factor(x)) as.character(x) else x)
}