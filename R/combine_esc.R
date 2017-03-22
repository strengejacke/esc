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
#' @seealso \code{\link{write_esc}}
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
#' @importFrom purrr map map_df
#' @importFrom sjmisc remove_empty_cols is_empty empty_cols
#' @export
combine_esc <- function(...) {
  # get input
  obj <- list(...)

  # make sure we have same columns for all effects. for correlation effect sizes,
  # we have columns zr, ci.lo.zr and ci.hi.zr
  obj <- lapply(obj, function(x) {
    if (!inherits(x, "esc_d2r")) {
      class(x) <- c(class(x), "esc_d2r")
      # fill non existing values
      x$zr <- NA
      x$ci.lo.zr <- NA
      x$ci.hi.zr <- NA
    }

    x
  })

  # remove empty columns if we did not have any correlations
  result <- do.call(rbind, purrr::map(obj, ~as.data.frame(.x)))

  # remove empty cols, if any
  if (!sjmisc::is_empty(sjmisc::empty_cols(result)))
    result <- sjmisc::remove_empty_cols(result)

  # make factor columns character
  purrr::map_df(result, function(x) if (is.factor(x)) as.character(x) else x)
}


#' @title Write one or more 'esc' objects into an Excel csv-file
#' @name write_esc
#'
#' @description This method is a small wrapper to write csv-files. It writes
#'              the results from \code{\link{combine_esc}} into an Excel csv-file.
#'
#' @param ... One or more objects of class \code{esc}
#' @param path Path to write to, or just file name (to write to working directory).
#' @param sep The field separator string. In some Western European locales, Excel
#'            uses a semicolon by default, while in other locales the field
#'            separator string in Excel is a comma.
#'
#' @note For Western European locales, the \code{sep}-argument probably needs to
#'       be set to semicolon (\code{sep = ";"}), so Excel reads the csv-file properly.
#'       If \code{sep = ";"}, \code{\link[utils]{write.csv2}} is used to write the
#'       file. Else, \code{\link[readr]{write_excel_csv}} is used.
#'
#' @return Invisibly returns the combined data frame that is written to
#'         the csv-file (see \code{\link{combine_esc}}).
#'
#' @seealso \code{\link{combine_esc}}
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
#' # write to current working directory,
#' # file extension ".csv" is automatically added
#' \dontrun{
#' write_esc(e1, e2, e3, e4, path = "EffSizes")}
#'
#' @importFrom readr write_excel_csv
#' @importFrom utils write.csv2
#' @export
write_esc <- function(..., path, sep = ",") {
  # check if file extension exists
  has.extension <- (regexpr("\\.[^\\.]*$", path) != -1)
  if (!has.extension) path <- paste0(path, ".csv")

  # check fir valid file extention
  dot.start <- regexpr("\\.[^\\.]*$", path) + 1
  ext <- tolower(substring(path, dot.start, nchar(path)))

  # correct wrong file extension
  if (ext != "csv") {
    warning("Data can only be written to csv-files. Changing file extension to `.csv`.", call. = F)
    path <- paste0(substring(path, 0, last = dot.start - 1), "csv")
  }

  # combine esc-results
  x <- combine_esc(...)

  # tell user what's going on...
  cat("Writing Excel csv-file to:\n")
  cat(normalizePath(path = path, winslash = "/", mustWork = F))

  # write to excel
  if (sep == ";")
    utils::write.csv2(x, file = path, fileEncoding = "UTF-8")
  else
    readr::write_excel_csv(x = x, path = path)

  # return data frame
  invisible(x)
}