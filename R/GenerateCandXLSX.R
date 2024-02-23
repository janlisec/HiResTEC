#' @title GenerateCandXLSX.
#'
#' @description
#' \code{GenerateCandXLSX} will produce a XLSX of a list containing test results objects.
#'
#' @details
#' Not yet.
#'
#' @param res_list A list of result objects (each testing an individual mz pair).
#' @param xlsx_file File name.
#' @param rejected Logical. Return rejected if TRUE.
#'
#' @return
#' Candidate table as data.frame.
#'
#' @examples
#' # load evaluation result of example data
#' data(res_list)
#' # generate table within R (use xlsx_file to write to file)
#' str(GenerateCandXLSX(res_list))
#' GenerateCandXLSX(res_list)[, 1:5]
#'
#' @export
#'
#' @importFrom openxlsx write.xlsx
#' @importFrom plyr ldply
#'
GenerateCandXLSX <- function(res_list = NULL, xlsx_file = NULL, rejected = FALSE) {
  # combine results without err_msg into a candidate list
  cand <- plyr::ldply(res_list, function(x) {
    # is this a candidate (==no err msg)
    test <- is.null(x[["err_msg"]])
    # if rejected are requested invert logic
    if (rejected) test <- !test
    if (test) {
      data.frame(
        "ID" = x[["cand_id"]],
        "RT" = x[["rt"]],
        "mz" = x[["mz1"]],
        "dE" = x[["dE"]],
        "P" = x[["P_raw"]],
        "Name" = x[["FluxLib"]],
        "row" = x[["row"]],
        "BaseFormula" = paste0("C", x[["ng"]]),
        "Spectrum" = paste(apply(x[["s"]], 1, function(y) {
          paste(round(y[1], 4), y[2], sep = ":")
        }), collapse = " "),
        "OverlappingLine" = x[["OverlappingLine"]],
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  }, .id = NULL)

  # write Excel output
  if (!is.null(xlsx_file) && nrow(cand) >= 1) {
    # tmp <- try(xlsx::write.xlsx(x=cand, file=xlsx_file, sheetName="cand", col.names=TRUE, row.names=FALSE))
    tmp <- try(openxlsx::write.xlsx(x = cand, file = xlsx_file, sheetName = "cand"))
    if (inherits(tmp, "try-error")) warning("Can't write to Excel-File. Seems to be opened.")
  }

  return(cand)
}
