#' @title Generate a table for the candidates obtained in \link{EvaluateCandidateListAgainstRawData}.
#'
#' @description \code{GenerateCandXLSX} will produce a XLSX of a list containing
#'     test results objects.
#'
#' @details Just a wrapper, to get the important information in a tabular layout.
#'
#' @param res_list A list of result objects (each testing an individual mz pair).
#' @param xlsx_file File name.
#' @param rejected Logical. Prepare table of rejected candidates if TRUE.
#'
#' @return Candidate table as data.frame.
#'
#' @examples
#' # load evaluation result of example data and
#' # generate table within R (use parameter xlsx_file to write to file)
#' x <- GenerateCandXLSX(HiResTEC::res_list)
#' str(x)
#' x[,1:5]
#'
#' @export
#'
GenerateCandXLSX <- function(res_list = NULL, xlsx_file = NULL, rejected = FALSE) {
  verify_suggested("openxlsx")
  # combine results without err_msg into a candidate list
  cand <- ldply_base(res_list, function(x) {
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
  })

  # write Excel output
  if (!is.null(xlsx_file) && nrow(cand) >= 1) {
    # tmp <- try(xlsx::write.xlsx(x=cand, file=xlsx_file, sheetName="cand", col.names=TRUE, row.names=FALSE))
    tmp <- try(openxlsx::write.xlsx(x = cand, file = xlsx_file, sheetName = "cand"))
    if (inherits(tmp, "try-error")) warning("Can't write to Excel-File. Seems to be opened.")
  }

  return(cand)
}
