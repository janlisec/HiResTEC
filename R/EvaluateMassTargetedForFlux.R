#' @title EvaluateMassTargetedForFlux.
#'
#' @description
#' \code{EvaluateMassTargetedForFlux} will work on mz, rt information and evaluate a peak targeted for flux.
#'
#' @details This function .
#'
#' @param xg xcmsSet object with group information.
#' @param dat Raw data. List of xcmsRaw's.
#' @param gmd GMD.
#' @param tp Timepoint (numeric information).
#' @param gr Line/Group (class information).
#' @param mz Mass of M+H.
#' @param rt Retention time.
#' @param dmz_M0 Allowed mz deviation in Dalton to find M0 within xcmsSet.
#' @param drt_M0 Allowed rt deviation in seconds to find M0 within xcmsSet.
#' @param dmz Allowed mz deviation in Dalton for a corresponding candidate relative to M0.
#' @param drt Allowed rt deviation in seconds for a corresponding candidate relative to M0.
#' @param mz_iso The isotopic mass shift under investigation (e.g. 1.003355 for 13C experiments)
#' @param ng Limit analysis to a selected isotope? Specify ng=3 to specifically check M0 vs M+3. Leave NULL to check all present.
#' @param mfrow Layout for flute plotting. NULL for automatic layout.
#' @param pdf_file PDF output.
#' @param xlsx_file XLSX output.
#' @param testing testing (cf \link{EvaluateCandidateListAgainstRawData}).
#'
#' @return
#' Nothing.
#'
#' @keywords internal
#' @noRd
EvaluateMassTargetedForFlux <- function(xg = NULL, dat = NULL, tp = NULL, gr = NULL, mz = NULL, rt = NULL, dmz_M0 = 0.01, drt_M0 = 10, dmz = 0.05, drt = 0.5, mz_iso = 1.003355, ng = NULL, mfrow = NULL, pdf_file = "SinglePeak.pdf", xlsx_file = "SinglePeak.xlsx", testing = FALSE) {
  specific_row <- which(abs(xg@groups[, "mzmed"] - mz) < dmz_M0 & abs(xg@groups[, "rtmed"] - rt) < drt_M0)
  res <- NULL
  if (length(specific_row) >= 1) {
    if (length(specific_row) > 1) {
      print("Several M0 peaks were found for this mz/rt pair.")
      gv <- groupval(xg = xg)
      print(cbind(xg@groups[specific_row, c("mzmed", "rtmed", "npeaks")], "medianInt" = apply(gv[specific_row, ], 1, median, na.rm = T)))
      specific_row <- specific_row[which.min(abs(xg@groups[specific_row, "mzmed"] - mz))]
      print("Keep best mz.")
      print(cbind(xg@groups[specific_row, c("mzmed", "rtmed", "npeaks"), drop = F], "medianInt" = apply(gv[specific_row, , drop = FALSE], 1, median, na.rm = T)))
    }
    out <- EvaluatePairsFromXCMSSet(xg = xg, tp = tp, gr = gr, drt = drt, dmz = dmz, mz_iso = mz_iso, n = ifelse(is.null(ng), 6, ng), specific_row = specific_row)
    if (nrow(out) >= 1) {
      if (is.null(ng)) {
        res <- EvaluateCandidateListAgainstRawData(x = out, tp = tp, gr = gr, dat = dat)
        GenerateQCPlots(res_list = res, mfrow = mfrow, pdf_file = pdf_file)
        GenerateCandXLSX(res_list = res, xlsx_file = xlsx_file)
      }
      if (is.numeric(ng) && ng %in% round(apply(out[, c("mz1", "mz2")], 1, diff))) {
        res <- EvaluateCandidateListAgainstRawData(x = out[which(round(apply(out[, c("mz1", "mz2")], 1, diff)) == ng), , drop = F], tp = tp, gr = gr, dat = dat)
        GenerateQCPlots(res, mfrow = mfrow)
      } else {
        print("A peak was found for this mz/rt pair but not the specified matching isotope peak.")
      }
    } else {
      print("A peak was found for this mz/rt pair but no corresponding/matching peak.")
    }
  } else {
    print("No peak was found for this mz/rt pair.")
  }
  invisible(res)
}
