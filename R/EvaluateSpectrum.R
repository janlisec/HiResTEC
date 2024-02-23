#' @title EvaluateSpectrum.
#'
#' @description \code{EvaluateSpectrum} will extend a res-object (results)
#'     by information regarding mass spectra extracted from raw data files.
#'
#' @details If still present, rows containing NA in mz2 column will be filtered additionally.
#'
#' @param res res-object.
#' @param dat List of raw data files.
#' @param drt Time deviation. Passed to \link{DeconvoluteSpectrum}.
#' @param dmz Mass deviation. Passed to \link{DeconvoluteSpectrum}.
#' @param ionization Ion source. Passed to \link{DeconvoluteSpectrum}.
#' @param smooth Smoothing parameter. Passed to \link{DeconvoluteSpectrum}.
#'
#' @return Extended res-object.
#'
#' @keywords internal
#' @noRd
EvaluateSpectrum <- function(res = NULL, dat = NULL, drt = 2, dmz = 0.005, ionization = c("APCI", "ESI")[1], smooth = 0) {
  # deconvolute up to 10 spectra of early and late time point
  tp <- res[["tp"]]

  use.decon <- which(tp == min(tp, na.rm = T))
  res[["s"]] <- try(DeconvoluteSpectrum(dat = dat[use.decon], rt = res[["rt"]], mz1 = res[["mz1"]], use.mz.adjust = FALSE, rt_dev = drt, mz_dev = dmz, ionization = ionization, smooth = smooth))
  use.decon <- which(tp == max(tp, na.rm = T))
  res[["s2"]] <- try(DeconvoluteSpectrum(dat = dat[use.decon], rt = res[["rt"]], mz1 = res[["mz1"]] + res[["ng"]] * res[["mz_iso"]], use.mz.adjust = FALSE, rt_dev = drt, mz_dev = dmz, ionization = ionization, smooth = smooth))

  # perform some tests on the deconvoluted t_min spectrum
  if (nrow(res[["s"]]) == 0) {
    res[["err_msg"]] <- c(res[["err_msg"]], "deconvolution failed")
  } else {
    # find line of mz1
    l <- which(abs(res[["s"]][, "mz"] - res[["mz1"]]) < dmz)
    if (length(l) >= 2) l <- l[which.max(res[["s"]][l, "int"])]
    if (length(l) < 1) {
      res[["err_msg"]] <- c(res[["err_msg"]], "no mz1 in spectrum within dmz")
    } else {
      if (res[["s"]][l, "int"] < 0.1 * max(res[["s"]][, "int"])) {
        res[["err_msg"]] <- c(res[["err_msg"]], "mz1 is <10% of base peak")
      }
      # cluster s according to mz isotopic groups
      mass_groups <- InterpretMSSpectrum::GetGroupFactor(res[["s"]][, "mz"], 2)
      if (!max(res[["s"]][mass_groups == mass_groups[l], "int"]) %in% res[["s"]][l, "int"]) {
        res[["err_msg"]] <- c(res[["err_msg"]], "mz1 not main peak in MID at t_initial")
      }
    }
  }

  # return test result
  invisible(res)
}
