#' @title CalcEnrichment.
#' @description \code{CalcEnrichment} will calculate from a vector of isotope intensities the enrichment of 13C.
#' @param x Vector of intensities.
#' @param sig Round enrichment to this precision.
#' @param robust Provide an intensity threshold that at least one ion has to exceed. Otherwise NA is returned.
#' @return The enrichment, defined as the ratio of 13C/total_C.
#' @keywords internal
#' @noRd
CalcEnrichment <- function(x = c(100, 10, 1), sig = 4, robust = 0) {
  if (any(is.finite(x))) {
    if (!any(x > robust, na.rm = T)) {
      return(NA)
    } else {
      n <- length(x) - 1
      return(100 * round(sum(seq(0, n) * (x / n), na.rm = T) / sum(x, na.rm = T), sig))
    }
  } else {
    return(NA)
  }
}

#' @title CandidateBoxplot.
#' @description \code{CandidateBoxplot} will plot enrichment of tracer incorporationn candidates.
#' @details Plot will be annotated with mz/rt information plus ANOVA(lm()) P-values.
#' @param res Candidate evaluation result.
#' @return A list with all important information including deconvoluted spectrum and linear model.
#' @keywords internal
#' @noRd
CandidateBoxplot <- function(res = NULL) {

  verify_suggested("beeswarm")

  # get tp info from res object
  tp <- res[["tp"]]
  tp_fac <- as.factor(tp)

  # requires sam$pchs and sam$cols
  if (!exists("sam") | (exists("sam") && nrow(sam) != length(tp)) | (exists("sam") && !all(c("pchs", "cols") %in% colnames(sam)))) {
    sam <- data.frame("cols" = grDevices::rainbow(length(levels(tp_fac)))[as.numeric(tp_fac)], "pchs" = rep(c(21, 22, 24, 25), length.out = length(levels(res[["gr"]])))[as.numeric(res[["gr"]])], stringsAsFactors = FALSE)
  }

  enrg <- attr(res[["enr"]], "Enrichment")
  # filter estimation by boxplot outlier calculation
  flt <- !(enrg %in% c(NA, graphics::boxplot(enrg ~ factor(tp), plot = FALSE)$out))
  ids <- which(flt)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(par(op))
  graphics::par(mfrow = c(1, 2))
  graphics::par(mar = c(5, 4, 4, 2) + 0.1)

  # beeswarm plot + annotation
  tmp.x <- beeswarm::beeswarm(enrg[flt] ~ tp[flt], method = "square", pwpch = sam[flt, "pchs"], pwbg = sam[flt, "cols"], main = paste("Candidate", res[["cand_id"]]), cex = 3, ylab = "Enrichment (calculated for all isotopes from mz1 to mz1+n)", xlab = "Timepoint", axes = F)[, c("x", "y")]
  graphics::text(x = tmp.x, labels = unlist(split(ids, tp[flt])))
  graphics::axis(1, at = 1:length(unique(tp[flt])), labels = sort(unique(tp[flt])))
  graphics::axis(2)
  graphics::box()
  # browser()
  graphics::mtext(text = paste("mz1 =", round(res[["mz1"]], 4)), side = 3, line = -1.2 * 1, adj = 0.02)
  graphics::mtext(text = paste("n =", res[["ng"]]), side = 3, line = -1.2 * 2, adj = 0.02)
  graphics::mtext(text = paste("RT =", round(res[["rt"]], 1)), side = 3, line = -1.2 * 3, adj = 0.02)
  graphics::mtext(text = paste("dE =", res[["dE"]]), side = 3, line = -1.2 * 5, adj = 0.02)
  graphics::mtext(text = paste("row =", res[["row"]]), side = 3, line = -1.2 * 6, adj = 0.02)
  y <- median(attr(res[["enr"]], "Enrichment")[tp == min(tp)], na.rm = T)
  if (is.finite(y)) graphics::axis(side = 2, at = y, labels = "", tcl = 1, line = NULL)
  y <- median(attr(res[["enr"]], "Enrichment")[tp == max(tp)], na.rm = T)
  if (is.finite(y)) graphics::axis(side = 4, at = y, labels = "", tcl = 1, line = NULL)

  # boxplot with anova results
  if (!is.null(res[["err_msg"]])) graphics::mtext(text = paste("Reject", ifelse(length(res[["err_msg"]]) > 1, paste0("(", length(res[["err_msg"]]), ")"), ""), "=", res[["err_msg"]][1]), side = 3, line = -1.2 * 8, adj = 0.02)
  if (all(res[["enr"]][1, flt] > 0, na.rm = T)) {
    plot(enrg[flt] ~ res[["inter"]][flt], col = sapply(split(sam[flt, "cols"], res[["inter"]][flt]), function(x) {
      ifelse(length(x) >= 1, unique(x), 0)
    }), ylab = "", main = res[["FluxLib"]], xlab = "Interaction Group*Time")
    if (inherits(res[["lm"]], "lm")) {
      anova_res <- stats::anova(res[["lm"]])
      for (k in 1:(nrow(anova_res) - 1)) {
        graphics::mtext(text = paste0("P_", rownames(anova_res)[k], " = ", formatC(anova_res[k, 5], format = "e", digits = 2)), side = 3, line = -1.2 * k, adj = 0.02)
      }
    }
  } else {
    plot(1, 1, ann = F)
  }
  invisible(NULL)
}

#' @title DetermineElementComposition.
#' @description \code{DetermineElementComposition} will compute the Mass
#'    Distribution Vectors of isotopologues based on measured intensities
#'    at [13C]=1\%.
#' @details Not yet.
#' @param int Measured intensities.
#' @param mz Mass ob peak.
#' @param ionization Currently only APCI is supported.
#' @return The most likely combinations of C and Si atoms for this formula.
#' @keywords internal
#' @noRd
DetermineElementComposition <- function(int = NULL, mz = 300, ionization = "APCI") {

  verify_suggested("CorMID")

  # set nC and nSi range based on mz
  # ToDo
  nC_rng <- 1:ceiling(mz / 12)

  # set reference MID
  ref <- c(100, rep(0, length(int) - 1))

  # find best fitting MID brute force
  if (ionization == "APCI") {
    test <- lapply(nC_rng, function(nCbio) {
      sapply(0:min(c(floor(nCbio / 3), 8)), function(nSi) {
        mid <- CorMID::CorMID(int = int, fml = paste0("C", nCbio + 3 * nSi, "Si", nSi))
        return(sqrt(sum((mid - ref[1:length(mid)])^2)))
      })
    })
  }
  # return respective nC and nSi
  nC <- which.min(sapply(test, min))
  nSi <- which.min(test[[nC]])
  return(paste0("C", nC, "Si", nSi))
}

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

#' @title ExportHeuristicsTable.
#' @description \code{ExportHeuristicsTable}.
#' @param res_list1 res_list1.
#' @param res_list2 res_list2.
#' @param out1 out1.
#' @param out2 out2.
#' @param p_thr p_thr.
#' @param xls_name xls_name.
#' @keywords internal
#' @noRd
ExportHeuristicsTable <- function(res_list1 = NULL, res_list2 = NULL, out1 = NULL, out2 = NULL, p_thr = 0.01, xls_name = "Tables.xlsx") {

  verify_suggested("openxlsx")

  # helper function
  MakeResListDF <- function(res_list1 = NULL, out1 = NULL, name = "res_list") {
    candidates1 <- sapply(res_list1, function(x) {
      x[["cand_id"]] != ""
    })
    table(candidates1)
    good_rows1 <- sapply(res_list1[candidates1], function(x) {
      y <- x[["row"]]
      names(y) <- x[["cand_id"]]
      return(y)
    })
    raw_all <- sum(out1$P <= p_thr, na.rm = T)
    raw_cnd <- sum(rownames(out1)[which(out1$P <= p_thr)] %in% good_rows1)
    raw_rej <- sum(rownames(out1)[which(out1$P <= p_thr)] %in% sapply(res_list1[!candidates1], function(x) {
      x[["row"]]
    }))
    tmp <- data.frame("x" = c(
      paste0(nrow(out1), " (", raw_all, ")"),
      paste0(length(res_list1), " (", raw_cnd + raw_rej, ")"),
      paste0(sum(candidates1), " (", raw_cnd, ")"),
      paste0(sum(!candidates1), " (", raw_rej, ")"),
      paste0(nrow(out1) - length(res_list1), " (", raw_all - raw_cnd - raw_rej, ")")
    ))
    colnames(tmp) <- name
    attr(tmp, "good_rows") <- good_rows1
    return(tmp)
  }

  # Tab.1
  if (!is.null(res_list1)) {
    nr1 <- ifelse(is.null(names(res_list1)), "res_list1", names(res_list1)[1])
    r1 <- MakeResListDF(res_list1 = res_list1, out1 = out1, name = nr1)
    total1 <- length(res_list1) - sum(sapply(res_list1, function(x) {
      is.null(x[["err_msg"]])
    }))
    tmp1 <- ldply_base(sort(table(unlist(sapply(res_list1, function(x) {
      x[["err_msg"]]
    })))), function(x) {
      data.frame("abs" = x, "rel" = round(100 * x / total1, 2))
    })
  } else {
    r1 <- NULL
    tmp1 <- NULL
  }

  if (!is.null(res_list2)) {
    if (is.null(out2)) out2 <- out1
    nr2 <- ifelse(is.null(names(res_list2)), "res_list2", names(res_list2)[1])
    r2 <- MakeResListDF(res_list1 = res_list2, out1 = out2, name = nr2)
    total2 <- length(res_list2) - sum(sapply(res_list2, function(x) {
      is.null(x[["err_msg"]])
    }))
    tmp2 <- ldply_base(sort(table(unlist(sapply(res_list2, function(x) {
      x[["err_msg"]]
    })))), function(x) {
      data.frame("abs" = x, "rel" = round(100 * x / total2, 2))
    })
  } else {
    r2 <- NULL
    tmp2 <- NULL
  }

  fract <- data.frame("Fraction" = c("all", "*tested", "**candidates", "**rejected", "*untested"))
  descr <- data.frame("Description" = c(
    "number of mass pairs detected in preCL (number of significant results at p>0.05 without raw data evaluation)",
    "number of evaluated mass pairs from preCL: evaCL (w/o QC)",
    "number of positively evaluated mass pairs (w/o QC)",
    "number of rejected mass pairs (w/o QC)",
    "number of untested mass pairs which correlate with a candidate (w/o QC)"
  ))
  openxlsx::write.xlsx(x = cbind(fract, r1, r2, descr), file = xls_name, sheetName = "Tab.1")

  # Tab.2 (Error message heuristic)
  tmp <- data.frame("QC test" = unique(tmp1[, 1], tmp2[, 1]), stringsAsFactors = F, check.names = F)
  if (!is.null(tmp1)) {
    tmp[tmp[, 1] %in% tmp1[, 1], "APCI method"] <- paste0(tmp1[, 2], " (", round(tmp1[, 3], 1), "%)")
    colnames(tmp)[ncol(tmp)] <- nr1
  }
  if (!is.null(tmp1)) {
    tmp[tmp[, 1] %in% tmp2[, 1], "Test method"] <- paste0(tmp2[, 2], " (", round(tmp2[, 3], 1), "%)")
    colnames(tmp)[ncol(tmp)] <- nr2
  }
  openxlsx::write.xlsx(x = tmp, file = xls_name, sheetName = "Tab.2", append = TRUE, row.names = F)

  invisible(list(attr(r1, "good_rows"), attr(r2, "good_rows")))
}

#' @title groupval.
#'
#' @description
#' \code{groupval} is a dirty rewrite of the same function available in xcms.
#'
#' @details Comment on 2019-05-14 by JL: 'xcms' needs to be put to suggest due to mzR
#'     on request by CRAN --> replaced groupval from xcms by own function.
#'
#' @param xg xcmsSet object.
#'
#' @keywords internal
#' @noRd
groupval <- function(xg = NULL) {
  gv <- matrix(NA, nrow = nrow(xg@groups), ncol = nrow(xg@phenoData))
  for (i in 1:nrow(gv)) {
    idx <- xg@peaks[xg@groupidx[[i]], "sample"]
    if (!any(duplicated(idx))) {
      gv[i, idx] <- xg@peaks[xg@groupidx[[i]], "maxo"]
    } else {
      flt <- !duplicated(idx)
      gv[i, idx[flt]] <- xg@peaks[xg@groupidx[[i]][flt], "maxo"]
    }
  }
  return(gv)
}

#' @title RankCandidateList.
#' @description \code{RankCandidateList} will reorder a candidate list of
#'     mz pairs according to maximum intensity.
#' @details If still present, rows containing NA in mz2 column will be filtered
#'     additionally.
#' @param x Candidate dataframe.
#' @return Reordered candidate data.frame.
#' @keywords internal
#' @noRd
RankCandidateList <- function(x) {
  x <- x[!is.na(x[, "mz2"]), , drop = F]
  if (nrow(x) >= 2) {
    int_cols <- grep("^I[12]", colnames(x))
    num_groups <- length(int_cols) / 2
    col_mz1 <- 1:num_groups
    col_mz2 <- (num_groups + 1):(2 * num_groups)
    tmp <- x[, int_cols, drop = FALSE]
    tmp[is.na(tmp)] <- 0

    # next line contains the simple version of just adding up intensities
    # tmp <- tmp[,col_mz1,drop=FALSE]+tmp[,col_mz2,drop=FALSE]

    # version including normalization of comparison groups to highest mass to avoid absolute increase/decrease effects over time
    # !! this would work within comparisons regarding the same mz1, but potentially some M+3/M+4 candidates might get higher ranking than the one containing M+0
    # tmp <- lapply(1:num_groups, function(i) {tmp[,c(i,i+num_groups)]})
    # tmp <- sapply(tmp, function(x) { apply(x/max(x), 1, sum) })

    # normalize peaks to median of mz1 at t=min
    # this will scale up all mz1 values to the same value and scale all mz2 values according to the correction factor obtained for mz1
    # thereby we take into account absolute increase/decrease effects over time
    tp <- as.numeric(sapply(strsplit(colnames(tmp), "_"), function(x) {
      x[2]
    }))
    max_ints_tmin <- apply(tmp[, intersect(which(tp == min(tp)), col_mz1), drop = F], 1, max, na.rm = T)
    cor_mat <- tmp[, col_mz1, drop = F] / max_ints_tmin
    tmp <- tmp[, col_mz1, drop = F] / cor_mat + tmp[, col_mz2, drop = F] / cor_mat # the first summation coefficient could be dropped (adds only a constant)

    # compute order out of this
    tmp_order <- order(1 - is.na(x[, "P"]), apply(tmp, 1, max), decreasing = TRUE, na.last = TRUE)
    x <- x[tmp_order, ]
  }
  return(x)
}

#' @title verify_suggested.
#' @description Check if packages are available and stop function otherwise.
#' @param pkg Package names to be checked.
#' @return NULL.
#' @keywords internal
#' @noRd
verify_suggested <- function(pkg) {
  # verify that suggested packages are available
  check_pkg <- sapply(pkg, requireNamespace, quietly = TRUE)
  if (!all(check_pkg)) {
    msg <- paste0(
      "The use of this function requires package", ifelse(sum(!check_pkg)>1, "s", ""),
      paste(names(check_pkg)[!check_pkg], collapse=", "),
      ". Please install."
    )
    stop(msg)
  }
  invisible(NULL)
}

#' @title ldply_base
#' @param .data list.
#' @param .fun fun.
#' @param .progress Show progress bar if 'text'.
#' @keywords internal
#' @noRd
ldply_base <- function(.data, .fun = identity, .progress = "none") {
  n <- length(.data)

  if (.progress == "text") {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  }

  result <- vector("list", n)
  for (i in seq_along(.data)) {
    if (.progress == "text") {
      utils::setTxtProgressBar(pb, i)
    }
    result[[i]] <- .fun(.data[[i]])
  }

  if (.progress == "text") {
    close(pb)
  }

  df <- do.call(rbind, result)
  df <- data.frame(df, row.names = NULL, check.names = FALSE)
  return(df)
}
