#' @title Identify and evaluate mz pairs from a peak list.
#'
#' @description \code{EvaluatePairsFromXCMSSet} will analyze an `xcmsSet` result
#'     or a generic peak list from a mass spectrometry experiment for mass
#'     pairs (mz1, mz2) with changes due to any tracer incorporation.
#'
#' @details Using 'APCI' as method assumes that (i) you analyze TMS-derivatized
#'     compounds and (ii) your MS resolution does not allow to separate Si and C
#'     isotopes but reports an intermediate mass as m/z. In this case you will
#'     find carbon isotopes below there expected masses, i.e. M+1 would be
#'     1.001 mDa apart from M+0 instead of 1.003.
#'     The effect is increased with isotope number, i.e. M+6 will be ~20 mDa
#'     below the expected value. Hence, selecting method 'APCI' will combine your
#'     selected dmz with a allowed deviation due to mass shifts caused by Si
#'     isotopes. Use 'ESI' if you are not sure if this effect takes place in your
#'     settings.
#'
#' @param xg xcmsSet object with group information. Alternatively, can be a numeric
#'     matrix containing `mz` and `rt` information in the first two columns
#'     followed by peak intensities of all samples in the same order as in
#'     parameters `tp` and `gr`.
#' @param tp Time point information for all samples (obviously required,
#'     internally converted to factor).
#' @param gr Group information for all samples, e.g. different genotypes or
#'     concentrations (optional, factor).
#' @param dmz Allowed mass deviation in Da.
#' @param drt Allowed rt deviation in time units of xcmsSet (usually seconds) to
#'     test for candidates.
#' @param mz_iso Mass defect of the isotope under investigation (use 1.00335 for ^13^C experiments.
#' @param n Number of maximal incorporated carbons to test.
#' @param method Currently APCI or ESI. If APCI, dmz will be modified depending
#'     on n (see details).
#' @param specific_row A single row of the peak list to process.
#' @param testing Stop in function using browser() if specific_row is specified;
#'     can be a isotope number, i.e. 3 will stop at third isotope.
#' @param silent Suppress warnings and console output if TRUE.
#'
#' @return A dataframe with all observable pairs within the provided xg object
#'     (usually an `xcmsSet` peak) list including mean group intensities and P values.
#'
#' @examples
#' # The example measurement data provided with HiResTEC contain a peak at 1026s
#' raw <- HiResTEC::raw
#' sam <- HiResTEC::sam
#' mz <- c(556.26, 561.26, 564.26)
#'
#' # extract the peak intensities for 3 m/z of this peak
#' int <- sapply(raw, function(x) {
#'   tmp <- getMultipleBPC(x = x, mz = mz, mz_dev = 0.04, rt = 1026)
#'   tmp[attr(tmp, "maxBPC"),]
#' })
#' colnames(int) <- sam$ID; rownames(int) <- NULL
#' xg <- data.frame(
#'  "mz" = mz,
#'  "rt" = rep(1026.5, 3),
#'  int
#' )
#'
#' # evaluate this peak list for interesting pairs
#' EvaluatePairsFromXCMSSet(xg=xg, tp=sam$TP, gr=sam$Group, silent=TRUE, n=8)
#'
#' @importFrom plyr ldply
#'
#' @export
#'
EvaluatePairsFromXCMSSet <- function(xg = NULL, tp = NULL, gr = NULL, drt = 1, dmz = 0.025, mz_iso = 1.00335, n = 6, method = c("APCI", "ESI")[1], specific_row = NULL, testing = FALSE, silent = FALSE) {
  # check if grouped xcmsSet provided and extract data from object
  # [20190618 JL needed to remove class information regarding xcms objects for CRAN]
  # stopifnot(class(xg)=="xcmsSet")


  if (inherits(xg, "data.frame")) {
    gv <- as.matrix(xg[,-c(1:2)])
    xg <- as.matrix(xg[,c(1:2)])
    colnames(xg) <- c("mzmed", "rtmed")
  } else {
    # [20190514 JL xcms needs to be put to suggest due to mzR on request by CRAN --> replaced xcms::groupval by quick and dirty own function]
    gv <- groupval(xg = xg)

    # if (requireNamespace("xcms", quietly = TRUE)) {
    #   gv <- xcms::groupval(xg, value = "maxo", method="maxint")
    # } else {
    #   warning("Package xcms not available.")
    # }

    xg <- xg@groups
  }


  utils::data("mz_shift_corrector", package = "HiResTEC")
  mz_shift_corrector <- mz_shift_corrector[[method]][1:(n + 1)]

  # check if time (and group) information provided and prepare
  if (is.numeric(tp)) tp <- factor(tp)
  stopifnot(!is.null(tp))
  stopifnot(length(levels(tp)) >= 2)
  if (is.null(gr)) gr <- gl(1, ncol(gv)) else gr <- factor(gr)
  inter <- interaction(tp, gr, drop = TRUE, sep = "_")
  tp <- as.numeric(as.character(tp))

  # check if processing is limited to specific row
  if (is.null(specific_row)) specific_row <- 1:nrow(xg)
  stopifnot(all(specific_row %in% 1:nrow(xg)))

  # find all pairs
  out <- plyr::ldply(specific_row, function(i) {
    #message(i)
    # filter for coeluting peaks
    rtok <- abs(xg[i, "rtmed"] - xg[, "rtmed"]) < drt # table(rtok)

    # specify potential mz2 candidate masss which may pair with mz1
    mz_cand <- xg[i, "mzmed"] + c(1:n * mz_iso)

    # check each candidate for being present
    plyr::ldply(mz_cand, function(x) {
      # filter for fitting masses within coeluting peaks
      # mzok <- which(abs(x-xg[which(rtok),"mzmed"]) < dmz)
      # browser()
      #message(x)
      ni <- round(x - xg[i, "mzmed"])
      md <- xg[which(rtok), "mzmed"] - x
      # upper boud defined by dmz; lower bound defined by dmz+correction for Si-Isotopes
      mzok <- which(md < dmz & md > -(mz_shift_corrector[ni] + dmz))
      if (length(mzok) >= 1) {
        if (length(mzok) >= 2) {
          if (!silent) warning(paste("\nFound more than one fitting mz2 for line", i, "and mz2", x, " and keep closest by dmz."))
          mzok <- mzok[which.min(abs(x - xg[which(rtok)[mzok], "mzmed"]))]
        }
        j <- which(rtok)[mzok]

        # get peaks and ratios
        pks <- gv[c(i, j), ]
        # ratios <- pks[2,]/pks[1,]
        # 'ratios' are PseudoEnrichments now!!!
        ratios <- apply(pks, 2, CalcEnrichment)

        # evaluate intensity group means
        int_mz1 <- data.frame(lapply(split(pks[1, ], inter), function(x) {
          ifelse(any(is.finite(x)), mean(x, na.rm = T), NA)
        }), check.names = FALSE)
        colnames(int_mz1) <- paste0("I1_", colnames(int_mz1))
        int_mz2 <- data.frame(lapply(split(pks[2, ], inter), function(x) {
          ifelse(any(is.finite(x)), mean(x, na.rm = T), NA)
        }), check.names = FALSE)
        colnames(int_mz2) <- paste0("I2_", colnames(int_mz2))

        # set a limit for peak ratio at t=0 (to avoid wrong candidates due to [M+] vs [M+H] ratios)
        if (any(tp == 0)) {
          test_int_iso <- int_mz1 / int_mz2
          flt <- grep("I1_0", names(test_int_iso))
          test_int_iso <- any(test_int_iso[flt] > 0.7, na.rm = TRUE) | all(is.na(test_int_iso[flt]))
          # old version
          # test_int_iso <- any((int_mz1/int_mz2)>0.7, na.rm=TRUE)
        } else {
          test_int_iso <- TRUE
        }

        # compute peak distance (dRT)
        dRT <- abs(xg[i, "rtmed"] - xg[j, "rtmed"])

        # stop process here if requested by user via option 'testing'
        if (is.numeric(testing) && testing == which(mz_cand == x)) testing <- TRUE
        browser(expr = !is.null(specific_row) && is.logical(testing) && testing)

        # compute ANOVA model for 'ratios ~ time * group'
        # if (test_int_iso & sum(sapply(split(ratios, inter), function(y) { any(is.finite(y)) }))>=0.5*length(levels(inter))) {
        test_all_finite <- sum(sapply(split(ratios, inter), function(y) {
          all(is.finite(y))
        })) >= floor(0.5 * length(levels(inter)))
        if (test_int_iso & test_all_finite) {
          # check if peak ratios change over time using a linear model with interaction
          if (length(levels(gr)) > 1) {
            # if groups are specified use:
            ratios_lm <- try(lm(ratios ~ tp * gr), silent = TRUE)
          } else {
            # if only timepoints are specified use:
            ratios_lm <- try(lm(ratios ~ tp), silent = TRUE)
          }
          if (!attr(ratios_lm, "class") == "lm") {
            if (!silent) cat(paste0("\ni=", i, "  x=", x, "  lm failed -> set P=1 and dR=0"))
            P <- 1
            dR <- 0
          } else {
            alm <- anova(ratios_lm)
            alm_tp <- grep("tp", rownames(alm))
            if (length(alm_tp) >= 1) {
              P <- min(alm[alm_tp, "Pr(>F)"], na.rm = T)
              # [ToDo] think about best option to calculate dR
              # dR <- 100*diff(range(sapply(split(ratios, inter), median,na.rm=T)))
              # dR <- coef(ratios_lm)["tp"]*diff(range(tp))*100
              # browser()
              # taking only the extreme time points avoids partially FalsePositives
              # filt <- tp %in% range(tp)
              # dR <- diff(sapply(split(ratios[filt], tp[filt]), median, na.rm=T))
              dR <- round(median(ratios[tp == max(tp)], na.rm = T) - median(ratios[tp == min(tp)], na.rm = T), 2)
            } else {
              # no time dependent P-value could be found using a linear model
              # a reason could be that mz2 was n.d. at t=0 and mz1 n.d. at t=1
              P <- 1
              dR <- 0
            }
          }
          return(data.frame("peak_idx" = i, "rt" = xg[i, "rtmed"], "mz1" = xg[i, "mzmed"], "mz2" = xg[j, "mzmed"], "P" = P, "dR" = dR, "dRT" = dRT, int_mz1, int_mz2))
        } else {
          return(data.frame("peak_idx" = i, "rt" = xg[i, "rtmed"], "mz1" = xg[i, "mzmed"], "mz2" = xg[j, "mzmed"], "P" = NA, "dR" = NA, "dRT" = dRT, int_mz1, int_mz2))
        }
      } else {
        dummy <- data.frame(as.list(rep(NA, 2 * length(levels(inter)))))
        colnames(dummy) <- paste0(rep(c("I1_", "I2_"), each = length(levels(inter))), levels(inter))
        return(data.frame("peak_idx" = i, "rt" = xg[i, "rtmed"], "mz1" = xg[i, "mzmed"], "mz2" = NA, "P" = NA, "dR" = NA, "dRT" = NA, dummy)[-1, , drop = F])
      }
    })
  }, .progress = ifelse(silent, "none", "text"))
  return(out)
}
