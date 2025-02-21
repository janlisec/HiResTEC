#' @title CandidateBoxplot.
#' @description \code{CandidateBoxplot} will plot enrichment of tracer incorporationn candidates.
#' @details Plot will be annotated with mz/rt information plus ANOVA(lm()) P-values.
#' @param res Candidate evaluation result.
#' @return A list with all important information including deconvoluted spectrum and linear model.
#' @importFrom beeswarm beeswarm
#' @keywords internal
#' @noRd
CandidateBoxplot <- function(res = NULL) {
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
