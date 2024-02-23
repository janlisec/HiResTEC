#' @title ExportHeuristicsTable.
#'
#' @description \code{ExportHeuristicsTable}.
#'
#' @param res_list1 res_list1.
#' @param res_list2 res_list2.
#' @param out1 out1.
#' @param out2 out2.
#' @param p_thr p_thr.
#' @param xls_name xls_name.
#'
#' @importFrom openxlsx write.xlsx
#'
#' @keywords internal
#' @noRd
ExportHeuristicsTable <- function(res_list1 = NULL, res_list2 = NULL, out1 = NULL, out2 = NULL, p_thr = 0.01, xls_name = "Tables.xlsx") {
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
    tmp1 <- plyr::ldply(sort(table(unlist(sapply(res_list1, function(x) {
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
    tmp2 <- plyr::ldply(sort(table(unlist(sapply(res_list2, function(x) {
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
