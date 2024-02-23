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
