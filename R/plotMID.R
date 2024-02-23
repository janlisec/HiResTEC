#' @title plotMID.
#'
#' @description
#' \code{plotMID} will plot a Mass Isotopomer Distribution (MID).
#'
#' @details
#' Not yet.
#'
#' @param mid Matrix of corrected mass isotopomer distributions. Samples in columns, MID values in rows.
#' @param gr Groups, a factor vector of length ncol(mid).
#' @param name Name of metabolite.
#' @param contr Contrasts. Not yet clear if useful.
#' @param stackedbars Alternative plotting layout using stacked bar plot.
#' @param subplot_ylim Calculate ylim individually per subplot if 0, show full range in all subplots if 100 and limit to the minimal specified number otherwise.
#' @param ... Further arguments to 'boxplot' or 'barplot' (depending on 'stackedbars').
#'
#' @return
#' A plot.
#'
#' @importFrom graphics barplot
#' @importFrom grDevices grey
#'
#' @export
#'
#' @examples
#' mid <- matrix(c(seq(0, 0.3, 0.1), seq(1, 0.7, -0.1)), byrow = TRUE, nrow = 2)
#' gr <- gl(2, 2, labels = letters[1:2])
#' plotMID(mid = mid, gr = gr, name = "Metabolite X")
#' plotMID(mid = mid, gr = gr, stackedbars = TRUE, las = 1, xlab = "MID", legend.text = c("x", "y"))
#' lt <- paste0("M", 0:1)
#' rownames(mid) <- lt
#' plotMID(mid = mid, gr = gr, stackedbars = TRUE, las = 1, xlab = "MID", legend.text = lt)
#' plotMID(mid = mid[, 2, drop = FALSE], stackedbars = TRUE, col = c(3, 4))
#' colnames(mid) <- paste0("S", 1:4)
#' gr2 <- gl(n = 1, k = 1, labels = "bla")
#' plotMID(mid = mid[, 2, drop = FALSE], gr = gr2, stackedbars = TRUE, name = NULL)
#' plotMID(mid = mid, gr = factor(colnames(mid)), stackedbars = TRUE, name = NULL)
plotMID <- function(mid = NULL, gr = NULL, name = "unknown", contr = NULL, stackedbars = FALSE, subplot_ylim = 100, ...) {
  stopifnot(is.matrix(mid) && is.numeric(mid) && nrow(mid) >= 2 && ncol(mid) >= 1)
  argg <- c(as.list(environment()), list(...))
  opar <- par(no.readonly = TRUE)
  if (is.null(gr)) gr <- gl(n = 1, k = ncol(mid), labels = "")
  if (stackedbars) {
    # get group medians
    tmp <- plyr::adply(unname(mid), 1, function(x) {
      sapply(split(x, gr), median, na.rm=TRUE)
    }, .id = NULL)
    rownames(tmp) <- rownames(mid)
    colnames(tmp) <- levels(gr)
    # readjust to sum=100
    tmp <- apply(tmp, 2, function(x) {
      100 * x / sum(x)
    })
    # readjust default colors
    col <- NULL
    if (!"col" %in% names(argg)) {
      # default colors for APCI fragments
      default_color_fragments <- unlist(list("M+H" = grey(0), "M+" = grey(0.6), "M-H" = grey(0.3), "M+H2O-CH4" = 2))
      if (!is.null(rownames(tmp)) && all(rownames(tmp) %in% names(default_color_fragments))) {
        col <- default_color_fragments[rownames(tmp)]
      }
      default_color_mid <- 1:7
      names(default_color_mid) <- paste0("M", 0:6)
      if (!is.null(rownames(tmp)) && all(rownames(tmp) %in% names(default_color_mid))) {
        col <- default_color_mid[rownames(tmp)]
      }
    } else {
      col <- argg$col
    }
    # set a number of parameters to default values and substitute from argument list if provided by user
    xlab <- ""
    if ("xlab" %in% names(argg) && !is.null(argg$xlab)) {
      xlab <- argg$xlab
    }
    ylab <- ifelse(is.character(name), name, "")
    if ("ylab" %in% names(argg) && !is.null(argg$ylab)) {
      ylab <- argg$ylab
    }
    legend.text <- NULL
    if ("legend.text" %in% names(argg)) {
      legend.text <- argg$legend.text
    }
    args.legend <- NULL
    if (!"args.legend" %in% names(argg)) {
      args.legend <- argg$args.legend
    }
    las <- 1
    if ("las" %in% names(argg)) {
      las <- argg$las
    }
    horiz <- TRUE
    if ("horiz" %in% names(argg)) {
      horiz <- argg$horiz
    }
    axes <- TRUE
    if ("axes" %in% names(argg)) {
      axes <- argg$axes
    }
    mar_b <- ifelse(nchar(xlab) >= 1, 5, 3)
    mar_l <- ifelse(nchar(ylab) >= 1, 4, 2)
    if ("mar" %in% names(argg)) {
      par(mar = argg$mar)
    } else {
      par(mar = c(mar_b, mar_l, 1, 0) + 0.5)
    }
    graphics::barplot(tmp, ylab = ylab, col = col, xlab = xlab, legend.text = legend.text, las = las, args.legend = args.legend, axes = axes, horiz = horiz)
    par(mar = opar$mar)
  } else {
    tmp <- apply(mid, 1, function(x) {
      split(x, gr)
    })
    par(mfrow = c(1, length(tmp)))
    par(mar = c(ifelse("xlab" %in% names(argg), 5, 3), 4, 1, 0) + 0.5)
    ylim <- range(mid, na.rm = T)
    for (k in 1:length(tmp)) {
      if (subplot_ylim < 100) {
        ylim <- range(tmp[[k]], na.rm = T)
        while (diff(ylim) < subplot_ylim) {
          ylim <- round(ylim + c(-1, 1))
          if (ylim[1] < 0) ylim[1] <- 0
          if (ylim[2] > 100) ylim[2] <- 100
        }
      }
      graphics::boxplot(tmp[[k]], main = "", ylab = "", ylim = ylim, ...)
      graphics::mtext(text = paste0("M", k - 1), side = 3, adj = 1)
      if (k == 1) {
        graphics::title(ylab = name, cex.lab = 2)
      }
      if (!is.null(contr) && all(contr %in% 1:length(tmp[[1]]))) {
        if (length(contr) == 1) {
          for (l in (1:length(tmp[[1]]))[-contr]) {
            p <- try(stats::t.test(tmp[[k]][[l]], tmp[[k]][[contr]])$p.value)
            if (is.numeric(p)) {
              pos <- ifelse(mean(tmp[[k]][[l]]) < 50, 3, 1)
              graphics::text(x = l + 0.5, y = mean(tmp[[k]][[l]]), labels = formatC(p, format = "e", digits = 1), pos = pos, col = ifelse(p < 0.01, 2, grDevices::grey(0.9)))
            }
          }
        } else {
          for (l in contr) {
            p <- try(stats::t.test(tmp[[k]][[l]], tmp[[k]][[l + 1]])$p.value)
            if (is.numeric(p)) {
              pos <- ifelse(mean(tmp[[k]][[l]]) < 50, 3, 1)
              graphics::text(x = l + 0.5, y = mean(tmp[[k]][[l]]), labels = formatC(p, format = "e", digits = 1), pos = pos, col = ifelse(p < 0.01, 2, grDevices::grey(0.9)))
            }
          }
        }
      }
    }
    par(mar = opar$mar)
    par(mfrow = opar$mfrow)
  }
  invisible(NULL)
}
