#' GSEA plot
#'
#' @param rnk ranking matrix
#' @param gset a list of gene set
#'
#' @noRd
#' @import graphics
#' @import stats
#' @importFrom gplots bluered
#'
gsea.enplot <- function(rnk, gset, names = NULL, main = NULL,
                        decreasing = TRUE, cex = 1, cex.main = 0.9, len.main = 40,
                        lab.line = c(0.8, 2), cex.lab = 0.8, main.line = 0.3,
                        xlab = "Rank in ordered dataset", res = 1200,
                        ylab = "Rank metric") {
  if (0) {
    names <- NULL
    main <- NULL
    decreasing <- TRUE
    cex.main <- 0.9
    len.main <- 40
    cex <- 1
    xlab <- "Rank in ordered dataset"
    res <- 1200
    ylab <- "Ranked list metric"
  }
  if (!is.null(names)) names(rnk) <- names
  rnk <- rnk[!is.na(rnk)]
  rnk <- rnk[order(rnk + 1e-8 * rnorm(length(rnk)), decreasing = decreasing)]

  ## ranked list metric
  ii <- (1:length(rnk))
  if (length(ii) > res) ii <- ii[seq(1, length(ii), length(ii) / res)]
  qq <- stats::quantile(rnk[ii], probs = c(0.01, 0.99), na.rm = TRUE)
  y1 <- qq[2]
  y0 <- qq[1]
  dy <- 0.25 * (y1 - y0)
  ## par(mgp=c(2.3,0.6,0))
  graphics::plot(ii, rnk[ii],
    type = "h", col = "grey", ylim = c(y0 - dy, y1),
    xlab = NA, ylab = NA, xaxt = "n"
  )
  graphics::mtext(xlab, 1, line = lab.line[1], cex = cex.lab)
  graphics::mtext(ylab, 2, line = lab.line[2], cex = cex.lab)
  graphics::abline(h = 0, lty = 2, lwd = 0.5)

  ## gene set barcode
  jj <- match(gset, names(rnk))

  w1 <- ifelse(length(jj) < 100, 0.6, 0.3)
  w1 <- ifelse(length(jj) < 50, 1, w1)
  arrows(jj, (y0 - dy), jj, y0, col = "grey10", lwd = w1 * cex, length = 0)

  ## red/blue bar at bottom
  kk <- c(seq(1, length(rnk) * 0.99, floor(length(rnk) / 20)), length(rnk))
  length(kk)
  i <- 1
  for (i in 1:(length(kk) - 1)) {
    r <- mean(rnk[kk[c(i, i + 1)]])
    r1 <- (r / max(abs(rnk), na.rm = TRUE))
    r1 <- abs(r1)**0.5 * sign(r1)
    irnk <- floor(31 * (1 + r1) / 2)
    cc <- gplots::bluered(32)[1 + irnk]
    graphics::rect(kk[i], y0 - 1.05 * dy, kk[i + 1], y0 - 0.65 * dy, col = cc, border = NA)
  }

  ## weighted cumulative random walk
  x0 <- 1 * (names(rnk) %in% gset)
  x0 <- x0 * abs(rnk)
  n0 <- sum(!(names(rnk) %in% gset))
  n1 <- sum(names(rnk) %in% gset)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))

  rnk.trace <- (r1 - r0)
  rnk.trace <- rnk.trace / max(abs(rnk.trace)) * 0.9
  if (max(rnk.trace) >= abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y1)
  if (max(rnk.trace) < abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y0)
  if (!decreasing) rnk.trace <- -1 * rnk.trace
  graphics::lines(ii, rnk.trace[ii], col = "green", type = "l", lwd = 2.4)

  if (is.null(main)) main <- "Enrichment plot"
  tt.main <- as.character(main)

  breakstring <- function(s, n, nmax = 999, force = FALSE, brk = "\n") {
    if (is.na(s)) {
      return(NA)
    }
    s <- substring(as.character(s), 1, nmax)
    if (nchar(s) < n) {
      return(s)
    }
    b <- substring(s, 1, n)
    n1 <- n + 1
    for (i in 1:10) {
      if (n1 > nchar(s)) break
      b1 <- substring(s, n1, n1 + n)
      b <- paste0(b, brk, b1)
      n1 <- n1 + n + 1
    }
    return(b)
  }

  if (nchar(tt.main) > len.main) {
    tt.main < breakstring(tt.main, len.main) ## pgx-funtions.R
  }
  graphics::title(main = tt.main, cex.main = cex.main, line = main.line)
}

