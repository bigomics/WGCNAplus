#' Plot consensus module overlap heatmap
#' @param net1 First network object.
#' @param net2 Second network object.
#' @param setLabels Labels for the two sets.
#' @param lab.line Margin line for axis labels.
#' @param plotDendro If TRUE, plot dendrograms above.
#' @param setpar If TRUE, set par layout.
#' @return NULL (invisible). Generates a plot.
#' @export
plotConsensusOverlapHeatmap <- function(net1, net2,
                                              setLabels = NULL,
                                              lab.line = c(8, 8),
                                              plotDendro = FALSE,
                                              setpar = TRUE) {

  if (is.null(setLabels)) setLabels <- c("Set1", "Set2")
  if (length(setLabels) == 1) setLabels <- paste0(setLabels, 1:2)

  if (plotDendro) {

    layout.matrix <- matrix(c(1, 2, 5, 3, 4, 5), nrow = 3, ncol = 2)
    layout(mat = layout.matrix, heights = c(0.8, 0.2, 2.5), widths = c(1, 1))

    WGCNA::plotDendroAndColors(
      dendro = net1$dendrograms[[1]],
      colors = net1$colors,
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = FALSE,
      guideHang = 0.05,
      setLayout = FALSE,
      main = setLabels[1]
    )

    WGCNA::plotDendroAndColors(
      dendro = net2$dendrograms[[1]],
      colors = net2$colors,
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = FALSE,
      guideHang = 0.05,
      setLayout = FALSE,
      main = setLabels[2]
    )

  }


  firstColors <- labels2colors(net1$colors)
  secondColors <- labels2colors(net2$colors)
  overlap <- overlapTable(firstColors, secondColors)

  firstModTotals <- rowSums(overlap$countTable)
  secondModTotals <- colSums(overlap$countTable)
  firstModules <- rownames(overlap$countTable)
  secondModules <- colnames(overlap$countTable)

  # Truncate p values smaller than 10^{-50} to 10^{-50}
  pTable <- -log10(overlap$pTable)
  pTable[is.infinite(pTable)] <- 1.3 * max(pTable[is.finite(pTable)])
  pTable[pTable > 50] <- 50

  if (setpar) {
    par(mfrow = c(1, 1))
    par(cex = 1.0)
    par(mar = c(10, 12.4, 2.7, 1) + 0.3)
  }

  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  ss <- paste("Correspondence of", setLabels[1], "and ", setLabels[2], "modules")
  WGCNA::labeledHeatmap(
    Matrix = t(pTable),
    xLabels = paste(" ", firstModules),
    yLabels = paste(" ", secondModules),
    colorLabels = TRUE,
    xSymbols = paste0(firstModules, " (", firstModTotals, ")"),
    ySymbols = paste0(secondModules, " (", secondModTotals, ")"),
    textMatrix = t(overlap$countTable),
    colors = WGCNA::blueWhiteRed(100)[50:100],
    main = ss,
    cex.text = 1.0,
    cex.lab = 1.0,
    setStdMargins = FALSE
  )
  mtext(toupper(setLabels[1]), side = 1, line = lab.line[1], cex = 1.1)
  mtext(toupper(setLabels[2]), side = 2, line = lab.line[2], cex = 1.1)
}


#' Plot module preservation summary statistics
#' @param pres Preservation results object.
#' @param setpar If TRUE, set par layout.
#' @return NULL (invisible). Generates a plot.
#' @export
plotPreservationSummaries <- function(pres, setpar = TRUE) {

  # Create a simple bar plot of Zsummary:
  Z <- pres$Zsummary
  ntest <- ncol(Z)

  if (setpar) par(mfrow = c(3, ntest), mar = c(5, 5, 4, 1))

  xylist <- list(
    c("moduleSize", "Zsummary.pres"),
    c("moduleSize", "medianRank.pres"),
    c("Zsummary.pres", "medianRank.pres")
  )

  for (xy in xylist) {
    for (k in colnames(Z)) {
      X <- data.frame(
        color = substring(names(pres$moduleSize), 3, 99),
        moduleSize = pres$moduleSize,
        Zsummary.pres = pres$Zsummary[, k],
        medianRank.pres = pres$medianRank[, k]
      )
      xvar <- xy[1]
      yvar <- xy[2]
      ylim <- c(0, max(X[, yvar]))
      if (yvar == "medianRank.pres") ylim <- rev(ylim)
      plot(
        X[, xvar],
        X[, yvar],
        pch = 21,
        cex = 2,
        bg = X$color,
        ylim = ylim,
        xlab = xvar,
        ylab = yvar
      )
      title(yvar, cex.main = 1.4, line = 2.2)
      sub <- paste(k, "vs.", "reference")
      title(sub, cex.main = 1, line = 0.9)
      if (yvar == "Zsummary.pres") abline(h = c(2, 10), lty = 2)
    }
  }
}

#' Plot preservation and module-trait heatmaps
#' @param pres Preservation results object.
#' @param subplots Subplots to show: "zsummary", "consmt", "wt.consmt".
#' @param order.by Order modules by "name", "zsummary", or "clust".
#' @param setpar If TRUE, set par layout.
#' @param rm.na If TRUE, remove all-NA columns.
#' @return NULL (invisible). Generates a plot.
#' @export
plotPreservationModuleTraits <- function(pres,
                                         subplots = c("zsummary", "consmt", "wt.consmt"),
                                         order.by = "name",
                                         setpar = TRUE,
                                         rm.na = FALSE) {

  if (all(is.numeric(subplots))) {
    subplots <- c("zsummary", "consmt", "wt.consmt")[subplots]
  }

  if (setpar) par(mfrow = c(2, 2), mar = c(14, 12, 4, 2))

  ## compute consensus
  Zsummary <- pres$Zsummary
  cR <- pres$modTraits
  ydim <- sapply(pres$layers, function(w) nrow(w$datTraits))
  consZ <- computeConsensusMatrix(cR, ydim, psig = 1, consfun = "gmean")

  ## match
  ii <- intersect(rownames(Zsummary), rownames(consZ))
  Zsummary <- Zsummary[ii, , drop = FALSE]
  consZ <- consZ[ii, , drop = FALSE]

  ## order
  if (order.by == "name") {
    ii <- order(rownames(Zsummary))
    Zsummary <- Zsummary[ii, , drop = FALSE]
    consZ <- consZ[ii, , drop = FALSE]
  }
  if (order.by == "zsummary") {
    ii <- order(-rowMeans(Zsummary**2))
    Zsummary <- Zsummary[ii, , drop = FALSE]
    consZ <- consZ[ii, , drop = FALSE]
  }
  if (order.by == "clust") {
    consZ1 <- consZ
    consZ1[is.na(consZ1)] <- 0
    ii <- hclust(dist(consZ1))$order
    jj <- hclust(dist(t(consZ1)))$order
    Zsummary <- Zsummary[ii, , drop = FALSE]
    consZ <- consZ[ii, jj, drop = FALSE]
  }

  ## Zsummary heatmap
  if ("zsummary" %in% subplots) {
    WGCNA::labeledHeatmap(
      Matrix = Zsummary,
      xLabels = colnames(Zsummary),
      yLabels = rownames(Zsummary),
      ySymbols = rownames(Zsummary),
      colors = tail(WGCNA::blueWhiteRed(100), 50),
      colorLabels = TRUE,
      setStdMargins = FALSE
    )
    title("Module preservation (Zsummary)", line = 1.2, cex.main = 1.2)
  }

  ## Consensus Module-Trait
  validcol <- function(R) {
    which(colMeans(is.na(R)) < 1 &
      matrixStats::colSds(R, na.rm = TRUE) > 0.01)
  }

  if ("consmt" %in% subplots) {
    clim <- max(abs(consZ), na.rm = TRUE)
    cval <- seq(-clim, clim, length.out = 201)
    ii <- which(cval >= min(consZ, na.rm = TRUE) & cval <= max(consZ, na.rm = TRUE))
    col2 <- WGCNA::blueWhiteRed(201)[ii]
    jj <- 1:ncol(consZ)
    if (rm.na) jj <- validcol(consZ)
    WGCNA::labeledHeatmap(
      Matrix = consZ[, jj, drop = FALSE],
      xLabels = colnames(consZ)[jj],
      yLabels = rownames(consZ),
      ySymbols = rownames(consZ),
      colors = col2,
      colorLabels = TRUE,
      setStdMargins = FALSE
    )
    title("Consensus Module-Traits", line = 1.2, cex.main = 1.2)
  }

  ## Preservation-weighted Consensus Module-Trait
  if ("wt.consmt" %in% subplots) {
    wz <- rowMeans(Zsummary**2, na.rm = TRUE)
    wz <- wz / max(wz)
    consW <- consZ * wz[rownames(consZ)]

    clim <- max(abs(consW), na.rm = TRUE)
    cval <- seq(-clim, clim, length.out = 201)
    ii <- which(cval >= min(consW, na.rm = TRUE) & cval <= max(consW, na.rm = TRUE))
    col2 <- WGCNA::blueWhiteRed(201)[ii]

    jj <- 1:ncol(consW)
    if (rm.na) jj <- validcol(consW)
    WGCNA::labeledHeatmap(
      Matrix = consW[, jj, drop = FALSE],
      xLabels = colnames(consW)[jj],
      yLabels = rownames(consW),
      ySymbols = rownames(consW),
      colors = col2,
      colorLabels = TRUE,
      setStdMargins = FALSE
    )
    title("Preservation-weighted Consensus\nModule-Traits", line = 1, cex.main = 1.2)
  }

}
