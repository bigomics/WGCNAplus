#' Plot TOM heatmap with module colors
#'
#' @param wgcna A WGCNA result object.
#' @param justdata If TRUE, return dissimilarity matrix only.
#' @param block Block number to plot.
#' @param legend If TRUE, add color legend.
#' @param downsample Downsample to this many genes.
#' @return NULL (invisible). Generates a plot.
#' @export
plotTOM <- function(wgcna, justdata = FALSE, block = NULL,
                          legend = TRUE, downsample = NULL) {
  datExpr <- wgcna$datExpr
  wTOM <- NULL
  if (!is.null(wgcna$TOM)) {
    wTOM <- wgcna$TOM
  }
  ## if SV of TOM is stored, reconstruct TOM
  if (is.null(wTOM) && !is.null(wgcna$svTOM)) {
    wTOM <- tcrossprod(wgcna$svTOM)
  }
  if (is.null(wTOM)) {
    message("[plotTOM] ERROR. no TOM matrix")
    return(NULL)
  }
  dissTOM <- 1 - wTOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)

  ## clustering wgcnaults
  moduleColors <- NULL
  if (is.null(block) && "merged_dendro" %in% names(wgcna$net)) {
    geneTree <- wgcna$net$merged_dendro
    gg <- geneTree$labels
    if (length(gg) > 0) {
      dissTOM <- dissTOM[gg, gg]
      moduleColors <- labels2colors(wgcna$net$colors[gg])
    }
  }

  if (is.null(moduleColors)) {
    if (is.null(block)) block <- 1
    geneTree <- wgcna$net$dendrograms[[block]]
    ii <- which(wgcna$net$blocks == block & wgcna$net$goodGenes == TRUE)
    gg <- names(wgcna$net$color)[ii]
    dissTOM <- dissTOM[gg, gg]
    moduleColors <- labels2colors(wgcna$net$colors[gg])
  }

  if (justdata) {
    return(dissTOM)
  }

  if (!is.null(downsample) && ncol(dissTOM) > downsample) {
    ii <- seq(1, ncol(dissTOM), length.out = downsample)
    dissTOM <- dissTOM[ii, ii]
    moduleColors <- moduleColors[ii]
    geneTree <- fastcluster::hclust(as.dist(dissTOM),
      method = "average"
    )
  }

  if (legend) {
    par(oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0))
    plotly::layout(
      matrix(c(
        0, 0, 5, 0,
        0, 0, 2, 0,
        4, 1, 3, 6
      ), nr = 3, byrow = TRUE),
      widths = c(2.3, 0.5, 10, 3),
      heights = c(2.3, 0.5, 10)
    )
    WGCNA::TOMplot(
      dissTOM^7,
      geneTree,
      moduleColors,
      setLayout = FALSE,
      main = NULL
    )

    ## add color legend
    frame()
    legend(
      -0.1, 1,
      fill = wgcna$me.colors,
      legend = names(wgcna$me.colors),
      cex = 1.2, bty = "n", x.intersp = 0.5
    )
  } else {
    WGCNA::TOMplot(
      dissTOM^7,
      geneTree,
      moduleColors,
      setLayout = TRUE,
      main = NULL
    )
  }
}

#' Plot module-trait correlation heatmap
#'
#' @param wgcna A WGCNA result object.
#' @param setpar If TRUE, set par margins.
#' @param cluster If TRUE, cluster rows and columns.
#' @param multi If TRUE, treat as multi-dataset.
#' @param main Plot title.
#' @param justdata If TRUE, return data only.
#' @param transpose If TRUE, transpose the matrix.
#' @param colorlabel If TRUE, use color labels.
#' @param show Show "both", "traits", or "contrasts".
#' @param nmax Max number of modules to show.
#' @param tmax Max number of traits to show.
#' @param text If TRUE, display correlation values.
#' @param pstar If TRUE, show significance stars.
#' @return NULL (invisible). Generates a plot.
#' @export
plotModuleTraitHeatmap <- function(wgcna, setpar = TRUE, cluster = FALSE,
                                   multi = FALSE, main = NULL, justdata = FALSE,
                                   transpose = FALSE, colorlabel = TRUE,
                                   show = c("both","traits","contrasts")[1],
                                   nmax = -1, tmax = -1,
                                   text = TRUE, pstar = TRUE) {
  
  if(!multi) layers <- list(gx=wgcna)
  if(multi && !is.null(wgcna$layers)) layers <- wgcna$layers
  if(multi && is.null(wgcna$layers)) layers <- wgcna    

  MEs <- lapply(layers, function(w) as.matrix(w$net$MEs))
  MEs <- mergeME(MEs)

  Y <- layers[[1]]$datTraits
  sel <- 1:ncol(Y)
  if(show=="traits") sel <- grep("_vs_",colnames(Y),invert=TRUE)
  if(show=="contrasts") sel <- grep("_vs_",colnames(Y))
  Y <- Y[,sel,drop=FALSE]

  moduleTraitCor <- cor(MEs, Y, use = "pairwise.complete")

  #nSamples <- nrow(layers[[1]]$datExpr)
  nSamples <- t(!is.na(MEs)) %*% (!is.na(Y))

  if(nmax > 0) {
    sel <- head(order(-apply(abs(moduleTraitCor), 1, max, na.rm=TRUE)),nmax)
    moduleTraitCor <- moduleTraitCor[sel,,drop=FALSE]
    nSamples <- nSamples[sel,,drop=FALSE]
  }
  if(tmax > 0) {
    sel <- head(order(-apply(abs(moduleTraitCor), 2, max, na.rm=TRUE)),tmax)
    moduleTraitCor <- moduleTraitCor[,sel,drop=FALSE]
    nSamples <- nSamples[,sel,drop=FALSE]
  }

  if (transpose) {
    moduleTraitCor <- t(moduleTraitCor)
    nSamples <- t(nSamples)
  }

  plotLabeledCorrelationHeatmap(
    R = moduleTraitCor,
    nSamples = nSamples,
    setpar = setpar,
    cluster = cluster,
    text = text,
    main = main,
    justdata = justdata,
    colorlabel = colorlabel,
    pstar = pstar
  )

}

#' Plot labeled correlation heatmap with p-values
#'
#' @param R Correlation matrix.
#' @param nSamples Number of samples or sample matrix.
#' @param cluster If TRUE, cluster rows and columns.
#' @param text If TRUE, display correlation values.
#' @param main Plot title.
#' @param justdata If TRUE, return data only.
#' @param colorlabel If TRUE, use color labels.
#' @param pstar If TRUE, show significance stars.
#' @param zlim Color scale limits.
#' @param colorpal Color palette function.
#' @param cex.text Text size for cell labels.
#' @param cex.lab Label text size.
#' @param setpar If TRUE, set par margins.
#' @param is.dist If TRUE, treat R as distance matrix.
#' @return NULL (invisible). Generates a plot.
#' @export
plotLabeledCorrelationHeatmap <- function(R, nSamples,
                                                cluster = TRUE, text = TRUE,
                                                main = NULL, justdata = FALSE,
                                                colorlabel = TRUE, pstar = TRUE,
                                                zlim = NULL, colorpal = NULL,
                                                cex.text = 0.7, cex.lab = NULL,
                                                setpar = TRUE, is.dist = FALSE) {
  ## Define numbers of genes and samples
  if (cluster && nrow(R) > 1 && ncol(R) > 1) {
    R0 <- R
    R0[is.na(R0)] <- 0
    is.sym <- nrow(R) == ncol(R) && all(rownames(R) == colnames(R))
    if (is.dist) {
      ii <- hclust(as.dist(abs(R0)))$order
      jj <- ii
    } else if (is.sym) {
      ii <- hclust(dist(R0), method = "average")$order
      jj <- ii
    } else {
      ii <- hclust(dist(R0), method = "average")$order
      jj <- hclust(dist(t(R0)), method = "average")$order
    }
    R <- R[ii, jj]
  }

  R0 <- pmax(pmin(R, 1, na.rm=TRUE), -1, na.rm=TRUE)
  ii <- which(nSamples < 3)
  nSamples <- pmax(nSamples, 3)
  Pvalue <- WGCNA::corPvalueStudent(R0, nSamples)
  if(is.matrix(nSamples) && length(ii)>0) {
    Pvalue[ii] <- NA
  }

  if (justdata) {
    return(R)
  }

  ## Will display correlations and their p-values
  if (pstar) {
    textPv <- cut(Pvalue,
      breaks = c(-1, 0.001, 0.01, 0.05, 99),
      # labels = c("★★★", "★★", "★", "")
      # labels = c("***", "**", "*", "")
      labels = c("+++", "++", "+", "")
    )
  } else {
    textPv <- paste0("(", signif(Pvalue, 1), ")")
  }

  textMatrix <- NULL
  if (text && pstar) textMatrix <- paste0(signif(R, 2), "\n", textPv)
  if (text && !pstar) textMatrix <- paste0(signif(R, 2))
  if (!text && pstar) textMatrix <- textPv
  if (!text && !pstar) textMatrix <- NULL

  if (!is.null(textMatrix)) {
    textMatrix[which(is.na(R))] <- ""
    dim(textMatrix) <- dim(R)
  }

  if (!colorlabel) {
    colnames(R) <- paste0(" ", colnames(R))
    rownames(R) <- paste0(" ", rownames(R))
  }

  if (setpar) par(mar = c(8, 8, 3, 3))
  if (is.null(main)) main <- "Correlation heatmap"

  ## set colorscale. make sure 0 is white if non-symmetric
  col1 <- "grey90"
  if (is.null(zlim)) {
    zlim <- c(min(R, na.rm = TRUE), max(R, na.rm = TRUE))
  }
  rlim <- max(abs(zlim), na.rm = TRUE)
  if (is.null(colorpal)) colorpal <- WGCNA::blueWhiteRed
  if (rlim > 0) {
    rval <- seq(-rlim, rlim, length.out = 201)
    ii <- which(rval >= zlim[1] & rval <= zlim[2])
    col1 <- colorpal(201)[ii]
  }

  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(
    Matrix = R,
    # xLabels = paste0(1:ncol(R),":",colnames(R)),
    xLabels = colnames(R),
    # xLabels = paste0(" ",colnames(R)),
    yLabels = rownames(R),
    xSymbols = colnames(R),
    ySymbols = rownames(R),
    colorLabels = TRUE,
    colors = col1,
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = cex.text,
    cex.lab = cex.lab,
    zlim = zlim,
    main = main
  )
}

#' Plot module expression or correlation heatmap
#'
#' @param wgcna A WGCNA result object.
#' @param module Module name to plot.
#' @param genes Gene names to include.
#' @param rgamma Gamma exponent for correlation.
#' @param min.rho Minimum correlation threshold.
#' @param cex Text size for labels.
#' @param nmax Max number of genes to show.
#' @param cluster If TRUE, cluster genes.
#' @param type Type: "expression" or "correlation".
#' @param heatmap.mar Heatmap margins as c(bottom, left).
#' @param main Plot title.
#' @return NULL (invisible). Generates a plot.
#' @export
plotModuleHeatmap <- function(wgcna,
                                    module,
                                    genes = NULL,
                                    rgamma = 4,
                                    min.rho = 0,
                                    cex = 0.8,
                                    nmax = -1,
                                    cluster = TRUE,
                                    type = c("expression", "correlation")[1],
                                    heatmap.mar = c(7, 7),
                                    main = NULL) {
  if (!is.null(module) && is.null(genes)) {
    genes <- wgcna$me.genes[[module]]
  }
  if (is.null(genes) || length(genes) == 0) {
    stop("must specify genes or module")
  }

  if (nmax > 0) {
    sdx <- matrixStats::colSds(wgcna$datExpr[, genes])
    genes <- head(genes[order(-sdx)], nmax)
  }

  if (type == "expression") {
    X <- t(wgcna$datExpr[, genes])
    annot <- wgcna$datTraits
    gx.heatmap(X,
      nmax = nmax,
      ## col.annot = annot,
      key = FALSE, keysize = 0.5, mar = heatmap.mar
    )
  }

  if (type == "correlation") {
    R <- cor(wgcna$datExpr[, genes])
    R <- sign(R) * abs(R)**rgamma
    if (cluster) {
      ii <- hclust(as.dist(1 - R), method = "average")$order
      R <- R[ii, ii]
    }
    R[abs(R) < min.rho] <- NA
    image(R)
    mtext(rownames(R),
      side = 4, at = seq(0, 1, 1 / (nrow(R) - 1)),
      las = 2, adj = 0, cex = cex, line = 0.5
    )
  }

  if (is.null(main) && !is.null(module)) main <- module
  if (is.null(main)) main <- "Module Heatmap"
  if (!is.null(main) && main != "") {
    title(main, line = 2, cex.main = 1.3)
  }
}
