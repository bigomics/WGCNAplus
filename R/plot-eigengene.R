#' Plot eigengene cluster dendrogram
#' @param wgcna A WGCNA result object.
#' @param ME Eigengene matrix to use directly.
#' @param add_traits Include traits in clustering.
#' @param horiz Plot dendrogram horizontally.
#' @param setMargins Adjust plot margins automatically.
#' @param method Clustering method to use.
#' @param showlabels Show dendrogram labels.
#' @param plot Whether to generate the plot.
#' @param multi Use multi-dataset mode.
#' @param main Plot title string.
#' @return Hierarchical clustering object (invisible).
#' @export
plotEigenGeneClusterDendrogram <- function(wgcna = NULL,
                                           ME = NULL,
                                           add_traits = TRUE,
                                           horiz = FALSE,
                                           setMargins = TRUE,
                                           method = "wgcna",
                                           showlabels = TRUE,
                                           plot = TRUE,
                                           multi = FALSE,
                                           main = NULL) {

  # Matrix with eigengenes and traits
  if (is.null(wgcna) && is.null(ME)) {
    stop("ERROR: wgcna or ME must be given")
  }

  if (is.null(ME)) {
    if (multi) {
      ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
      ME <- mergeME(ME)
      Y <- wgcna[[1]]$datTraits
    } else {
      ME <- wgcna$net$MEs
      Y <- wgcna$datTraits
    }

    if (length(add_traits)==1 && is.logical(add_traits) && add_traits==TRUE) {
      ME <- mergeME(ME, Y)
    } else if (length(add_traits) > 0 && !is.logical(add_traits)) {
      sel <- intersect(add_traits, colnames(Y))
      if(length(sel)) ME <- mergeME(ME, Y[,sel])
    }
  }

  impME <- svdImpute2(as.matrix(ME))
  ME <- WGCNA::orderMEs(impME)

  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!
  if (is.null(main)) main <- "Eigengene Dendrogram"

  hc <- NULL
  if (method == "wgcna") {
    ## plot dendrogram with WGCNA function
    WGCNA::plotEigengeneNetworks(
      ME, main,
      setMargins = setMargins,
      marDendro = c(0, 4, 2, 0),
      plotHeatmaps = FALSE
    )
  } else {
    ## plot dendrogram with hclust function
    if (setMargins && horiz) par(mar = c(4, 4, 4, 8))
    if (setMargins && !horiz) par(mar = c(8, 4, 4, 1))
    hc <- hclust(as.dist(1 - cor(ME)), method = "average")
    if (plot) {
      save.labels <- hc$labels
      if (!showlabels) hc$labels <- rep("", ncol(ME))
      plot(as.dendrogram(hc), horiz = horiz, main = main)
      hc$labels <- save.labels
    }
  }

  invisible(hc)

}


#' Plot adjacency correlation heatmap matrix of eigengenes with or
#' without traits. This can show how traits cluster with the eigengenes.
#' @param wgcna A WGCNA result object.
#' @param add_traits Include traits in heatmap.
#' @param traits Specific traits to include.
#' @param add_me Include module eigengenes.
#' @param marx Margin expansion factor.
#' @param main Plot title string.
#' @param multi Use multi-dataset mode.
#' @param phenotype Phenotype for conditioning.
#' @param colorlabel Color-code row/column labels.
#' @param text Show correlation values as text.
#' @param pstar Show significance stars.
#' @param power Adjacency power exponent.
#' @param setMargins Adjust plot margins automatically.
#' @param mar1 Margins for dendrogram panel.
#' @param mar2 Margins for heatmap panel.
#' @param cex.lab Label character expansion factor.
#' @param cex.text Text character expansion factor.
#' @param plotDendro Plot the dendrogram panel.
#' @param plotHeatmap Plot the heatmap panel.
#' @param dendro.horiz Plot dendrogram horizontally.
#' @param dendro.width Relative width of dendrogram.
#' @param dendro.labels Show dendrogram labels.
#' @param nmax Maximum number of eigengenes shown.
#' @param fixclust Fix cluster ordering.
#' @param mask.intra Mask intra-module correlations.
#' @param justdata Return data without plotting.
#' @return Correlation matrix (invisible), or matrix if justdata is TRUE.
#' @export
plotEigenGeneAdjacencyHeatmap <- function(wgcna,
                                          add_traits = TRUE,
                                          traits = NULL,
                                          add_me = TRUE,
                                          marx = 1, main = NULL,
                                          multi = FALSE,
                                          phenotype = NULL,
                                          colorlabel = TRUE,
                                          text = FALSE,
                                          pstar = TRUE,
                                          power = 1,
                                          setMargins = TRUE,
                                          mar1 = c(5.6, 4.5, 1.8, 0),
                                          mar2 = c(8, 10, 4, 2),
                                          cex.lab = 0.8,
                                          cex.text = 0.7,
                                          plotDendro = TRUE,
                                          plotHeatmap = TRUE,
                                          dendro.horiz = TRUE,
                                          dendro.width = 0.3,
                                          dendro.labels = TRUE,
                                          nmax = -1,
                                          fixclust = FALSE,
                                          mask.intra = FALSE,
                                          justdata = FALSE) {

  if(!multi) wgcna <- list(gx=wgcna)

  # Matrix with eigengenes and traits
  ME <- NULL
  if (add_me) {
    ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
    ME <- mergeME(ME)
  }
  Y <- wgcna[[1]]$datTraits

  if (add_traits) {
    sel <- colnames(Y)
    if (!is.null(traits)) {
      sel <- intersect(traits, sel)
    }
    if (is.null(ME)) {
      ME <- Y[,sel,drop=FALSE]
    } else {
      ME <- mergeME(ME, Y[,sel,drop=FALSE])
    }
  }

  if (!add_traits && !is.null(phenotype)) {
    if (is.null(ME)) {
      ME <- Y[,phenotype,drop=FALSE]
    } else {
      ME <- mergeME(ME, Y[,phenotype,drop=FALSE])
    }
  }

  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!

  ## Compute eigengene correlation matrix.
  ## Repeat 'power' times for higher order adjacency.
  power <- round(power)
  R <- ME
  for (i in 1:power) {
    tt <- cortest(R, R)
    R <- tt$rho
    nSamples <- tt$n
  }

  ## If phenotype is given, we condition the heatmap
  ## using the correlation to the phenotype.
  if (!is.null(phenotype)) {
    ## proper sign in case of inhibitor layer (like miRNA)
    layersign <- rep(1, length(wgcna))
    names(layersign) <- names(wgcna)
    layersign[grep("^mi", names(wgcna), ignore.case = TRUE)] <- -1
    ff <- list()
    for (k in names(wgcna)) {
      rho <- cor(ME, Y[,phenotype], use="pairwise")[,1]
      ff[[k]] <- layersign[k] * rho
    }
    names(ff) <- NULL
    ff <- unlist(ff)
    ff <- 0.5 * (1 + ff)
    ff <- ff[match(rownames(R), names(ff))]
    names(ff) <- rownames(R)
    ff[is.na(ff)] <- 1
    ww <- outer(ff, ff)
    ww <- ww / max(ww, na.rm = TRUE)
    R <- R * ww
  }

  if (nmax>0) {
    if (!is.null(phenotype)) {
      y1 <- Y[,phenotype]
      y1 <- y1[match(rownames(ME),names(y1))]
      rho <- cor(ME, y1, use="pairwise")[,1]
      ii <- head(order(-abs(rho)), nmax)
    } else {
      ii <- head(order(-Matrix::rowMeans(R**2)), nmax)
    }
    R <- R[ii, ii]
  }

  if (justdata) return(R)

  # Plot correlation heatmap matrix.
  # ps: this plot will overwrite he dendrogram plot
  if (is.null(main)) main <- "Eigengene Adjacency Heatmap"

  if (plotDendro && plotHeatmap) {
    layout.matrix <- matrix(1:2, nrow = 1, ncol = 2)
    layout(mat = layout.matrix, heights = 1, widths = c(dendro.width, 1))
    if (dendro.horiz && dendro.labels) {
      mar1[4] <- mar2[2] ## copy left margin
    }
  }
  if (plotDendro) par(mar = mar1)

  R0 <- R
  R0[is.na(R0)] <- 0

  if (fixclust) {
    ii <- rownames(R)
    hc <- hclust(as.dist(1 - R0[ii, ii]), method = "average")
  } else {
    hc <- hclust(as.dist(1 - R0), method="average")
  }
  if (plotDendro) {
    par(cex = cex.lab)
    plot(as.dendrogram(hc), horiz = TRUE, ylab = "Eigengene dendrogram")
    par(cex = 1)
  }

  if (plotHeatmap) {
    ii <- hc$labels[hc$order]
    ii <- intersect(ii, rownames(R))
    R1 <- R[rev(ii), ii]
    nsamples <- nSamples[rownames(R1),colnames(R1)]
    par(mar=mar2)
    plotLabeledCorrelationHeatmap(
      R1,
      nSamples = nsamples,
      text = text,
      pstar = pstar,
      colorlabel = colorlabel,
      cluster = FALSE,
      setpar = FALSE,
      main = main,
      cex.lab = cex.lab,
      cex.text = cex.text
    )
  }

  invisible(R)

}

#' Plot inter-dataset eigengene correlation heatmaps
#' @param wgcna Named list of WGCNA objects.
#' @param addtraits Include traits in correlation.
#' @param phenotype Phenotype for conditioning.
#' @param nmax Maximum number of eigengenes shown.
#' @param main Plot title string.
#' @param showvalues Show correlation values.
#' @param showsig Show significance stars.
#' @param cex.text Text character expansion factor.
#' @param cex.lab Label character expansion factor.
#' @param fixcluster Fix cluster ordering.
#' @param setpar Set plotting parameters.
#' @return NULL (invisible). Generates a plot.
#' @export
plotMultiEigengeneCorrelation <- function(wgcna,
                                          addtraits = TRUE,
                                          phenotype = NULL,
                                          nmax = -1,
                                          main = NULL,
                                          showvalues = FALSE,
                                          showsig = TRUE,
                                          cex.text = 0.7,
                                          cex.lab = 0.8,
                                          fixcluster = TRUE,
                                          setpar = TRUE) {

  ## Show inter-correlation of modules
  me <- lapply(wgcna, function(w) w$net$MEs)
  if (length(me) == 1) {
    me <- list(me[[1]], me[[1]])
  }

  comb <- combn(length(me), 2)
  ncomb <- ncol(comb)
  nsamples <- nrow(wgcna[[1]]$datExpr)
  Y <- wgcna[[1]]$datTraits

  ## for miRNA we need to flip sign
  msign <- c(1, -1)[1 + 1 * (names(wgcna) %in% c("mi", "mir"))]

  if (setpar) {
    nc <- ceiling(sqrt(ncomb))
    nr <- ceiling(ncomb / nc)
    par(mfrow = c(nr, nc), mar = c(8, 10, 3, 1))
  }

  for (k in 1:ncol(comb)) {
    i <- comb[1, k]
    j <- comb[2, k]
    M1 <- me[[i]]
    M2 <- me[[j]]

    if (addtraits) {
      M1 <- cbind(M1, Y)
      M2 <- cbind(M2, Y)
    }

    if (FALSE && !addtraits && !is.null(phenotype)) {
      y <- Y[, phenotype, drop = FALSE]
      M1 <- cbind(M1, y)
      M2 <- cbind(M2, y)
    }

    R1 <- cor(M1, M2, use = "pairwise.complete")

    if (nmax > 0) {
      ii <- head(order(-apply(abs(R1), 1, max)), nmax)
      jj <- head(order(-apply(abs(R1), 2, max)), nmax)
      R1 <- R1[ii, jj]
    }

    ## cluster unweighted matrix
    ii <- hclust(dist(R1), method = "average")$order
    jj <- hclust(dist(t(R1)), method = "average")$order
    R1 <- R1[ii, jj]

    ## This conditions the correlation on phenotype. Important.
    do.condition <- !is.null(phenotype)
    if (do.condition) {
      y <- Y[, phenotype]
      w1 <- cor(M1[, rownames(R1)], y, use = "pairwise")[, 1]
      w2 <- cor(M2[, colnames(R1)], y, use = "pairwise")[, 1]
      if (msign[i] != 0) w1 <- msign[i] * w1
      if (msign[j] != 0) w2 <- msign[j] * w2
      w1 <- pmax(w1, 0)
      w2 <- pmax(w2, 0)
      ww <- outer(w1, w2)
      ww <- ww / max(ww, na.rm = TRUE)
      R1 <- R1 * ww
    }

    main <- paste(names(me)[i], "vs.", names(me)[j])
    if (do.condition) main <- paste(main, "(conditioned)")

    plotLabeledCorrelationHeatmap(
      R1,
      nsamples,
      text = showvalues,
      pstar = showsig,
      is.dist = FALSE,
      cluster = !fixcluster,
      setpar = FALSE,
      main = main,
      cex.text = cex.text,
      cex.lab = cex.lab
    )
  }

}


#' Plot eigengene network as graph
#' @param wgcna A WGCNA result object.
#' @param add_traits Include traits in graph.
#' @param main Plot title string.
#' @param multi Use multi-dataset mode.
#' @param vcex Vertex size scaling factor.
#' @param labcex Label size scaling factor.
#' @return NULL (invisible). Generates a plot.
#' @export
plotEigenGeneGraph <- function(wgcna,
                               add_traits = TRUE,
                               main = NULL,
                               multi = FALSE,
                               vcex = 1,
                               labcex = 1) {

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required for eigengene graph plotting")
  }

  if (multi) {
    ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
    ME <- mergeME(ME)
    if (add_traits)  ME <- cbind(ME, wgcna[[1]]$datTraits)
  } else {
    ME <- wgcna$net$MEs
    if (add_traits) ME <- cbind(ME, wgcna$datTraits)
  }

  if (NCOL(ME) <= 2) ME <- cbind(ME, ME)

  sdx <- matrixStats::colSds(as.matrix(ME * 1), na.rm = TRUE)
  if (any(sdx == 0)) ME <- ME + runif(length(ME), 0, 1e-5)

  ## Recalculate MEs with color as labels
  corx <- cor(ME, use = "pairwise")
  corx[is.na(corx)] <- 0
  clust <- hclust(as.dist(1 - corx))
  phylo <- ape::as.phylo(clust)
  gr <- igraph::as.igraph(phylo, directed = FALSE)

  is.node <- grepl("Node", igraph::V(gr)$name)
  module.name <- igraph::V(gr)$name
  if (multi) {
    module.size <- lapply(wgcna, function(w) table(w$net$labels))
    names(module.size) <- NULL
    module.size <- unlist(module.size)
    module.colors <- sapply(wgcna, function(w) w$me.colors)
    names(module.colors) <- NULL
    module.colors <- unlist(module.colors)
  } else {
    module.size <- table(wgcna$net$labels)
    module.colors <- wgcna$me.colors
  }
  module.size <- module.size / mean(module.size)
  igraph::V(gr)$label <- igraph::V(gr)$name
  igraph::V(gr)$label[is.node] <- NA
  igraph::V(gr)$color <- module.colors[module.name]
  igraph::V(gr)$size <- vcex * 18 * (module.size[module.name])**0.4
  igraph::V(gr)$size[is.na(igraph::V(gr)$size)] <- 0

  ## par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0)
  igraph::plot.igraph(
    gr,
    layout = igraph::layout.kamada.kawai,
    vertex.label.cex = 0.85 * labcex,
    edge.width = 3
  )
  if (!is.null(main)) title(main, line = -1.5)

}


#' Plot consensus sample dendrogram with colors
#' @param cons A consensus WGCNA result object.
#' @param i Index of consensus layer.
#' @param what What to show: "both", "me", or "traits".
#' @param show.me Show module eigengenes.
#' @param show.traits Show trait values.
#' @param show.contrasts Show contrast values.
#' @param clust.expr Cluster by expression similarity.
#' @param setLayout Whether to set layout.
#' @param marAll Margin sizes vector.
#' @param colorHeightMax Maximum color row height.
#' @param main Plot title string.
#' @return NULL (invisible). Generates a plot.
#' @export
plotConsensusSampleDendroAndColors <- function(cons, i,
                                               what = c("both", "me", "traits")[1],
                                               show.me = TRUE,
                                               show.traits = TRUE,
                                               show.contrasts = TRUE,
                                               clust.expr = TRUE,
                                               setLayout = TRUE,
                                               marAll = c(0.2, 7, 1.5, 0.5),
                                               colorHeightMax = 0.6,
                                               main = NULL) {

  plotSampleDendroAndColors(
    wgcna = cons$layers[[i]],
    main = main,
    datExpr = cons$datExpr[[i]],
    datTraits = cons$datTraits,
    datME = cons$net$multiME[[i]]$data,
    what = what,
    show.me = show.me,
    show.traits = show.traits,
    show.contrasts = show.contrasts,
    marAll = marAll,
    clust.expr = clust.expr,
    setLayout = setLayout,
    colorHeightMax = colorHeightMax
  )

}

#' Plot sample dendrogram with trait colors
#' @param wgcna A WGCNA result object.
#' @param input.type Input type: "wgcna" or "net".
#' @param what What to show: "me", "traits", or "both".
#' @param show.me Show module eigengenes.
#' @param show.traits Show trait values.
#' @param show.contrasts Show contrast values.
#' @param datTraits Trait data matrix override.
#' @param datExpr Expression data matrix override.
#' @param datME Eigengene data matrix override.
#' @param clust.expr Cluster by expression similarity.
#' @param setLayout Whether to set layout.
#' @param marAll Margin sizes vector.
#' @param colorHeightMax Maximum color row height.
#' @param main Plot title string.
#' @param justdata Return data without plotting.
#' @return NULL (invisible). Generates a plot. If justdata, returns eigengene matrix.
#' @export
plotSampleDendroAndColors <- function(wgcna,
                                      input.type = "wgcna",
                                      what = c("me", "traits", "both")[3],
                                      show.me = TRUE,
                                      show.traits = TRUE,
                                      show.contrasts = TRUE,
                                      datTraits = NULL,
                                      datExpr = NULL,
                                      datME = NULL,
                                      clust.expr = TRUE,
                                      setLayout = TRUE,
                                      marAll = c(0.2, 7, 1.5, 0.5),
                                      colorHeightMax = 0.6,
                                      main = NULL,
                                      justdata = FALSE) {

  if (input.type == "net") {
    ME0 <- wgcna$MEs
    if (is.null(datExpr)) stop("must supply datExpr")
    if (is.null(datTraits)) stop("must supply datTraits")
  } else {
    ME0 <- wgcna$net$MEs
    datTraits <- 1 * wgcna$datTraits
    datExpr <- wgcna$datExpr
  }

  if (!is.null(datME)) ME0 <- datME

  ME <- ME0[, 0]
  samples <- rownames(ME)
  if (show.me) ME <- cbind(ME, ME0)

  if (show.traits) {
    sel <- grep("_vs_", colnames(datTraits), invert = TRUE)
    ME <- cbind(ME, datTraits[samples, sel, drop = FALSE])
  }

  if (show.contrasts) {
    sel <- grep("_vs_", colnames(datTraits))
    ME <- cbind(ME, datTraits[samples, sel, drop = FALSE])
  }

  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!
  sdx <- matrixStats::colSds(as.matrix(ME * 1), na.rm = TRUE)
  ME <- ME[, which(sdx > 0), drop = FALSE]

  ## Recalculate MEs with color as labels
  if (clust.expr) {
    corx <- cor(t(datExpr), use = "pairwise")
  } else {
    corx <- cor(t(ME0), use = "pairwise")
  }
  corx[is.na(corx)] <- 0
  sampleTree <- hclust(as.dist(1 - corx), method = "average")

  corx <- cor(ME, use = "pairwise")
  corx[is.na(corx)] <- 0
  jj <- hclust(as.dist(1 - corx))$order
  colors <- WGCNA::numbers2colors(ME[, jj])

  if (justdata) return(ME)

  if (is.null(main)) {
    if (what == "me") main <- "Sample dendrogram and module heatmap"
    if (what == "traits") main <- "Sample dendrogram and trait heatmap"
    if (what == "both") main <- "Sample dendrogram and module+traits heatmap"
  }

  WGCNA::plotDendroAndColors(
    dendro = sampleTree,
    colors = colors,
    groupLabels = colnames(ME)[jj],
    dendroLabels = rownames(ME),
    hang = 0.03,
    addGuide = FALSE,
    guideHang = 0.05,
    setLayout = setLayout,
    marAll = marAll,
    main = main,
    colorHeightMax = colorHeightMax
  )

}
