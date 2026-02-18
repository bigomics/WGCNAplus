#' Plot gene dendrogram with module colors
#'
#' @param wgcna A WGCNA result object.
#' @param main Plot title string.
#' @param extra.colors Additional color rows to display.
#' @param show.kme Show module eigengene correlations.
#' @param show.traits Show trait correlations.
#' @param show.contrasts Show contrast correlations.
#' @param show.mat Custom matrix to show.
#' @param clust Cluster columns by correlation.
#' @param use.tree Tree index for consensus.
#' @param block Block number to plot.
#' @param rm.na Remove NA-only columns.
#' @param sd.wt Standard deviation weighting factor.
#' @param nmax Maximum number of columns shown.
#' @param marAll Margin sizes vector.
#' @param setLayout Whether to set layout.
#' @param ... Additional arguments to WGCNA::plotDendroAndColors.
#' @return NULL (invisible). Generates a plot.
#' @export
plotDendroAndColors <- function(wgcna, main=NULL,
                                extra.colors=NULL,
                                show.kme = FALSE,
                                show.traits = FALSE,
                                show.contrasts = FALSE,
                                show.mat = NULL,
                                clust = TRUE,
                                use.tree = 0,
                                block = 1,
                                rm.na = TRUE,
                                sd.wt = 0,
                                nmax = -1,
                                marAll = c(0.4, 5, 1, 0.2),
                                setLayout = TRUE, ...) {

  if ("net" %in% names(wgcna)) {
    net <- wgcna$net
    if ("colors" %in% names(wgcna)) {
      net$netcolors <- net$colors
      net$colors <- wgcna$colors
    }
  } else {
    net <- wgcna
  }

  ## select dendrogram
  dendro <- net$dendrograms
  if ("layers" %in% names(wgcna) && use.tree > 0) {
    dendro <- wgcna$layers[[use.tree]]$net$dendrograms
  }
  if (length(dendro) > 1) {
    message("warning: this wgcna has multiple blocks")
  }
  geneTree <- dendro[[block]]

  colors <- cbind(labels2colors(net$colors))
  if (NCOL(colors) == 1) colnames(colors)[1] <- "Module colors"

  gg <- geneTree$labels
  if (is.null(gg) && !is.null(block)) {
    ii <- which(net$blocks == block & net$goodGenes == TRUE)
    gg <- names(net$color)[ii]
  }
  if (is.null(gg) && is.null(block)) {
    ii <- which(net$goodGenes == TRUE)
    gg <- names(net$color)[ii]
  }
  ## colors <- colors[gg,,drop=FALSE]
  colors <- colors[match(gg, rownames(colors)), , drop = FALSE]
  if (!is.null(extra.colors)) {
    jj <- match(gg, rownames(extra.colors))
    colors <- cbind(colors, 0, extra.colors[jj, ])
  }

  calcKMEcolors <- function(X, Y) {
    kme <- cor(X, Y, use="pairwise")
    kme[is.na(kme)] <- 0
    sdx <- matrixStats::colSds(X,na.rm=TRUE)
    if(sd.wt>0) kme <- kme * (sdx / max(abs(sdx),na.rm=TRUE))**sd.wt
    if(nmax>0 && nmax<ncol(kme)) {
      sel <- head(order(-colMeans(kme**2)),nmax)
      kme <- kme[,sel,drop=FALSE]
    }
    kmeColors <- rho2bluered(kme, a = 0.8)
    kmeColors <- kmeColors[gg,,drop=FALSE]
    if(clust && ncol(kme)>2) {
      cor.kme <- cor(kme,use="pairwise")
      cor.kme[is.na(cor.kme)] <- 0
      ii <- hclust(as.dist(1 - cor.kme))$order
      kmeColors <- kmeColors[,ii,drop=FALSE]
    }
    kmeColors
  }

  is.multi <- is.list(wgcna$datExpr)
  is.multi
  if(!is.multi) {
    if(show.kme) {
      X <- wgcna$datExpr
      Y <- net$MEs
      kmeColors <- calcKMEcolors(X, Y)
      colors <- cbind(colors, 0, kmeColors)
    }
    if (show.traits) {
      X <- wgcna$datExpr
      Y <- wgcna$datTraits
      Y <- Y[,grep("_vs_",colnames(Y),invert=TRUE),drop=FALSE]
      if(NCOL(Y)>0) {
        kmeColors <- calcKMEcolors(X, Y)
        colors <- cbind(colors, 0, kmeColors)
      }
    }
    if (show.contrasts) {
      X <- wgcna$datExpr
      Y <- wgcna$datTraits
      Y <- Y[,grep("_vs_",colnames(Y)),drop=FALSE]
      if(NCOL(Y)>0) {
        kmeColors <- calcKMEcolors(X, Y)
        colors <- cbind(colors, 0, kmeColors)
      }
    }
    if(!is.null(show.mat)) {
      X <- wgcna$datExpr
      Y <- show.mat
      if(NCOL(Y)>0) {
        kmeColors <- calcKMEcolors(X, Y)
        colors <- cbind(colors, 0, kmeColors)
      }
    }
  }

  if (is.multi) {
    if (show.kme) {
      X <- wgcna$datExpr
      Y <- wgcna$net$multiMEs
      for(i in 1:length(X)) {
        kmeColors <- calcKMEcolors(X[[i]], Y[[i]]$data)
        colnames(kmeColors) <- paste0(names(X)[i],":",colnames(kmeColors))
        colors <- cbind(colors, 0, kmeColors)
      }
    }
    if (show.traits) {
      X <- wgcna$datExpr
      Y <- wgcna$datTraits
      for(i in 1:length(X)) {
        Y1 <- Y[rownames(X[[i]]),]
        Y1 <- Y1[,grep("_vs_",colnames(Y1),invert=TRUE),drop=FALSE]
        if(NCOL(Y1)) {
          kmeColors <- calcKMEcolors(X[[i]], Y1)
          colnames(kmeColors) <- paste0(names(X)[i],":",colnames(kmeColors))
          colors <- cbind(colors, 0, kmeColors)
        }
      }
    }
    if (show.contrasts) {
      X <- wgcna$datExpr
      Y <- wgcna$datTraits
      for(i in 1:length(X)) {
        Y1 <- Y[rownames(X[[i]]),]
        Y1 <- Y1[,grep("_vs_",colnames(Y1)),drop=FALSE]
        if(NCOL(Y1)) {
          kmeColors <- calcKMEcolors(X[[i]], Y1)
          colnames(kmeColors) <- paste0(names(X)[i],":",colnames(kmeColors))
          colors <- cbind(colors, 0, kmeColors)
        }
      }
    }
    if(!is.null(show.mat)) {
      X <- wgcna$datExpr
      Y <- show.mat
      for(i in 1:length(X)) {
        Y1 <- Y[rownames(X[[i]]),]
        if(NCOL(Y1)) {
          kmeColors <- calcKMEcolors(X[[i]], Y1)
          colnames(kmeColors) <- paste0(names(X)[i],":",colnames(kmeColors))
          colors <- cbind(colors, 0, kmeColors)
        }
      }
    }
  }

  if (rm.na) {
    all.eq <- rowMeans(t(colors) == colors[1, ]) == 1
    sel <- colMeans(is.na(colors)) < 1 & !all.eq
    sel <- sel | colnames(colors) %in% c("", NA)
    colors <- colors[, sel, drop = FALSE]
  }

  if (is.null(main)) main <- "Gene dendrogram and module colors"
  ## Plot the dendrogram and the module colors underneath
  WGCNA::plotDendroAndColors(
    dendro = geneTree,
    colors = colors,
    groupLabels = colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = FALSE,
    guideHang = 0.05,
    marAll = marAll,
    setLayout = setLayout,
    main = main,
    ...
  )

}


#' Plot dendrograms for multiple WGCNA objects
#'
#' @param multi_wgcna Named list of WGCNA objects.
#' @param block Block number to plot.
#' @param extra.colors Additional color rows to display.
#' @param show.kme Show module eigengene correlations.
#' @param show.traits Show trait correlations.
#' @param show.contrasts Show contrast correlations.
#' @param show.mat Custom matrix to show.
#' @param clust Cluster columns by correlation.
#' @param use.tree Tree index for consensus.
#' @param rm.na Remove NA-only columns.
#' @param sd.wt Standard deviation weighting factor.
#' @param nmax Maximum number of columns shown.
#' @param main Plot title string.
#' @param colorHeight Relative height of color rows.
#' @param marAll Margin sizes vector.
#' @param cex Character expansion factor.
#' @return NULL (invisible). Generates a plot.
#' @export
plotMultiDendroAndColors <- function(wgcna,
                                     block = 1,
                                     extra.colors = NULL,
                                     show.kme = FALSE,
                                     show.traits = FALSE,
                                     show.contrasts = FALSE,
                                     show.mat = NULL,
                                     clust = TRUE,
                                     use.tree = 0,
                                     rm.na = TRUE,
                                     sd.wt = 0,
                                     nmax = -1,
                                     main = 'Dendro and colors',
                                     colorHeight = 0.5,
                                     marAll = c(0.4,5,1,0.2),
                                     cex = 1
                                     )
{

  layers <- wgcna
  if(!is.null(wgcna$layers)) {
    layers <- wgcna$layers
  } else {
    layers <- list( ' '=wgcna )
  }
  
  nw <- length(layers)
  nc <- ceiling(sqrt(nw))
  nr <- ceiling(nw / nc)
  hh <- rep(c((1 - colorHeight), colorHeight), nr)
  hh
  nf <- layout(matrix(1:(2*nr*nc), nrow=2*nr, ncol=nc,
    byrow=FALSE), heights = hh )

  tt <- paste0(main, " ",toupper(names(layers)))

  ##layout.show(nf)
  par(cex=cex)
  for (k in 1:nw) {
    plotDendroAndColors(
      layers[[k]],
      marAll = marAll,
      show.traits = show.traits,
      show.contrasts = show.contrasts,
      show.kme = show.kme,
      show.mat = show.mat,
      use.tree = use.tree,
      clust = clust,
      sd.wt = sd.wt,
      nmax = nmax,
      setLayout = FALSE,
      main = tt[k]
    )
  }
}


#' plotDendroAndTraits for Consensus output
#'
#' @param cons A consensus WGCNA result object.
#' @param show.traits Show trait correlations.
#' @param show.contrasts Show contrast correlations.
#' @param traits Traits to include.
#' @param main Plot title string.
#' @param rm.na Remove NA-only columns.
#' @param use.tree Tree index for consensus.
#' @param marAll Margin sizes vector.
#' @param setLayout Whether to set layout.
#' @param ... Additional arguments passed on.
#' @return NULL (invisible). Generates a plot.
plotDendroAndColors_cons <- function(cons,
                                     show.traits = TRUE,
                                     show.contrasts = TRUE,
                                     traits = NULL,
                                     main = NULL,
                                     rm.na = TRUE,
                                     use.tree = 0,
                                     marAll = c(0.2, 8, 2, 0.2),
                                     setLayout = TRUE,
                                     ...) {
  message("DEPRECATED: please use plotDendroAndColors")
  
  if (0) {
    show.traits <- TRUE
    traits <- NULL
    main <- NULL
    rm.na <- TRUE
    use.tree <- 0
    marAll <- c(0.2, 8, 2, 0.2)
    setLayout <- TRUE
  }

  ## quick hack to use plotDendroAndColors_multi()
  multi <- c(list(Consensus = cons), cons$layers)
  use.tree0 <- use.tree
  if (use.tree %in% 0:99) use.tree <- as.integer(use.tree)
  if (is.character(use.tree)) {
    use.tree <- match(use.tree, names(multi))
  } else {
    use.tree <- as.integer(use.tree) + 1
  }
  if (is.na(use.tree)) {
    message("ERROR: invalid class(use.tree) = ", class(use.tree0))
    message("ERROR: invalid use.tree = ", use.tree0)
    return(NULL)
  }

  plotDendroAndColors_multi(
    multi,
    show.traits = show.traits,
    show.contrasts = show.contrasts,
    traits = traits,
    main = main,
    rm.na = rm.na,
    use.tree = use.tree,
    marAll = marAll,
    setLayout = setLayout,
    ...
  )
}


#' Plot dendrogram with traits for multiple datasets
#'
#' @param multi Named list of WGCNA objects.
#' @param show.traits Show trait correlations.
#' @param show.contrasts Show contrast correlations.
#' @param traits Traits to include.
#' @param main Plot title string.
#' @param rm.na Remove NA-only columns.
#' @param use.tree Which tree layer to use.
#' @param marAll Margin sizes vector.
#' @param setLayout Whether to set layout.
#' @param ... Additional arguments passed on.
#' @return NULL (invisible). Generates a plot.
#' @export
plotDendroAndColors_multi <- function(multi,
                                      show.traits = TRUE,
                                      show.contrasts = TRUE,
                                      traits = NULL,
                                      main = NULL,
                                      rm.na = TRUE,
                                      use.tree = 1,
                                      marAll = c(0.2, 8, 2, 0.2),
                                      setLayout = TRUE,
                                      ...) {
  message("DEPRECATED: please use plotDendroAndColors")

  ## module colors
  colors <- sapply(multi, function(m) m$net$colors)

  if (show.traits || show.contrasts) {
    traitSig <- list()
    nsets <- length(multi)
    i <- 1
    for (k in names(multi)) {
      if (k == "Consensus") next
      w <- multi[[k]]
      Y <- w$datTraits
      sel1 <- sel2 <- NULL
      if (show.traits) sel1 <- grep("_vs_", colnames(Y), invert = TRUE)
      if (show.contrasts) sel2 <- grep("_vs_", colnames(Y))
      sel <- c(sel1, sel2)
      if (!is.null(traits)) sel <- intersect(sel, traits)
      X <- w$datExpr
      kk <- intersect(rownames(X), rownames(Y))
      traitSig[[k]] <- cor(X[kk, ], Y[kk, sel], use = "pairwise")
    }

    if (rm.na) {
      for (i in 1:length(traitSig)) {
        sel <- colMeans(is.na(traitSig[[i]])) < 1
        traitSig[[i]] <- traitSig[[i]][, sel, drop = FALSE]
      }
    }

    ## prepend datatype/set name
    for (k in names(traitSig)) {
      colnames(traitSig[[k]]) <- paste0(k, ":", colnames(traitSig[[k]]))
    }

    traitSig2 <- c()
    for (i in 1:length(traitSig)) {
      traitSig2 <- cbind(traitSig2, traitSig[[i]])
      if (i < length(traitSig)) traitSig2 <- cbind(traitSig2, 0)
    }
    traitColors <- rho2bluered(traitSig2, f = 0.95, a = 0.8)
    ii <- which(colnames(traitColors) == "")
    if (length(ii)) traitColors[, ii] <- "#FFFFFF"
    if (is.null(colors)) {
      colors <- traitColors
    } else {
      colors <- cbind(colors, 0, traitColors)
    }
  }

  message("using tree of layer: ", names(multi)[use.tree])
  geneTree <- multi[[use.tree]]$net$dendrograms[[1]]

  if (is.null(main)) main <- "Gene Dendrogram, Modules and Trait Correlation"

  WGCNA::plotDendroAndColors(
    geneTree,
    colors = colors,
    colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    marAll = marAll,
    main = main,
    ...
  )
}
