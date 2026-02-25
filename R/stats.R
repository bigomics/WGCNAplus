## ---------------------------------------------------------------------
## Gene statistics
## ---------------------------------------------------------------------

#' Pairwise correlation test with p-values
#' @param X Numeric matrix of variables.
#' @param Y Numeric matrix of variables.
#' @return List with rho, pvalue, and n.
#' @keywords internal
cortest <- function(X, Y) {

  rho <- cor(X, Y, use = "pairwise.complete")
  nSamples <- t(!is.na(X)) %*% (!is.na(Y))
  ii <- which(nSamples < 3)
  nSamples0 <- pmax(nSamples, 3)
  Pvalue <- WGCNA::corPvalueStudent(rho, nSamples0)
  if (length(ii)) Pvalue[ii] <- NA
  return(list(rho = rho, pvalue = Pvalue, n = nSamples))

}

#' Compute general feature statistics after WGCNA results.
#' @param net WGCNA network object with MEs and colors.
#' @param datExpr Expression data matrix.
#' @param datTraits Trait data matrix.
#' @param TOM Topological overlap matrix or NULL.
#' @return List of gene statistic matrices.
computeGeneStats <- function(net,
                             datExpr,
                             datTraits,
                             TOM) {

  kk <- intersect(rownames(datExpr), rownames(datTraits))
  datExpr <- datExpr[kk, , drop = FALSE]
  datTraits <- datTraits[kk, , drop = FALSE]

  ## Define numbers of genes and samples
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)

  ## Recalculate MEs with color labels
  tt <- cortest(net$MEs, datTraits)
  moduleTraitCor <- tt$rho
  moduleTraitPvalue <- tt$pvalue

  ## Module membership correlation (with p-values)
  mm <- cortest(datExpr, net$MEs)
  moduleMembership <- mm$rho
  MMPvalue <- mm$pvalue

  ## Gene-trait significance (trait correlation) (with p-values)
  ts <- cortest(datExpr, datTraits)
  traitSignificance <- ts$rho
  TSPvalue <- ts$pvalue

  ## Fold-change
  foldChange <- NULL
  foldChangePvalue <- NULL
  is.binary <- apply(datTraits, 2, function(a) length(unique(a[!is.na(a)])) == 2)
  binY <- NULL
  if (any(is.binary)) {
    binY <- datTraits[, which(is.binary), drop = FALSE]
    nmin <- apply(binY, 2, function(x) min(table(x)))
    binY <- binY[, which(nmin >= 2), drop = FALSE]
  }
  if (!is.null(binY) && NCOL(binY) > 0) {
    lm <- list()
    for (i in 1:ncol(binY)) {
      y <- 1 * binY[, i]
      X <- t(datExpr)
      suppressWarnings(suppressMessages(
        res <- try(gx.limma(X, y, lfc = 0, fdr = 1, sort.by = "none", verbose = 0, max.na = 1))
      ))
      if (!"try-error" %in% class(res)) {
        k <- colnames(binY)[i]
        lm[[k]] <- res
      }
    }
    lm <- lm[!sapply(lm, is.null)]
    foldChange <- sapply(lm, function(m) m$logFC)
    foldChangePvalue <- sapply(lm, function(m) m$P.Value)
    if (length(lm) == 1) {
      foldChange <- cbind(foldChange)
      foldChangePvalue <- cbind(foldChangePvalue)
    }
    rownames(foldChange) <- rownames(lm[[1]])
    rownames(foldChangePvalue) <- rownames(lm[[1]])
  }

  # Continuous traits (not always present)
  contY <- datTraits[, which(!is.binary), drop = FALSE]
  foldChange.cont <- NULL
  foldChangePvalue.cont <- NULL
  if (NCOL(contY) > 0) {
    contlm <- apply(contY, 2, function(y) {
      tt <- cortest(datExpr, y)
      rho <- tt$rho[,1]
      P.Value <- tt$pvalue[,1]
      data.frame(rho, P.Value)
    })
    contlm <- contlm[!sapply(contlm, is.null)]
    foldChange.cont <- sapply(contlm, function(m) m$rho)
    foldChangePvalue.cont <- sapply(contlm, function(m) m$P.Value)
    if (length(contlm) == 1) {
      foldChange.cont <- cbind(foldChange.cont)
      foldChangePvalue.cont <- cbind(foldChangePvalue.cont)
    }
    rownames(foldChange.cont) <- rownames(contlm[[1]])
    rownames(foldChangePvalue.cont) <- rownames(contlm[[1]])
    colnames(foldChange.cont) <- colnames(contY)
    colnames(foldChangePvalue.cont) <- colnames(contY)
  }

  # Merge
  foldChange <- cbind(foldChange, foldChange.cont)
  foldChangePvalue <- cbind(foldChangePvalue, foldChangePvalue.cont)

  ## Gene Centrality. Compute centrality of gene in Module subgraph
  ## using TOM matrix. WARNING: this can create large TOM matrix
  geneCentrality <- NULL
  if (!is.null(TOM)) {
    if (nrow(TOM) != ncol(TOM)) TOM <- TOM %*% Matrix::t(TOM)
    if (is.null(dimnames(TOM))) {
      dimnames(TOM) <- list(colnames(datExpr), colnames(datExpr))
    }
    diag(TOM) <- 0
    TOM[which(abs(TOM) < 0.01)] <- 0
    gr <- igraph::graph_from_adjacency_matrix(TOM, mode = "undirected",
      weighted = TRUE, diag = FALSE)
    geneCentrality <- rep(NA, nrow(TOM))
    names(geneCentrality) <- rownames(TOM)
    me.genes <- tapply(names(net$colors), net$colors, list)
    for (gg in me.genes) {
      gr1 <- igraph::subgraph(gr, gg)
      ct <- igraph::page_rank(gr1, weights = NULL)$vector
      ct <- ct / mean(ct, na.rm = TRUE)
      geneCentrality[gg] <- ct
    }
  }

  ## force align. Sometime they are shorter for some reason...
  gg <- rownames(moduleMembership)
  matMatch <- function(m, gg) {
    if (is.null(m)) return(NULL)
    m <- m[match(gg, rownames(m)), , drop = FALSE]
    rownames(m) <- gg
    return(m)
  }

  MMPvalue <- matMatch(MMPvalue, gg)
  traitSignificance <- matMatch(traitSignificance, gg)
  TSPvalue <- matMatch(TSPvalue, gg)
  foldChange <- matMatch(foldChange, gg)
  foldChangePvalue <- matMatch(foldChangePvalue, gg)

  stats <- list(
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    moduleMembership = moduleMembership,
    MMPvalue = MMPvalue,
    traitSignificance = traitSignificance,
    TSPvalue = TSPvalue,
    foldChange = foldChange,
    foldChangePvalue = foldChangePvalue,
    geneCentrality = geneCentrality
  )

  return(stats)

}

#' Get gene statistics for a trait and module
#' @param wgcna WGCNA result object or NULL.
#' @param trait Trait name or names.
#' @param module Module name to filter or NULL.
#' @param plot Logical; create scatterplot pairs.
#' @param stats Precomputed stats list or NULL.
#' @param labels Module label vector or NULL.
#' @param showlogFC Logical; show fold-change in plot.
#' @param col Color vector for plot points.
#' @param main Plot title string or NULL.
#' @return Data frame of gene statistics.
#' @export
getGeneStats <- function(wgcna,
                         trait,
                         module = NULL,
                         plot = TRUE,
                         stats = NULL,
                         labels = NULL,
                         showlogFC = TRUE,
                         col = NULL,
                         main = NULL) {

  if (!is.null(stats)) {
    features <- rownames(stats[["moduleMembership"]])
    if (is.null(stats)) stop("must give stats")
    if (is.null(labels)) stop("must give labels")
  } else if (!is.null(wgcna)) {
    labels <- wgcna$net$labels
    features <- names(wgcna$net$colors)
    stats <- wgcna$stats
  } else {
    stop("must supply wgcna or trait")
  }

  is.color <- mean(labels %in% c("grey", WGCNA::standardColors())) > 0.9
  if (is.color) {
    prefix <- substring(rownames(stats[["moduleTraitCor"]]), 1, 2)[1]
    labels <- paste0(prefix, labels)
  }

  p1 <- c("moduleMembership", "MMPvalue")
  p2 <- c("traitSignificance", "TSPvalue", "foldChange", "foldChangePvalue")
  p3 <- c("geneCentrality")

  df <- data.frame(feature = features, module = labels)

  ## get moduleMembership
  mm.stats <- stats[p1]
  mm.label <- colnames(mm.stats[[1]])
  idx <- cbind(1:length(labels), match(labels, mm.label))
  A1 <- sapply(mm.stats, function(x) x[idx])
  rownames(A1) <- labels
  df <- cbind(df, A1)

  ## get traitSig columns for trait
  tt.cols <- colnames(stats[[p2[1]]])
  if (is.null(trait)) trait <- tt.cols
  trait <- intersect(trait, tt.cols)

  if (length(trait) > 1) {
    A2 <- lapply(stats[p2], function(x) x[, trait])
    for (i in 1:length(A2)) colnames(A2[[i]]) <- paste0(names(A2)[i], ".", colnames(A2[[i]]))
    A2 <- do.call(cbind, A2)
    df <- cbind(df, A2)
  } else if (length(trait) == 1) {
    A2 <- lapply(stats[p2], function(x) x[, trait])
    A2 <- do.call(cbind, A2)
    df <- cbind(df, A2)
  } else {
    message("[getGeneStats] ERROR: trait not in stats object")
    return(NULL)
  }

  A3 <- stats[[p3]]
  if (!is.null(A3)) {
    df <- cbind(df, centrality = A3)
  }

  ## calculate score
  sel <- c("moduleMembership", "traitSignificance", "foldChange", "centrality")
  sel <- intersect(sel, colnames(df))
  df1 <- as.matrix(abs(df[, sel]))
  score <- exp(rowMeans(log(1e-8 + df1))) * sign(df[, "foldChange"])
  df <- data.frame(df[, 1:2], score = score, df[, -c(1, 2)])

  if (!is.null(module)) {
    sel <- grep(paste0(module, "$"), df$module)
    df <- df[sel, , drop = FALSE]
  }

  ## reorder, sort on score
  score.sign <- sign(median(df$score, na.rm = TRUE))
  df <- df[order(-df$score * score.sign), ]

  if (plot) {
    cols <- c("moduleMembership", "traitSignificance")
    if (showlogFC) cols <- c(cols, "foldChange", "centrality")
    cols <- intersect(cols, colnames(df))
    df1 <- df[, cols]
    col1 <- labels2colors(labels[rownames(df1)])
    if (!is.null(col)) col1 <- col
    pairs(df1, col = col1, oma = c(1, 1, 1, 1) * 1.8)
    if (is.null(main)) {
      main <- paste("Gene significance for module", module, "and trait", trait)
    }
    title(main, line = 3, cex.main = 1.15)
  }

  rownames(df) <- df$feature

  return(df)

}

#' Calculate compound significance scores per gene
#' @param wgcna WGCNA result object with stats.
#' @return Data frame of compound significance scores.
#' @export
calculateCompoundSignificance <- function(wgcna) {

  Q <- list()
  ww <- list(gx = wgcna)
  if (!is.null(wgcna$layers)) ww <- wgcna$layers    
  for (k in names(ww)) {
    stats <- ww[[k]]$stats
    m1 <- stats$moduleMembership
    t1 <- stats$traitSignificance
    f1 <- stats$foldChange
    c1 <- ww[[k]]$net$colors[rownames(m1)]
    rxs <- function(x,k=2) apply(x**k,1,max,na.rm=TRUE)^(1/k)
    Q1 <- data.frame(c1, rxs(m1,k=1), rxs(t1), rxs(f1))
    colnames(Q1) <- c("color","MM","TS","FC")
    Q1$score1 <- apply(Q1[,c(2,3)],1,prod)
    Q1$score2 <- apply(Q1[,c(2,4)],1,prod)    
    Q1 <- Q1[order(-Q1$score2),]
    Q[[k]] <- Q1
  }

  if (length(ww)==1) Q <- Q[[1]]

  return(Q)

}

#' Compute gene statistics with consensus colors per layer
#' @param cons Consensus WGCNA result object.
#' @return List of gene statistics per layer.
#' @export
computeConsensusGeneStats <- function(cons) {

  k <- names(cons$layers)[1]
  stats <- list()
  for(k in names(cons$layers)) {
    w <- cons$layers[[k]]
    colors <- cons$net$colors
    wMEs <- cons$net$multiMEs[[k]]$data
    wnet <- list( MEs = wMEs, colors = colors)
    stats[[k]] <- computeGeneStats(
      wnet, w$datExpr, w$datTraits, TOM=NULL)
  }
  
  return(stats)

}


#' Get consensus gene statistics across layers
#' @param cons Consensus WGCNA result object.
#' @param stats Precomputed consensus stats list.
#' @param trait Trait name to extract.
#' @param module Module name to filter or NULL.
#' @return List with consensus and full data frames.
#' @export
getConsensusGeneStats <- function(cons,
                                  stats,
                                  trait,
                                  module = NULL) {

  ## create extended color vector
  labels = paste0("ME",cons$net$colors)
  gstats <- list()
  for (k in names(stats)) {
    gstats[[k]] <- getGeneStats(
      wgcna = NULL,
      stats = stats[[k]],
      labels = labels,
      trait = trait,
      plot = FALSE,
      module = module,A
      col = NULL,
      main = NULL
    )
  }

  ## Align rows
  ff <- gstats[[1]]$feature
  for (k in names(gstats)) {
    ii <- match(ff, gstats[[k]]$feature)
    gstats[[k]] <- gstats[[k]][ii,]
  }

  ## Compute consensus statistics. Consensus statistics are computed
  ## as geometric mean of score variables, and/or maximum pvalue for
  ## p.value columns.
  xcols <- c(3,4,6,8)
  pcols <- c(10,5,7,9)
  pcols1 <- c(5,7,9)
  xcols <- c("score","moduleMembership","traitSignificance","foldChange")
  pcols <- c("scorePvalue","MMPvalue","TSPvalue","foldChangePvalue")
  pcols1 <- pcols[-1]

  for (i in 1:length(gstats)) {
    gstats[[i]][,'scorePvalue'] <- apply(gstats[[i]][,pcols1],1,max,na.rm=TRUE)
  }

  xc <- lapply(gstats, function(x) log(abs(x[,xcols])))
  xc <- exp(Reduce('+', xc) / length(xc))
  xp <- Reduce(pmax, lapply(gstats, function(x) x[,pcols]))
  df3 <- data.frame( gstats[[1]][,1:2], xc, xp)
  df3 <- df3[,colnames(gstats[[1]])]

  ## Determine consensus status. Feature is 'C' (concordant) if sign
  ## in all layers are equal and significant. 'D' (discordant) if sign
  ## if not equal in all layers but significant. 'N' is any is non-significant.
  sign.pos <- Reduce('*',lapply(gstats,function(g) sign(g$score) == 1))
  sign.neg <- Reduce('*',lapply(gstats,function(g) sign(g$score) == -1))
  allsig   <- Reduce('*',lapply(gstats,function(g) (g$scorePvalue) < 0.05))
  consensus <- c("D","C")[ 1 + 1*(sign.pos | sign.neg)]
  consensus[which(allsig==0)] <- 'N'
  cons.df <- data.frame(df3[,1:2], consensus, df3[,-c(1,2)])

  ## This creates the full stats matrix (all subgroups)
  df1 <- gstats[[1]][,c("feature","module")]
  df2 <- gstats[[1]][,0]
  cols <- colnames(gstats[[1]])[-c(1:2)]
  for (k in cols ) {
    xx <- sapply(gstats, function(g) g[,k])
    df2[[k]] <- I(xx)
  }
  df2 <- do.call(cbind, lapply(df2,unclass))
  newcols <- unlist(lapply(cols, function(k) paste0(k,'.',names(gstats))))
  colnames(df2) <- newcols
  full.df <- data.frame(df1, consensus=cons.df$consensus, df2)

  ii <- order(-cons.df$score * sign(mean(cons.df$score,na.rm=TRUE)))
  cons.df <- cons.df[ii, ]
  full.df <- full.df[ii, ]
  
  return(list(consensus = cons.df, full = full.df))
  
}
