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
      module = module,
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


## ========================================================================
## ======================== Fisher test based =============================
## ========================================================================
#' Perform Fisher's exact test on gene sets
#' This function performs Fisher's exact test on a set of genes within a given set of gene sets.
#' It returns a data frame containing the results of the test, including p-values, q-values,
#' odds ratios, and gene set overlaps.
#' @param genes A character vector containing the names of genes.
#' @param genesets A list of gene sets, where each element is a character vector representing a gene set.
#' @param background A character vector containing the names of genes in the background set.
#'                   Defaults to `NULL`, which means all genes are considered.
#' @param fdr The false discovery rate (FDR) threshold for multiple testing adjustment. Defaults to 0.05.
#' @param mc A logical value indicating whether to perform multiple testing adjustment
#'           using Monte Carlo simulation. Defaults to `TRUE`.
#' @param sort.by The statistic used to sort the results. Defaults to "p.value".
#' @param nmin The minimum number of genes required in a gene set for the test to be performed.
#' @param min.genes The minimum number of genes in a gene set to be considered. Defaults to 15.
#' @param max.genes The maximum number of genes in a gene set to be considered. Defaults to 500.
#' @param method The method used for computing p-values. Defaults to "fast.fisher".
#' @param check.background A logical value indicating whether to check the presence of genes in
#'                         the background set. Defaults to `TRUE`.
#' @param report.genes A logical value indicating whether to report genes in gene sets. Defaults to `TRUE`.
#' @param verbose A numeric value indicating the level of verbosity. Defaults to 1.
#' @export
#' @return A data frame containing the results of the Fisher's exact test. The data frame
#'         includes columns such as "p.value" (p-value), "q.value" (adjusted p-value),
#'         "odd.ratio" (odds ratio), "overlap" (gene set overlap), and optionally "genes" (common genes).
#' @export
gset.fisher <- function(genes,
                        genesets,
                        background = NULL,
                        fdr = 0.25,
                        mc = TRUE,
                        sort.by = "p.value",
                        nmin = 3,
                        min.genes = 15,
                        max.genes = 500,
                        method = "fast.fisher",
                        check.background = TRUE,
                        report.genes = FALSE,
                        no.pass = NA,
                        verbose = 1) {

  if (class(genesets) == "list") {
    gsnames <- names(genesets)
    res <- gset.fisherLIST(
      genes = genes,
      genesets = genesets,
      background = background,
      fdr = fdr,
      mc = mc,
      sort.by = sort.by,
      nmin = nmin,
      min.genes = min.genes,
      max.genes = max.genes,
      method = method,
      check.background = check.background,
      report.genes = report.genes,
      no.pass = no.pass,
      verbose = verbose
    )
  } else if (inherits(genesets, "Matrix")) {
    gsnames <- colnames(genesets)
    res <- gset.fastFET(
      genes,
      G = genesets,
      bg = background,
      report.genes = report.genes
    )
  } else {
    stop("[WGCNAplus::gset.fisher] FATAL ERROR")
  }

  if (nrow(res) > 0) {
    size.ok <- res$size >= min.genes & res$size <= max.genes
    jj <- which(res$q.value <= fdr & res$N >= nmin & size.ok)
    res <- res[jj, ]
  }

  if (nrow(res) > 0) {
    if (sort.by %in% colnames(res)) {
      order.sign <- ifelse(sort.by %in% c("p.value", "q.value"), +1, -1)
      res <- res[order(order.sign * res[, sort.by]), ]
    } else {
      gsnames <- intersect(gsnames, rownames(res))
      res <- res[gsnames, ]
    }
  }

  return(res)

}

#' @export
gset.fisherLIST <- function(genes,
                            genesets,
                            background = NULL,
                            fdr = 0.05,
                            mc = TRUE,
                            sort.by = "p.value",
                            nmin = 3,
                            min.genes = 15,
                            max.genes = 500,
                            method = "fast.fisher",
                            check.background = TRUE,
                            report.genes = FALSE,
                            no.pass = NA,
                            verbose = 1) {

  bgNULL <- FALSE
  if (is.null(background)) {
    message("[WGCNAplus::gset.fisherLIST] note: it is recommended to specify background")
    background <- unique(unlist(genesets))
    if (verbose > 0) {
      cat("setting background to ", length(background), "genes covered in genesets\n")
    }
    bgNULL <- TRUE
  }

  if (check.background && !bgNULL) {
    genes <- intersect(genes, background)
    genesets <- lapply(genesets, function(s) intersect(s, background))
  }

  if (!is.null(min.genes) && min.genes > 0) {
    genesets.len <- sapply(genesets, length)
    genesets <- genesets[order(-genesets.len)]
    if (sum(duplicated(names(genesets))) > 0) {
      if (verbose > 0) cat("warning: duplicated gene set names. taking largest.\n")
      genesets <- genesets[which(!duplicated(names(genesets)))]
    }
    genesets.len <- sapply(genesets, length)
    genesets <- genesets[genesets.len >= min.genes & genesets.len <= max.genes]
    if (verbose > 0) {
      cat("testing", length(genesets), "genesets with", length(genes),
        "genes (background", length(background), "genes)\n")
    }
    if (length(genesets) == 0) {
      cat("warning: no gene sets passed size filter\n")
      rr <- data.frame(p.value = NA, q.value = NA, odd.ratio = NA,
        N = NA, size = NA, n.overlap = NA, genes = NA)
      rownames(rr) <- NULL
      return(rr[0, ])
    }
  }

  ## odd ratio
  ## see http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf
  n.size <- sapply(genesets, length)
  bg0 <- setdiff(background, genes)
  nbackground0 <- length(background)
  nbackground1 <- length(bg0)

  ## this can become slow
  a <- unlist(lapply(genesets, function(x) sum(x %in% genes)))
  b <- (n.size - a)
  c <- unlist(lapply(genesets, function(x) sum(!(genes %in% x))))
  d <- (nbackground1 - b)
  odd.ratio <- (a / c) / (b / d) ## note: not exactly same as from fishertest

  ## report intersection genes (slow..) like leading edge list.
  gsgenes <- NULL
  if (report.genes) {
    gsgenes <- unlist(lapply(genesets, function(x) paste(sort(intersect(genes, x)), collapse = "|")))
  }

  ## compute fisher-test (should be one-sided?)
  test.fisher <- function(gs) {
    a0 <- table(background %in% gs, background %in% genes)
    if (NCOL(a0) == 1 || colSums(a0)[2] == 0) return(NA)
    stats::fisher.test(a0, alternative = "greater")$p.value
  }

  test.chisq <- function(gs) {
    a0 <- table(background %in% gs, background %in% genes)
    if (NCOL(a0) == 1 || colSums(a0)[2] == 0) return(NA)
    stats::chisq.test(a0)$p.value
  }

  pv <- rep(NA, length(genesets))
  names(pv) <- names(genesets)

  if (method == "fast.fisher") {
    pv <- rep(NA, length(a))
    ii <- 1:length(a)
    ii <- which((a + c) > 0)
    d1 <- d + 1 * (d == 0) ## hack to avoid crash...
    b1 <- b + 1 * (b == 0) ## hack to avoid crash...
    pv1 <- try(corpora::fisher.pval(a[ii], (a + b1)[ii], c[ii], (c + d1)[ii], alternative = "greater"),
      silent = TRUE)
    if (class(pv1) != "try-error") {
      pv[ii] <- pv1
    } else {
      message("[playbase::gset.fisher] fast.fisher failed. Testing with standard fisher.")
      method <- "fisher"
    }
  }

  if (method == "fisher") {
    if (mc) {
      pv <- unlist(lapply(genesets, test.fisher))
    } else {
      i <- 1
      for (i in 1:length(genesets)) {
        pv[i] <- test.fisher(genesets[[i]])
      }
    }
  } else if (method == "chisq") {
    if (mc) {
      pv <- unlist(lapply(genesets, test.chisq))
    } else {
      for (i in 1:length(genesets)) {
        pv[i] <- test.chisq(genesets[[i]])
      }
    }
  } else if (method == "fast.fisher") {
  } else {
    stop("unknown method")
  }

  if (any(is.na(pv))) pv[is.na(pv)] <- no.pass

  qv <- rep(NA, length(pv))
  qv <- stats::p.adjust(pv, method = "fdr")

  v1 <- as.character(paste0(a, "/", n.size))
  rr <- data.frame(p.value = pv, q.value = qv,
    odd.ratio = odd.ratio, N = a, size = n.size, overlap = v1)

  if (!is.null(gsgenes)) rr <- cbind(rr, genes = gsgenes)

  rownames(rr) <- names(genesets)

  if (nrow(rr) > 0) {
    jj <- which(rr$q.value <= fdr & n.size >= nmin)
    rr <- rr[jj, ]
    if (sort.by %in% c("pvalue", "p.value", "p")) {
      rr <- rr[order(rr$p.value), ]
    } else {
      rr <- rr[order(rr$odd.ratio, decreasing = TRUE), ]
    }
  }

  return(rr)
  
}


#' Calculate fast Fisher exact test.
#' @param genes Vector of significant genes
#' @param G    Sparse matrix containing gene sets
#' @param bg   Vector of genes as background
#' @param report.genes  Logical to report gene set genes in output
#' @export
gset.fastFET <- function(genes,
                         G,
                         bg,
                         report.genes = FALSE) {

  bgnull <- FALSE

  if (is.null(bg)) {
    message("[WGCNAplus::gset.fastFET] note: it is recommended to specify background")
    bg <- unique(c(genes, rownames(G)))
    bgnull <- TRUE
  }

  if (length(bg) > 1 && !bgnull) {
    genes <- intersect(genes, bg)
    G <- G[intersect(bg, rownames(G)), , drop = FALSE]
  }

  length.bg <- NULL

  if (length(bg) == 1 && is.integer(bg)) {
    length.bg <- as.integer(bg)
  } else if (length(bg) > 1) {
    length.bg <- length(bg)
  }

  if (is.null(length.bg)) stop("[WGCNAplus::gset.fastFET] error: invalid background. bg:", head(bg))
  if (length(genes) == 0) stop("[WGCNAplus::gset.fastFET] error: zero genes length")
  if (nrow(G) == 0) stop("[WGCNAplus::gset.fastFET] error: empty gene set matrix G")

  genes <- intersect(genes, rownames(G))
  gsize <- Matrix::colSums(G != 0)
  genes <- intersect(genes, rownames(G))
  a <- Matrix::colSums(G[genes, ] != 0)
  b <- length(genes) - a
  c <- gsize - a
  d <- length.bg - (a + b + c)
  pv <- corpora.fastFET(a, b, c, d)

  names(pv) <- colnames(G)
  odd.ratio <- (a / b) / (c / d)
  qv <- p.adjust(pv, method = "fdr")
  overlap <- paste0(a, "/", gsize)

  gsgenes <- NULL
  if (report.genes) {
    gsgenes <- apply(G[genes, ], 2, function(x) paste(sort(genes[which(x != 0)]), collapse = "|"))
  }

  df <- data.frame(p.value = pv, q.value = qv,
    odd.ratio = odd.ratio, N = a, size = gsize, overlap = overlap)

  if (!is.null(gsgenes)) df <- cbind(df, genes = gsgenes)

  return(df)

}


#' Wrapper superfast version of Fisher Exact Test from 'corpora' R package.
#' This is the fastest implementation currently available. Uses phyper inside.
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#' @export
corpora.fastFET <- function(a, b, c, d,
                            alternative = c("two.sided", "less", "greater")[3],
                            log.p = FALSE) {

  pv <- rep(NA, length(a))
  ii <- 1:length(a)
  ii <- which((a + b) > 0)
  d1 <- d + 1 * (d == 0) ## avoid crash...
  c1 <- c + 1 * (c == 0) ## avoid crash...

  k1 <- a[ii]
  n1 <- (a + c1)[ii]
  k2 <- b[ii]
  n2 <- (b + d1)[ii]

  .match.len <- function(vars, len = NULL, adjust = FALSE, check.numeric = TRUE,
                         envir = parent.frame()) {
    vecs <- setNames(lapply(vars, get, envir = envir), vars)
    ok <- sapply(vecs, is.numeric)
    if (check.numeric && any(!ok)) {
      stop(
        "argument(s) ", paste(vars[!ok], collapse = ", "),
        " must be numeric vector(s)"
      )
    }

    if (is.null(len)) len <- max(sapply(vecs, length))

    for (v in vars) {
      if (length(vecs[[v]]) == 1) {
        if (adjust) assign(v, rep(vecs[[v]], len), envir = envir)
      } else if (length(vecs[[v]]) != len) {
        stop(sprintf(
          "argument %s should be of length %d or a scalar (%s must have same length)",
          v, len, paste(vars, collapse = ", ")
        ))
      }
    }
    invisible(len)
  }

  alternative <- match.arg(alternative)
  l <- .match.len(c("k1", "n1", "k2", "n2"), adjust = TRUE)

  if (any(k1 < 0) || any(k1 > n1) || any(n1 <= 0)) {
    stop("k1 and n1 must be integers with 0 <= k1 <= n1")
  }

  if (any(k2 < 0) || any(k2 > n2) || any(n2 <= 0)) {
    stop("k2 and n2 must be integers with 0 <= k2 <= n2")
  }

  if (any(k1 + k2 <= 0)) {
    stop("either k1 or k2 must be non-zero")
  }

  k <- k1 + k2
  if (alternative == "two.sided") {
    if (log.p) {
      pval <- pmin(
        phyper(k1 - 1, n1, n2, k, lower.tail = FALSE, log.p = TRUE),
        phyper(k1, n1, n2, k, lower.tail = TRUE, log.p = TRUE)
      ) + log(2)
      pval <- pmin(pval, 0)
    } else {
      pval <- 2 * pmin(
        phyper(k1 - 1, n1, n2, k, lower.tail = FALSE),
        phyper(k1, n1, n2, k, lower.tail = TRUE)
      )
      pval <- pmax(0, pmin(1, pval))
    }
  } else if (alternative == "greater") {
    pval <- phyper(k1 - 1, n1, n2, k, lower.tail = FALSE, log.p = log.p)
  } else if (alternative == "less") {
    pval <- phyper(k1, n1, n2, k, lower.tail = TRUE, log.p = log.p)
  }

  return(pval)
  
}
