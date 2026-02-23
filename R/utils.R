#' Purple-grey-yellow color palette
#' @param n Number of colors to generate.
#' @return Character vector of hex colors.
#' @keywords internal
purpleGreyYellow <- function(n) {

  colorRampPalette(c("purple", "grey65", "yellow"))(n)

}

#' Convert correlation to blue-red colors
#' Converts correlation values [-1;1] to blue-white-red colors. Good
#' for creating color labels for labeledHeatmaps that expect colors.
#' NOTE: use WGCNA::numbers2colors???
#' @param R Numeric correlation matrix or vector.
#' @param a Exponent for nonlinear scaling.
#' @param f Color attenuation factor.
#' @return Character matrix or vector of colors.
#' @keywords internal
rho2bluered <- function(R, a = 1, f = 0.95) {

  BLUERED <- WGCNA::blueWhiteRed(100)

  if (a != 1) R <- sign(R) * abs(R)**a

  if (is.null(ncol(R))) {
    col <- BLUERED[1 + round(99 * (1 + R) / 2)]
  } else {
    col <- apply(R, 2, function(x) BLUERED[1 + round(99 * (1 + x) / 2)])
    dimnames(col) <- dimnames(R)
  }

  if (f < 1) {
    col <- apply(col, 2, adjustcolor, red.f = f, green.f = f, blue.f = f)
  }

  if (NCOL(col) == 1) col <- cbind(col)
  dimnames(col) <- dimnames(R)

  return(col)

}


#' Convert numeric labels to colors
#' Converts WGCNA labels (numeric or color) to colors.
#' @param colors Numeric or character label vector.
#' @param ... Additional arguments passed to WGCNA::labels2colors.
#' @return Character vector of color names.
#' @export
labels2colors <- function(colors, ...) {

  if (all(is.numeric(colors))) {
    colors <- WGCNA::labels2colors(colors, ...)
    return(colors)
  }

  stdColors <- c("grey", WGCNA::standardColors())

  if (all(colors %in% stdColors)) return(colors)

  icolors <- as.integer(factor(as.character(colors)))
  colors <- WGCNA::standardColors()[icolors]

  return(colors)

}


#' Filter module colors by KME
#' Filter color vector by minimum KME and mergeCutHeight. Set color of
#' features with KME smaller than minKME to grey (or 0) group. Merge
#' similar modules with (module) correlation larger than
#' (1-mergeCutHeight) together.
#' @param X Numeric expression matrix.
#' @param colors Module color assignments.
#' @param minKME Minimum KME threshold.
#' @param mergeCutHeight Cut height for merging modules.
#' @param minmodsize Minimum module size.
#' @param ntop Maximum features per module.
#' @return Filtered color assignment vector.
#' @export
filterColors <- function(X,
                         colors,
                         minKME = 0.3,
                         mergeCutHeight = 0.15,
                         minmodsize = 20,
                         ntop = -1) {

  sX <- X + 1e-8 * matrix(rnorm(length(X)), nrow(X), ncol(X))
  sX <- t(scale(t(sX)))

  ## get singular vectors and correct sign
  vv <- tapply(1:nrow(sX), colors, function(i) svd(sX[i, ], nv = 1)$v[, 1])
  mm <- tapply(1:nrow(sX), colors, function(i) colMeans(sX[i, ]))
  vv.sign <- mapply(function(a, b) sign(cor(a, b)), mm, vv)
  vv <- mapply(function(a, b) a * b, vv, vv.sign, SIMPLIFY = FALSE)

  kme <- rep(NA, nrow(X))
  names(kme) <- rownames(X)
  names(colors) <- rownames(X)

  grey.val <- NULL
  is.color <- mean(colors %in% WGCNA::standardColors(435)) > 0.8
  if (is.numeric(colors)) {
    colors <- as.integer(colors)
    grey.val <- 0
  } else {
    colors <- as.character(colors)
    grey.val <- "---"
    if (is.color) grey.val <- "grey"
  }

  names(colors) <- rownames(X)
  new.colors <- colors

  if (minKME > 0) {
    for (i in 1:length(vv)) {
      ii <- which(colors == names(vv)[i])
      r <- cor(t(X[ii, ]), vv[[i]])[, 1]
      jj <- ii[which(r < minKME)]
      if (length(jj)) new.colors[jj] <- NA
      kme[ii] <- r
    }
    new.colors[is.na(new.colors)] <- grey.val
  }

  ## merge groups
  if (mergeCutHeight > 0) {
    mx <- rowmean(X, new.colors)
    rr <- cor(t(mx))
    diag(rr) <- 0
    merge.idx <- which(rr > (1 - mergeCutHeight), arr.ind = TRUE)
    if (nrow(merge.idx) > 0) {
      for (i in 1:nrow(merge.idx)) {
        aa <- rownames(rr)[merge.idx[i, ]]
        jj <- which(new.colors %in% aa)
        max.color <- names(which.max(table(new.colors[jj])))
        new.colors[jj] <- max.color
      }
    }
  }

  ## remove small groups
  modsize <- table(new.colors)
  if (min(modsize) < minmodsize) {
    small.mod <- names(which(modsize < minmodsize))
    sel <- which(new.colors %in% small.mod)
    new.colors[sel] <- NA
  }

  ## Filter by KME score
  if (ntop > 0) {
    keep <- tapply(names(kme), new.colors, function(i) head(names(sort(-kme[i])), ntop))
    keep <- unlist(keep)
    not.keep <- setdiff(names(kme), keep)
    if (length(not.keep)) new.colors[not.keep] <- NA
  }

  new.colors[which(is.na(new.colors))] <- grey.val

  return(new.colors)

}

#' TOM-based hierarchical clustering
#' Wrapper to hclust from matrix using default WGCNA parameters.
#' @param X Numeric expression matrix.
#' @param power Soft-thresholding power.
#' @return An hclust dendrogram object.
#' @keywords internal
tomclust <- function(X, power = 6) {

  A <- WGCNA::adjacency(t(X), power = power, type = "signed")
  TOM <- fastTOMsimilarity(A, tomtype = "signed", lowrank = 40)
  hc <- hclust(as.dist(1 - TOM), method = "average")
  return(hc)

}


#' Validate dendrogram heights across powers
#' @param datExpr Numeric expression data matrix.
#' @param n Number of top-variance features to sample.
#' @param powers Numeric vector of powers to test.
#' @param maxpower Maximum power to evaluate.
#' @return List with quantiles, IQR, and optimal power.
#' @keywords internal
checkDendroHeights <- function(datExpr,
                               n = 200,
                               powers = NULL,
                               maxpower = 20) {

  ii <- 1:ncol(datExpr)
  if (n < ncol(datExpr)) {
    ii <- head(order(-matrixStats::colSds(datExpr)), n)
  }

  tX <- datExpr[, ii]
  ht <- list()
  p <- 9
  p <- 24

  if (is.null(powers)) {
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    if (maxpower > 20) {
      powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
    }
  }

  for (i in 1:length(powers)) {
    A <- WGCNA::adjacency(tX, power = powers[i], type = "signed")
    TOM <- fastTOMsimilarity(A, tomtype = "signed", lowrank = 40)
    hc <- hclust(as.dist(1 - TOM), method = "average")
    ht[[i]] <- hc$height
  }

  names(ht) <- paste0("p=", powers)
  S <- sapply(ht, quantile, probs = c(0.25, 0.5, 0.75))
  iqr <- (S[3, ] - S[1, ])
  optK <- powers[which.max(iqr)]

  return(list(quantiles = S, IQR = iqr, optK = optK))
  
}

#' Plot soft threshold power analysis
#' @param datExpr Numeric expression data matrix.
#' @param networktype Network type (e.g. "signed").
#' @param cex Text size for plot labels.
#' @param maxpower Maximum power to evaluate.
#' @param nmax Maximum features to subsample.
#' @param plots Character vector of plot types.
#' @param main Plot main title.
#' @param RsquaredCut R-squared cutoff for fit.
#' @param setPar Whether to set par layout.
#' @return Invisibly NULL. Plots are drawn.
#' @export
plotPowerAnalysis <- function(datExpr,
                              networktype = "signed",
                              cex = 1,
                              maxpower = 20,
                              nmax = 2000,
                              plots = c(
                                "sft.modelfit", "mean.k",
                                "dendro.IQR"
                              ),
                              main = NULL,
                              RsquaredCut = 0.85,
                              setPar = TRUE) {

  RsquaredCut <- RsquaredCut[1]

  ## Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  if (maxpower > 20) {
    powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
  }

  ## subsample for speed
  if (ncol(datExpr) > nmax && nmax > 0) {
    ii <- sample(1:ncol(datExpr), nmax)
    datExpr <- datExpr[, ii]
  }

  ## Call the network topology analysis function
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers,
    RsquaredCut = RsquaredCut, networkType = networktype, verbose = 0)

  if (setPar) {
    np <- length(plots)
    nc <- ceiling(sqrt(np))
    par(mfrow = c(nc, nc), mar = c(3.3, 3.5, 1, 1), mgp = c(2, 0.9, 0))
    par(mfrow = c(1, np), mar = c(3.8, 3.8, 1, 1), mgp = c(2.4, 0.95, 0))
  }

  ## Plot results
  if ("sft.modelfit" %in% plots) {
    ## Scale-free topology fit index as function of soft-thresholding power
    y <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    base::plot(
      x = sft$fitIndices[, 1],
      y = y,
      ylim = c(min(y), 1),
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "SFT model fit (signed R^2)",
      main = main
    )
    abline(h = 0, col = "black", lty = 3)
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      labels = powers, cex = cex, col = "red"
    )
    ## this line corresponds to using an R^2 cut-off of h
    abline(h = RsquaredCut, col = "red", lty = 2)
  }

  ## Mean connectivity as a function of the soft-thresholding power
  if ("mean.k" %in% plots) {
    base::plot(sft$fitIndices[, "Power"], sft$fitIndices[, "mean.k."],
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "Mean connectivity",
      main = main
    )
    text(sft$fitIndices[, "Power"], sft$fitIndices[, "mean.k."],
      labels = powers, cex = cex, col = "red")
  }

  ht <- NULL
  if ("dendro.IQR" %in% plots) {
    ht <- checkDendroHeights(datExpr, n = 200, powers = powers)
    base::plot(
      sft$fitIndices[, 1], ht$IQR,
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "Dendrogram height IQR",
      main = main
    )
    text(sft$fitIndices[, 1], ht$IQR, labels = powers, cex = cex, col = "red")
  }

}

#' Multi-dataset power analysis plot
#' @param exprList Named list of expression matrices.
#' @param cex Text size for plot labels.
#' @param maxpower Maximum power to evaluate.
#' @param nmax Maximum features to subsample.
#' @param networktype Network type (e.g. "signed").
#' @param plots Character vector of plot types.
#' @param main Plot main title.
#' @param cex.legend Legend text size.
#' @param RsquaredCut R-squared cutoff for fit.
#' @param setPar Whether to set par layout.
#' @return Invisibly NULL. Plots are drawn.
#' @export
plotPowerAnalysis_multi <- function(exprList,
                                    cex = 1, maxpower = 20,
                                    nmax = 2000,
                                    networktype = "signed",
                                    plots = c(
                                      "sft.modelfit", "mean.k",
                                      "dendro.IQR"
                                    ),
                                    main = NULL,
                                    cex.legend = 1,
                                    RsquaredCut = 0.85,
                                    setPar = TRUE) {

  RsquaredCut <- RsquaredCut[1]

  ## Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  if (maxpower > 20) {
    powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
  }

  ## process each list
  sft <- list()
  for (i in 1:length(exprList)) {
    datExpr <- Matrix::t(exprList[[i]])

    ## subsample for speed
    if (nmax > 0 && nrow(exprList[[i]]) > nmax) {
      ii <- sample(1:ncol(datExpr), nmax)
      datExpr <- datExpr[, ii]
    }

    ## Call the network topology analysis function
    k <- names(exprList)[i]
    sft[[k]] <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      RsquaredCut = RsquaredCut,
      networkType = networktype,
      verbose = 0
    )
  }

  if (setPar) {
    np <- length(plots)
    par(mfrow = c(1, np), mar = c(3.8, 4.5, 3, 1), mgp = c(2.6, 0.95, 0))
  }

  ## Plot results:
  if ("sft.modelfit" %in% plots) {
    ## Scale-free topology fit index as function of soft-thresholding power
    Y <- c()
    for (i in 1:length(sft)) {
      y1 <- -sign(sft[[i]]$fitIndices[, "slope"]) * sft[[i]]$fitIndices[, "SFT.R.sq"]
      Y <- cbind(Y, y1)
    }
    colnames(Y) <- names(sft)
    x <- sft[[1]]$fitIndices[, "Power"]
    Y <- pmax(Y, 0)
    matplot(
      x = x,
      y = Y,
      ylim = c(0, 1),
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "SFT model fit (signed R^2)",
      main = main
    )

    for (i in 1:ncol(Y)) {
      text(powers, Y[, i], labels = "\u2588", cex = cex, col = "white")
      text(powers, Y[, i], labels = powers, cex = cex, col = 1 + i)
    }

    ## this line corresponds to using an R^2 cut-off of h
    abline(h = RsquaredCut, col = "grey10", lty = 2)
    legend("bottomright", legend = colnames(Y), fill = 2:10,
      cex = cex.legend, y.intersp = 0.9)
    title("SFT model fit")
  }

  ## Mean connectivity as a function of the soft-thresholding power
  if ("mean.k" %in% plots) {
    Y <- sapply(sft, function(s) s$fitIndices[, "mean.k."])
    matplot(
      powers,
      Y,
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "Mean connectivity",
      main = main
    )
    for (i in 1:ncol(Y)) {
      text(powers, Y[, i], labels = "\u2588", cex = cex, col = "white")
      text(powers, Y[, i], labels = powers, cex = cex, col = 1 + i)
    }
    legend("topright",
      legend = colnames(Y), fill = 2:10,
      cex = cex.legend, y.intersp = 0.9
    )
    title("Mean connectivity")
  }

  ht <- NULL
  if ("dendro.IQR" %in% plots) {
    ht <- list()
    for (i in 1:length(exprList)) {
      ht[[i]] <- checkDendroHeights(
        Matrix::t(exprList[[i]]),
        n = 200, powers = powers
      )
    }
    Y <- sapply(ht, function(h) h$IQR)
    matplot(
      powers,
      Y,
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "Dendrogram height IQR",
      main = main
    )

    for (i in 1:ncol(Y)) {
      text(powers, Y[, i], labels = "\u2588", cex = cex, col = "white")
      text(powers, Y[, i], labels = powers, cex = cex, col = 1 + i)
    }

    legend("bottomright", legend = names(exprList),
      fill = 2:10, cex = cex.legend, y.intersp = 0.9)
    title("Dendrogram IQR")
  }

}


#' Pick soft thresholding power
#' Better (?) method to pick soft threshold (aka power).
#' @param datExpr Numeric expression data matrix.
#' @param sft Pre-computed soft threshold result.
#' @param rcut R-squared cutoff for model fit.
#' @param method Selection method: "sft" or "iqr".
#' @param nmax Maximum features to subsample.
#' @param powers Numeric vector of powers to test.
#' @param verbose Verbosity level.
#' @return Integer optimal soft-thresholding power.
#' @export
pickSoftThreshold <- function(datExpr,
                              sft = NULL,
                              rcut = 0.85,
                              method = c("sft", "iqr")[1],
                              nmax = -1,
                              powers = NULL,
                              verbose = 1) {

  if (is.null(powers)) {
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  }

  ## subsample for speed
  if (ncol(datExpr) > nmax && nmax > 0) {
    ii <- sample(1:ncol(datExpr), nmax)
    datExpr <- datExpr[, ii]
  }

  if (is.null(sft)) {
    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers,
      networkType = "signed", verbose = verbose)
  }

  optPower <- NULL
  if (method == "sft") {
    ## Pick power according to scale-free (SFT) parameter
    sqr <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    if (max(sqr, na.rm = TRUE) >= rcut) {
      optPower <- min(powers[which(sqr >= rcut)])
    } else {
      ## remove initial value that are possible negative
      if (sqr[1] < 0.05) {
        for (i in 1:length(sqr)) sqr[i] <- ifelse(sqr[i] < 0.05, NA, sqr[i])
      }
      ds <- 0.5 * median(abs(diff(sqr)), na.rm = TRUE) ## small step
      if (any(diff(sqr) < -ds, na.rm = TRUE)) {
        i <- min(which(diff(sqr) < -ds)) + 1
        sqr[i:length(sqr)] <- NA
      }
      optPower <- powers[which.max(sqr)]
    }
  } else if (method == "iqr") {
    ht <- checkDendroHeights(datExpr, n = 200, powers = powers)
    optPower <- powers[which.max(ht$IQR)]
  } else {
    stop("[pickSoftThreshold] invalid method = ", method)
  }

  if (verbose > 0) {
    message("[pickSoftThreshold] sft$powerEstimate = ", sft$powerEstimate)
    message("[pickSoftThreshold] optPower = ", optPower)
  }

  return(optPower)

}


#' Scale TOM matrices to equal quantiles
#' Scale a list of TOM matrices so that the quantiles (default p=0.95)
#' are equal after scaling with respect to the first TOM matrix.
#' @param TOMs List of TOM matrices.
#' @param scaleP Reference quantile for scaling.
#' @return List of scaled TOM matrices.
#' @keywords internal
scaleTOMs <- function(TOMs, scaleP = 0.95) {

  nGenes <- nrow(TOMs[[1]])

  nSets <- length(TOMs)

  # Sample sufficiently large number of TOM entries
  nSamples <- as.integer(1 / (1 - scaleP) * 1000)

  # Choose the sampled TOM entries
  scaleSample <- sample(nGenes * (nGenes - 1) / 2, size = nSamples)

  TOMScalingSamples <- list()

  # These are TOM values at reference percentile
  scaleQuant <- rep(1, nSets)

  # Scaling powers to equalize reference TOM values
  scalePowers <- rep(1, nSets)

  # Loop over sets
  set <- 1
  for (set in 1:nSets) {
    # Select the sampled TOM entries
    tval <- as.dist(TOMs[[set]])[scaleSample]
    # Calculate the 95th percentile
    scaleQuant[set] <- quantile(tval, probs = scaleP, type = 8)
    TOMScalingSamples[[set]] <- tval
    # Scale the TOM
    if (set > 1) {
      scalePowers[set] <- log(scaleQuant[1]) / log(scaleQuant[set])
      TOMs[[set]] <- TOMs[[set]]^scalePowers[set]
    }
  }

  return(TOMs)

}

#' Get module-trait correlation data
#' @param wgcna A WGCNA result object.
#' @return Numeric module-trait correlation matrix.
#' @keywords internal
get_modTraits <- function(wgcna) {

  if(!is.null(wgcna$modTraits)) {
    M <- wgcna$modTraits
  } else {
    M <- cor( wgcna$net$MEs, wgcna$datTraits, use="pairwise")
  }

  M[is.na(M)] <- 0

  return(M)

}


## =========================================================================
## Functions inlined from playbase to minimize external dependencies
## =========================================================================

#' Check if value is a Date
#' @param x Value to test.
#' @return Logical TRUE if Date parseable.
#' @keywords internal
is.Date <- function(x) {

  if (!all(is.na(as.Date(
    as.character(x),
    format = c("%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d", "%Y-%m-%d")
  )))) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}

#' Collapse expanded trait matrix
#'
#' @param Y Expanded trait matrix with "=" columns.
#' @return Collapsed matrix with category columns.
#' @keywords internal
collapseTraitMatrix <- function(Y) {

  if (sum(grepl("=", colnames(Y))) < 2) {
    return(Y)
  }

  is.cat <- grepl("=", colnames(Y))
  M <- Y[, which(!is.cat), drop = FALSE]
  categories <- unique(sub("=.*", "", colnames(Y)[which(is.cat)]))
  y <- categories[1]

  for (y in categories) {
    ii <- which(sub("=.*", "", colnames(Y)) == y)
    Y1 <- Y[, ii]
    colnames(Y1) <- sub(".*=", "", colnames(Y1))
    m1 <- colnames(Y1)[max.col(Y1)]
    M <- cbind(M, m1)
    colnames(M)[ncol(M)] <- y
  }

  return(M)

}

#' Log counts per million normalization
#' @param counts Count matrix (dense or sparse).
#' @param total Library size normalization target.
#' @param prior Pseudocount before log transform.
#' @param log Whether to log2-transform.
#' @return Normalized expression matrix.
#' @keywords internal
logCPM <- function(counts,
                   total = 1e6,
                   prior = 1,
                   log = TRUE) {

  if (is.null(total)) {
    total0 <- mean(Matrix::colSums(counts, na.rm = TRUE))
    total <- ifelse(total0 < 1e6, total0, 1e6)
    message("[logCPM] setting column sums to = ", round(total, 2))
  }

  if (any(class(counts) == "dgCMatrix")) {
    cpm <- counts
    cpm[is.na(cpm)] <- 0
    cpm@x <- total * cpm@x / rep.int(Matrix::colSums(cpm), diff(cpm@p))
    if (log) cpm@x <- log2(prior + cpm@x)
    return(cpm)
  } else {
    totcounts <- Matrix::colSums(counts, na.rm = TRUE)
    cpm <- sweep(counts, 2, totcounts, FUN = "/") * total
    if (log) cpm <- log2(prior + cpm)
    return(cpm)
  }

}

#' Make contrast matrix from label matrix
#' @param lab.matrix Label matrix with "_vs_" column names.
#' @return Numeric contrast matrix.
#' @keywords internal
makeContrastsFromLabelMatrix <- function(lab.matrix) {

  if (!all(grepl("_vs_", colnames(lab.matrix)))) {
    stop("[makeContrastsFromLabelMatrix] FATAL:: all contrast names must include _vs_")
  }

  ct.names <- colnames(lab.matrix)
  main.grp <- sapply(strsplit(ct.names, split = "_vs_"), "[", 1)
  ctrl.grp <- sapply(strsplit(ct.names, split = "_vs_"), "[", 2)
  main.grp <- sub(".*:", "", main.grp)
  ctrl.grp <- sub("@.*", "", ctrl.grp)

  contr.mat <- matrix(0, nrow(lab.matrix), ncol(lab.matrix))
  rownames(contr.mat) <- rownames(lab.matrix)
  colnames(contr.mat) <- colnames(lab.matrix)
  for (i in 1:ncol(lab.matrix)) {
    lab1 <- trimws(lab.matrix[, i])
    lab1x <- setdiff(lab1, c(NA, ""))
    grps <- c(main.grp[i], ctrl.grp[i])
    if (all(lab1x %in% grps)) {
      j1 <- which(lab1 == main.grp[i])
      j0 <- which(lab1 == ctrl.grp[i])
    } else {
      j1 <- grep(paste0("^", toupper(main.grp[i])), toupper(lab1))
      j0 <- grep(paste0("^", toupper(ctrl.grp[i])), toupper(lab1))
    }
    contr.mat[j1, i] <- +1 / length(j1)
    contr.mat[j0, i] <- -1 / length(j0)
  }

  return(contr.mat)

}

#' Column-bind sparse matrices
#' @param m1 First sparse matrix.
#' @param m2 Second sparse matrix.
#' @return Combined sparse matrix.
#' @keywords internal
cbind_sparse_matrix <- function(m1, m2) {

  gene_vector <- unique(c(rownames(m1), rownames(m2)))

  if (!all(gene_vector %in% rownames(m1))) {
    missing_genes_m1 <- setdiff(gene_vector, rownames(m1))
    zero_rows_m1 <- Matrix::Matrix(0, nrow = length(missing_genes_m1), ncol = ncol(m1), sparse = TRUE)
    rownames(zero_rows_m1) <- missing_genes_m1
    m1 <- rbind(m1, zero_rows_m1)
  }

  if (!all(gene_vector %in% rownames(m2))) {
    missing_genes_m2 <- setdiff(gene_vector, rownames(m2))
    zero_rows_m2 <- Matrix::Matrix(0, nrow = length(missing_genes_m2), ncol = ncol(m2), sparse = TRUE)
    rownames(zero_rows_m2) <- missing_genes_m2
    m2 <- rbind(m2, zero_rows_m2)
  }

  m1 <- m1[gene_vector, , drop = FALSE]
  m2 <- m2[gene_vector, , drop = FALSE]
  combined_gmt <- cbind(m1, m2)
  combined_gmt <- combined_gmt[, order(-Matrix::colSums(combined_gmt != 0))]
  combined_gmt <- combined_gmt[, !duplicated(colnames(combined_gmt))]

  return(combined_gmt)

}

#' Merge two sparse matrices
#' @param m1 First sparse matrix (or NULL).
#' @param m2 Second sparse matrix (or NULL).
#' @param margin Unused, reserved for future.
#' @param verbose Verbosity level.
#' @return Combined sparse matrix.
#' @keywords internal
merge_sparse_matrix <- function(m1,
                                m2,
                                margin = NULL,
                                verbose = 1) {

  if (is.null(m1)) return(m2)

  if (is.null(m2)) return(m1)

  return(cbind_sparse_matrix(m1 = m1, m2 = m2))

}

#' Strip MOFA feature prefix
#' @param xx Character vector, matrix, or list.
#' @return Input with prefixes removed.
#' @keywords internal
mofa.strip_prefix <- function(xx) {

  if (class(xx) == "character") {
    xx <- sub("[A-Za-z0-9]+:", "", xx)
    return(xx)
  }

  if (class(xx) == "matrix") {
    rownames(xx) <- sub("[A-Za-z0-9]+:", "", rownames(xx))
    return(xx)
  }

  if (class(xx) %in% c("list", "array") || is.list(xx)) {
    for (i in 1:length(xx)) {
      dt <- paste0("^", names(xx)[i], ":")
      if (is.null(dim(xx[[i]]))) {
        names(xx[[i]]) <- sub(dt, "", names(xx[[i]]))
      } else {
        rownames(xx[[i]]) <- sub(dt, "", rownames(xx[[i]]))
      }
    }
    return(xx)
  }

  return(xx)

}

#' Add MOFA feature prefix
#' @param xx Named list of matrices or vectors.
#' @return Input with layer-name prefixes added.
#' @keywords internal
mofa.prefix <- function(xx) {

  xx <- mofa.strip_prefix(xx)

  for (i in 1:length(xx)) {
    dt <- paste0(names(xx)[i], ":")
    if (is.null(dim(xx[[i]]))) {
      names(xx[[i]]) <- paste0(dt, names(xx[[i]]))
    } else {
      rownames(xx[[i]]) <- paste0(dt, rownames(xx[[i]]))
    }
  }

  return(xx)

}

#' Merge MOFA data layers by row
#' @param xx Named list of data matrices.
#' @return Combined matrix with prefixed rownames.
#' @keywords internal
mofa.merge_data <- function(xx) { do.call(rbind, mofa.prefix(xx)) }

#' Split MOFA data by prefix
#' @param X Matrix with prefixed rownames.
#' @param keep.prefix Whether to keep prefix in names.
#' @return Named list of matrices per layer.
#' @keywords internal
mofa.split_data <- function(X, keep.prefix = FALSE) {

  if (!all(grepl("[:]|SOURCE|SINK", rownames(X)))) {
    rownames(X) <- paste0("x:", rownames(X))
  }

  dtype <- sub(":.*", "", rownames(X))
  xx <- tapply(1:nrow(X), dtype, function(i) X[i, , drop = FALSE])
  if (!keep.prefix) xx <- mofa.strip_prefix(xx)

  return(xx)

}

#' Select top SD features per MOFA layer
#' @param xdata Matrix or list of matrices.
#' @param ntop Number of top features to keep.
#' @return Filtered data with top-variance features.
#' @keywords internal
mofa.topSD <- function(xdata, ntop) {

  if (is.list(xdata)) {
    res <- lapply(xdata, function(x) {
      sdx <- matrixStats::rowSds(x, na.rm = TRUE)
      head(x[order(-sdx), , drop = FALSE], ntop)
    })
  } else if (is.matrix(xdata)) {
    if (all(grepl(":", rownames(xdata)))) {
      xdata <- mofa.split_data(xdata)
      res <- lapply(xdata, function(x) {
        sdx <- matrixStats::rowSds(x, na.rm = TRUE)
        head(x[order(-sdx), , drop = FALSE], ntop)
      })
      res <- mofa.merge_data(res)
    } else {
      sdx <- matrixStats::rowSds(xdata, na.rm = TRUE)
      res <- head(xdata[order(-sdx), , drop = FALSE], ntop)
    }
  } else {
    message("[mofa.topSD] WARNING: could not detect type")
    res <- xdata
  }

  return(res)

}

#' Compute row means by group
#' @param X Numeric matrix (dense or sparse).
#' @param group Group labels for rows.
#' @param reorder Whether to preserve group order.
#' @return Matrix of group-level row means.
#' @keywords internal
rowmean <- function(X, group = rownames(X), reorder = TRUE) {

  if (nrow(X) == 1) return(X)

  ngroup <- length(unique(group))
  if (ngroup == 1) {
    newX <- matrix(Matrix::colMeans(X, na.rm = TRUE), nrow = 1, ncol = ncol(X))
    dimnames(newX) <- list(group[1], colnames(X))
    return(newX)
  }

  if (is.matrix(X) || any(class(X) %in% c("matrix"))) {
    sumX <- base::rowsum(as.matrix(X), group, na.rm = TRUE)
    nX <- base::rowsum(1 * (!is.na(as.matrix(X))), group)
    newX <- sumX / nX
  } else if (sum(is.na(X)) == 0) {
    group_mat <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + group))
    rownames(group_mat) <- sub("^group", "", rownames(group_mat))
    group_mat <- group_mat / Matrix::rowSums(group_mat)
    newX <- group_mat %*% X
  } else {
    group_mat <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + group))
    rownames(group_mat) <- sub("^group", "", rownames(group_mat))
    X0 <- X
    X0[is.na(X0)] <- 0
    nc <- group_mat %*% (!is.na(X))
    newX <- (group_mat %*% X0) / nc
    newX <- Matrix::Matrix(newX, sparse = TRUE)
  }

  if (reorder) {
    ii <- match(unique(group), rownames(newX))
    newX <- newX[ii, , drop = FALSE]
  }

  return(newX)

}

#' Extract MOFA prefix from names
#' @param x Character vector, matrix, or data.frame.
#' @return Character vector of extracted prefixes.
#' @keywords internal
mofa.get_prefix <- function(x) {

  if (class(x) %in% c("matrix", "data.frame") || !is.null(dim(x))) {
    x <- rownames(x)
  }

  ifelse(grepl(":", x), sub(":.*", "", x), "")

}

#' Merge MOFA data with flexible options
#' @param xdata Named list of data matrices.
#' @param merge.rows Row merge strategy: "prefix", "union", or "intersect".
#' @param merge.cols Column merge strategy: "prefix", "union", or "intersect".
#' @return Merged numeric matrix.
#' @keywords internal
mofa.merge_data2 <- function(xdata, merge.rows="prefix", merge.cols="union") {

  n1 <- length(Reduce(intersect,lapply(xdata,rownames)))
  n2 <- length(Reduce(intersect,lapply(xdata,colnames)))  
  rdim <- sapply(xdata,nrow)
  cdim <- sapply(xdata,ncol)

  if (n1 < min(rdim) && merge.rows!="prefix") {
    message("WARNING: rows do not match")
  }

  if (n2 < min(cdim) && merge.cols!="prefix") {
    message("WARNING: columns do not match")
  }

  prefix.rows <- (merge.rows=="prefix")
  prefix.cols <- (merge.cols=="prefix")

  if (prefix.cols) {
    ## prefix the column names. i.e. different datasets.
    for(i in 1:length(xdata)) {
      nn <- sub("^[A-Za-z]+:","",colnames(xdata[[i]]))
      colnames(xdata[[i]]) <- paste0(names(xdata)[i],":",nn)
    }
    merge.cols <- "union"    
  }

  if (prefix.rows) {
    ## if columns overlap (i.e. same samples), prefix feature names.
    for (i in 1:length(xdata)) {
      nn <- sub("^[A-Za-z]+:","",rownames(xdata[[i]]))
      rownames(xdata[[i]]) <- paste0(names(xdata)[i],":",nn)
    }
    merge.rows <- "union"
  }

  if (merge.rows == "intersect") {
    allfeatures <- Reduce(intersect,lapply(xdata,rownames))
  } else {
    allfeatures <- unique(unlist(lapply(xdata, rownames)))
  }

  if (merge.cols == "intersect") {
    allsamples  <- Reduce(intersect,lapply(xdata,colnames))
  } else {
    allsamples  <- unique(unlist(lapply(xdata, colnames)))
  }

  D  <- matrix(0, length(allfeatures), length(allsamples))
  nn <- matrix(0, length(allfeatures), length(allsamples))
  rownames(D) <- allfeatures
  colnames(D) <- allsamples
  for (i in 1:length(xdata)) {
    A <- xdata[[i]]
    ii <- match(rownames(D), rownames(A))
    jj <- match(colnames(D), colnames(A))
    A1 <- A[ii, jj]
    nn <- nn + !is.na(A1) * 1
    A1[is.na(A1)] <- 0
    D <- D + A1
  }

  D <- D / nn
  D[which(nn == 0)] <- NA
  rownames(D) <- allfeatures
  colnames(D) <- allsamples
  return(D)

}

#' Convert probe IDs to gene symbols
#' @param probes Character vector of probe IDs.
#' @param annot_table Annotation data.frame with mappings.
#' @param query Target column name to return.
#' @param key Column name to match probes.
#' @param fill_na Fill missing with probe IDs.
#' @param add_datatype Prepend data type prefix.
#' @return Character vector of mapped symbols.
#' @keywords internal
probe2symbol <- function(probes,
                         annot_table,
                         query = "symbol",
                         key = NULL,
                         fill_na = FALSE,
                         add_datatype = FALSE) {

  annot_table <- cbind(rownames = rownames(annot_table), annot_table)
  id.cols <- intersect(c("feature", "gene_name", "symbol"), colnames(annot_table))

  if (length(id.cols) > 0) {
    stripped_annot <- apply(annot_table[, id.cols, drop = FALSE], 2, function(a) sub("^[A-Za-z]+:", "", a))
    annot_table <- cbind(annot_table, stripped_annot)
  }

  probes1 <- setdiff(probes, c(NA, ""))

  if (is.null(key) || !key %in% colnames(annot_table)) {
    key <- which.max(apply(annot_table, 2, function(a) sum(probes1 %in% a)))
  }

  if (is.null(key)) {
    message("[probe2symbol] FATAL. could not get key column.")
    return(NULL)
  }

  query <- head(intersect(query, colnames(annot_table)), 1)

  if (length(query) == 0) {
    message("ERROR. no symbol column.")
    return(NULL)
  }

  if (query == "symbol" && !"symbol" %in% colnames(annot_table) &&
    "gene_name" %in% colnames(annot_table)) {
    query <- "gene_name"
  }

  ii <- match(probes, annot_table[, key])
  query_col <- annot_table[ii, query]

  if (fill_na) {
    query_col <- ifelse(query_col == "" | is.na(query_col),
      yes = probes, no = query_col)
  }

  if (add_datatype && "data_type" %in% colnames(annot_table)) {
    datatype_col <- annot_table[ii, "data_type"]
    has_datatype <- !is.na(datatype_col) & datatype_col != ""
    already_has_prefix <- startsWith(query_col, paste0(datatype_col, ":"))
    should_add <- has_datatype & !already_has_prefix
    query_col <- ifelse(should_add, paste0(datatype_col, ":", query_col), query_col)
  }

  return(query_col)

}

#' Rename matrix rows by annotation mapping
#' @param counts Matrix, data.frame, or named vector.
#' @param annot_table Annotation data.frame with mappings.
#' @param new_id Target identifier column name.
#' @param na.rm Remove rows with NA names.
#' @param unique Average duplicate row names.
#' @param keep.prefix Preserve MOFA-style prefix.
#' @return Input with renamed rows.
#' @keywords internal
rename_by2 <- function(counts,
                       annot_table,
                       new_id = "symbol",
                       na.rm = TRUE,
                       unique = TRUE,
                       keep.prefix = FALSE) {

  annot_table$rownames <- rownames(annot_table)
  annot_table$rownames2 <- sub("^[A-Za-z]+:", "", rownames(annot_table))

  if (is.matrix(counts) || inherits(counts, "Matrix") ||
    is.data.frame(counts) || !is.null(dim(counts))) {
    type <- "matrix"
    probes <- rownames(counts)
  } else {
    type <- "vector"
    probes <- names(counts)
  }

  probe_match <- apply(annot_table, 2, function(x) sum(probes %in% x))

  if (max(probe_match, na.rm = TRUE) == 0) return(counts)

  if (type == "vector") counts <- cbind(counts)

  from_id <- names(which.max(probe_match))
  if (new_id == "symbol" && !"symbol" %in% colnames(annot_table) &&
    "gene_name" %in% colnames(annot_table)) {
    new_id <- "gene_name"
  }

  if (new_id == from_id) {
    sel <- which(probes %in% annot_table[, from_id])
    counts <- counts[sel, , drop = FALSE]
    if (type == "vector") counts <- counts[, 1]
    return(counts)
  }

  if (type == "vector") counts <- cbind(counts)

  keep.prefix <- (keep.prefix && all(grepl(":", probes)))

  from <- annot_table[, from_id]
  ii <- match(probes, from)
  if (keep.prefix) {
    dt <- mofa.get_prefix(probes)
    new.name <- annot_table[ii, new_id]
    new.name <- paste0(dt, ":", new.name)
  } else {
    new.name <- annot_table[ii, new_id]
  }
  rownames(counts) <- new.name

  if (na.rm) {
    counts <- counts[!rownames(counts) %in% c("", "NA", NA), , drop = FALSE]
  }

  ndup <- sum(duplicated(rownames(counts)))
  if (unique && ndup > 0) {
    rowdup <- rownames(counts)[which(duplicated(rownames(counts)))]
    ii <- which(rownames(counts) %in% rowdup)
    nodup.counts <- rowmean(counts[ii, , drop = FALSE], rownames(counts)[ii])
    rown <- unique(rownames(counts))
    counts <- rbind(counts[-ii, , drop = FALSE], nodup.counts)
    counts <- counts[rown, , drop = FALSE]
  }

  if (type == "vector") counts <- counts[, 1]

  return(counts)

}

#' Tidy a dataframe with type inference
#' @param Y Data.frame or matrix to tidy.
#' @return Data.frame with inferred column types.
#' @keywords internal
tidy.dataframe <- function(Y) {

  Y <- Y[, which(colMeans(is.na(Y)) < 1), drop = FALSE]
  Y <- apply(Y, 2, function(x) sub("^NA$", NA, x))
  Y <- Y[, which(colMeans(is.na(Y)) < 1), drop = FALSE]
  Y <- apply(Y, 2, function(x) gsub("^[ ]*|[ ]*$", "", x))

  suppressWarnings(num.Y <- apply(Y, 2, function(x) as.numeric(as.character(x))))

  is.numeric <- (0.8 * colMeans(is.na(num.Y)) <= colMeans(is.na(Y)))
  nlevel <- apply(Y, 2, function(x) length(unique(x)))

  is.factor <- (!is.numeric | (is.numeric & nlevel <= 3))
  is.factor <- (is.factor | grepl("batch|replicat|type|clust|group", colnames(Y)))

  new.Y <- data.frame(Y, check.names = FALSE)
  new.Y[, which(is.numeric)] <- num.Y[, which(is.numeric), drop = FALSE]

  for (i in which(is.numeric)) new.Y[[i]] <- num.Y[, i]

  for (i in which(is.factor)) new.Y[[i]] <- factor(as.character(new.Y[, i]))

  new.Y <- data.frame(new.Y, check.names = FALSE)

  return(new.Y)

}

#' Expand phenotype matrix to binary design
#' @param M Phenotype data.frame or matrix.
#' @param drop.ref Drop reference level columns.
#' @param keep.numeric Keep numeric columns as-is.
#' @param check Perform type checking.
#' @return Expanded binary design matrix.
#' @keywords internal
expandPhenoMatrix <- function(M,
                              drop.ref = TRUE,
                              keep.numeric = FALSE,
                              check = TRUE) {

  a1 <- tidy.dataframe(M)
  nlevel <- apply(a1, 2, function(x) length(setdiff(unique(x), NA)))
  nterms <- colSums(!is.na(a1))
  nratio <- nlevel / nterms

  if (inherits(a1, "data.frame")) {
    a1.typed <- utils::type.convert(a1, as.is = TRUE)
    y.class <- sapply(a1.typed, function(a) class(a)[1])
  } else {
    a1.typed <- utils::type.convert(a1, as.is = TRUE)
    y.class <- apply(a1.typed, 2, function(a) class(a)[1])
  }

  is.fac <- rep(FALSE, ncol(a1))
  is.int <- (y.class == "integer")
  ii <- which(is.int)
  is.fac[ii] <- apply(a1[, ii, drop = FALSE], 2, function(x) {
    nlev <- length(unique(x[!is.na(x)]))
    max(x, na.rm = TRUE) %in% c(nlev, nlev - 1)
  })
  is.fac2 <- (y.class == "integer" & nlevel <= 3 & nratio < 0.66)
  y.class[is.fac | is.fac2] <- "character"

  y.isnum <- (y.class %in% c("numeric", "integer"))
  kk <- which(y.isnum | (!y.isnum & nlevel > 1 & nratio < 0.66))

  if (length(kk) == 0) {
    kk <- which(y.isnum | (!y.isnum & nlevel > 1))
  }

  if (length(kk) == 0) return(NULL)

  a1 <- a1[, kk, drop = FALSE]
  a1.isnum <- y.isnum[kk]

  m1 <- list()
  for (i in 1:ncol(a1)) {
    if (a1.isnum[i]) {
      suppressWarnings(x <- as.numeric(a1[, i]))
      if (keep.numeric) {
        m0 <- matrix(x, ncol = 1)
        colnames(m0) <- "#"
      } else {
        if (drop.ref) {
          m0 <- matrix((x > stats::median(x, na.rm = TRUE)), ncol = 1)
          colnames(m0) <- "high"
        } else {
          mx <- stats::median(x, na.rm = TRUE)
          m0 <- matrix(cbind(x <= mx, x > mx), ncol = 2)
          colnames(m0) <- c("low", "high")
        }
      }
    } else if (drop.ref && nlevel[i] == 2) {
      x <- as.character(a1[, i])
      x1 <- tail(sort(x), 1)
      m0 <- matrix(x == x1, ncol = 1)
      colnames(m0) <- x1
    } else {
      x <- as.character(a1[, i])
      x[is.na(x) | x == "NA" | x == " "] <- "_"
      m0 <- stats::model.matrix(~ 0 + x)
      colnames(m0) <- sub("^x", "", colnames(m0))
    }
    rownames(m0) <- rownames(a1)
    if ("_" %in% colnames(m0)) {
      m0 <- m0[, -which(colnames(m0) == "_")]
    }
    m1[[i]] <- m0
  }

  names(m1) <- colnames(a1)

  for (i in 1:length(m1)) {
    colnames(m1[[i]]) <- paste0(names(m1)[i], "=", colnames(m1[[i]]))
  }

  m1 <- do.call(cbind, m1)
  colnames(m1) <- sub("=#", "", colnames(m1))
  rownames(m1) <- rownames(M)

  return(m1)

}

#' Compute supercells for single-cell data
#' @param counts Count matrix (genes x cells).
#' @param meta Cell metadata data.frame.
#' @param group Grouping variable name or vector.
#' @param gamma Graining level for SuperCell.
#' @param nvargenes Number of variable genes to use.
#' @param log.transform Whether to log-CPM transform.
#' @return List with counts, meta, and membership.
#' @keywords internal
pgx.supercell <- function(counts,
                          meta,
                          group = NULL,
                          gamma = 20,
                          nvargenes = 1000,
                          log.transform = TRUE) {

  if (!requireNamespace("SuperCell", quietly = TRUE)) {
    stop("Package 'SuperCell' is required for pgx.supercell(). Please install it.")
  }

  if (log.transform) {
    X <- logCPM(counts, total = 1e4)
  } else {
    X <- counts
  }

  if (is.null(group) && "group" %in% colnames(meta)) {
    message("using group column detected in meta\n")
    group <- meta[, "group"]
  }

  if (!is.null(group) && any(group %in% colnames(meta))) {
    group <- intersect(group, colnames(meta))
    message("using groups: ", paste(group, collapse = "."))
    group <- meta[, group]
    if (NCOL(group) > 1) group <- apply(group, 1, paste, collapse = ".")
  }

  SC <- SuperCell::SCimplify(X, gamma = gamma,
    n.var.genes = nvargenes, cell.split.condition = group)
  message("[pgx.supercell] SuperCell::SCimplify completed")

  meta <- as.data.frame(meta)
  dsel <- which(sapply(meta, class) %in% c("factor", "character", "logical"))
  group.argmax <- function(x) tapply(x, SC$membership, function(x) names(which.max(table(x))))
  dmeta <- apply(meta[, dsel, drop = FALSE], 2, function(x) as.character(group.argmax(x)))
  rownames(dmeta) <- sort(unique(SC$membership))
  csel <- which(sapply(meta, class) %in% c("numeric", "integer"))
  group.mean <- function(x) tapply(x, SC$membership, function(x) mean(x, na.rm = TRUE))
  cmeta <- apply(meta[, csel, drop = FALSE], 2, function(x) group.mean(x))

  sc.meta <- data.frame(dmeta)
  if (length(csel) > 0) sc.meta <- cbind(sc.meta, cmeta)
  ii <- setdiff(match(colnames(meta), colnames(sc.meta)), NA)
  sc.meta <- sc.meta[, ii, drop = FALSE]

  counts <- as.matrix(counts)
  if (log.transform) {
    sc.counts <- SuperCell::supercell_GE(counts, mode = "sum", groups = SC$membership)
  } else {
    sc.counts <- SuperCell::supercell_GE(counts, mode = "average", groups = SC$membership)
  }

  message("[pgx.supercell] SuperCell::supercell_GE completed")
  sc.membership <- paste0("mc", SC$membership)
  colnames(sc.counts) <- paste0("mc", 1:ncol(sc.counts))
  rownames(sc.meta) <- colnames(sc.counts)

  return(list(counts = sc.counts, meta = sc.meta, membership = sc.membership))

}


#' Calculate sparse correlation matrix handling missing values
#' @param G Sparse matrix containing gene sets
#' @param mat Matrix of values
#' @return Correlation matrix between G and mat
#' @details If mat has no missing values, calculates corr. using corSparse.
#' Otherwise computes column-wise correlations only using non-missing values.
#' @export
cor_sparse_matrix <- function(G, mat) {

  if (sum(is.na(mat)) == 0) {
    cor_matrix <- qlcMatrix::corSparse(G, mat)
  } else {
    message("matrix has missing values: computing column-wise reduced cor")
    corSparse.vec <- function(X, y) {
      jj <- which(!is.na(y))
      qlcMatrix::corSparse(X[jj, , drop = FALSE], cbind(y[jj]))
    }
    cor_matrix <- lapply(1:ncol(mat), function(i) corSparse.vec(G, mat[, i]))
    cor_matrix <- do.call(cbind, cor_matrix)
  }

  return(cor_matrix)

}

#' Calculate gene set rank correlation
#' Compute rank correlation between a gene rank vector/matrix and gene sets
#' @param rnk Numeric vector or matrix of gene ranks, with genes as row names
#' @param gset Numeric matrix of gene sets, with genes as row/column names
#' @param compute.p Logical indicating whether to compute p-values
#' @param use.rank Logical indicating whether to rank transform rnk before correlation
#' @return Named list with components:
#' \itemize{
#'  \item rho - Matrix of correlation coefficients between rnk and gset
#'  \item p.value - Matrix of p-values for correlation (if compute.p = TRUE)
#'  \item q.value - Matrix of FDR adjusted p-values (if compute.p = TRUE)
#' }
#' @details This function calculates sparse rank correlation between rnk and each
#' column of gset using \code{qlcMatrix::corSparse()}. It handles missing values in
#' rnk by computing column-wise correlations.
#' P-values are computed from statistical distribution
#' @examples
#' \dontrun{
#' librart(playbase)
#' ranks <- sample(1:10000, 1000, replace = TRUE)
#' names(ranks) <- replicate(1000, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
#' genesets <- matrix(rnorm(1000 * 20), ncol = 20)
#' rownames(genesets) <- names(ranks)
#' gset.rankcor(ranks, genesets, compute.p = TRUE)
#' }
#' @export
gset.rankcor <- function(rnk, gset, compute.p = FALSE, use.rank = TRUE) {

  if (ncol(gset) == 0 || NCOL(rnk) == 0) {
    if (ncol(gset) == 0) message("ERROR. gset has zero columns")
    if (NCOL(rnk) == 0) message("ERROR. rnk has zero columns")
    return(NULL)
  }

  if (!inherits(gset, "Matrix")) stop("gset must be a matrix")

  is.vec <- (NCOL(rnk) == 1 && !any(class(rnk) %in% c("matrix", "Matrix")))
  if (is.vec && is.null(names(rnk))) stop("rank vector must be named")
  if (!is.vec && is.null(rownames(rnk))) stop("rank matrix must have rownames")
  if (is.vec) rnk <- matrix(rnk, ncol = 1, dimnames = list(names(rnk), "rnk"))

  n1 <- sum(rownames(rnk) %in% colnames(gset), na.rm = TRUE)
  n2 <- sum(rownames(rnk) %in% rownames(gset), na.rm = TRUE)

  if (n1 > n2) gset <- Matrix::t(gset)

  gg <- intersect(rownames(gset), rownames(rnk))
  rnk1 <- rnk[gg, , drop = FALSE]
  gset <- gset[gg, , drop = FALSE]

  if (use.rank) {
    if (inherits(rnk1, "dgCMatrix")) {
      rnk1 <- sparseMatrixStats::colRanks(rnk1, na.last = "keep", ties.method = "random", preserveShape = TRUE)
    } else {
      rnk1 <- matrixStats::colRanks(rnk1, na.last = "keep", ties.method = "random", preserveShape = TRUE)
    }
  }

  ## (1) If no missing values: use corSparse on whole matrix.
  ## (2) If rnk matrix has missing values: proceed 1-column at time
  ## and do reduced corSparse on intersection of genes.
  rho1 <- cor_sparse_matrix(gset, rnk1)

  rownames(rho1) <- colnames(gset)
  colnames(rho1) <- colnames(rnk1)
  rho1[is.nan(rho1)] <- NA

  .cor.pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))

  df <- list(rho = rho1, p.value = NA, q.value = NA)
  if (compute.p) {
    pv <- apply(rho1, 2, function(x) .cor.pvalue(x, n = nrow(rnk1)))
    pv[is.nan(pv)] <- NA
    qv <- apply(pv, 2, stats::p.adjust, method = "fdr")
    df[["p.value"]] <- pv
    df[["q.value"]] <- qv
  }
  
  return(df)

}


#' Convert binary matrix to GMT list
#' This function converts binary matrix to a GMT (Gene Matrix
#' Transposed) list, The binary matrix indicates presence or absence of
#' genes in each gene set: rows are genes and columns are gene sets.
#' @param mat Matrix with non-zero entries representing genes in
#'   each gene set where rows are genes and columns are  gene sets.
#' @return List of vector representing each gene set. Each list
#'   element correspond to a gene set and is a vector of genes
#'@export
mat2gmt <- function(mat) {

  idx <- Matrix::which(mat != 0, arr.ind = TRUE)
  gmt <- tapply(as.character(rownames(idx)), idx[, 2], list)
  names(gmt) <- colnames(mat)[as.integer(names(gmt))]

  return(gmt)

}

#' @export
ai.ask <- function(question,
                   model,
                   engine = c("ellmer", "tidyprompt")[2]) {

  if (model == "ellmer" && grepl("grok", model)) model <- "tidyprompt"

  if (engine == "ellmer") {
    resp <- ai.ask_ellmer(question = question, model = model, prompt = NULL) 
  }

  if (engine == "tidyprompt") {
    resp <- ai.ask_tidyprompt(question = question, model = model) 
  }

  return(resp)

}

#' @export
ai.ask_ellmer <- function(question,
                          model = DEFAULT_LLM,
                          prompt = NULL) {

  chat <- NULL

  if (inherits(model, "Chat")) {
    chat <- model
  } else if (is.character(model)) {
    if (model %in% OLLAMA_MODELS || grepl("^ollama:", model) ) {
      model1 <- sub("^ollama:", "", model)
      chat <- ellmer::chat_ollama(model = model1, system_prompt = prompt)
    } else if (grepl("^gpt|^openai:",model) && Sys.getenv("OPENAI_API_KEY") != "") {
      message("warning: using remote GPT model:", model)
      model1 <- sub("^openai:", "", model)
      key <- Sys.getenv("OPENAI_API_KEY")
      chat <- ellmer::chat_openai(model = model1, system_prompt = prompt, api_key = key)
    } else if (grepl("^grok|^xai:",model) && Sys.getenv("XAI_API_KEY") != "") {
      model1 <- sub("^xai:","",model)
      key <- Sys.getenv("XAI_API_KEY")
      chat <- ellmer::chat_openai(model = model1, system_prompt = prompt,
        api_key = key, base_url = "https://api.x.ai/v1/")
    } else if (grepl("^groq:",model) && Sys.getenv("GROQ_API_KEY") != "") {
      model1 <- sub("groq:", "", model)
      key <- key Sys.getenv("GROQ_API_KEY")
      chat <- ellmer::chat_groq(model = model1, system_prompt = prompt,api_key = key)
    } else if (grepl("^gemini|^google:",model) && Sys.getenv("GEMINI_API_KEY")!="") {
      model1 <- sub("^google:","",model)
      key <- Sys.getenv("GEMINI_API_KEY")
      chat <- ellmer::chat_google_gemini(model = model1, system_prompt = prompt, api_key = key)
    }
  }

  if (is.null(chat)) {
    message("ERROR. could not create model ", model)
    return(NULL)
  }

  . <- chat$chat(question, echo = FALSE)

  chat$last_turn()@text

}

ai.ask_tidyprompt <- function(question,
                              model,
                              verbose = 0) {

  llm <- NULL
  if (model %in% OLLAMA_MODELS || grepl("^ollama:", model) ) {
    model1 <- sub("^ollama:", "", model)
    prms <- list(model = model1)
    llm <- tidyprompt::llm_provider_ollama(parameters = prms)
  } else if (grepl("^remote:", model) ) {
    remotesrv <- Sys.getenv("OLLAMA_REMOTE")
    if (remotesrv == "") message("error: please set OLLAMA_REMOTE")
    if (remotesrv != "") {
      model1 <- sub("^remote:", "", model)    
      if (verbose > 0) {
        message("connecting to remote ollama server = ", remotesrv)
        message("remote model = ", model1)        
      }
      prms <- list(model = model1)
      url <- paste0("http://", remotesrv, "/api/chat")
      llm <- tidyprompt::llm_provider_ollama(parameters = prms, url = url)
    }
  } else if (grepl("^groq:", model)) {
    model2 <- sub("groq:", "", model)
    prms <- list(model = model2)
    llm <- tidyprompt::llm_provider_groq(parameters = prms)
  } else if (grepl("^grok|^xai:", model)) {
    model2 <- sub("^xai:", "", model)
    prms <- list(model = model2)
    llm <- tidyprompt::llm_provider_xai(parameters = prms)
  } else if (grepl("^gpt-|^openai:", model)) {
    model2 <- sub("^openai:", "", model)
    prms <- list(model = model2)
    llm <- tidyprompt::llm_provider_openai(parameters = prms)
  } else if (grepl("^gemini-|^google:", model)) {
    model2 <- sub("^google:", "", model)
    prms <- list(model = model2)
    key <- Sys.getenv("GEMINI_API_KEY")
    llm <- tidyprompt::llm_provider_google_gemini(parameters = prms, api_key = key)
  }

  if (is.null(llm)) {
    message("warning. unsupported model: ", model)
    return(NULL)
  }
  
  if (verbose > 0) {
    message("model = ", model)
    message("question = ", question)    
  }

  resp <- NULL

  resp <- question |>
    tidyprompt::send_prompt(
      llm_provider = llm,
      clean_chat_history = TRUE,
      verbose = FALSE,
      return_mode = "only_response"
    )

  resp <- sub("<think>.*</think>", "", resp)

  return(resp)

}
