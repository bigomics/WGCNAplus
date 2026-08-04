#' Purple-grey-yellow color palette
#' @param n Number of colors to generate.
#' @return Character vector of hex colors.
#' @keywords internal
#' @export
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
#' @export
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

#' Validate dendrogram heights across powers
#' @param datExpr Numeric expression data matrix.
#' @param n Number of top-variance features to sample.
#' @param powers Numeric vector of powers to test.
#' @param maxpower Maximum power to evaluate.
#' @return List with quantiles, IQR, and optimal power.
#' @keywords internal
#' @export
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
    TOM <- fastTOMsimilarity(A, lowrank = 40)
    hc <- fastcluster::hclust(as.dist(1 - TOM), method = "average")
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

#' Get module-trait correlation data
#' @param wgcna A WGCNA result object.
#' @return Numeric module-trait correlation matrix.
#' @keywords internal
#' @export
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
#' @export
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
#' @param Y Expanded trait matrix with "=" columns.
#' @return Collapsed matrix with category columns.
#' @keywords internal
#' @export
collapseTraitMatrix <- function(Y) {

  if (sum(grepl("=", colnames(Y))) < 2) return(Y)

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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
mofa.merge_data <- function(xx) { do.call(rbind, mofa.prefix(xx)) }

#' Split MOFA data by prefix
#' @param X Matrix with prefixed rownames.
#' @param keep.prefix Whether to keep prefix in names.
#' @return Named list of matrices per layer.
#' @keywords internal
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
mofa.merge_data2 <- function(xdata,
                             merge.rows = "prefix",
                             merge.cols = "union") {

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
    for (i in 1:length(xdata)) {
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
ai.get_ollama_models <- function(models = NULL,
                                 size = NULL) {

  available.models <- system("ollama list | tail -n +2 | cut -d' ' -f 1", intern=TRUE)

  models.sizes <- system("ollama list | tail -n +2 | tr -s ' ' | cut -d' ' -f 3", intern=TRUE)
  models.sizes <- as.numeric(models.sizes)
  models.sizes <- ifelse(models.sizes < 100, models.sizes, models.sizes/1000)
  names(models.sizes) <- available.models

  if (!is.null(models) && !any(models == "OLLAMA_MODELS")) {
    available.models <- intersect(models,available.models)
  }

  msize <- models.sizes[available.models] 
  if (!is.null(size) && size == "S") {
    available.models <- available.models[which(msize <= 3)]
  }

  if (!is.null(size) && size == "M") {
    sel <- which( msize > 3 & msize <= 6)
    available.models <- available.models[sel]
  }

  if (!is.null(size) && size == "L") {
    msize <- models.sizes[available.models] 
    available.models <- available.models[which(msize > 6)]
  }
  
  return(available.models)

}
OLLAMA_MODELS <- ai.get_ollama_models()

#' @export
ai.ask <- function(question,
                   model,
                   engine = c("ellmer", "tidyprompt")[2]) {

  if (model == "ellmer" && grepl("grok", model)) model <- "tidyprompt"

  if (! engine %in% c("ellmer", "tidyprompt"))
    stop("[WGCNAplus::ai.ask] Error. 'engine' must be 'ellmer' or 'tidyprompt'")
  
  if (engine == "ellmer")
    resp <- ai.ask_ellmer(question = question, model = model, prompt = NULL) 

  if (engine == "tidyprompt")
    resp <- ai.ask_tidyprompt(question = question, model = model) 

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
      key <- Sys.getenv("GROQ_API_KEY")
      chat <- ellmer::chat_groq(model = model1, system_prompt = prompt, api_key = key)
    } else if (grepl("^gemini|^google:",model) && Sys.getenv("GEMINI_API_KEY") != "") {
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

#' @export
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

#' Generate image with Gemini (aka Nano Banana).
#' Note this model handles very large prompts correctly.
#' @export
ai.create_image_gemini <- function(prompt,
                                   model = "gemini-2.5-flash-image",
                                   api_key = Sys.getenv("GEMINI_API_KEY"),
                                   format = c("file","base64","raw")[1], 
                                   filename = "image.png",
                                   aspectRatio = "16:9",
                                   imageSize = "1K",
                                   base_url = "https://generativelanguage.googleapis.com/v1beta") {

  assertthat::assert_that(assertthat::is.string(prompt), assertthat::noNA(prompt))
  assertthat::assert_that(assertthat::is.string(model), assertthat::noNA(model))
  assertthat::assert_that(assertthat::is.string(api_key), assertthat::noNA(api_key))
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for ai.create_image_gemini()")
  }
  `%>%` <- dplyr::`%>%`

  if (nchar(api_key) == 0) {
    stop("GEMINI_API_KEY environment variable is not set", call. = FALSE)
  }

  message("calling gemini image generation...")
  model <- sub("^google:", "", model)
  url <- glue::glue("{base_url}/models/{model}:generateContent")

  headers <- c(`x-goog-api-key` = api_key, `Content-Type` = "application/json")

  body <- list(
    contents = list(list(parts = list(list(text = prompt)))),
    generationConfig = list(
      responseModalities = list("TEXT", "IMAGE"),
      imageConfig =  list(aspectRatio = aspectRatio, imageSize = imageSize)      
    )
  )

  if (grepl("gemini-2.5",model)) {
    body$generationConfig$imageConfig <- list(aspectRatio = aspectRatio)
  }
  
  response <- httr::POST(
    url = url,
    httr::add_headers(.headers = headers),
    body = jsonlite::toJSON(body, auto_unbox = TRUE),
    encode = "raw"
  )

  httr::http_type(response)
  if (httr::http_type(response) != "application/json") {
    stop("Gemini API returned unexpected content type", call. = FALSE)
  }

  parsed <- response %>%
    httr::content(as = "text", encoding = "UTF-8") %>%
    jsonlite::fromJSON(flatten = TRUE)

  httr::http_error(response)
  if (httr::http_error(response)) {
    error_msg <- if (!is.null(parsed$error$message)) parsed$error$message else "Unknown error"
    stop(paste0("Gemini API request failed [", httr::status_code(response), "]: ", error_msg), call. = FALSE)
  }
  
  parts <- parsed$candidates$content.parts  
  b64 <- NULL
  mimetype <- NULL
  for (part in parts) {
    if (!is.null(part$inlineData.data)) {
      b64 <- part$inlineData.data
      b64 <- head(b64[!is.na(b64)],1)
      mimetype <- part$inlineData.mimeType
      mimetype <- head(mimetype[!is.na(mimetype)],1)
      break()
    }
  }
  
  if (is.null(b64) || length(b64)==0) stop("No image data found in response")

  if (format == "file") {
    raw_image <- base64enc::base64decode(b64)    
    filetype <- sub("jpeg", "jpg", sub("image/","",mimetype))
    filename2 <- paste0(sub("[.].*$", "", filename), ".", filetype)
    writeBin(raw_image, filename2)
    message("Saved image to: ", filename2)
    return(invisible(filename2))
  }

  if (format == "raw") {
    raw_image <- base64enc::base64decode(b64)    
    return(invisible(raw_image))
  }

  if (format == "base64") return(invisible(b64))

  stop("return error")

}

#' Generate image with OpenAI's dallE.
#' Note: limitation of the prompt is about 1000 characters. 
#' @export
ai.create_image_openai <- function (prompt,
                                    model = NULL, 
                                    size = c("512x512","1024x1024","256x256")[1], 
                                    format = c("file","base64","raw"),
                                    filename = "image.png",
                                    api_key = Sys.getenv("OPENAI_API_KEY"),
                                    base_url = "https://api.openai.com/v1",
                                    user = NULL,
                                    organization = NULL)  {
  
  format <- match.arg(format)
  assertthat::assert_that(assertthat::is.string(prompt), assertthat::noNA(prompt))
  assertthat::assert_that(assertthat::is.string(size), assertthat::noNA(size))
  assertthat::assert_that(assertthat::is.string(format), assertthat::noNA(format))

  if (!is.null(user)) {
    assertthat::assert_that(assertthat::is.string(user), assertthat::noNA(user))
  }

  assertthat::assert_that(assertthat::is.string(api_key), assertthat::noNA(api_key))

  if (!is.null(organization)) {
    assertthat::assert_that(assertthat::is.string(organization), assertthat::noNA(organization))
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for ai.create_image_openai()")
  }
  `%>%` <- dplyr::`%>%`

  model <- sub("^openai:", "", model)
  
  if (grepl("api.x.ai", base_url, fixed = TRUE)) {
    message("calling grok ($0.07 per image)")    
    if (is.null(model)) model <- "grok-2-image-1212"
    size <- NULL
  } else if (grepl("api.openai.com", base_url, fixed = TRUE)) {
    message("calling openai ($0.05 per image)")    
  } else {
    stop("invalid base_url =", base_url)
  }
  
  url <- glue::glue("{base_url}/images/generations")
  headers <- c(Authorization = paste("Bearer", api_key), `Content-Type` = "application/json")

  body <- list()
  body[["model"]] <- model
  body[["prompt"]] <- prompt
  body[["n"]] <- 1
  body[["response_format"]] <- "b64_json"
  body[["size"]] <- size
  body[["user"]] <- user
  response <- httr::POST(url = url, httr::add_headers(.headers = headers), body = body, encode = "json")

  httr::http_type(response)
  if (httr::http_type(response) != "application/json") {
    paste("OpenAI API probably has been changed. Please check online documentation.") %>% stop()
  }

  parsed <- response %>% httr::content(as = "text", encoding = "UTF-8") %>% jsonlite::fromJSON(flatten = TRUE)
  if (httr::http_error(response)) {
    error_msg <- parsed$error
    if(is.list(error_msg)) error_msg <- parsed$error$message
    paste0("OpenAI API request failed [", httr::status_code(response), 
      "]:\n\n", error_msg) %>% stop(call. = FALSE)
  }

  b64 <- parsed$data[['b64_json']]
  if (is.null(b64) || length(b64) == 0) stop("No image data found in parsed response")

  if (format == "file") {
    raw_image <- base64enc::base64decode(b64)    
    writeBin(raw_image, filename)
    message("Saved image to: ", filename)
    return(invisible(filename))
  }

  if (format == "raw") {
    raw_image <- base64enc::base64decode(b64)    
    return(invisible(raw_image))
  }

  if (format == "base64") return(invisible(b64))

  stop("return error")
  
}

#' Generate image with Grok (which uses Flux).
#' Note: limitation of the prompt is about 1000 characters.
#' @export
ai.create_image_grok <- function(prompt,
                                 model = "grok-2-image-1212", 
                                 format = c("file","base64","raw")[1], 
                                 api_key = Sys.getenv("XAI_API_KEY"),
                                 base_url = "https://api.x.ai/v1",
                                 filename = "image.png",
                                 user = NULL,
                                 organization = NULL) {

  model <- sub("^grok:", "", model)
  ai.create_image_openai(
    prompt,
    size = "default",
    format = format,
    filename = filename,
    model = model,
    base_url = base_url,
    api_key = api_key,
    user = user,
    organization = organization
  )

}


## =========================================================================
## Statistics helpers inlined from playbase (to drop the playbase dependency)
## =========================================================================

#' @title SVD-based missing value imputation
#' @description Iterative SVD imputation. Inlined verbatim from playbase.
#' @param X Numeric matrix with NAs to impute.
#' @param nv Number of singular vectors.
#' @param init Initialisation: NULL (row/col medians), "<q>%" quantile, or an
#'   imputeMissing method name.
#' @return Matrix with imputed values.
#' @keywords internal
#' @export
svdImpute2 <- function(X, nv = 10, threshold = 0.001, init = NULL,
                       maxSteps = 100, fill.empty = "median",
                       infinite.na = TRUE,
                       verbose = FALSE) {
  if (infinite.na) X[is.infinite(X)] <- NA
  ind.missing <- which(is.na(X), arr.ind = TRUE)
  empty.rows <- which(rowMeans(is.na(X)) == 1)
  empty.cols <- which(colMeans(is.na(X)) == 1)

  if (is.null(nv)) {
    nv <- max(1, round(mean(is.na(X)) * min(dim(X)))) ## heuristic..
    message(paste0("setting nv to ", nv))
  }
  nv <- min(nv, round(min(dim(X)) / 3))

  init.methods <- c("MinDet", "MinProb", "QRILC", "min")

  if (!is.null(init) && is.character(init) && grepl("%", init)) {
    ## initialize missing values with quantile fixed value
    q <- as.numeric(sub("%", "", init))
    init <- quantile(X[!is.na(X)], probs = q * 0.01)[1]
    message(paste0("setting initial values to ", q, "%. init=", round(init, 4)))
    X[ind.missing] <- init
  } else if (!is.null(init) && is.character(init) &&
    init %in% init.methods) {
    ## initialize missing values with other impute method
    message(paste("setting initial values using", init))
    initX <- imputeMissing(X, method = init)
    X[ind.missing] <- initX[ind.missing]
  } else {
    ## initialize missing values with col/row medians
    row.mx <- apply(X, 1, median, na.rm = TRUE)
    col.mx <- apply(X, 2, median, na.rm = TRUE)
    row.mx[is.nan(row.mx)] <- NA
    col.mx[is.nan(col.mx)] <- NA
    X[ind.missing] <- row.mx[ind.missing[, 1]]
    ind.missing2 <- which(is.na(X), arr.ind = TRUE)
    X[ind.missing2] <- col.mx[ind.missing2[, 2]]
  }

  ## do SVD iterations
  count <- 0
  error <- Inf
  Xold <- X
  nv <- min(nv, dim(X) - 1)
  while ((error > threshold) && (count < maxSteps)) {
    if (nv < min(dim(X)) / 5) {
      res <- irlba::irlba(X, nv = nv, nu = nv)
    } else {
      res <- svd(X, nv = nv, nu = nv)
      res$d <- res$d[1:nv]
    }
    if (nv == 1) {
      imx <- res$d * (res$u %*% t(res$v))
    } else {
      imx <- res$u %*% (diag(res$d) %*% t(res$v))
    }
    X[ind.missing] <- imx[ind.missing]
    count <- count + 1
    if (count > 0) {
      error <- sqrt(sum((Xold - X)^2, na.rm = TRUE) / sum(Xold^2, na.rm = TRUE))
      if (verbose) {
        cat(count, ": change in estimate: ", error, "\n")
      }
    }
    Xold <- X
  }

  ## extra corrections (refill empty columns or rows)
  has.empty <- (length(empty.rows) > 0 || length(empty.cols) > 0)
  if (has.empty && fill.empty == "NA") {
    if (length(empty.rows)) X[empty.rows, ] <- NA
    if (length(empty.cols)) X[, empty.cols] <- NA
  }
  if (has.empty && fill.empty == "sample") {
    ii <- which(
      (!ind.missing[, 1] %in% empty.rows) &
        (!ind.missing[, 2] %in% empty.cols)
    )
    mm <- X[ind.missing[ii, ]]
    if (length(empty.rows) && length(mm)) {
      n1 <- length(X[empty.rows, ])
      message("[svdImpute2] warning: empty rows : n1 = ", n1)
      X[empty.rows, ] <- sample(mm, n1, replace = TRUE)
    }
    if (length(empty.cols) && length(mm)) {
      n2 <- length(X[, empty.cols])
      message("[svdImpute2] warning: empty cols : n2 = ", n2)
      X[, empty.cols] <- sample(mm, n2, replace = TRUE)
    }
  }

  return(X)
}

#' @title Impute missing values
#' @description Inlined from playbase, trimmed to the pure-R methods WGCNAplus
#'   uses. Only "SVD2" and "rowmeans" are kept; the original's LLS / bpca /
#'   msImpute / NMF / RF / MSnbase branches pulled pcaMethods, MSnbase,
#'   missForest and msImpute and are dropped.
#' @param X Numeric matrix to impute.
#' @param method Imputation method(s); one or both of "SVD2", "rowmeans".
#' @return Matrix with imputed values, or NULL if no method matched.
#' @keywords internal
#' @export
imputeMissing <- function(X, method = "SVD2", nv = 5,
                          keep.limits = FALSE, infinite.na = TRUE) {
  ## ponytail: trimmed to pure-R methods; add a branch back if a caller needs
  ## LLS/bpca/etc (each drags its own Bioc/CRAN impute dependency).
  if (infinite.na) X[is.infinite(X)] <- NA
  impX <- list()

  if ("rowmeans" %in% method) {
    cx <- X
    ii <- which(is.na(cx), arr.ind = TRUE)
    cx[ii] <- rowMeans(cx, na.rm = TRUE)[ii[, 1]]
    ii <- which(is.na(cx), arr.ind = TRUE)
    cx[ii] <- colMeans(cx, na.rm = TRUE)[ii[, 2]]
    ii <- which(is.na(cx))
    cx[ii] <- mean(cx, na.rm = TRUE)
    impX[["rowmeans"]] <- cx
  }

  if ("SVD2" %in% method) {
    impX[["SVD2"]] <- svdImpute2(X, nv = nv, init = "5%")
  }

  if (length(impX) == 0) {
    return(NULL)
  }

  if (keep.limits) {
    minx <- min(X, na.rm = TRUE)
    maxx <- max(X, na.rm = TRUE)
    impX <- lapply(impX, function(x) pmin(pmax(x, minx), maxx))
  }

  if (length(impX) > 1) {
    metaX <- do.call(cbind, lapply(impX, as.vector))
    metaX <- matrix(rowMeans(metaX, na.rm = TRUE),
      nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X)
    )
  } else {
    metaX <- impX[[1]]
  }

  ## any remaining NA we fill with col/row medians
  if (any(is.na(metaX))) {
    missing <- which(is.na(metaX), arr.ind = TRUE)
    row.mx <- apply(metaX, 1, median, na.rm = TRUE)
    col.mx <- apply(metaX, 2, median, na.rm = TRUE)
    metaX[missing] <- rowMeans(cbind(row.mx[missing[, 1]], col.mx[missing[, 2]]), na.rm = TRUE)
  }

  metaX
}

#' @title Differential expression analysis using limma
#' @description Inlined from playbase. Only the default \code{method = 1} path
#'   is kept (the only one WGCNAplus calls); the method 2/3 branches needed
#'   playbase's makeDirectContrasts/makeFullContrasts and are dropped.
#' @param X Expression matrix, genes in rows, samples in columns.
#' @param pheno Phenotype vector/factor.
#' @return Data frame with limma results.
#' @keywords internal
#' @export
gx.limma <- function(X,
                     pheno,
                     B = NULL,
                     remove.na = TRUE,
                     fdr = 0.05,
                     compute.means = TRUE,
                     lfc = 0.20,
                     max.na = 0.20,
                     sort.by = "FC",
                     ref = c(
                       "ctrl", "ctr", "control", "ct", "dmso", "nt", "0", "0h", "0hr",
                       "non", "no", "not", "neg", "negative", "ref", "veh", "vehicle",
                       "wt", "wildtype", "untreated", "normal", "false", "healthy"
                     ),
                     trend = FALSE,
                     robust = FALSE,
                     method = 1,
                     f.test = FALSE,
                     verbose = 1) {
  if (method != 1) {
    stop("[gx.limma] vendored copy supports method = 1 only")
  }
  if (!is.null(B) && NCOL(B) == 1) {
    B <- matrix(B, ncol = 1)
    rownames(B) <- rownames(pheno)
    colnames(B) <- "batch"
  }

  ## detect single sample case
  is.single <- (max(table(pheno), na.rm = TRUE) == 1)
  if (is.single) {
    message("[gx.limma] WARNING: no replicates, duplicating samples...")
    X <- cbind(X, X)
    pheno <- c(pheno, pheno)
    if (!is.null(B)) B <- rbind(B, B)
  }

  ## filter probes and samples
  ii <- which(rowMeans(is.na(X)) <= max.na)
  jj <- 1:ncol(X)
  if (remove.na && any(is.na(pheno))) {
    jj <- which(!is.na(pheno))
    if (verbose > 0) message("[gx.limma] ", sum(is.na(pheno) > 0), " with missing phenotype")
  }
  X0 <- X[ii, jj, drop = FALSE]
  pheno0 <- as.character(pheno[jj])
  X0 <- X0[!(rownames(X0) %in% c(NA, "", "NA")), ]
  B0 <- NULL
  if (!is.null(B)) B0 <- B[jj, , drop = FALSE]

  if (verbose > 0) {
    message("[gx.limma] X contains ", sum(duplicated(rownames(X))), "duplicated rownames")
    message("[gx.limma] analyzing ", ncol(X0), " samples")
    message("[gx.limma] table.pheno: ", table(pheno), "samples")
    message("[gx.limma] testing ", nrow(X0), " features")
    message("[gx.limma] lfc = ", lfc)
    message("[gx.limma] fdr = ", fdr)
    message("[gx.limma] max.na = ", max.na)
    if (!is.null(B0)) message("[gx.limma] including ", ncol(B0), " batch covariates")
  }

  ## auto-detect reference
  pheno.ref <- c()
  ref.detected <- FALSE
  ref <- toupper(ref)
  is.ref <- (toupper(pheno0) %in% toupper(ref))
  is.ref2 <- grepl(paste(paste0("^", ref), collapse = "|"), pheno0, ignore.case = TRUE)
  if (!any(is.ref) && !all(is.ref2)) is.ref <- is.ref2
  ref.detected <- (sum(is.ref) > 0 && sum(!is.ref) > 0)

  if (ref.detected) {
    pheno.ref <- unique(pheno0[is.ref])
    if (verbose > 0) message("[gx.limma] setting reference to y=", pheno.ref, "\n")
    levels <- c(pheno.ref, sort(setdiff(unique(pheno0), pheno.ref)))
  } else {
    if (verbose > 0) message("[gx.limma] WARNING: could not auto-detect reference\n")
    levels <- as.character(sort(unique(pheno0)))
    if (verbose > 0) message("[gx.limma] setting reference to first class:", levels[1], "\n")
  }

  if (length(levels) != 2 & !f.test) {
    stop("[gx.limma] ERROR: only two class comparisons. Please activate F test using f.test=TRUE\n\n")
    return(NULL)
  }

  pheno1 <- stats::relevel(factor(pheno0), ref = levels[1])

  ## setup model and perform LIMMA (method = 1, no contrast matrix)
  design <- cbind(1, pheno0 == levels[2])
  colnames(design) <- c("intercept", "main_vs_ref")
  if (f.test) {
    design <- stats::model.matrix(~pheno1)
    colnames(design)[2:ncol(design)] <- paste0(levels(pheno1)[-1], "_vs_", levels(pheno1)[1])
  }
  if (!is.null(B0)) {
    if (verbose > 0) {
      message("[gx.limma] augmenting design matrix with: ", paste(colnames(B0)), "\n")
    }
    sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
    design <- cbind(design, B0[, sel, drop = FALSE])
  }
  fit <- limma::lmFit(X0, design)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  coef <- if (f.test) NULL else "main_vs_ref"
  top <- suppressMessages(limma::topTable(fit, coef = coef, number = nrow(X0), sort.by = "none"))

  if (f.test) {
    cols <- c("logFC", "AveExpr", "F", "P.Value", "adj.P.Val")
    kk <- intersect(colnames(top), cols)
    top <- top[, kk]
    top$B <- NULL
  }

  if ("ID" %in% colnames(top)) {
    rownames(top) <- top$ID
    top$ID <- NULL
  }

  if (f.test) {
    kk <- setdiff(colnames(top), colnames(design))
    top <- top[, kk, drop = FALSE]
  }

  top <- top[rownames(X0), , drop = FALSE]

  if (f.test) {
    ## compute averages
    avg <- do.call(cbind, tapply(1:ncol(X0), pheno1, function(i) {
      rowMeans(X0[, i, drop = FALSE], na.rm = TRUE)
    }))
    if (!"logFC" %in% colnames(top)) {
      maxFC <- apply(avg, 1, max, na.rm = TRUE) - apply(avg, 1, min, na.rm = TRUE)
      top$logFC <- NULL
      top <- cbind(logFC = maxFC, top)
      rownames(top) <- rownames(X0)
    }
  }

  if (!is.null(fdr) && !is.null(lfc)) {
    ii <- which(top$adj.P.Val <= fdr & abs(top$logFC) >= lfc)
    top <- top[ii, ]
    if (verbose > 0) {
      message(
        "[gx.limma] Found ", nrow(top), " significant at fdr = ",
        fdr, " and minimal FC = ", lfc, "\n"
      )
    }
  }

  if (compute.means && nrow(top) > 0) {
    if (f.test) {
      avg <- avg[rownames(top), ]
    } else {
      avg <- t(apply(X0[rownames(top), ], 1, function(x) tapply(x, pheno0, mean, na.rm = TRUE)))
      avg <- avg[, as.character(levels), drop = FALSE]
    }
    colnames(avg) <- paste0("AveExpr.", colnames(avg))
    top <- cbind(top, avg)
  }
  top$B <- NULL

  if (is.single) {
    top$P.Value <- NA
    top$adj.P.Val <- NA
    top$t <- NA
    if (f.test) top$F <- NA
  }

  if (sort.by == "FC") top <- top[order(abs(top$logFC), decreasing = TRUE), ]
  if (sort.by == "p") top <- top[order(top$P.Value), ]

  return(top)
}

#' Merge multiple ME matrices into one. Allow different dimensions.
#' @param mlist List of ME matrices.
#' @param me2 Optional second ME matrix.
#' @param prefix Prefix columns with list names.
#' @return Merged eigengene matrix.
#' @keywords internal
mergeME <- function(mlist,
                    me2 = NULL,
                    prefix = FALSE) {

  if (!is.null(me2) && !inherits(mlist,"list")) mlist <- list(mlist, me2)

  all.samples <- unique(unlist(sapply(mlist, rownames, simplify=FALSE)))

  if (prefix) {
    for (i in 1:length(mlist)) {
      colnames(mlist[[i]]) <- paste0(names(mlist)[i],":",colnames(mlist[[i]]))
    }
  }

  is.mat <- all(sapply(mlist, inherits, what="matrix"))
  all.me <- unique(unlist(sapply(mlist, colnames, simplify=FALSE)))

  M <- as.data.frame(matrix(NA, nrow = length(all.samples), ncol = length(all.me)))
  rownames(M) <- all.samples
  colnames(M) <- all.me

  for (i in 1:length(mlist)) {
    ii <- match(rownames(mlist[[i]]), rownames(M))
    jj <- match(colnames(mlist[[i]]), colnames(M))
    M[ii, jj] <- mlist[[i]]
  }

  if (is.mat) M <- as.matrix(M)

  return(M)

}

#' Get top correlated modules
#' @param wgcna A WGCNA result object.
#' @param topratio Ratio threshold for top selection.
#' @param kx Power exponent for ranking.
#' @param rm.grey Remove grey module.
#' @param multi Use multi-omics mode.
#' @return Character vector of top module names.
#' @export
getTopModules <- function(wgcna,
                          topratio = 0.85,
                          kx = 4,
                          rm.grey = TRUE,
                          multi = FALSE) {

  if (!multi) {
    ww <- list(gx = wgcna)  ## single-omics wgcna object
  } else if(!is.null(wgcna$layers)) {
    ww <- wgcna$layers
  } else {
    ww <- wgcna
  }

  M <- list()
  for (i in 1:length(ww)) {
    me <- ww[[i]]$net$MEs
    dt <- ww[[i]]$datTraits
    M[[i]] <- cor(me, dt, use="pairwise")
  }

  top.modules <- c()
  for (i in 1:length(M)) {
    mx <- rowMeans(abs(M[[i]]**kx),na.rm=TRUE)**(1/kx)
    tt <- names(which( mx > topratio * max(mx)))
    top.modules <- c(top.modules, tt)
  }

  if (rm.grey) {
    sel.grey <- grepl("[A-Z]{2}grey$",top.modules)
    top.modules <- top.modules[!sel.grey]
  }
  
  return(top.modules)

}

#' Get multi-dataset top genes and sets tables
#' @param multi_wgcna Multi-omics WGCNA object.
#' @param annot Annotation table or NULL.
#' @param module Module names to select.
#' @param psig P-value significance threshold.
#' @param ntop Number of top entries.
#' @param level Feature level or NULL.
#' @param rename Column name for renaming.
#' @return List with top sets, genes, and pheno.
#' @keywords internal
#' @export
getTopTables <- function(wgcna,
                         annot = NULL,
                         module = NULL,
                         psig = 0.05,
                         ntop = 40,
                         level = NULL,
                         rename = "symbol") {
  
  if ("layers" %in% names(wgcna)) {
    layers <- wgcna$layers
  } else if (all(c("datExpr","datTraits") %in% names(wgcna))) {
    layers <- list(gx = wgcna)
  } else {
    layers <- wgcna
  }

  ## set level
  nw <- length(layers)
  if (!is.null(level)) {
    level <- head(rep(level, nw),nw)
  } else {
    level <- c("gene","geneset")[1 + 1*grepl("^gs|^gset|geneset",names(layers))]
  }
  names(level) <- names(layers)

  toplist <- list()
  for (k in names(layers)) {
    topk <- getTopGenesAndSets(layers[[k]], module = module,  annot = annot,
      ntop = ntop, psig = psig, level = level[[k]], rename = rename)
    if (!is.null(module)) {
      topk <- lapply( topk, function(s) s[which(names(s) %in% module)] )
    }
    toplist[[k]] <- topk
  }

  top <- list()
  top$genes <- lapply(toplist, function(t) t[['genes']])
  names(top$genes) <- NULL
  top$genes <- unlist(top$genes, recursive = FALSE)

  top$sets <- lapply(toplist, function(t) t[["sets"]])
  names(top$sets) <- NULL
  top$sets <- unlist(top$sets, recursive = FALSE)

  top$pheno <- lapply(toplist, function(t) t[["pheno"]])
  names(top$pheno) <- NULL
  top$pheno <- unlist(top$pheno, recursive = FALSE)

  return(top)

}


#' Get top genes and gene sets per module
#' @param wgcna A WGCNA result object.
#' @param annot Annotation table or NULL.
#' @param module Module names to select.
#' @param ntop Number of top entries.
#' @param psig P-value significance threshold.
#' @param level Feature level: "gene" or "geneset".
#' @param rename Column name for renaming.
#' @return List with top sets, genes, and pheno.
#' @keywords internal
#' @export
getTopGenesAndSets <- function(wgcna,
                               annot = NULL,
                               module = NULL,
                               ntop = 40,
                               psig = 0.05,
                               level = "gene",
                               rename = "symbol") {

  
  if ("layers" %in% names(wgcna) && class(wgcna$datExpr) == "list") {
    message("[getTopGenesAndSets] multilayer object...")
    cons <- .getConsensusTopGenesAndSets(wgcna, annot=annot, module=module,  ntop=ntop, rename=rename)
    return(cons)
  }

  stats <- NULL
  if (!"stats" %in% names(wgcna) || is.null(wgcna$stats) ) {
    stats <- computeGeneStats(wgcna$net, wgcna$datExpr, wgcna$datTraits, wgcna$svTOM)
  } else {
    stats <- wgcna$stats
  }

  if (!any(c("gse","gsea") %in% names(wgcna))) {
    warning("object has no enrichment results (gsea)")
  }

  ## get top genes by centrality-weighted-meanFC2
  mm <- stats$moduleMembership
  mm.sig <- 1*(stats$MMPvalue <= psig)
  ff <- sqrt(rowMeans(stats$foldChange**2, na.rm=TRUE))
  mm <- mm * mm.sig * ff
  if (!is.null(annot)) {
    annot$gene_title <- paste0(annot$gene_title," (",annot$symbol,")")
    mm <- rename_by2(mm, annot, new_id=rename)
  }
  gg <- rownames(mm)
  mm <- as.list(data.frame(mm))
  if (!is.null(module)) mm <- mm[which(names(mm) %in% module)]
  for (i in 1:length(mm)) names(mm[[i]]) <- gg
  mm <- lapply(mm, function(x) x[x!=0])
  topgenes <- lapply(mm, function(x) names(head(sort(-x),ntop)))

  ## top genesets
  topsets <- NULL
  if (any(c("gse","gsea") %in% names(wgcna))) {
    if (!is.null(wgcna$gsea)) ee <- wgcna$gsea
    if (!is.null(wgcna$gse)) ee <- wgcna$gse
    if (!is.null(module)) ee <- ee[which(names(ee) %in% module)]
    topsets <- lapply(ee,function(x) head(rownames(x),ntop))
  }

  ## top correlated phenotypes
  M <- get_modTraits(wgcna)
  toppheno <- apply(M, 1, function(x) names(which(x > 0.8*max(x, na.rm=TRUE))))

  if (level == "geneset") {
    topsets <- topgenes
    topgenes <- NULL
  }

  LL <- list(sets = topsets, genes = topgenes, pheno = toppheno)

  return(LL)
  
}


#' Get consensus top genes and sets
#' @param cons Consensus WGCNA object.
#' @param annot Annotation table or NULL.
#' @param module Module names to select.
#' @param ntop Number of top entries.
#' @param level Feature level: "gene" or "geneset".
#' @param rename Column name for renaming.
#' @return List with top sets, genes, and pheno.
#' @keywords internal
.getConsensusTopGenesAndSets <- function(cons,
                                         annot = NULL,
                                         module = NULL,
                                         ntop = 40,
                                         level = c("gene","geneset")[1],
                                         rename = "symbol" ) {

  if (!"stats" %in% names(cons)) stop("object has no stats")
  if (!any(c("gse","gsea") %in% names(cons))) {
    warning("object has no enrichment results (gsea)")
  }

  if (!is.null(annot)) {
    annot$gene_title <- paste0(annot$gene_title," (",annot$symbol,")")
  }

  ## get top genes (highest kME)
  topgenesx <- list()
  for (i in 1:length(cons$stats)) {
    mm <- cons$stats[[i]]$moduleMembership
    if (!is.null(annot)) mm <- rename_by2(mm, annot, rename)
    gg <- rownames(mm)
    mm <- as.list(data.frame(mm))
    if (!is.null(module)) mm <- mm[module]
    sel.topgenes <- lapply(mm, function(x) head(order(-x), 3 * ntop))
    topgenesx[[i]] <- lapply(sel.topgenes, function(i) gg[i])
  }

  ## intersect topgenes across all datatypes
  topgenes <- topgenesx[[1]]
  k <- 2
  for (k in 2:length(topgenesx)) {
    topgenes <- mapply(intersect, topgenes, topgenesx[[k]], SIMPLIFY = FALSE)
  }
  topgenes <- lapply(topgenes, head, ntop)

  if (!is.null(module)) {
    sel <- intersect(names(topgenes),module)
    topgenes <- topgenes[sel]
  }

  ## top genesets (as symbol!)
  topsets <- NULL
  if (any(c("gse","gsea") %in% names(cons))) {
    if (!is.null(cons$gsea)) ee <- cons$gsea
    if (!is.null(cons$gse)) ee <- cons$gse
    ee <- ee[match(names(topgenes),names(ee))]
    names(ee) <- names(topgenes)
    topsets <- lapply(ee,function(x) head(rownames(x),ntop))
  }

  ## module traits
  M <- lapply(cons$net$multiMEs, function(x) as.matrix(x$data))
  Y <- lapply(M, function(m) cons$datTraits[rownames(m),])
  R <- mapply( function(x,y) abs(cor(x,y,use="pairwise")), M, Y, SIMPLIFY=FALSE)
  R <- Reduce('+', R)
  toppheno <- apply(R, 1, function(x) names(which(x > 0.9*max(x,na.rm=TRUE))),
    simplify = FALSE)
  toppheno

  if (level == "geneset") {
    topsets <- topgenes
    topgenes <- NULL
  }

  return(list(sets = topsets, genes = topgenes, pheno = toppheno))

}
