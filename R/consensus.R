## =========================================================================
## CONSENSUS WGCNA
## =========================================================================

#' Run consensus WGCNA across multiple datasets
#' @param exprList List of expression matrices.
#' @param phenoData Sample phenotype data frame.
#' @param contrasts Optional contrast matrix.
#' @param GMT Gene-set matrix for enrichment.
#' @param annot Annotation table for gene mapping.
#' @param ngenes Maximum number of genes to use.
#' @param power Soft-thresholding power.
#' @param minModuleSize Minimum module size.
#' @param minKME Minimum KME to stay in module.
#' @param mergeCutHeight Cut height for merging modules.
#' @param deepSplit Deep split sensitivity level.
#' @param maxBlockSize Maximum block size for computation.
#' @param addCombined Add combined dataset layer.
#' @param calcMethod TOM calculation method.
#' @param drop.ref Drop reference level in traits.
#' @param cons.psig P-value threshold for consensus.
#' @param compute.stats Compute gene statistics.
#' @param compute.enrichment Compute module enrichment.
#' @param summary Generate AI module summaries.
#' @param ai_model LLM model for annotation.
#' @param experiment Experiment description string.
#' @param gsea.mingenes Minimum genes for enrichment.
#' @param gsea.ntop Top genes for enrichment.
#' @param gset.methods Enrichment methods to use.
#' @param verbose Verbosity level.
#' @param progress Optional progress callback.
#' @return List with consensus WGCNA results.
#' @export
runConsensusWGCNA <- function(exprList,
                              phenoData,
                              contrasts = NULL,
                              GMT = NULL,
                              annot = NULL,
                              ngenes = 2000,
                              power = 12,
                              minModuleSize = 20,
                              minKME = 0.3,
                              mergeCutHeight = 0.15,
                              deepSplit = 2,
                              maxBlockSize = 9999,
                              addCombined = FALSE,
                              calcMethod = "fast",
                              drop.ref = FALSE,
                              cons.psig = 0.05,
                              compute.stats = TRUE,
                              compute.enrichment = TRUE,
                              summary = TRUE,
                              ai_model = getOption("WGCNAplus.default_llm"),
                              experiment = "",
                              gsea.mingenes = 10,
                              gsea.ntop = 1000,
                              gset.methods = c("fisher", "gsetcor", "xcor"),
                              verbose = 1,
                              progress = NULL) {
  
  colors <- NULL
  
  ## Align and reduce matrices if needed
  gg <- Reduce(intersect, lapply(exprList, rownames))
  exprList <- lapply(exprList, function(x) x[gg, , drop = FALSE])
  if (length(gg) > ngenes) {
    sdx <- Reduce("*", lapply(exprList, function(x) matrixStats::rowSds(x)))
    ii <- head(order(-sdx), ngenes)
    exprList <- lapply(exprList, function(x) x[ii, , drop = FALSE])
  }

  if (addCombined) {
    exprList[["Combined"]] <- do.call(cbind, exprList)
  }

  exprsamples <- unlist(lapply(exprList, colnames))
  if (!all(exprsamples %in% rownames(phenoData))) {
    stop("samples mismatch for exprList and phenoData")
  }

  multiExpr <- WGCNA::list2multiData(lapply(exprList, Matrix::t))
  cor <- WGCNA::cor

  if (!is.null(power) && length(power) == 1) {
    power <- rep(power, length(multiExpr))
  }

  # module detection procedure
  layers <- list()
  for (i in 1:length(multiExpr)) {
    X <- Matrix::t(multiExpr[[i]]$data)
    message("[runConsensusWGCNA] Computing WGCNA for ", names(multiExpr)[i])
    layers[[names(multiExpr)[i]]] <- computeWGCNA(
      X = X,
      samples = phenoData,
      contrasts = contrasts,
      ngenes = ngenes,
      power = power[i],
      minmodsize = minModuleSize,
      calcMethod = calcMethod,
      deepsplit = deepSplit,
      mergeCutHeight = mergeCutHeight,
      numericlabels = FALSE,
      minKME = minKME,
      maxBlockSize = maxBlockSize,
      compute.stats = compute.stats,
      sv.tom = 40,
      verbose = verbose
    )
  }

  # Consensus module detection
  message("[runConsensusWGCNA] >>> computing CONSENSUS modules...")
  consPower <- unlist(sapply(layers, function(w) w$net$power))
  if (is.null(consPower) && !is.null(power)) consPower <- power
  if (is.null(consPower)) consPower <- rep(12, length(layers))

  sel <- setdiff(names(multiExpr), c("Combined"))

  cons <- WGCNA::blockwiseConsensusModules(
    multiExpr[sel],
    power = as.numeric(consPower),
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = as.integer(minModuleSize),
    deepSplit = as.integer(deepSplit),
    mergeCutHeight = as.numeric(mergeCutHeight),
    numericLabels = FALSE,
    minKMEtoStay = as.numeric(minKME),
    maxBlockSize = as.integer(maxBlockSize),
    saveTOMs = FALSE,
    useDiskCache = FALSE,
    verbose = verbose
  )

  cons$power <- consPower

  ## create and match colors
  for (i in 1:length(layers)) {
    layers[[i]] <- matchColors(layers[[i]], cons$colors)
  }
  layers.colors <- sapply(layers, function(r) r$net$colors)
  colors <- cbind(Consensus = cons$colors, layers.colors)

  ## add labels to dendrogram
  for (i in 1:length(cons$dendrograms)) {
    ii <- which(cons$goodGenes & cons$blocks == i)
    cons$dendrograms[[i]]$labels <- names(cons$colors)[ii]
  }

  ## merge dendrograms
  message("[runConsensusWGCNA] merge_block_dendrograms...")
  multiX <- Matrix::t(do.call(rbind, lapply(exprList, function(x) scale(t(x)))))
  merged <- try(merge_block_dendrograms(cons, multiX))
  if (!inherits(merged, "try-error")) {
    cons$merged_dendro <- merged
  } else {
    cons$merged_dendro <- NULL
  }

  ## create module-trait matrices for each set
  message("[runConsensusWGCNA] >>> computing module-traits matrices...")
  epm <- expandPhenoMatrix(phenoData, drop.ref = drop.ref, keep.numeric = TRUE)
  datTraits <- 1 * epm  

  if (!is.null(contrasts)) {
    message("[runConsensusWGCNA] adding contrasts to datTraits")
    ctx <- makeContrastsFromLabelMatrix(contrasts)
    ctx <- sign(ctx)
    ctx[ctx == 0] <- NA
    ctx[ctx == -1] <- 0
    datTraits <- cbind(datTraits, ctx)
  }

  zlist <- list()
  for (k in names(cons$multiME)) {
    M <- (cons$multiME[[k]][[1]])
    Z <- datTraits
    kk <- intersect(rownames(M), rownames(Z))
    zrho <- cor(M[kk, ], Z[kk, ], use = "pairwise")
    zrho[is.na(zrho)] <- 0
    zlist[[k]] <- zrho
  }

  ## create consensus module-trait matrix
  ydim <- sapply(exprList, ncol)
  consZ <- computeConsensusMatrix(zlist, ydim = ydim, psig = cons.psig)
  avgZ <- Reduce("+", zlist) / length(zlist)
  
  ## add slots
  datExpr <- lapply(exprList, Matrix::t)

  res <- list(
    net = cons,
    layers = layers,
    datExpr = datExpr,
    datTraits = datTraits,
    modTraits = avgZ,
    consModTraits = consZ,
    dendro = cons$merged_dendro,
    colors = colors,
    zlist = zlist,
    ydim = ydim,
    class = "consensus"
  )

  ## run stats
  if (compute.stats) {
    message("[runConsensusWGCNA] >>> computing consensus stats...")
    res$stats <- computeConsensusGeneStats(res)
  }

  ## run enrichment
  if (compute.enrichment && !is.null(GMT)) {
    message("[runConsensusWGCNA] >>> computing module enrichment...")
    if (!is.null(annot)) GMT <- rename_by2(GMT, annot, "symbol")
    res$gsea <- computeConsensusModuleEnrichment(
      res,
      GMT = GMT,
      method = gset.methods,
      annot = annot,
      min.genes = gsea.mingenes,
      ntop = gsea.ntop
    )
  }

  if (summary) {
    message("Annotating modules using ", ai_model)
    ai <- describeModules(
      res,
      multi = FALSE,
      ntop = 50,
      model = ai_model,
      annot = annot,
      experiment = experiment,
      verbose = FALSE
    )
    res$summary <- ai$answers
    res$prompts <- ai$questions
  }

  return(res)

}


#' Create aligned consensus WGCNA layers
#' @param exprList List of expression matrices.
#' @param samples Sample metadata data frame.
#' @param contrasts Optional contrast matrix.
#' @param ngenes Maximum number of genes to use.
#' @param power Soft-thresholding power.
#' @param minModuleSize Minimum module size.
#' @param deepSplit Deep split sensitivity level.
#' @param mergeCutHeight Cut height for merging modules.
#' @param minKME Minimum KME to stay in module.
#' @param maxBlockSize Maximum block size for computation.
#' @param prefix Layer name prefixes.
#' @param verbose Verbosity level.
#' @return List of aligned WGCNA layer results.
#' @export
createConsensusLayers <- function(exprList,
                                        samples,
                                        contrasts = NULL,
                                        ngenes = 2000,
                                        power = 12,
                                        minModuleSize = 20,
                                        deepSplit = 2,
                                        mergeCutHeight = 0.15,
                                        minKME = 0.3,
                                        maxBlockSize = 9999,
                                        prefix = NULL,
                                        verbose = 1) {

  
  if (is.null(prefix)) prefix <- names(exprList)
  nx <- length(exprList)
  prefix <- head(rep(prefix, nx), nx)

  message("[computeConsensusLayers] Aligning matrices...")
  gg <- Reduce(intersect, lapply(exprList, rownames))
  exprList <- lapply(exprList, function(x) x[gg, ])

  if (length(gg) > ngenes) {
    message("[computeConsensusLayers] Reducing to ", ngenes, " genes")
    sdx <- Reduce("*", lapply(exprList, function(x) matrixStats::rowSds(x)))
    ii <- head(order(-sdx), ngenes)
    exprList <- lapply(exprList, function(x) x[ii, ])
  }

  multiExpr <- WGCNA::list2multiData(lapply(exprList, Matrix::t))

  ## determine power vector
  if (is.null(power) || any(is.na(power))) power <- "sft"
  if (as.character(power[1]) %in% c("sft", "iqr")) {
    ## Estimate best power
    power <- power[1]
    message("[createConsensusLayers] optimal power method = ", power)
    est.power <- rep(NA, length(exprList))
    for (i in 1:length(exprList)) {
      p <- pickSoftThreshold(
        Matrix::t(exprList[[i]]),
        sft = NULL, rcut = 0.85, powers = NULL,
        method = power, nmax = 1000, verbose = 0
      )
      if (length(p) == 0 || is.null(p)) p <- NA
      est.power[i] <- p
    }
    power <- ifelse(is.na(est.power), 12, est.power)
  } else {
    power <- as.numeric(power)
  }
  nw <- length(exprList)
  power <- head(rep(power, nw), nw)
  names(power) <- names(exprList)

  message("[computeConsensusLayers] Computing consensus modules...")
  cons <- WGCNA::blockwiseConsensusModules(
    multiExpr,
    power = as.numeric(power),
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = as.integer(minModuleSize),
    deepSplit = as.integer(deepSplit),
    mergeCutHeight = as.numeric(mergeCutHeight),
    numericLabels = FALSE,
    minKMEtoStay = as.numeric(minKME),
    maxBlockSize = as.integer(maxBlockSize),
    saveTOMs = FALSE,
    useDiskCache = FALSE,
    verbose = verbose
  )

  message("[computeConsensusLayers] Creating consensus layers...")
  aligned <- list()
  for (i in 1:length(exprList)) {
    k <- names(exprList)[i]
    sel <- c(
      "colors", "unmergedColors", "goodSamples", "goodGenes",
      "dendrograms", "blockGenes", "blocks"
    )
    net <- cons[sel]
    net$power <- power[i]
    X <- exprList[[i]]
    w <- computeWGCNA(
      X = exprList[[i]],
      samples = samples,
      contrasts = contrasts,
      prefix = prefix[i],
      ngenes = -1,
      net = net,
      calcMethod = "fast",
      sv.tom = 0
    )
    aligned[[k]] <- w
  }

  return(aligned)

}

#' Match module colors to reference labels
#' @param wgcna WGCNA result object.
#' @param refcolors Reference color assignments.
#' @return Updated WGCNA object with matched colors.
#' @export
#' @keywords internal
matchColors <- function(wgcna, refcolors) {

  oldcolors <- wgcna$net$colors
  newcolors <- WGCNA::matchLabels(oldcolors, refcolors)
  lut <- table(oldcolors, newcolors)
  old2new <- colnames(lut)[max.col(lut)]
  names(old2new) <- rownames(lut)
  prefix <- substring(names(wgcna$me.colors),1,2)[1]
  old2newME <- paste0(prefix,old2new)
  names(old2newME) <- paste0(prefix,names(old2new))
  old2new <- c(old2new, old2newME)

  newcol <- function(x) {
    array(old2new[x], dimnames = list(names(x)))
  }

  ## rename everything in net object
  wgcna$net$colors <- newcol(wgcna$net$colors)
  if ("labels" %in% names(wgcna$net)) {
    wgcna$net$labels <- newcol(wgcna$net$labels)
  }

  ## rename unmergedColors
  if ("unmergedColors" %in% names(wgcna$net)) {
    wgcna$net$unmergedColors <- newcol(wgcna$net$unmergedColors)
  }
  names(wgcna$net$MEs) <- newcol(names(wgcna$net$MEs))

  ## rename everything in wgcna object
  names(wgcna$me.genes) <- newcol(names(wgcna$me.genes))
  wgcna$me.colors <- newcol(wgcna$me.colors)
  names(wgcna$me.colors) <- newcol(names(wgcna$me.colors))
  colnames(wgcna$W) <- newcol(colnames(wgcna$W))
  if (!is.null(wgcna$modTraits)) {
    rownames(wgcna$modTraits) <- newcol(rownames(wgcna$modTraits))
  }

  ## rename everything in stats object
  if ("stats" %in% names(wgcna) && !is.null(wgcna$stats)) {
    rownames(wgcna$stats[['moduleTraitCor']]) <- newcol(rownames(wgcna$stats[['moduleTraitCor']]))
    rownames(wgcna$stats[['moduleTraitPvalue']]) <- newcol(rownames(wgcna$stats[['moduleTraitPvalue']]))
    colnames(wgcna$stats[['moduleMembership']]) <- newcol(colnames(wgcna$stats[['moduleMembership']]))
    colnames(wgcna$stats[['MMPvalue']]) <- newcol(colnames(wgcna$stats[['MMPvalue']]))
  }

  return(wgcna)

}

#' Get top correlated genes across modules
#' @param wgcna WGCNA result object or list.
#' @param ref Reference layer name.
#' @param ngenes Number of top genes per module.
#' @param multi Use multi-layer mode.
#' @param modules Subset of modules to query.
#' @return List of data frames with cross-module genes.
#' @export
#' @keywords internal
getModuleCrossGenes <-  function(wgcna,
                                 ref = NULL,
                                 ngenes = 100,
                                 multi = TRUE,
                                 modules = NULL) {


  if (!multi) {
    wgcna <- list(gx = wgcna)
    ref <- 'gx'
  }

  if (is.null(ref)) ref <- head(intersect(names(wgcna),c("gx","px")),1)
  if (is.null(ref) || !ref %in% names(wgcna)) ref <- names(wgcna)[1]

  W <- wgcna[[ref]]
  geneX <- W$datExpr

  MEx <- sapply(wgcna, function(w) as.matrix(w$net$MEs))
  MEx <- do.call(cbind, MEx)

  if (!is.null(modules)) {
    modules <- intersect(modules, colnames(MEx))
    MEx <- MEx[,modules,drop=FALSE]
  }

  nbx.cor <- cor(geneX, MEx)

  nbx.list <- list()

  for (k in colnames(nbx.cor)) {
    ii <- head(order(-nbx.cor[,k]), ngenes)
    rho <- nbx.cor[ii,k]
    gene <- rownames(nbx.cor)[ii]
    me <- W$net$labels[gene]
    nbx.list[[k]] <- data.frame( gene = gene, rho = rho, module = me)
  }

  return(nbx.list)

}

#' Compute consensus matrix from list of matrices. The consensus
#' matrix checks for consistent sign and minimal threshold for each
#' matrix. Optionally filters on consistent p-value.
#' @param ydim original dimension of data
#' @export
#' @keywords internal
computeConsensusMatrix <- function(matlist,
                                   ydim,
                                   psig = 0.05,
                                   consfun = "min") {

  if (length(ydim) == 1) ydim <- rep(ydim[1], length(matlist))

  pv <- mapply(function(z, n) {
    WGCNA::corPvalueStudent(z, n)
  }, matlist, ydim, SIMPLIFY = FALSE)

  for (i in 1:length(pv)) pv[[i]][is.na(pv[[i]])] <- 1 ## missing???

  ## create consensus module-trait matrix
  matsign <- list()
  for (i in 1:length(matlist)) {
    matsign[[i]] <- sign(matlist[[i]]) * (pv[[i]] <= psig)
  }
  matsign <- lapply(matsign, function(x) { x[is.na(x)] <- 0; x})

  all.pos <- Reduce("*", lapply(matsign, function(z) (z >= 0)))
  all.neg <- Reduce("*", lapply(matsign, function(z) (z <= 0)))
  concordant <- (all.pos | all.neg)

  matlistN <- Reduce("+", lapply(matlist, function(x) !is.na(x)))
  matlist0 <- lapply(matlist, function(x) { x[is.na(x)] <- 0; x})

  zsign <- sign(Reduce("+", matsign)) ## mean sign??
  if (consfun == "min") {
    pminFUN <- function(...) pmin(..., na.rm = TRUE)
    consZ <- Reduce(pminFUN, lapply(matlist, abs)) * zsign
  } else if (consfun == "gmean") {
    ## geometric mean
    matlistG <- lapply(matlist, function(x) {
      x <- log(abs(x))
      x[is.na(x)] <- 0
      x
    })
    consZ <- exp(Reduce("+", matlistG) / matlistN)
    consZ <- consZ * zsign
  } else {
    ## mean
    consZ <- Reduce("+", matlist0) / matlistN
  }
  consZ[!concordant] <- NA

  if (psig < 1) {
    ## enforce strong consensus. All layers must be strictly significant.
    all.sig <- Reduce("*", lapply(pv, function(p) 1 * (p <= psig)))
    consZ[!all.sig] <- NA
  }

  return(consZ)

}

#' Compute consensus enrichment by calculating overlapping enriched terms.
computeConsensusModuleEnrichment <- function(cons,
                                             GMT,
                                             annot,
                                             methods = c("fisher","gsetcor","xcor"),
                                             min.genes = 3,
                                             ntop = 400) {

  if (is.null(GMT)) {
    message("ERROR: must provide GMT")
    return(NULL)
  }

  gseaX <- list()
  for (i in 1:length(cons$datExpr)) {

    geneX <- t(cons$datExpr[[i]])

    if (!is.null(annot)) {
      geneX <- rename_by2(geneX, annot, "symbol")
      GMT   <- rename_by2(GMT, annot, "symbol")
    }

    symbols <- intersect(rownames(GMT), rownames(geneX))
    if (length(symbols) == 0) {
      message("[computeConsensusModuleEnrichment] ERROR. No symbol overlap.")
      return(NULL)
    }
    geneX <- geneX[symbols, ]
    GMT <- GMT[symbols, ]

    ## select on minimum gene sets size
    sel <- which(Matrix::colSums(GMT!=0) >= min.genes)
    GMT <- GMT[, sel]

    ## Create extended Eigengene matrix (ME). ME should be nicely
    ## normalized/scaled so we just rbind across datasets
    ME <- cons$net$multiMEs[[i]]$data

    ## get genes in modules
    me.genes <- tapply(names(cons$net$colors), cons$net$colors, list)
    names(me.genes) <- paste0("ME",names(me.genes))

    if (!is.null(annot)) {
      me.genes <- lapply(me.genes, function(gg) probe2symbol(gg, annot))
    }
    me.genes <- lapply(me.genes, function(g) intersect(g, symbols))
    colnames(geneX) <- rownames(ME)

    gseaX[[names(cons$datExpr)[i]]] <- run_enrichment_methods(
      ME,
      me.genes = me.genes,
      GMT= GMT,
      geneX = geneX,
      methods = methods,
      min.genes = min.genes,
      ntop = ntop
    )
    
  }

  cons.gsea <- list()
  for (m in names(gseaX[[1]])) {
    xx <- lapply( gseaX, function(g) g[[m]] )
    sel <- Reduce(intersect, lapply(xx, rownames))
    if (length(sel) > 0) {
      if (length(sel)==1) sel <- c(sel,sel) ## length==1 crashes...
      xx <- lapply(xx, function(x) x[sel,,drop=FALSE] )
      xx.score <- sapply(xx, function(x) x[,"score"])
      colnames(xx.score) <- paste0("score.",colnames(xx.score))

      xx.pvalue <- lapply(xx, function(x) x[,grep("^p",colnames(x))])
      xx.pvalue <- do.call(cbind, xx.pvalue)

      m.score <- rowMeans(xx.score,na.rm=TRUE)
      m.pvalue  <- apply(sapply(xx, function(x) x[,"p.value"]), 1, max, na.rm=TRUE)
      m.qvalue  <- p.adjust( m.pvalue )
      df <- data.frame(
        module = xx[[1]]$module,
        geneset = xx[[1]]$geneset,
        score = m.score,
        xx.score,
        p.value = m.pvalue,
        q.value = m.qvalue,
        overlap = xx[[1]]$overlap,
        genes = xx[[1]]$genes,
        xx.pvalue
      )
      df <- df[order(df$p.value),]
      cons.gsea[[m]] <- df
    }
  }

  return(cons.gsea)

}
