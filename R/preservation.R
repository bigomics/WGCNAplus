#' Run module preservation analysis across datasets
#' @param exprList Named list of expression matrices.
#' @param phenoData Sample phenotype data frame.
#' @param contrasts Contrast definitions or NULL.
#' @param power Soft-thresholding power.
#' @param reference Reference network index or name.
#' @param add.merged Logical; add merged expression layer.
#' @param ngenes Number of top variable genes.
#' @param minModuleSize Minimum module size.
#' @param deepSplit Tree cut sensitivity parameter.
#' @param annot Gene annotation data frame or NULL.
#' @param compute.stats Logical; compute gene statistics.
#' @param compute.enrichment Logical; compute gene-set enrichment.
#' @param GMT Gene-set matrix or NULL.
#' @param gset.methods Enrichment methods to use.
#' @return Preservation result object with Z-summary.
runPreservationWGCNA <- function(exprList,
                                 phenoData,
                                 contrasts = NULL,
                                 power = 12,
                                 reference = 1,
                                 add.merged = FALSE,
                                 ngenes = 2000,
                                 minModuleSize = 20,
                                 deepSplit = 2,
                                 annot = NULL,
                                 compute.stats = TRUE,
                                 compute.enrichment = TRUE,
                                 GMT = NULL,
                                 gset.methods = c("fisher", "gsetcor", "xcor"),
                                 ai_model = "",
                                 summary = FALSE) {

  if (is.character(reference)) {
    reference <- match(reference, names(exprList))
  }

  if (reference > 0) {
    reference.name <- names(exprList)[reference]
  } else {
    reference.name <- "Consensus"
  }

  ## multiset WGCNA
  pres <- runConsensusWGCNA(
    exprList,
    phenoData = phenoData,
    contrasts = contrasts,
    GMT = NULL, ## no enrichment now
    annot = annot, 
    ngenes = ngenes,
    power = power,
    minModuleSize = minModuleSize,
    minKME = 0.3,
    mergeCutHeight = 0.15,
    deepSplit = deepSplit,
    maxBlockSize = 9999,
    addCombined = FALSE,
    calcMethod = "fast",
    drop.ref = FALSE,
    compute.stats = FALSE,
    compute.enrichment = FALSE,
    gsea.mingenes = 10,
    gset.methods = NULL,
    ai_model = ai_model,
    summary = summary
  )

  colorList <- lapply(pres$layers, function(w) w$net$colors)
  names(colorList) <- names(pres$layers)
  exprList <- lapply(pres$layers, function(w) t(w$datExpr))

  if (add.merged || reference == 0) {
    message("[runPreservationWGCNA] adding merged layer...")
    cX <- lapply(exprList, function(x) x - rowMeans(x))
    merged <- do.call(cbind, cX)
    exprList$Merged <- NULL
    exprList <- c(list(Merged = merged), exprList)
    cons.colors <- pres$net$colors
    colorList <- c(list(Consensus = cons.colors), colorList)
    reference <- reference + 1
  }

  message("[runPreservationWGCNA] running WGCNA::modulePreservation...")
  multiExpr <- WGCNA::list2multiData(lapply(exprList, Matrix::t))
  mp <- WGCNA::modulePreservation(
    multiExpr,
    colorList,
    referenceNetworks = reference,
    nPermutations = 10,
    networkType = "signed",
    quickCor = 0,
    verbose = 2,
    indent = 0
  )

  ## Zsummary tables
  mp.tables <- mp$preservation$Z[[1]][-reference]
  Z <- sapply(mp.tables, function(mat) mat[, "Zsummary.pres"])
  rownames(Z) <- rownames(mp.tables[[1]])
  rownames(Z) <- paste0("ME", rownames(Z))
  colnames(Z) <- names(multiExpr)[-reference]

  ## median rank
  mp.tables <- mp$preservation$observed[[1]][-reference]
  M <- sapply(mp.tables, function(mat) mat[, "medianRank.pres"])
  rownames(M) <- rownames(mp.tables[[1]])
  rownames(M) <- paste0("ME", rownames(M))
  colnames(M) <- names(multiExpr)[-reference]

  ## module size
  moduleSize <- mp.tables[[1]][, "moduleSize"]
  names(moduleSize) <- rownames(Z)

  ## module-traits. We need to recompute the MEs (module eigengenes)
  ## using the color coding of the reference set.
  refColors <- colorList[[1]]
  MEx <- lapply(exprList, function(x) {
    WGCNA::moduleEigengenes(t(x), colors = refColors)$eigengenes
  })

  ## Compute module-trait correlation matrices
  Y <- lapply(pres$layers, function(w) w$datTraits)
  if ("Merged" %in% names(MEx) && !"Merged" %in% names(Y)) {
    kk <- rownames(MEx[["Merged"]])
    Y[["Merged"]] <- pres$datTraits[kk, ]
    Y <- Y[names(MEx)]
  }

  kk <- Reduce(union, lapply(Y, colnames))
  Y <- lapply(Y, function(y) y[, match(kk, colnames(y)), drop = FALSE])
  for (i in 1:length(Y)) colnames(Y[[i]]) <- kk

  R <- mapply(cor, MEx, Y, use = "pairwise", SIMPLIFY = FALSE)

  ## gene statistics of reference layer
  if (compute.stats) {
    message("[runPreservationWGCNA] computing gene statistics...")
    ref <- reference.name
    wnet <- list(MEs = MEx[[ref]], colors = pres$colors[, ref])
    pres$stats <- computeGeneStats(wnet, pres$datExpr[[ref]], pres$datTraits, TOM = NULL)
  }

  ## geneset enrichment of reference layer
  if (compute.enrichment && !is.null(GMT)) {
    message("[runPreservationWGCNA] computing geneset enrichment...")
    if(!is.null(annot)) GMT <- rename_by2(GMT, annot, "human_ortholog")
    pres$gsea <- computeModuleEnrichment(
      pres$layers[[ref]],
      GMT = GMT,
      annot = annot,
      methods = gset.methods,
      ntop = 1000,
      xtop = 100,
      filter = NULL,
      add.wgcna = FALSE
    )
  }

  pres$modulePreservation <- mp
  pres$Zsummary <- Z
  pres$medianRank <- M
  pres$moduleSize <- moduleSize
  pres$modTraits <- R
  pres$MEs <- MEx

  return(pres)

}
