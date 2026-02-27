#' Run multi-omics WGCNA analysis
#' @param dataX Named list of data matrices.
#' @param samples Sample metadata data frame.
#' @param contrasts Contrast matrix or NULL.
#' @param power Soft-thresholding power.
#' @param ngenes Number of top variable genes.
#' @param datanames Names for data layers.
#' @param clustMethod Clustering method name.
#' @param cutMethod Tree cutting method name.
#' @param minmodsize Minimum module size.
#' @param minKME Minimum module membership.
#' @param deepsplit Deep split sensitivity level.
#' @param mergeCutHeight Module merge cut height.
#' @param compute.enrichment Compute gene set enrichment.
#' @param xref Cross-reference layer names.
#' @param annot Annotation table or NULL.
#' @param GMT Gene set matrix or NULL.
#' @param drop.ref Drop reference level.
#' @param add.pheno Add phenotype matrix.
#' @param add.gsets Add gene set scores.
#' @param do.consensus Use consensus WGCNA.
#' @param gset.methods Enrichment methods vector.
#' @param gset.ntop Top gene sets count.
#' @param gset.xtop Top cross-omics count.
#' @param report Generate report.
#' @param ai_model LLM model name or NULL.
#' @param experiment Experiment description string.
#' @param verbose Verbosity level.
#' @param progress Optional Shiny progress object.
#' @return List with layers, report, and settings.
#' @export
computeWGCNA_multiomics <- function(dataX,
                                    samples,
                                    contrasts = NULL,
                                    power = 12,
                                    ngenes = 2000,
                                    datanames = NULL,
                                    clustMethod = "average",
                                    cutMethod = "hybrid",
                                    minmodsize = 10,
                                    minKME = 0.3,
                                    deepsplit = 2,
                                    mergeCutHeight = 0.15,
                                    compute.enrichment = TRUE,
                                    xref = NULL,
                                    annot = NULL,
                                    GMT = NULL,
                                    drop.ref = FALSE,
                                    add.pheno = FALSE,
                                    add.gsets = FALSE,
                                    do.consensus = FALSE,
                                    gset.methods = c("fisher", "gsetcor", "xcor"),
                                    gset.ntop = 1000,
                                    gset.xtop = 100,
                                    report = TRUE,
                                    ai_model = getOption("WGCNAplus.default_llm"),
                                    experiment = "",
                                    verbose = 1,
                                    progress = NULL
                                    ) {


  if (inherits(dataX,"matrix")) {
    dataX <- mofa.split_data(dataX, keep.prefix=FALSE)
  }

  if (!is.null(annot)) {
    dataX <- lapply(dataX, function(x) rename_by2(x, annot, "symbol"))
  }

  compute.enrichment <- (compute.enrichment && !is.null(GMT))
  if (!is.null(annot) && !is.null(GMT)) {
    GMT <- rename_by2(GMT, annot, "symbol")
  }


  ## add pheno matrix??
  if (add.gsets && !is.null(GMT)) {
    X <- mofa.merge_data2(dataX, merge.rows="union")
    if (!is.null(annot)) X <- rename_by2(X, annot, "symbol")
    kk <- intersect(rownames(X), rownames(GMT))
    if (length(kk) == 0) {
      message("Error: X and GMT do not share features")
    }
    if (length(kk)) {
      if (!requireNamespace("plaid", quietly = TRUE)) {
        stop("Package 'plaid' is required for gene set integration")
      }
      gsetX <- plaid::plaid(X[kk,], GMT[kk,])
      dataX$gs <- gsetX
    }
  }

  ## add phenomatrix??
  if (add.pheno) {
    phenoX <- expandPhenoMatrix(samples, keep.numeric = TRUE, drop.ref = drop.ref)
    dataX$ph <- t(phenoX)
  }

  dt.na <- which(unlist(lapply(dataX, function(x) sum(is.na(x)))) > 0)
  if (any(dt.na)) {
    dataX[dt.na] <- lapply(dataX[dt.na], imputeMissing, method = "SVD2")
  }

  if (is.null(power)) power <- NA
  nw <- length(dataX)
  if (length(power)<nw) power <- head(rep(power, nw),nw)
  names(power) <- names(dataX)
  if (any(is.na(power))) power[is.na(power)] <- "sft"

  if (length(minKME)<nw) minKME <- head(rep(minKME, nw),nw)
  if (length(deepsplit)<nw) deepsplit <- head(rep(deepsplit, nw),nw)
  names(minKME) <- names(dataX)
  names(deepsplit) <- names(dataX)

  if (any(as.character(power) %in% c("sft","iqr"))) {
    ii <- which(as.character(power) %in% c("sft","iqr"))
    message("[compute_multiomics] estimating power with method = ", power[ii])
    for (i in ii) {
      p <- pickSoftThreshold(
        Matrix::t(dataX[[i]]), sft=NULL, rcut=0.85, powers = NULL,
        method=power[i], nmax=1000, verbose=1)
      if (length(p)==0 || is.null(p) ) p <- NA
      power[i] <- p
    }
    power <- ifelse (is.na(power), 12, power)
    message("[compute_multiomics] estimated powers = ", power)
  }
  power <- as.numeric(power)
  names(power) <- names(dataX)
  
  ## This runs WGCNA on an expression list.
  layers <- list()
  has.gxpx <- all(c("gx","px") %in% names(dataX))
  if (do.consensus && has.gxpx) {
    cat("[compute_multiomics] computing WGCNA consensus layers for GX+PX \n")
    nn <- mean(rownames(dataX[["gx"]]) %in% rownames(dataX[["px"]]))
    if (nn < 0.10) {
      message("[compute_multiomics] ERROR: gx and px features do not overlap")
    } else {
      layers <- createConsensusLayers(
        dataX[c('gx','px')],
        samples = samples,
        contrasts = contrasts,
        prefix = c("GX", "PX"),
        ngenes = ngenes,
        power = power[1],
        minModuleSize = minmodsize,
        deepSplit = deepsplit[1],
        mergeCutHeight = mergeCutHeight,
        minKME = minKME[c('gx','px')],
        maxBlockSize = 9999,
        verbose = 1
      )
    }
  }

  dtlist <- setdiff(names(dataX), names(layers))
  for (dt in dtlist) {
    cat("[compute_multiomics] computing WGCNA for", dt, "-------------\n")
    minkme1 <- ifelse(dt=='ph', 0, minKME[dt])
    minmodsize <- ifelse(dt=='ph', 1, minmodsize)
    layers[[dt]] <- computeWGCNA(
      X = dataX[[dt]],
      samples = samples,
      contrasts = contrasts,
      ngenes = ngenes,
      calcMethod = "fast",
      power = power[dt],
      lowrank = 40,
      clustMethod = clustMethod,
      cutMethod = cutMethod,
      deepsplit = deepsplit[dt],
      minKME = minkme1,
      minmodsize = minmodsize,
      mergeCutHeight = mergeCutHeight,
      compute.stats = TRUE,
      sv.tom = 40,
      prefix = toupper(dt),
      drop.ref = drop.ref,
      is.multiomics = FALSE,
      verbose = verbose
    )
  }

  layers <- layers[names(dataX)]

  ## get members
  me.genes <- lapply(layers, function(m) m$me.genes)
  names(me.genes) <- NULL
  me.genes <- unlist(me.genes, recursive=FALSE)

  ## get colors
  me.colors <- lapply(names(me.genes), function(m) rep(m,length(me.genes[[m]])))
  me.colors <- unlist(me.colors)
  names(me.colors) <- unlist(me.genes,use.names=FALSE)
  mm <- lapply(layers,function(m) unique(m$net$labels))
  dt <- unlist(lapply(names(mm),function(i) rep(i,length(mm[[i]]))))
  names(dt) <- unlist(mm)
  names(me.colors) <- paste0(dt[me.colors],":",names(me.colors))
  
  ## Compute enrichment
  gsea <- NULL
  if (compute.enrichment) {
    message("[compute_multiomics] computing module enrichment...")
    gsea <- computeModuleEnrichment(
      wgcna = layers,
      multi = TRUE,
      methods = gset.methods,
      ntop = gset.ntop,
      xtop = gset.xtop,
      xref = xref,
      annot = annot,
      GMT = GMT,
      filter = NULL,
      add.wgcna = FALSE
    )

    ## split up results?? still needed in old formats
    for (k in names(layers)) {
      mm <- names(layers[[k]]$me.genes)
      gg <- gsea[mm]
      names(gg) <- mm
      layers[[k]]$gsea <- gg
    }

  }

  lasagna.model <- NULL
  lasagna.graph <- NULL
  do.lasagna = TRUE
  if (do.lasagna) {
    dbg("[compute_multiomics] >>> creating lasagna ")

    ## Get eigengene matrices, remove grey modules
    ww <- lapply(layers, function(w) t(w$net$MEs))
    ww <- lapply(ww, function(w) w[!grepl("[A-Z]{2}grey$", rownames(w)), , drop=FALSE])
    ww <- ww[which(sapply(ww,nrow)>0)]

    datTraits <- layers[[1]]$datTraits
    gdata <- list(X = ww, samples = datTraits)

    ## Create lasagna model
    lasagna.model <- lasagna.create_model(
      gdata,
      pheno = "expanded",
      ntop = 2000,
      nc = 20,
      add.sink = FALSE,
      intra = TRUE,
      fully_connect = FALSE,
      add.revpheno = TRUE
    )

    dbg("[compute_multiomics] conditioning... ")
    ## Multi-condition edge weighting
    lasagna.graph <- lasagna.multisolve(
      lasagna.model,
      min_rho = 0.1,
      max_edges = 1000,
      fc.weight = TRUE,
      sp.weight = FALSE,
      prune = FALSE
    )

  }

  report.out <- NULL
  if (report) {
    dbg("[compute_multiomics] >>> creating report")
    ## Create summaries of each module.
    if (!is.null(ai_model)) message("Creating report using ", ai_model)
    if (is.null(ai_model)||ai_model=="") message("Creating dummy report")
    report.out <- create_report(
      layers, ai_model,
      annot = annot,
      multi = TRUE,
      graph = lasagna.graph,
      topratio = 0.85,
      psig = 0.05,
      verbose = 1
    )
  }

  ## get some settings
  power <- sapply(layers, function(a) a$net$power, USE.NAMES=FALSE)
  names(power) <- names(layers)

  dbg("[compute_multiomics] copying settings...")

  ## Like consensus: put datatype wgcna in layers slot.
  settings <- list(
    minmodsize = minmodsize,
    power = power,
    mergeCutHeight = mergeCutHeight,
    deepsplit = deepsplit,
    minKME = minKME,
    networktype = "signed",
    tomtype = "signed",
    gset.methods = gset.methods,
    NULL
  )

  dbg("[compute_multiomics] creating out list")

  out <- list(
    layers = layers,
    me.genes = me.genes,
    me.colors = me.colors,
    gsea = gsea,
    report = report.out,
    datanames = datanames,
    lasagna = lasagna.model,
    graph = lasagna.graph,
    experiment = experiment,
    settings = settings,
    class = "multiomics"
  )

  return(out)

}


## ----------------------------------------------------
## Perform geneset analysis on modules
## ----------------------------------------------------
#' Merge multiple ME matrices into one. Allow different dimensions.
#' @param mlist List of ME matrices.
#' @param me2 Optional second ME matrix.
#' @param prefix Prefix columns with list names.
#' @return Merged eigengene matrix.
#' @keywords internal
mergeME <- function(mlist,
                    me2 = NULL,
                    prefix = FALSE) {

  if (!is.null(me2) && !inherits(mlist,"list")) {
    mlist <- list(mlist, me2)
  }

  all.samples <- unique(unlist(sapply(mlist, rownames, simplify=FALSE)))

  if (prefix) {
    for (i in 1:length(mlist)) {
      colnames(mlist[[i]]) <- paste0(names(mlist)[i],":",colnames(mlist[[i]]))
    }
  }

  is.mat <- all(sapply(mlist, inherits, what="matrix"))
  all.me <- unique(unlist(sapply(mlist, colnames, simplify=FALSE)))

  M <- as.data.frame(matrix(NA, nrow=length(all.samples), ncol=length(all.me)))
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
