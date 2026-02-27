#' Compute enrichment of each WGCNA module using various
#' methods. Handles single-type and multi-omics WGCNA objects.
#' @export
computeModuleEnrichment <- function(wgcna,
                                    GMT,
                                    multi = FALSE,
                                    methods = c("fisher", "gsetcor", "xcor"),
                                    ntop = 200,
                                    xtop = 100,
                                    annot = NULL,
                                    xref = NULL,
                                    min.genes = 3,
                                    min.rho = 0.8,
                                    filter = NULL,
                                    add.wgcna = TRUE) {
  
  if (!multi) {
    layers <- list(gx = wgcna)
    if (!is.null(annot)) rownames(annot) <- paste0("gx:",rownames(annot))
  } else if (!is.null(wgcna$layers)) {
    layers <- wgcna$layers
  } else {
    layers <- wgcna
  }

  if (is.null(GMT)) {
    message("[computeModuleEnrichment] ERROR: must provide GMT")
    return(NULL)
  }

  ## create fake annotation table if no annotation table is given.
  if (is.null(annot)) {
    gg <- lapply(layers, function(w) colnames(w$datExpr))
    ff <- list()
    for (i in 1:length(gg)) ff[[i]] <- paste0(names(layers)[i],":",gg[[i]])
    annot <- data.frame(feature = unlist(ff), symbol = unlist(gg))
    rownames(annot) <- NULL
  }

  ## make sure GMT features  are in symbols
  symbol.col <- intersect(c("symbol","gene_name"),colnames(annot))[1]
  GMT <- rename_by2(GMT, annot, symbol.col)

  ## add cross-referencing data??
  xref <- intersect(xref, names(layers))
  if (length(xref)) {
    message("[computeModuleEnrichment] NOTE. Adding cross-correlation datatype: ", xref)
  }

  gsea <- list()
  dtype = names(layers)[1]
  for (dtype in names(layers)) {

    ## collapse features to symbol
    sel <- unique(c(dtype, xref))
    datExpr <- lapply(layers[sel], function(w) t(as.matrix(w$datExpr)))
    geneX <- mofa.merge_data2(datExpr, merge.rows="prefix", merge.cols="union")
    geneX <- rename_by2(geneX, annot, symbol.col)

    ## check if overlap exists
    bg <- intersect(rownames(geneX), rownames(GMT))
    if (length(bg) == 0) {
      message("[computeModuleEnrichment] WARNING. no overlapping genes for ", dtype)
      next()
    }
    G1 <- GMT[bg, , drop = FALSE]

    if (!is.null(filter)) {
      sel <- grep(filter, colnames(G1))
      if (length(sel)) G1 <- G1[, sel, drop = FALSE]
    }

    G1 <- G1[, which(Matrix::colSums(G1 != 0) >= min.genes), drop = FALSE]

    if (nrow(G1) >= 3  &&  ncol(G1) >= 3) {
      ## get eigengene members. convert to symbols.
      me.genes <- layers[[dtype]]$me.genes
      me.genes <- lapply(me.genes, function(g) probe2symbol(g, annot, query = symbol.col))
      dt.gsea <- run_enrichment_methods(
        ME = as.matrix(layers[[dtype]]$net$MEs), ## eigengene matrix,
        me.genes = me.genes,
        GMT = G1,
        geneX = geneX,
        methods = methods,
        ntop = ntop,
        xtop = xtop,
        min.rho = min.rho
      )
      
      gsea <- c(gsea, dt.gsea)
    }
  }

  if (add.wgcna) {
    wgcna$gsea <- gsea
    out <- wgcna
  } else {
    out <- gsea
  }

  message("[computeModuleEnrichment] done!")

  return(out)

}


#' Run enrichment analysis methods on modules
#' @param ME Module eigengene matrix.
#' @param me.genes List of gene vectors per module.
#' @param GMT Gene-set membership sparse matrix.
#' @param geneX Gene expression matrix (genes x samples).
#' @param methods Character vector of methods to run.
#' @param ntop Maximum gene sets to return.
#' @param xtop Maximum cross-correlated features to add.
#' @param min.genes Minimum genes per gene set.
#' @param min.rho Minimum correlation for cross-features.
#' @return List of enrichment data frames per module.
#' @keywords internal
run_enrichment_methods <- function(ME,
                                   me.genes,
                                   GMT,
                                   geneX,
                                   methods = c("fisher","gsetcor","xcor"),
                                   ntop = 400,
                                   xtop = 100,
                                   min.genes = 3,
                                   min.rho = 0.8) {

  rho.list <- list()
  pval.list <- list()

  bg <- intersect(rownames(GMT), rownames(geneX))
  GMT <- GMT[bg, ]
  geneX <- geneX[bg, ]

  ## select on minimum genes
  sel <- which(Matrix::colSums(GMT!=0) >= min.genes)
  GMT <- GMT[, sel]

  gsetX <- plaid::plaid(geneX, matG=GMT)
  message("Computing enrichment for ", nrow(gsetX), " genesets")

  ## Add highly cross-correlated genes. limit xtop if geneX is too small.
  if (xtop > 0) {
    message("Adding high cross-correlated features. min.rho = ", min.rho)
    xtop <- min(xtop, round(nrow(geneX)/4))
    nbx.genes <- list()
    for (k in colnames(ME)) {
      ss <- intersect(colnames(geneX),rownames(ME))
      if (length(ss)<3) next()
      gx <- geneX[,ss,drop=FALSE]
      mx <- ME[ss,k,drop=FALSE]
      cx <- cor(t(gx), mx, use="pairwise")[,1]
      cx <- cx[!is.na(cx) & abs(cx) > min.rho]
      if (length(cx)) {
        nbx.genes[[k]] <- head(names(sort(-cx)), xtop)
      } else {
        nbx.genes[[k]] <- c()
      }
    }

    for (k in names(me.genes)) {
      me.genes[[k]] <- unique(c(me.genes[[k]], nbx.genes[[k]]))
    }
  }

  ## Correlate geneset score (averageCLR) with module eigengene (ME).
  ## The aim is to select genesets correlated with the ME.
  if ("gsetcor" %in% methods) {
    message("[run_enrichment_methods] calculating single-sample geneset correlation...")
    rc.rho <- matrix(NA, ncol(GMT), ncol(ME))
    rc.pvalue <- matrix(NA, ncol(GMT), ncol(ME))
    dimnames(rc.rho) <- list(colnames(GMT), colnames(ME))
    dimnames(rc.pvalue) <- list(colnames(GMT), colnames(ME))
    jj <- which(rownames(gsetX) %in% colnames(GMT))
    kk <- intersect(colnames(gsetX), rownames(ME)) ## common samples
    tt <- cortest(t(gsetX[jj,kk]), ME[kk,])
    rho.jj <- tt$rho
    pvalue.jj <- tt$pvalue
    ii <- match(rownames(gsetX)[jj], rownames(rc.rho))
    rc.rho[ii, ] <- rho.jj
    rc.pvalue[ii, ] <- pvalue.jj
    rho.list[["gsetcor"]] <- rc.rho
    pval.list[["gsetcor"]] <- rc.pvalue
  }

  ## Correlate module eigengene (ME) with genes and
  ## then do a gset.rankcor() on the ME correlation.
  if ("xcor" %in% methods) {
    message("[run_enrichment_methods] calculating eigengene correlation...")
    ss <- intersect(colnames(geneX),rownames(ME))
    rho <- cor(t(geneX[,ss,drop=FALSE]), ME[ss,,drop=FALSE], use="pairwise")
    rho[is.na(rho)] <- 0
    rc <- gset.rankcor(rho, GMT, compute.p = TRUE) ## NEEDS CHECK!!!
    rho.list[["xcor"]] <- rc$rho
    pval.list[["xcor"]] <- rc$p.value
  }

  gmt <- mat2gmt(GMT)

  if (1) {
    Pmin <- sapply(pval.list, function(P) apply(P, 1, min))
    sel <- head(order(rowMeans(apply(Pmin, 2, rank))), 5 * ntop)
    message("[run_enrichment_methods] preselecting ", length(sel), " sets for fgsea/Fisher test")
    sel <- rownames(Pmin)[sel]
    gmt <- gmt[sel]
  }

  ## Compute gene correlation to eigengenes ME,
  ## then do pre-ranked enrichment on rho value.
  if ("fgsea" %in% methods) {
    message("[run_enrichment_methods] calculating module fgsea...")
    ss <- intersect(colnames(geneX),rownames(ME))
    xrho <- cor(t(geneX[,ss,drop=FALSE]), ME[ss,,drop=FALSE], use = "pairwise")
    xrho[is.na(xrho)] <- 0
    res <- list()
    for (i in 1:ncol(xrho)) {
      k <- colnames(xrho)[i]
      res[[k]] <- fgsea::fgsea(gmt, xrho[, i]) ## NEEDS CHECK!!!
    }
    pw <- res[[1]]$pathway
    res <- lapply(res, function(r) r[match(pw, r$pathway), ])
    nes <- sapply(res, function(r) r$NES)
    pval <- sapply(res, function(r) r$pval)
    rownames(nes) <- rownames(pval) <- pw
    colnames(nes) <- colnames(pval) <- names(res)
    rho.list[["fgsea"]] <- nes
    pval.list[["fgsea"]] <- pval
  }

  ## Perform fisher-test on (extended) ME genes. The ME genes might
  ## have been extended with most correlated genes.
  if ("fisher" %in% methods) {
    message("[run_enrichment_methods] calculating Fisher tests...")
    rho <- matrix(NA, length(gmt), ncol(ME))
    pval <- matrix(NA, length(gmt), ncol(ME))
    dimnames(rho) <- list(names(gmt), colnames(ME))
    dimnames(pval) <- list(names(gmt), colnames(ME))

    ## perform Fisher test for all modules using the module genes
    for (i in 1:ncol(rho)) {
      k <- colnames(rho)[i]
      gg <- me.genes[[k]]
      rr <- try(gset.fisher(gg, GMT, background = bg, fdr = 1,
        min.genes = -1, verbose = 0, sort.by = "none", no.pass = 1))
      if (!"try-error" %in% class(rr)) {
        rr <- rr[match(rownames(rho), rownames(rr)), ]
        rho[, i] <- rr$odd.ratio
        pval[, i] <- rr$p.value
      }
    }

    ## handle infinite or NA
    rho[is.infinite(rho)] <- 2 * max(rho, na.rm = TRUE) ## Inf odd.ratio
    pval[is.na(pval)] <- 1
    rho[is.na(rho)] <- 0
    
    rho.list[["fisher"]] <- rho
    pval.list[["fisher"]] <- pval
  }

  gsets <- Reduce(intersect, lapply(rho.list, rownames))
  modules <- Reduce(intersect, lapply(rho.list, colnames))
  rho.list <- lapply(rho.list, function(x) x[gsets, modules, drop = FALSE])
  pval.list <- lapply(pval.list, function(x) x[gsets, modules, drop = FALSE])

  ## Compute meta rank and pval. Handle NA for failing methods.
  pvalNA <- lapply(pval.list, function(x) { x[is.na(x)]=0; x })
  meta.p <- Reduce(pmax, pvalNA)
  meta.q <- apply(meta.p, 2, p.adjust, method = "fdr")

  ## NEED RETHINK: how about negative FC???
  rnk.list <- lapply(rho.list, function(x) apply(x, 2, rank, na.last = "keep") / nrow(x))
  meta.rnk <- Reduce("+", rnk.list) / length(rnk.list)
  rnk.NAZERO <- lapply(rnk.list, function(x) { x[is.na(x)]=0; x })
  rnk.NSUM <- Reduce("+", lapply(rnk.list, function(x) !is.na(x)))
  meta.rnk <- Reduce("+", rnk.NAZERO) / rnk.NSUM

  ## create dataframe by module
  message("[run_enrichment_methods] creating dataframes...")
  gse.list <- list()
  for (i in 1:ncol(meta.p)) {
    k <- colnames(meta.p)[i]
    pv <- sapply(pval.list, function(x) x[, i])
    colnames(pv) <- paste0("p.", colnames(pv))
    df <- data.frame(
      module = k,
      geneset = rownames(meta.p),
      score = meta.rnk[, i],
      p.value = meta.p[, i],
      q.value = meta.q[, i],
      pv
    )
    df <- df[order(-abs(df$score)), ]
    df <- head(df, ntop)
    gse.list[[k]] <- df
  }

  ## add genes
  k <- names(gse.list)[1]
  for (k in names(gse.list)) {
    gset <- rownames(gse.list[[k]])
    gg <- me.genes[[k]]
    set.genes <- lapply(gmt[gset], function(s) intersect(s, gg))
    n0 <- sapply(gmt[gset], length)
    n1 <- sapply(set.genes, length)
    set.genes <- sapply(set.genes, function(g) paste(sort(g), collapse = "|"))
    gse.list[[k]]$overlap <- paste0(n1, "/", n0)
    gse.list[[k]]$genes <- set.genes
  }

  message("[run_enrichment_methods] done!")

  return(gse.list)

}
