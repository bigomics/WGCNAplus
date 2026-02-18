#' @title WGCNAplus: Extended WGCNA for Multi-Omics Data Analysis
#'
#' @description Extended Weighted Gene Coexpression Network Analysis (WGCNA)
#' with support for multi-omics, consensus analysis, module preservation,
#' AI-powered module annotation, and comprehensive visualization.
#'
#' @name WGCNAplus-package
#' @aliases WGCNAplus
#' @docType package
#'
#' @import WGCNA
#'
#' @importFrom playbase gx.limma gx.heatmap gset.fisher ai.ask
#'   pgx.scatterPlotXY imputeMissing svdImpute2
#'   normalizeMultiOmics pgx.clusterBigMatrix
#'   lasagna.create_model lasagna.multisolve
#'
#' @importFrom Matrix t crossprod
#' @importFrom matrixStats rowSds colSds rowMaxs colMaxs
#' @importFrom irlba irlba
#' @importFrom igraph graph_from_adjacency_matrix V E
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom fgsea fgsea
#' @importFrom fastcluster hclust
#' @importFrom ggplot2 ggplot aes geom_bar geom_text theme_minimal coord_flip
#'   labs scale_fill_manual element_text element_blank theme
#' @importFrom stats as.dist dist cutree p.adjust cor.test
#'   median quantile setNames
#' @importFrom grDevices colorRampPalette hcl.colors
#' @importFrom graphics barplot par plot image text mtext axis rect abline
#'   layout lines legend box title segments
#' @importFrom utils head tail type.convert
#' @importFrom methods is
#'
"_PACKAGE"

## ---------------------------------------------------------------------------
## Internal utility functions (copied from playbase::utils.R, not exported)
## ---------------------------------------------------------------------------

#' Get VM RSS memory usage
#'
#' @param digits Number of decimal places.
#' @return Memory usage string in MB.
#' @keywords internal
mem.vmrss <- function(digits = 0) {
  mem <- "[? MB]"
  if (Sys.info()["sysname"] %in% c("Linux")) {
    proc <- paste("/proc", Sys.getpid(), "status", sep = "/")
    rss <- gsub("VmRSS:[\t ]+| kB", "", system(paste("grep -i vmrss", proc), intern = TRUE))
    rss <- as.numeric(rss) / (1024) ## MB
    mem <- paste0(round(rss, digits), "MB")
  }
  mem
}

#' Get process virtual memory usage
#'
#' @param digits Number of decimal places.
#' @return Memory usage string in MB.
#' @keywords internal
mem.proc <- function(digits = 0) {
  mem <- "[? MB]"
  if (Sys.info()["sysname"] %in% c("Linux")) {
    file <- paste("/proc", Sys.getpid(), "stat", sep = "/")
    what <- vector("character", 52)
    vsz <- as.numeric(scan(file, what = what, quiet = TRUE)[23])
    vsz <- vsz / (1024**2) ## MB
    mem <- paste0(round(vsz, digits), "MB")
  }
  mem
}

#' Log info message with timestamp
#'
#' @param ... Message parts to concatenate.
#' @param type Log level label string.
#' @return NULL (called for side effect).
#' @keywords internal
info <- function(..., type = "INFO") {
  dd <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- "some message"
  msg <- sapply(list(...), paste, collapse = " ")
  dd <- paste0("[", dd, "]")
  mm <- paste0("[", mem.proc(), "/", mem.vmrss(), "]")
  type <- paste0("[", type, "]")
  message(paste0(type, dd, mm, " --- ", sub("\n$", "", paste(msg, collapse = " "))))
}

#' Log debug message with timestamp
#'
#' @param ... Message parts to concatenate.
#' @return NULL (called for side effect).
#' @keywords internal
dbg <- function(...) info(..., type = "DBUG")
