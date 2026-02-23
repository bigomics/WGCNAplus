#' Get merged multiomics geneset x feature matrix.
#' Combines geneset and metabolite sets.
#' @keywords internal
getPlaydataGMT <- function() {

  if (!requireNamespace("playdata", quietly = TRUE)) {
    stop("Package 'playdata' required. Install from GitHub: bigomics/playdata")
  }
  
  G1 <- Matrix::t(playdata::GSETxGENE)
  G2 <- Matrix::t(playdata::MSETxMETABOLITE)
  rownames(G2) <- sub(".*:","",rownames(G2))

  return(merge_sparse_matrix(G1, G2))

}
