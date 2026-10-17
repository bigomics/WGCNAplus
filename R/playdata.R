#' Get merged multiomics geneset x feature matrix from lib-x.
#' Combines the legacy gene-set and metabolite-set matrices of the release.
#' @keywords internal
#' @export
getLibxGMT <- function() {

  if (!requireNamespace("omicspgx", quietly = TRUE) || is.null(omicspgx::libx_dir())) {
    stop("Package 'omicspgx' with a configured lib-x directory is required")
  }

  G1 <- Matrix::t(omicspgx::libx_read("featuresets/legacy/gsetxgene"))
  G2 <- Matrix::t(omicspgx::libx_read("featuresets/legacy/msetxmetabolite"))
  rownames(G2) <- sub(".*:","",rownames(G2))

  return(merge_sparse_matrix(G1, G2))

}

#' @rdname getLibxGMT
#' @details `getPlaydataGMT()` is the former name and stays as an alias.
#' @keywords internal
#' @export
getPlaydataGMT <- function() getLibxGMT()
