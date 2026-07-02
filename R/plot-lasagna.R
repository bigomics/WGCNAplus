
#' Plot Lasagna graph with nodes colered by WGCNA ME colors
#' 
#' @export
wgcna.plot_lasagna <- function(wgcna, pheno, layout='tsne',
                               edge_color = c("blue", "magenta")) {

  xx <- lapply(wgcna$layers, function(s) t(s$datExpr))
  Y <- wgcna$layers[[1]]$datTraits
  
  obj <- lasagna::create_model(
    data = NULL,
    X = xx,
    meta = Y,
    meta.type = "traits",
    ntop = 1000,
    nc = 10,
    add.sink = TRUE,
    intra = FALSE,
    fully_connect = FALSE,
    add.revpheno = TRUE,
    condition.edges = 1
  )

  ## color by WGCNA clustering
  ## wgcna$me.genes
  ## table(wgcna$me.colors)
  ## head(wgcna$me.colors)
  V(obj$graph)$color <- wgcna$me.colors[V(obj$graph)$name]
  V(obj$graph)$color[is.na(V(obj$graph)$color)] <- "red"

  ## solve the graph for requested phenotype
  if(!pheno %in% colnames(obj$Y)) {
    stop("pheno not in traits matrix")
  }
  graph <- lasagna::solve(
    obj,
    pheno,
    min_rho = 0.01,
    max_edges = 1000,
    value = "rho",
    sp.weight = 1,
    prune = FALSE
  ) 

  ## 3D layout
  xpos <- layout_multipartite_3d(graph, obj$X, clust=layout)

  ## plot lasagna
  p <- plot_3d(graph, layout=xpos, draw_edges=TRUE,
    color.by="color", min_rho=0.0, sign_rho="both",
    edge_color = edge_color,
    cex=0.9, cex.gamma=0.5, num_edges=200, znames=NULL) 

  return(p)
}

