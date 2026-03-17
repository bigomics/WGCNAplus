#=========================================================================
#================== WGCNA PRESERVATION ANALYSIS ==========================
#=========================================================================
library(devtools)
load_all()
setwd("./vignettes/")

X <- read.csv("./data/liver/expression.csv", row.names = 1)
samples <- read.csv("./data/liver/samples.csv", row.names = 1)
contrasts <- read.csv("./data/liver/contrasts.csv", row.names = 1)
annot <- read.csv("./data/liver/annot.csv", row.names = 1)
#GMT <- readRDS("./data/liver/gmt.RDS")

group <- samples$sex
# group <- samples$group
# group <- samples$tumor_site

exprList <- tapply(1:ncol(X), group, function(i) X[,i])
exprList[[1]] <- as.matrix(exprList[[1]])
exprList[[2]] <- as.matrix(exprList[[2]])
lapply(exprList,dim)

head(rownames(exprList[[1]]))
for(i in 1:length(exprList)) {
exprList[[i]] <- rename_by2(exprList[[i]], annot, "human_ortholog")
}
head(rownames(exprList[[1]]))

res <- WGCNAplus::runPreservationWGCNA(
  exprList,
  phenoData = samples,
  contrasts = contrasts,
  #GMT = GMT,
  annot = annot,
  reference = 1, ## ref female
  ngenes = 2000,
  compute.stats = TRUE,
  compute.enrichment = FALSE,#TRUE,
  add.merged = FALSE,
  ai_model = "groq:openai/gpt-oss-20b",
  summary = TRUE
)

names(res)
names(res$gsea)
names(res$summary)
cat(res$summary[[1]])
cat(res$summary[[2]])

## Figures
par(mar = c(4, 10, 6, 2))
WGCNAplus::plotDendroAndColors(
  res,
  use.tree = 0,
  show.traits = 0,
  marAll = c(1, 15, 3, 0.2),
  show.kme = 0,
  show.contrasts = 1,
  rm.na = 1
)

## Zsummary and median.rank statistics plots
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2), cex = 1, las = 1)
WGCNAplus::plotPreservationSummaries(res, setpar = FALSE)

## Cross tabulation
pdf("./paper_figures/preservation.pdf", height = 8, width = 8)
par(mar = c(10, 15, 5, 2), cex = 1.5)
WGCNAplus::plotConsensusOverlapHeatmap(
  res$net,
  res$layer[[2]]$net,
  plotDendro = TRUE,
  lab.line = c(8, 12),
  setpar = FALSE
)
dev.off()

pdf("./paper_figures/preservation-nodendro.pdf", height = 10, width = 12)
par(mar = c(12, 10, 5, 2), cex = 1.2)
WGCNAplus::plotConsensusOverlapHeatmap(
  res$net,
  res$layer[[2]]$net,
  lab.line = c(8, 12),
  setpar = FALSE
)
dev.off()

## Dendro + heatmap (MEs preservation)
par(mar = c(2, 5, 4, 2), cex = 1.2)
WGCNA::plotEigengeneNetworks(
  res$net$multiMEs,
  setLabels = names(res$net$multiMEs),
  plotDendrograms = FALSE,
  marHeatmap = c(5, 5, 2, 2),
  printPreservation = TRUE
)

res$modulePreservation

## Module-trait heatmaps
par(mfrow = c(1, 3), mar = c(5, 9, 4, 6), cex = 1.15)
WGCNAplus::plotPreservationModuleTraits(res, setpar = FALSE, rm.na = TRUE, order.by = "zsummary")

## Module scores plots
colnames(res$layers[[1]]$datTrait)
trait = "BMI=obese"
trait = "BMI=overweight"

WGCNAplus::wgcna.plotTopModules_multi(res, trait, collapse = FALSE)
WGCNAplus::wgcna.plotModuleScores(res, trait, nmax = 16, collapse = FALSE, multi = TRUE)
WGCNAplus::plotModuleScores(res, trait, nmax = 16, collapse = FALSE)

WGCNAplus::wgcna.plotModuleScores(res, trait, nmax = 16, collapse = TRUE, multi = TRUE)
WGCNAplus::plotModuleScores(res, trait, nmax = 16, collapse = TRUE)

trait <- head(colnames(res$modTrait[[1]]), 4)
par(mfrow = c(2, 2), mar = c(10, 4, 4, 2), cex = 1.3)
WGCNAplus::plotTraitCorrelationBarPlots(
  res,
  trait[1:4],
  multi = TRUE,
  colored = TRUE,
  beside = TRUE,
  setpar = FALSE
)

WGCNAplus::plotTraitCorrelationBarPlots(
  res,
  trait[1:4],
  multi = TRUE,
  colored = FALSE,
  beside = TRUE,
  setpar = FALSE
)

par(mfrow = c(2, 2), mar = c(10, 4, 4, 2), cex = 1.3)
WGCNAplus::wgcna.plotTraitCorrelationBarPlots(
  res,
  trait[1:2],
  multi = TRUE,
  colored = TRUE,
  beside = FALSE,
  setpar = FALSE
)

WGCNAplus::wgcna.plotTraitCorrelationBarPlots(
  res,
  trait[1:2],
  multi = TRUE,
  colored = FALSE,
  beside = FALSE,
  setpar = FALSE
)



## Gset and pathway enrichment analysis of top modules
topmodules <- c("green", "blue", "turquoise", "brown")
cons_colors <- res$colors[, "Consensus"]
i=1; module_genes=list()
for(i in 1:length(topmodules)) {
  jj <- which(cons_colors == topmodules[i])
  module_genes[[topmodules[i]]] <- names(cons_colors)[jj]
}
enr <- gprofiler2::gost(
  query      = module_genes,
  organism   = "mmusculus",
  user_threshold = 0.01,
  sources    = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
  significant = TRUE,
  multi_query = FALSE
)
#gprofiler2::gostplot(enr, capped = TRUE, interactive = TRUE)
enr <- enr$result[, c("query", "source", "term_name", "p_value", "intersection_size")]
kk <- which(enr$p_value <= 0.01 & enr$intersection_size >= 20)
enr <- enr[kk, ]
enr <- enr[order(enr$p_value, decreasing = FALSE), ]
enr[1:10,-5]
dim(enr)
immune <- enr[grep("immune",enr$term_name), ]; dim(immune)
metabo <- enr[grep("taboli",enr$term_name), ]; nrow(mm)
homeo <- enr[grep("homeo",enr$term_name), ]; dim(homeo)

## ##------------------- Module subnetworks ---------------------NOT RUN
## net <- res$layers[[1]]
## names(net$me.genes)
## module = 'MEblack'
## module = 'MEgreenyellow'
## module = 'MEturquoise'

## genes <- net$me.genes[[module]]
## genes <- genes[order(-net$stats[['moduleMembership']][genes,module])]
## genes <- head(genes, 30)

## par(mfrow = c(2, 2), mar = c(0, 0, 3, 0))
## for(k in 1:length(res$layers)) {
##   wgcna.plotGeneNetwork(
##     res$layers[[k]],
##     sort(genes),
##     min.rho = 0.5,
##     edge.alpha = 0.3
##   )
##   title(names(res$layers)[k], line = 1, cex.main = 1.3)
##   title(paste(module,"module"), line = 0, cex.main = 1, font.main = 1)  
## }


## R <- cor(res$layers[[1]]$datExpr[,genes])
## ii <- hclust(as.dist(1-R), method="average")$order
## genes <- genes[ii]

## par(mar = c(3, 4, 4, 6))
## for(k in 1:length(res$layers)) {
##   w <- res$layers[[k]]
##   wgcna.plotModuleHeatmap(
##     w,
##     module = module,
##     genes = genes,
##     min.rho = 0.6,
##     cluster = FALSE,
##     main = ""
##   )
##   title(names(res$layers)[k], line = 1, cex.main = 1.3)  
## }
