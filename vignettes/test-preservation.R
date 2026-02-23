#=========================================================================
#================== WGCNA PRESERVATION ANALYSIS ==========================
#=========================================================================

library(devtools)
library(playbase)
library(WGCNA)
library(igraph)

pgx <- pgx.load("~/Playground/omicsplayground/data/multi-gbm.region.pgx")
##pgx <- pgx.load("~/Playground/omicsplayground/data/multi-liver.pgx")

samples <- pgx$samples
head(samples)
group <- pgx$samples[,"group"]
group <- pgx$samples[,"sex"]
group <- pgx$samples[,"tumor_site"]

table(group)
exprList <- tapply( 1:ncol(pgx$X), group, function(i) pgx$X[,i])
lapply(exprList,dim)

head(rownames(exprList[[1]]))
for(i in 1:length(exprList)) {
  exprList[[i]] <- rename_by2(exprList[[i]], pgx$genes, "human_ortholog")
}

##--------------------------------------------------------
##--------------------------------------------------------
##--------------------------------------------------------
##source("~/Playground/playbase/dev/include.R", chdir=TRUE)
load_all()

names(exprList)
lapply(exprList,dim)

res <- runPreservationWGCNA(
  exprList,
  phenoData = samples,
  contrasts = contrasts,
  GMT = pgx$GMT,
  annot = pgx$genes,
  reference = 1,
  ngenes = 2000,
  compute.stats = TRUE,
  compute.enrichment = TRUE,
  add.merged = FALSE,
  ai_model = "groq:openai/gpt-oss-20b",
  summary = TRUE
)

names(res)
names(res$gsea)
names(res$summary)

cat(res$summary[[1]])
cat(res$summary[[2]])

##----------------------------------------------------------
## Figures
##----------------------------------------------------------
##source("~/Playground/playbase/dev/include.R", chdir=TRUE)

par(mar=c(4,10,6,2))
plotDendroAndColors(
  res, use.tree=0, show.traits=1, marAll = c(1, 15, 3, 0.2),
  show.kme=0, show.contrasts=1, rm.na=1)

## Zsummary and median.rank statistics plots
par(mfrow=c(2,2),mar=c(5,5,4,2), cex=1.4)
plotPreservationSummaries(res, setpar=FALSE)

## Cross tabulation
par(mar=c(10,15,5,2), cex=1.7)
plotConsensusOverlapHeatmap(
  res$net, res$layer[[2]]$net,
  ##setLabels = names(res$layers)[1:2],
  plotDendro = TRUE,
  lab.line = c(8,12), setpar=FALSE)

## Dendro + heatmap (MEs preservation)
par(mar=c(2,5,4,2),cex=1.2)
WGCNA::plotEigengeneNetworks(
  res$net$multiMEs,
  setLabels = names(res$net$multiMEs),
  marHeatmap = c(3,4,2,2)
)

##----------------------------------------------------------
## Module-trait heatmaps
##----------------------------------------------------------
plotPreservationModuleTraits(res, rm.na=TRUE, order.by="name")
plotPreservationModuleTraits(res, rm.na=TRUE, order.by="zsummary")
plotPreservationModuleTraits(res, rm.na=TRUE, order.by="clust")

##----------------------------------------------------------
## Module scores plots
##----------------------------------------------------------

colnames(res$layers[[1]]$datTrait)
trait = "BMI=obese"
trait = "BMI=overweight"
trait

wgcna.plotTopModules_multi(res, trait, collapse=FALSE)
wgcna.plotModuleScores(res, trait, nmax=16, collapse=FALSE, multi=TRUE)
plotModuleScores(res, trait, nmax=16, collapse=FALSE)

wgcna.plotModuleScores(res, trait, nmax=16, collapse=TRUE, multi=TRUE)
plotModuleScores(res, trait, nmax=16, collapse=TRUE)

trait <- head(colnames(res$modTrait[[1]]),4)
trait
par(mfrow=c(2,2), mar=c(10,4,4,2),cex=1.3)
plotTraitCorrelationBarPlots(res, trait[1:4], multi=TRUE, colored=TRUE,
  beside=TRUE, setpar=FALSE)
plotTraitCorrelationBarPlots(res, trait[1:4], multi=TRUE, colored=FALSE,
  beside=TRUE, setpar=FALSE)

par(mfrow=c(2,2), mar=c(10,4,4,2),cex=1.3)
wgcna.plotTraitCorrelationBarPlots(res, trait[1:2], multi=TRUE,
  colored=TRUE, beside=FALSE, setpar=FALSE)
wgcna.plotTraitCorrelationBarPlots(res, trait[1:2], multi=TRUE,
  colored=FALSE, beside=FALSE, setpar=FALSE)

##------------------------------------------------------------
##------------------- Module subnetworks ---------------------
##------------------------------------------------------------
##source("~/Playground/playbase/dev/include.R", chdir=TRUE)

net <- res$layers[[1]]
names(net$me.genes)
module = 'MEblack'
module = 'MEgreenyellow'
module = 'MEturquoise'
module

genes <- net$me.genes[[module]]
genes <- genes[order(-net$stats[['moduleMembership']][genes,module])]
##genes <- sample(genes, 35)
genes <- head(genes, 30)
genes

par(mfrow=c(2,2), mar=c(0,0,3,0))
k=1
for(k in 1:length(res$layers)) {
  wgcna.plotGeneNetwork(
    res$layers[[k]],
    sort(genes),
    #edge.gamma = 8,
    min.rho = 0.5,
    edge.alpha=0.3
  )
  title(names(res$layers)[k], line=1, cex.main=1.3)
  title(paste(module,"module"), line=0, cex.main=1, font.main=1)  
}


R <- cor(res$layers[[1]]$datExpr[,genes])
ii <- hclust(as.dist(1-R), method="average")$order
genes <- genes[ii]

##par(mfrow=c(2,2))
par(mar=c(3,4,4,6))
k=1
for(k in 1:length(res$layers)) {
  w <- res$layers[[k]]
  wgcna.plotModuleHeatmap(
    w, module=module, genes=genes,
    min.rho = 0.6,
    cluster = FALSE,
    main = ""
  )
  title( names(res$layers)[k], line=1, cex.main=1.3)  
}


##------------------------------------------------------------
##------------------------------------------------------------
##------------------------------------------------------------

