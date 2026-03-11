##-----------------------
## Consensus WGCNA
##-----------------------
library(devtools)
load_all()
setwd("vignettes/")

X <- read.csv("./data/liver/expression.csv", row.names = 1)
samples <- read.csv("./data/liver/samples.csv", row.names = 1)
contrasts <- read.csv("./data/liver/contrasts.csv", row.names = 1)
annot <- read.csv("./data/liver/annot.csv", row.names = 1)
GMT <- readRDS("./data/liver/gmt.RDS")

dim(X)
group <- samples$sex
xx <- tapply(1:ncol(pgx$X), group, function(ii) pgx$X[,ii])

par(mfrow = c(2, 3), mar = c(5, 5, 3, 1), cex = 1.4)
plotPowerAnalysis(t(xx[[1]]), setPar = FALSE)
plotPowerAnalysis(t(xx[[2]]), setPar = FALSE)

## Run consensus WGCNA on an expression list
Y <- samples[, 3:6]
cons <- WGCNAplus::runConsensusWGCNA(
  xx,
  phenoData = Y,
  power = 6,
  annot = pgx$genes,
  compute.enrichment = 1,
  ngenes = 2000,
  minModuleSize = 40,
  maxBlockSize = 9999,
  minKME = 0.3,
  mergeCutHeight = 0.15,
  deepSplit = 2,
  calcMethod = "fast",
  drop.ref = FALSE,
  addCombined = FALSE,
  gsea.mingenes = 10,
  summary = TRUE,
  verbose = 1
)

names(cons)
cons$net$power

x11()
WGCNAplus::plotDendroAndColors(
  cons,
  marAll = c(2, 10, 3, 1),
  show.traits = FALSE, ## TRUE
  show.kme = 0,
  use.tree = 0,
  colorHeight = 0.2,
  colorHeightMax = 0.7,
  setLayout = 1
)

head(cons$datTraits)
WGCNAplus::plotModuleScores(cons, trait = "ab_fat", nmax = 9)
WGCNAplus::plotModuleScores(cons$layers[[1]], trait = "ab_fat", nmax = 9)

WGCNAplus::plotConsensusTraitCorrelation(cons, traits = NULL) 

## gene statistics
top <- WGCNAplus::getTopGenesAndSets(cons, module = NULL, ntop = 10) 
names(top)
lapply(top$genes,head)
lapply(top$sets,head)

stats <- WGCNAplus::computeConsensusGeneStats(cons)
names(stats)

head(stats[[1]][[1]])
module = "MEblue"
trait = "length_cm"
trait = "other_fat"
trait = "total_fat"

stats2 <- WGCNAplus::getConsensusGeneStats(cons, stats = stats, trait = trait, module = module)
names(stats2)
head(stats2[['full']])
head(stats2[['consensus']])

cons.modules <- stats2[["consensus"]]$module
cons.modules <- unique(cons.modules[which(stats2[["consensus"]]$consensus=="C")])
cons.modules

## enrichment
# GMT = Matrix::t(playdata::GSETxGENE)
lapply(cons$datExpr,dim)
cons$gsea <- WGCNAplus::computeConsensusModuleEnrichment(
  cons,
  GMT = GMT,
  annot = annot,
  methods = c("fisher","gsetcor","xcor"),
  min.genes = 5,
  ntop = 1000
)

names(cons$gsea)

top.gs <- head(cons$gsea[['MEred']]$geneset,100) 
head(top.gs,20)

top.gs <- head(cons$gsea[['MEblue']]$geneset,100) 
head(top.gs,20)


ai_model = "qwen3:1.7b"
ai_model = "groq:openai/gpt-oss-20b"
playbase::ai.genesets_summary(top.gs, model=ai_model)
playbase::ai.genesets_keywords(top.gs, model=ai_model) 

## Sample clustering dendrograms (useless?)
lapply(xx, dim)
dim(cons$modTraits)
head(cons$datTraits)[,1:4]

nsets <- length(xx)
layout.matrix <- matrix( 1:(2*nsets), nrow = 2, ncol = nsets)
layout(layout.matrix, heights=c(1,1), widths=rep(1,nsets))
for(i in 1:nsets) {
  dt <- toupper(names(cons$datExpr)[i])
  plotConsensusSampleDendroAndColors(
    cons,
    i,
    main = paste("sample tree and traits heatmap for", dt),
    what = c("me", "traits", "both")[2],
    marAll = c(0.2, 8, 1, 0.2),
    clust.expr = TRUE,
    setLayout = FALSE,
    colorHeightMax = 0.6
  ) 
}

## Module-trait heatmaps (important plots)
lapply(cons$zlist, dim)
ii <- hclust(dist(cons$zlist[[1]]))$order
jj <- hclust(dist(t(cons$zlist[[1]])))$order

## fix ordering of heatmaps
Z <- Reduce('+', lapply(cons$zlist,dist))
tZ <- Reduce('+', lapply(lapply(cons$zlist,t),dist))
ii <- hclust(Z)$order
jj <- hclust(tZ)$order

par(mfrow = c(2, 3), mar = c(8, 12, 3, 1), cex = 0.8)
for(i in 1:length(cons$zlist)) {
  k <- names(cons$zlist)[i]
  plotLabeledCorrelationHeatmap(
    cons$zlist[[i]][ii,jj],
    cons$ydim[i],
    cex.lab = 1.2,
    pstar = 1,
    text = FALSE,
    cluster = FALSE,
    setpar = FALSE,
    main = paste("module-trait for",toupper(k))
  )
}

## Consensus Module-Trait
matlist = cons$zlist
ydim <- sapply(cons$datExpr, nrow)
consZ <- computeConsensusMatrix(cons$zlist, ydim=ydim, psig=0.05) 
nsamples <- ncol(xx[[1]])
plotLabeledCorrelationHeatmap(
  consZ[ii,jj], nsamples, setpar=FALSE,
  text=FALSE, pstar=1, cluster=FALSE, cex.lab=1.2,
  main = "consensus Module-Trait"
)

## Distinct Module-Trait
diffZ <- computeDistinctMatrix(matlist, ydim = ydim, psig = 0.05, min.diff = 0.1) 
names(diffZ)
for(set in names(diffZ)) {
  plotLabeledCorrelationHeatmap(
    diffZ[[set]][ii,jj],
    nsamples,
    setpar = FALSE,
    text = FALSE,
    pstar = 1,
    cluster = FALSE,
    cex.lab = 1.2,
    main = paste("unique Module-Traits for", toupper(set))
  )
}

dim(consZ)
z0 <- consZ
z0[is.na(z0)] <- 0
z1 <- diffZ
for(i in 1:length(z1)) z1[[i]][is.na(z1[[i]])] <- 0
effZ <- abs(z0) - Reduce('+', lapply(z1,abs))
wgcna.plotLabeledCorrelationHeatmap(
  effZ[ii,jj],
  nsamples,
  setpar = FALSE,
  colorpal = purpleGreyYellow,
  text = FALSE,
  pstar = 1,
  cluster = FALSE,
  cex.lab = 1.2,
  main = "effective consensus"
)

#----------------------------------------------------------------
#  Conservation/Preservation Summary plot (nice plot!)
#----------------------------------------------------------------

MEs <- cons$net$multiMEs
lapply(MEs, function(m) dim(m$data))
modules <- colnames(MEs[[1]]$data)

par(mfrow = c(2, 1), mar = c(12, 4, 2, 2))
for(k in names(MEs)) {
  ME <- MEs[[k]]$data
  plotEigenGeneClusterDendrogram(
    wgcna = NULL,
    ME = ME,
    main = k,
    setMargins = FALSE,
    method = "hclust"
  )
}

par(cex = 0.9)
WGCNA::plotEigengeneNetworks(
  MEs,
  setLabels = names(xx),
  plotDendrograms = FALSE,
  marDendro = c(0,2,2,1)*2,
  marHeatmap = c(3,3,2,1)*2,
  zlimPreservation = c(0.5, 1),
  xLabelsAngle = 90
)


# Create correspondence table (not so usefull...)
net0 <- cons$net
names(cons$netList)
net1 <- cons$layers[[1]]$net
net2 <- cons$layers[[2]]$net 
names(cons$layers)
par(cex = 1.3)
plotConsensusOverlapHeatmap(net1, net2, setLabels=names(cons$netList))

