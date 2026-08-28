##-----------------------
## Consensus WGCNA
##-----------------------
library(devtools)
load_all()
setwd("vignettes/")

X <- read.csv("./data/liver/expression.csv", row.names = 1); dim(X)
samples <- read.csv("./data/liver/samples.csv", row.names = 1)
contrasts <- read.csv("./data/liver/contrasts.csv", row.names = 1)
annot <- read.csv("./data/liver/annot.csv", row.names = 1)
GMT <- readRDS("./data/liver/gmt.RDS")

X <- as.matrix(X)
all.equal(rownames(samples), colnames(X))
group <- samples$sex
exprList <- list()
exprList[["Female"]] <- X[, which(group == "Female")]; 
exprList[["Male"]] <- X[, which(group == "Male")]; 

## Run consensus WGCNA on an expression list
Y <- samples[, 3:6]
cons <- WGCNAplus::runConsensusWGCNA(
  exprList,
  phenoData = Y,
  power = 6,
  annot = annot,
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
cons$consModTraits
names(cons$layers)
mm <- cons$layers[["Male"]]
length(unique(mm$me.colors))
ff <- cons$layers[["Female"]]
length(unique(ff$me.colors))
cc <- cons$colors
head(cc)
unique(cc[,1])
table(cc[,1])
grey <- which(cc[,1] == "grey")
length(unique(rownames(cc)[-grey]))

ff <- "~/Desktop/WGCNAplus/WGCNAplus/vignettes/paper_figures_tables/Fig5A1.pdf"
pdf(file = ff, width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(6, 8, 1, 1))
net1 = cons$net
net2 = cons$layers[["Female"]]$net
lb = c("Consensus", "Female")
WGCNAplus::plotConsensusOverlapHeatmap(net1 = net1, net2 = net2, setLabels = lb)
dev.off()

ff <- "~/Desktop/WGCNAplus/WGCNAplus/vignettes/paper_figures_tables/Fig5A2.pdf"
pdf(file = ff, width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(6, 8, 1, 1))
net1 = cons$net
net2 = cons$layers[["Male"]]$net
lb = c("Consensus", "Male")
WGCNAplus::plotConsensusOverlapHeatmap(net1 = net1, net2 = net2, setLabels = lb)
dev.off()

lb1 <- cons$layers[["Female"]]$net$colors
lb2 <- cons$layers[["Male"]]$net$colors
ovl <- WGCNA::overlapTable(labels1 = lb1, labels2 = lb2)
ct <- ovl$countTable; dim(ct)
kk1 <- which(rownames(ct) != "grey")
kk2 <- which(colnames(ct) != "grey")
ct <- ct[kk1, kk2]
pv <- ovl$pTable
kk1 <- which(rownames(pv) != "grey")
kk2 <- which(colnames(pv) != "grey")
pv <- pv[kk1,kk2]
qv <- matrix(p.adjust(pv, method = "fdr"), nrow(pv), ncol(pv))
rownames(qv) <- rownames(pv)
colnames(qv) <- colnames(pv)
kk <- which(qv <= 0.05, arr.ind = TRUE)
cons.mods <- data.frame(
  m1 = rownames(pv)[kk[, 1]],
  m2 = colnames(pv)[kk[, 2]],
  nshared = ct[kk],
  p.value = pv[kk],
  q.value = qv[kk]
)
cons.mods <- cons.mods[order(cons.mods$q.value), ]
rownames(cons.mods) <- 1:nrow(cons.mods)
cons.mods


ff <- "~/Desktop/WGCNAplus/WGCNAplus/vignettes/paper_figures_tables/Fig60.pdf"
pdf(file = ff, width = 10, height = 6)
par(mfrow = c(1, 1), mar = c(6, 10, 1, 4))
nsamples <- min(ncol(exprList[["Female"]]), ncol(exprList[["Male"]]))
WGCNAplus::plotLabeledCorrelationHeatmap(cons$consModTraits, nsamples,
  cluster = TRUE, text = FALSE, pstar = TRUE, cex.lab=1.2, setpar = FALSE, main = "")
dev.off()


## ## gene statistics
## top <- WGCNAplus::getTopGenesAndSets(cons, module = NULL, ntop = 10) 
## names(top)
## lapply(top$genes,head)
## lapply(top$sets,head)

## stats <- WGCNAplus::computeConsensusGeneStats(cons)
## names(stats)

## head(stats[[1]][[1]])
## module = "MEblue"
## trait = "length_cm"
## trait = "other_fat"
## trait = "total_fat"

## stats2 <- WGCNAplus::getConsensusGeneStats(cons, stats = stats, trait = trait, module = module)
## names(stats2)
## head(stats2[['full']])
## head(stats2[['consensus']])

## cons.modules <- stats2[["consensus"]]$module
## cons.modules <- unique(cons.modules[which(stats2[["consensus"]]$consensus=="C")])
## cons.modules

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

## names(cons$zlist)
## par(mfrow = c(2, 3), mar = c(8, 12, 3, 1), cex = 0.8)
## for(i in 1:length(cons$zlist)) {
##   k <- names(cons$zlist)[i]
##   plotLabeledCorrelationHeatmap(
##     cons$zlist[[i]][ii,jj],
##     cons$ydim[i],
##     cex.lab = 1.2,
##     pstar = 1,
##     text = FALSE,
##     cluster = FALSE,
##     setpar = FALSE,
##     main = paste("module-trait for",toupper(k))
##   )
## }

## Consensus Module-Trait
matlist = cons$zlist
ydim <- sapply(cons$datExpr, nrow)
consZ <- computeConsensusMatrix(cons$zlist, ydim=ydim, psig=0.05) 
nsamples <- ncol(exprList[[2]])
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

## par(cex = 0.9)
## WGCNA::plotEigengeneNetworks(
##   MEs,
##   setLabels = names(exprList),
##   plotDendrograms = FALSE,
##   marDendro = c(0,2,2,1)*2,
##   marHeatmap = c(3,3,2,1)*2,
##   zlimPreservation = c(0.5, 1),
##   xLabelsAngle = 90
## )
