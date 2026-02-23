library(devtools)
library(igraph)

load_all("../../lasagna")
load_all()

## Imports
gset.rankcor <- playbase::gset.rankcor
mat2gmt <- playbase::mat2gmt
##mofa.merge_data2 <- playbase::mofa.merge_data2
ai.ask <- playbase::ai.ask
ai.create_image_gemini <- playbase::ai.create_image_gemini
lasagna.multisolve <- playbase::lasagna.multisolve

##--------------------------------------------------------------------
## Compute multi-omics
##--------------------------------------------------------------------

pgx <- playbase::pgx.load("~/ClientData/bruker/klidel-mox/klidel-multiomics.pgx")
pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/mox-brca.pgx")
pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/mox-geiger.pgx")

write.csv(pgx$X, file="brca/expression.csv")
write.csv(pgx$samples, file="brca/samples.csv")
write.csv(pgx$contrasts, file="brca/contrasts.csv")
save(pgx$GMT, file="brca/gmt.rda")
pgx <- playbase::pgx.initialize(pgx)

annot <- pgx$genes
pgx$X <- rename_by2(pgx$X, annot, "symbol", keep.prefix=TRUE)
dataX <- mofa.split_data(pgx$X)
names(dataX)
dataX <- dataX[c("mir","gx","px")]
lapply(dataX,dim)

annot <- pgx$genes
GMT0 <- rename_by2(getPlaydataGMT(), annot, "symbol")
GMT1 <- rename_by2(pgx$GMT, annot, "symbol")
GMT <-  merge_sparse_matrix(GMT0, GMT1)
dim(GMT)
head(rownames(GMT))
tail(rownames(GMT))

## full WGCNA on features
wgcna <- computeWGCNA_multiomics(
  dataX, pgx$samples, GMT=GMT, 
  power="iqr", minmodsize=3, minKME=0.3,
  mergeCutHeight=0.3  ## seems better than 0.15
)
wgcna$me.genes

pdf("geiger/wgcna-dendro.pdf", w=14, h=7)
par(mfrow=c(1,1),cex=1.4)
plotMultiDendroAndColors(
  wgcna, marAll=c(2,7,3,1),  
  show.traits=1, show.contrasts=1,
  show.kme=0, use.tree=0,
  colorHeight = 0.5
)
dev.off()

pdf("geiger/wgcna-moduletrait.pdf", w=14, h=7)
par(mfrow=c(1,1),mar=c(8,11,4,2), cex=1.2)
plotModuleTraitHeatmap(wgcna, cluster=TRUE,
  main="Module-Trait Heatmap", setpar=FALSE,
  transpose=TRUE, text=FALSE, multi=TRUE)
dev.off()

names(wgcna)

##--------------------------------------------------------------------
## LASAGNA (full)
##--------------------------------------------------------------------

data <- list(
  X = dataX,
  samples = pgx$samples
)

obj <- NULL
obj <- lasagna::create_model(data, pheno="pheno", ntop=1000, nc=10,
  add.sink=TRUE, intra=FALSE, fully_connect=FALSE, add.revpheno=TRUE,
  condition.edges=1)
names(obj)

## color by WGCNA clustering
wgcna$me.genes
table(wgcna$me.colors)
head(wgcna$me.colors)
V(obj$graph)$color <- wgcna$me.colors[V(obj$graph)$name]
##V(obj$graph)$color <- substring(wgcna$me.colors[V(obj$graph)$name],3,99)
V(obj$graph)$color[is.na(V(obj$graph)$color)] <- "red"
table(V(obj$graph)$color)

## solve the graph for a certain phenotype
colnames(obj$Y)
pheno = colnames(obj$Y)[1]
pheno = colnames(obj$Y)[3]
pheno = colnames(obj$Y)[13]
graph <- lasagna::solve(obj, pheno, min_rho=0.01, max_edges=1000,
  value="rho", sp.weight=1, prune=FALSE) 
graph

##tt <- V(graph)$name


## prune graph for cleaner plotting
pdf("geiger/lasagna-full.pdf", w=14, h=8)
par(mfrow=c(1,1), mar=c(1,1,1,1)*0, cex=0.9)
mp <- plot_multipartite(
  graph,
  min.rho = 0.9,
  ntop = 40,
  xdist = 2,
##  label = tt,
  color.var = "color",
##  labpos = c(2,2,4,4),
  cex.label = 0.85,
  vx.cex = 1.8,
  edge.cex = 1,
  edge.alpha = 0.25,
  edge.sign = "both",
  edge.type = "inter",
  edge.gamma = 4,
  yheight = 0.9,
  normalize.edges = 1,
  strip.prefix = TRUE,
  prune = 0
) 
dev.off()

##--------------------------------------------------------------------
## LASAGNA (modules)
##--------------------------------------------------------------------
  
library(igraph)
lasagna.create_model <- playbase::lasagna.create_model
lasagna.solve <- playbase::lasagna.solve
lasagna.multisolve <- playbase::lasagna.multisolve  
plotMultiPartiteGraph2 <- playbase::plotMultiPartiteGraph2

names(wgcna)
(wgcna$me.genes)

ww <- wgcna$layers
names(ww)
##ww <- ww[c("mir","gx","px")]
xx <- lapply(ww, function(m) t(m$net$MEs))
lapply(xx,dim)

data <- list(
  X = xx,
  samples = ww[[1]]$datTraits
)

lasagna <- lasagna::create_model(
  data,
  pheno = "expanded",
  ntop = 2000,
  nc = 40,
  add.sink = FALSE,
  intra = TRUE,
  fully_connect = FALSE,
  add.revpheno = TRUE,
  condition.edges = 0
)

names(lasagna)
class(lasagna$graph)
wgcna$graph <- lasagna$graph  
head(lasagna$Y)

solved <- lasagna::solve(
  lasagna,
  pheno = pheno,
  min_rho = 0.01,
  max_edges = 100,
  value.type = "rho",
  fc.weights = TRUE,
  sp.weight = FALSE,
  prune = FALSE
)

## transfer colors
V(solved)$color <- 'red'
ii <- grep("PHENO|SINK|SOURCE",V(solved)$name,invert=TRUE)
V(solved)$color[ii] <- sub(".*:","",V(solved)$name)[ii]

pdf("geiger/lasagna-wgcnaNC.pdf",w=12,h=6)
par(mfrow=c(1,1), mar = c(0,0,0,0), cex=1)
plot_multipartite(
  solved,
  min.rho = 0.2,
  ntop = -1,
  cex.label = 1.0,
  vx.cex = 3.4,
  color.var = "color",
  edge.cex = 1.7,
  edge.alpha = 0.3,
  edge.gamma = 2,
  edge.sign = "both",
  edge.type = "inter",
  yheight = 0.8,
  xdist = 3,
  normalize.edges = 1,
  strip.prefix2 = TRUE,
  prune = FALSE
)
title("WGCNA-Lasagna Colored by Module", line=-2, cex.main=1.6)
dev.off()


##--------------------------------------------------------------------
## Biomarkers
##--------------------------------------------------------------------
load_all()

QQ <- calculateCompoundSignificance(wgcna) 
names(QQ)
k=2
Q <- QQ[[k]]

pdf("geiger/biomarker-significance.pdf",w=12,h=8)
par(mfrow=c(2,1), mar=c(8,6,3,0), cex=0.9)
ntop=60
Q <- Q[order(-Q$score1),]
Q1 <- head(Q,ntop)
dt <- names(QQ)[k]
barplot( Q1$score1, col=Q1$color,
  ylab = "significance (MM*TS)", main=paste("biomarker significance",dt),
  cex.names=1, names.arg = rownames(Q1), las=3 )
Q <- Q[order(-Q$score2),]
Q1 <- head(Q,ntop)
barplot( Q1$score2, col=Q1$color,
  ylab = "significance (MM*FC)", main=paste("biomarker significance",dt),
  cex.names=1, names.arg = rownames(Q1), las=3 )
dev.off()

##--------------------------------------------------------------------
## TEST
##--------------------------------------------------------------------
  
group=cc
matrix_group_stats <- function(X, group, FUN="mean", probs=0.5) {
  nm <- length(unique(group))
  M <- matrix(NA, nm, nm)
  colnames(M) = rownames(M) = unique(group)  
  for(i in unique(group)) {
    for(j in unique(group)) {
      ii <- which(group == i)
      jj <- which(group == j)
      if(FUN == "mean") {
        M[i,j] = M[j,i] = mean( X[ii,jj], na.rm=TRUE )
      }
      if(FUN == "quantile") {
        probs2 <- c(1-probs,probs)
        pp <- quantile(X[ii,jj], probs=probs2)
        pp <- ifelse(-pp[1] > pp[2], pp[1], pp[2])
        M[i,j] = M[j,i] = pp
      }
    }
  }
  M
}

X <- t(pgx$X)
dim(X)
ii <- paste0("gx:",wgcna$me.genes[['GXblue']])
jj <- paste0("px:",wgcna$me.genes[['PXblack']])

Y <- wgcna$layers[[1]]$datTraits
dim(Y)
head(Y)

cc <- wgcna$me.colors
cc <- cc[colnames(X)]
length(cc)

rho <- cor(X, Y[,2])[,1]
head(rho)
length(rho)
dim(X)

rho <- stats::cor(X, Y, use = "pairwise.complete.obs")
maxrho <- apply(abs(rho), 1, max, na.rm = TRUE)
ii <- grep("SINK|SOURCE", names(maxrho))
if (length(ii)) maxrho[ii] <- 1
outer.rho <- outer(maxrho, maxrho)
dim(outer.rho)

outer.rho <- abs(rho %*% t(rho))**0.5
outer.rho <- abs(rho %*% t(rho))**1
outer.rho <- pmax(rho %*% t(rho),0)**0.5
outer.rho <- abs(outer(rho,rho))**0.5
dim(outer.rho)

rhoX <- cor(X) 
rhoX <- cor(X) * outer.rho
round(rhoX[ii,jj],4)
dim(rhoX)

R1 <- matrix_group_stats(rhoX, group=cc, FUN="mean", probs=0.5)
R1 <- matrix_group_stats(rhoX, group=cc, FUN="quantile", probs=0.5)
R1 <- matrix_group_stats(rhoX, group=cc, FUN="quantile", probs=0.99)
R1 <- matrix_group_stats(rhoX, group=cc, FUN="quantile", probs=1)
R1 <- matrix_group_stats(rhoX, group=cc, FUN="quantile", probs=0.8) 
sort(R1["GXblue",])

gx.heatmap(R1, scale="none", mar=c(1,1)*15,
  keysize=0.8, cexRow=2, cexCol=2)

d1 <- diag(R1)
r1.names <- colnames(R1)
R1 <- diag(1/sqrt(d1)) %*% R1 %*% diag(1/sqrt(d1))
rownames(R1) = colnames(R1) = r1.names
gx.heatmap(R1, scale="none", mar=c(1,1)*15,
  keysize=0.8, cexRow=2, cexCol=2)

me <- lapply(wgcna$layers,function(w) w$net$MEs)
names(me) <- NULL
M <- do.call(cbind, me)
dim(M)

MY <- cbind(M,Y)
dim(MY)

Rm <- cor(MY)
sort(Rm["PXblack",])
sort(Rm["GXblue",])

R1 <- cor(X[,ii],MY)
R1
colMeans(R1)
apply(R1,2,quantile,probs=c(0,0.01,0.5,0.99,1))

M <- M[,grep("gx.",colnames(M))]
colnames(M) <- paste("***",colnames(M),"***")

xx <- cbind( X[,ii], M)
Rx <- cor(t(X))
dim(Rx)

nm <- length(unique(cc))
M <- matrix(NA, nm, nm)
colnames(M) = rownames(M) = unique(cc)  
for(i in unique(cc)) {
  for(j in unique(cc)) {
    ii <- which(cc == i)
    jj <- which(cc == j)      
    M[i,j] = M[j,i] = mean( Rx[ii,jj] )
  }
}

R <- cor(M)
R <- cor(t(xx), t(X[jj,]))

pdf("moxbrca-module-corr.pdf",w=12,h=12)
gx.heatmap(R, scale='none', mar=c(25,20),
  keysize=0.6, cexRow=1.6, cexCol=1.6)
dev.off()

dim(R)  
mean(R)
quantile(R, probs=0.9)
quantile(R, probs=0.99)
max(R)  
