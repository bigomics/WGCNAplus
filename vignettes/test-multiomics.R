library(devtools)
load_all()

## Imports
## gset.rankcor <- playbase::gset.rankcor
## mat2gmt <- playbase::mat2gmt
## mofa.merge_data2 <- playbase::mofa.merge_data2
## ai.ask <- playbase::ai.ask
## ai.create_image_gemini <- playbase::ai.create_image_gemini
lasagna.multisolve <- playbase::lasagna.multisolve

##--------------------------------------------------------------------
## Compute multi-omics
##--------------------------------------------------------------------

pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/mox-brca.pgx")
pgx <- playbase::pgx.initialize(pgx)

dataX <- mofa.split_data(pgx$X)
names(dataX)
dataX <- dataX[c("mir","gx","px")]

annot <- pgx$genes
GMT0 <- rename_by2(getPlaydataGMT(), annot, "symbol")
GMT1 <- rename_by2(pgx$GMT, annot, "symbol")
GMT <-  merge_sparse_matrix(GMT0, GMT1)
dim(GMT)
head(rownames(GMT))
tail(rownames(GMT))

wgcna <- computeWGCNA_multiomics(pgx$X, pgx$samples, annot=pgx$genes, report=FALSE)
wgcna <- computeWGCNA_multiomics(dataX, pgx$samples, annot=pgx$genes,
  GMT=GMT, minmodsize = 5, mergeCutHeight = 0.6)
wgcna <- computeWGCNA_multiomics(dataX, pgx$samples, GMT=NULL)

wgcna <- computeWGCNA_multiomics(dataX, pgx$samples, GMT=GMT, 
  power=12, minmodsize = 3, minKME=0.1, mergeCutHeight = 0.5)

wgcna$me.genes

names(wgcna)
names(wgcna$layers)
names(wgcna$layers[[1]])

head(wgcna$layers[[1]]$gsea[[1]])
names(wgcna$layers[[1]]$gsea)
names(wgcna$layers[[2]]$gsea)
names(wgcna$layers[[3]]$gsea)

##--------------------------------------------------------------------
## Figures
##--------------------------------------------------------------------
load_all()

par(mfrow=c(1,1),cex=1.4)
plotMultiDendroAndColors(
  wgcna, marAll=c(2,7,3,1),  
  show.traits=1, show.contrasts=1,
  show.kme=0, use.tree=0,
  colorHeight = 0.5
)

par(mfrow=c(1,1),mar=c(6.5,9,4,2), cex=1.4)
plotModuleTraitHeatmap(
  wgcna, cluster=TRUE, multi=TRUE,
  main="Module-Trait Heatmap", setpar=FALSE,
  transpose=TRUE, text=FALSE)

## multi-omics ME clustering
par(mfrow=c(1,1),mar=c(6.5,9,4,2), cex=1.4)
plotEigenGeneAdjacencyHeatmap(
  wgcna$layers, multi=TRUE,
  plotDendro=0, plotHeatmap=1,
  add_traits = 1, add_me = TRUE,
  nmax = 30, text = FALSE, pstar = TRUE,
  mar1 = c(6, 4.5, 1.8, 0),
  mar2 = c(9, 12, 4, 2),
  cex.lab=1.2, colorlabel=TRUE
)

##--------------------------------------------------------------------
## Biomarkers
##--------------------------------------------------------------------
load_all()

topmodules <- getTopModules(wgcna, topratio=0.9, multi=TRUE)
topmodules

annot <- pgx$genes
names(wgcna)
annot=NULL; module=NULL,
psig=0.05; ntop=40; level=NULL;
rename="symbol"

top <- getTopTables(wgcna, psig=0.05, annot = annot, module=topmodules,
  rename="gene_title")

names(top)
names(top$genes)
lapply(top$genes, head)
lapply(top$sets, head)

QQ <- calculateCompoundSignificance(wgcna) 
names(QQ)
head(QQ)
k=3
Q <- QQ[[k]]

par(mfrow=c(2,1), mar=c(8,6,3,0), cex=1)
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


##--------------------------------------------------------------------
## LASAGNA
##--------------------------------------------------------------------

if(require("lasagna")) {
  
  library(igraph)
  lasagna.create_model <- playbase::lasagna.create_model
  lasagna.solve <- playbase::lasagna.solve
  lasagna.multisolve <- playbase::lasagna.multisolve  
  plotMultiPartiteGraph2 <- playbase::plotMultiPartiteGraph2

  names(wgcna)
  (wgcna$me.genes)
  
  ww <- wgcna$layers
  names(ww)
  ww <- ww[c("mir","gx","px")]
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
  
  eplot <- function(g,w=5,...) plot(g, edge.width=w*abs(E(g)$weight, ...))
  eplot(lasagna$graph)
  colnames(lasagna$Y)
  
  solved <- lasagna::solve(
    lasagna,
    #pheno = "daf2_vs_wildtype",
    #pheno = "activated=act",
    pheno = "condition=Her2",    
    min_rho = 0.01,
    max_edges = 100,
    value.type = "rho",
    fc.weights = TRUE,
    sp.weight = FALSE,
    prune = FALSE
  )
  eplot(solved,w=8)
  
  par(mfrow=c(1,1), mar = c(0,0,0,0), cex=1.4)
  lasagna::plot_multipartite(
    solved,
    min.rho = 0.01,
    ntop = -1,
    cex.label = 1.0,
    vx.cex = 3.4,
    edge.cex = 1.6,
    edge.alpha = 0.5,
    edge.gamma = 0.66,
    edge.sign = "both",
    edge.type = "inter",
    yheight = 0.8,
    xdist = 3,
    normalize.edges = 1,
    strip.prefix2 = TRUE,
    prune = FALSE
  )  
}


##--------------------------------------------------------------------
## LLM
##--------------------------------------------------------------------
load_all()

Sys.getenv("GEMINI_API_KEY")
Sys.getenv("GROQ_API_KEY")
Sys.getenv("XAI_API_KEY")

llm="gpt-5-nano"
llm="xai:grok-4-1-fast-non-reasoning"
llm="groq:llama-3.1-8b-instant"
llm="groq:openai/gpt-oss-120b"
llm="groq:openai/gpt-oss-20b"

ai.ask("hello", llm)

ai_model=llm
names(wgcna)
graph = wgcna$graph
annot = pgx$genes

cat(">>>> creating report with",llm,"...\n")
rpt <- WGCNAplus::create_report(
  wgcna, ai_model=llm, 
  ##graph = wgcna$graph,
  topratio=0.85, psig=0.05, annot=annot,
  multi=TRUE, verbose=1) 

names(rpt)
names(rpt$descriptions)
names(rpt$summaries)
cat(rpt$report)
cat(rpt$diagram)  

DiagrammeR::grViz(rpt$diagram)  
  
WGCNAplus::create_infographic(
  rpt$report, diagram=rpt$diagram,
  model = "gemini-3-pro-image-preview",
  api_key = Sys.getenv("GEMINI_API_KEY"),
  filename = "moxbrca-infographic.png"
)

modules <- names(rpt$summaries)
modules
m=modules[1]
for(m in modules) {
  fn <- paste0("moxbrca-module-",m,".png")
  message("saving to ", fn)
  create_module_infographic(
    rpt, module = m,
    model = "gemini-3-pro-image-preview",
    api_key = Sys.getenv("GEMINI_API_KEY"),
    filename = fn
  )
}



