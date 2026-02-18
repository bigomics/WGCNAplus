library(devtools)
load_all()

## Imports
gset.rankcor <- playbase::gset.rankcor
mat2gmt <- playbase::mat2gmt
##mofa.merge_data2 <- playbase::mofa.merge_data2
ai.ask <- playbase::ai.ask
ai.create_image_gemini <- playbase::ai.create_image_gemini

##--------------------------------------------------------------------
## Compute single-omics
##--------------------------------------------------------------------

pgx <- playdata::GEIGER_PGX
pgx <- playbase::pgx.initialize(pgx)

wgcna <- computeWGCNA(pgx$X, pgx$samples)

wgcna <- computeModuleEnrichment(
  wgcna,
  GMT = pgx$GMT,
  annot = pgx$genes
)

names(wgcna)
names(wgcna$gsea)
head(wgcna$gsea[[1]])

##--------------------------------------------------------------------
## Figures
##--------------------------------------------------------------------

par(mfrow=c(1,1),cex=1)
plotDendroAndColors(
  wgcna, marAll=c(2,7,3,1),  
  show.traits=1, show.contrasts=1,
  show.kme=0, use.tree=0,
  colorHeight = 0.5, main=""
)

par(mfrow=c(1,1),mar=c(6.5,9,4,2))
plotModuleTraitHeatmap(wgcna, cluster=TRUE,
  main="Module-Trait Heatmap", setpar=FALSE,
  transpose=TRUE, text=FALSE)

## multi-omics ME clustering
plotEigenGeneAdjacencyHeatmap(
  wgcna, multi=FALSE,
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

topmodules <- getTopModules(wgcna, topratio=0.9, multi=FALSE)
topmodules

annot <- pgx$genes
names(wgcna)

top <- getTopTables(wgcna,
  psig=0.05, annot = annot, module=topmodules,
  rename="gene_title")
names(top)
names(top$genes)
lapply(top$genes, head)
lapply(top$sets, head)

names(wgcna)
names(wgcna$layers)
names(wgcna$layers[[1]])

Q <- calculateCompoundSignificance(wgcna) 
head(Q)

par(mfrow=c(2,1), mar=c(8,6,3,0), cex=1.4)
ntop=60
Q <- Q[order(-Q$score1),]
Q1 <- head(Q,ntop)
barplot( Q1$score1, col=Q1$color,
  ylab = "significance (score1)", main="biomarker significance (MM*TS)",
  cex.names=0.85, names.arg = rownames(Q1), las=3 )
Q <- Q[order(-Q$score2),]
Q1 <- head(Q,ntop)
barplot( Q1$score2, col=Q1$color,
  ylab = "significance (score2)", main="biomarker significance (MM*FC)",
  cex.names=0.85, names.arg = rownames(Q1), las=3 )

##--------------------------------------------------------------------
## LASAGNA
##--------------------------------------------------------------------

if(require("lasagna")) {
  
  library(igraph)
  names(wgcna)
  
  xdata <- list(
    X = list(gx = t(wgcna$net$MEs)),
    samples = wgcna$datTraits
  )
  
  lasagna.create_model <- playbase::lasagna.create_model
  lasagna.solve <- playbase::lasagna.solve
  plotMultiPartiteGraph2 <- playbase::plotMultiPartiteGraph2
  
  lasagna <- lasagna.create_model(
    xdata,
    pheno = "expanded",
    ntop = 2000,
    nc = 20,
    add.sink = FALSE,
    intra = TRUE,
    fully_connect = FALSE,
    add.revpheno = TRUE
  )
  
  eplot <- function(g,w,...) plot(g, edge.width=w*abs(E(g)$weight, ...))
  eplot(lasagna$graph,w=5)
  colnames(lasagna$Y)
  
  solved <- lasagna.solve(
    lasagna,
    #pheno = "daf2_vs_wildtype",
    pheno = "activated=act",
    min_rho = 0.3,
    max_edges = 100,
    value.type = "rho",
    fc.weights = TRUE,
    sp.weight = FALSE,
    prune = TRUE
  )
  
  par(mfrow=c(1,1), mar = c(0,0,0,0), cex=1.4)
  plotMultiPartiteGraph2(
    solved,
    min.rho = 0.6,
    ntop = -1,
    ##labels = labels,
    cex.label = 1.0,
    vx.cex = 3.4,
    edge.cex = 1.1,
    edge.alpha = 0.3,
    edge.sign = "both",
    edge.type = "inter",
    yheight = 3,
    xdist = 3,
    normalize.edges = 1,
    strip.prefix2 = TRUE,
    prune = FALSE
  )

  names(lasagna)
  class(lasagna$graph)
  wgcna$graph <- lasagna$graph
  
}


##--------------------------------------------------------------------
## LLM
##--------------------------------------------------------------------
load_all("~/Playground/playbase")

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
  wgcna, ai_model=llm, topratio=0.85,
  ##graph = wgcna$graph,
  psig=0.05, annot=annot,
  multi=FALSE, verbose=1) 

names(rpt)
names(rpt$descriptions)
names(rpt$summaries)
cat(rpt$report)
cat(rpt$diagram)  

DiagrammeR::grViz(rpt$diagram)  
  
WGCNAplus::create_infographic(
  rpt$report, diagram=rpt$diagram,
  model = "gemini-3-pro-image-preview",
  filename = "infographic.png"
)

modules <- names(rpt$summaries)
modules
m=modules[1]
for(m in modules) {
  fn <- paste0("module-",m,".png")
  message("saving to ", fn)
  create_module_infographic(
    rpt, module = m,
    model = "gemini-3-pro-image-preview",
    api_key = Sys.getenv("GEMINI_API_KEY"),
    filename = fn
  )
}


prompt = "astronaut on horse on the moon"

prompt <- paste("Create a graphical abstract according to the given diagram and information in the WGCNA report. Use scientific visual style like Nature journals. Illustrate biological concepts with small graphics. \n\n", rpt$report, "\n---------------\n\n", rpt$diagram)
nchar(prompt)
cat(prompt)

load_all("~/Playground/playbase")

ai.create_image_gemini(
  prompt,
  model = "gemini-3-pro-image-preview",
  #aspectRatio = "16:9",  
  filename = "image-gemini1.png"
)

ai.create_image_grok(
  prompt,
  model = "grok-imagine-image",   
  aspect_ratio = "16:9",
  api_key = Sys.getenv("XAI_API_KEY"),
  base_url = "https://api.x.ai/v1",
  filename = "image-grok-imagine5.png",
  user = NULL, organization = NULL
)

ai.create_image_openai(
  prompt,
  model = "gpt-image-1.5",   
#  format = c("file","base64","raw")[1],
#  api_key = Sys.getenv("OPENAI_API_KEY"),
#  base_url = "https://api.x.ai/v1",
  filename = "image-openai-v15.png",
  user = NULL, organization = NULL
)

