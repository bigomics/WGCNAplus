##-----------------------
## Single-omics WGCNA
##-----------------------
library(devtools)
load_all()
setwd("vignettes/")

X <- read.csv("./data/Geiger/counts_so.csv", row.names = 1)
samples <- read.csv("./data/Geiger/samples_so.csv", row.names = 1)
contrasts <- read.csv("./data/Geiger/contrasts_so.csv", row.names = 1)
GMT <- readRDS("./data/Geiger/gmt.RDS")
annot <- read.csv("./data/Geiger/annot.csv", row.names = 1)

wgcna <- WGCNAplus::computeWGCNA(X, samples)
wgcna <- WGCNAplus::computeModuleEnrichment(wgcna, GMT = GMT, annot = annot)
names(wgcna)

par(mfrow = c(1, 2), cex = 2)
WGCNAplus::plotDendroAndColors(wgcna, marAll = c(2, 7, 3, 1), show.traits = 1,
  show.contrasts = 1, show.kme = 0, use.tree = 0, colorHeight = 0.5, main = "")

par(mfrow = c(1, 1), mar = c(6.5, 9, 4, 2))
WGCNAplus::plotModuleTraitHeatmap(wgcna, main = "Module-Trait Heatmap",
  cluster = TRUE, setpar = FALSE, transpose = TRUE, text = FALSE)

WGCNAplus::plotEigenGeneAdjacencyHeatmap(wgcna, multi = FALSE, ## multi=TRUE for mox
  plotDendro = 0, plotHeatmap = 1, add_traits = 1, add_me = TRUE,
  nmax = 30, text = FALSE, pstar = TRUE,
  mar1 = c(6, 4.5, 1.8, 0), mar2 = c(9, 12, 4, 2),
  cex.lab = 1.2, colorlabel = TRUE)


##-----------------------
## Biomarkers
##-----------------------
topmodules <- WGCNAplus::getTopModules(wgcna, topratio = 0.9, multi = FALSE)
top <- WGCNAplus::getTopTables(wgcna, annot = annot, module = topmodules, rename = "gene_title")
Q <- WGCNAplus::calculateCompoundSignificance(wgcna)

names(top)
names(top$genes)
lapply(top$genes, head)
lapply(top$sets, head)
lapply(top$pheno, head)
head(Q)

par(mfrow = c(2, 1), mar = c(8, 6, 3, 0), cex = 1.4)

Q <- Q[order(-Q$score1), ]
Q1 <- head(Q, 60)
barplot(Q1$score1, col = Q1$color,
  ylab = "significance (score1)",
  main="biomarker significance (MM*TS)",
  cex.names = 0.85, names.arg = rownames(Q1), las = 3)

Q <- Q[order(-Q$score2), ]
Q1 <- head(Q, 60)
barplot(Q1$score2, col = Q1$color,
  ylab = "significance (score2)",
  main = "biomarker significance (MM*FC)",
  cex.names = 0.85, names.arg = rownames(Q1), las = 3)


##----------------
## LASAGNA
##----------------
library(igraph)

ll <- list(gx = t(wgcna$net$MEs))
xdata <- list(X = ll, samples = wgcna$datTraits)

lasagna <- lasagna::create_model(
  xdata,
  meta.type = "expanded",
  ntop = 2000,
  nc = 20,
  add.sink = FALSE,
  intra = TRUE,
  fully_connect = FALSE,
  add.revpheno = TRUE
)

solved <- lasagna::solve(
  lasagna,
  pheno = "activated=act",
  min_rho = 0.3,
  max_edges = 100,
  value.type = "rho",
  fc.weights = TRUE,
  sp.weight = FALSE,
  prune = TRUE
)

wgcna$graph <- lasagna$graph


#eplot <- function(g, w,...) plot(g, edge.width=w*abs(E(g)$weight, ...))
#eplot(lasagna$graph,w=5)
#colnames(lasagna$Y)

par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), cex = 1.4)
lasagna::plot_multipartite(
  solved,
  min.rho = 0.6,
  ntop = -1,
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



##---------------
## LLM
##---------------
#Sys.getenv("GEMINI_API_KEY")
#Sys.getenv("GROQ_API_KEY")
#Sys.getenv("XAI_API_KEY")

llm = "xai:grok-4-1-fast-non-reasoning"
WGCNAplus::ai.ask(question = "hello", model = llm, engine = "ellmer")
WGCNAplus::ai.ask(question = "hello", model = llm, engine = "tidyprompt")

llm = "groq:llama-3.1-8b-instant"
WGCNAplus::ai.ask(question = "hello", model = llm, engine = "ellmer")
WGCNAplus::ai.ask(question = "hello", model = llm, engine = "tidyprompt")

llm = "groq:openai/gpt-oss-120b"
WGCNAplus::ai.ask(question = "hello", model = llm, engine = "ellmer")
WGCNAplus::ai.ask(question = "hello", model = llm, engine = "tidyprompt")

message(">>>> creating report with", llm, "...\n")

rpt <- WGCNAplus::create_report(
  wgcna,
  ai_model = llm,
  topratio = 0.85,
  graph = NULL,
  annot = annot,
  multi = FALSE
) 

rpt <- WGCNAplus::create_report(
  wgcna,
  ai_model = llm,
  topratio = 0.85,
  graph = wgcna$graph,
  annot = annot,
  multi = FALSE
) 

names(rpt)
names(rpt$descriptions)
names(rpt$summaries)
cat(rpt$report)
cat(rpt$diagram)  

DiagrammeR::grViz(rpt$diagram)  

model <- "gemini-3.1-flash-image-preview"
api_key <- Sys.getenv("GEMINI_API_KEY")
ff <- "infographic.png"

WGCNAplus::create_infographic(
  report = rpt$report,
  diagram = rpt$diagram,
  model = model,
  filename = ff
)

modules <- names(rpt$summaries); modules

for (m in modules) {
  fn <- paste0("module-", m, ".png")
  WGCNAplus::create_module_infographic(
    rpt = rpt,
    module = m,
    model = model,
    api_key = api_key,
    filename = fn
  )
}


prompt <- paste("Create a graphical abstract according to the given diagram and information in the WGCNA report. Use scientific visual style like Nature journals. Illustrate biological concepts with small graphics. \n\n", rpt$report, "\n---------------\n\n", rpt$diagram)

nchar(prompt)
cat(prompt)

ff <- "image-gemini1.png"
WGCNAplus::ai.create_image_gemini(prompt, model = model, filename = ff)

model <- "grok-imagine-image"
api_key <- Sys.getenv("XAI_API_KEY")
ff <- "image-grok-imagine5.png"
WGCNAplus::ai.create_image_grok(
  prompt,
  model = model,
  aspect_ratio = "16:9",
  api_key = api_key,
  base_url = "https://api.x.ai/v1",
  filename = ff,
  user = NULL,
  organization = NULL
)

model <- "gpt-image-1.5"
ff <- "image-openai-v15.png"
WGCNAplus::ai.create_image_openai(
  prompt,
  model = model,
  filename = ff,
  user = NULL,
  organization = NULL
)

