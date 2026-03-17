##-------------------------------------------------
## Compare WGCNA+ vs. WGCNA FastTom Calculation
##-------------------------------------------------
library(devtools)
load_all()
setwd("vignettes/")

X <- as.matrix(read.csv("./data/Geiger/counts_so.csv", row.names = 1))
samples <- read.csv("./data/Geiger/samples_so.csv", row.names = 1)

minmodsize = 20
power = 6
networkType = "signed"
TOMType = "signed"
cor <- WGCNA::cor

pdf("~/Desktop/WGCNAplus/WGCNAplus/vignettes/paper_figures/Fig1.pdf")
par(mfrow = c(3,2), mgp = c(3.3, 0.5, 0), tcl = -0.1, cex.lab = 1.7,
  cex.axis = 1.4, las = 1, cex.main = 2.2)
i=1; ss=c(1500, 3500, nrow(X))
for(i in 1:length(ss)) {

  X1 <- X
  sdx <- matrixStats::rowSds(X1, na.rm = TRUE)
  X1 <- X1[sdx > 0.1 * mean(sdx, na.rm = TRUE), ] ## filter low SD??
  X1 <- X1[order(-matrixStats::rowSds(X1, na.rm = TRUE)), ]
  X1 <- utils::head(X1, ss[i])

  datExpr <- t(X1)
  adjacency <- WGCNA::adjacency(datExpr, power = power, type = networkType)
  adjacency[is.na(adjacency)] <- 0

  st1 <- system.time(TOM <- WGCNA::TOMsimilarity(adjacency, TOMType = TOMType))['elapsed']
  st2 <- system.time(TOM <- WGCNAplus::fastTOMsimilarity(adjacency, lowrank = 20))['elapsed']
  
  MEM1 <- peakRAM::peakRAM(TOM <- WGCNA::TOMsimilarity(adjacency, TOMType = TOMType))$Peak_RAM_Used_MiB
  MEM2 <- peakRAM::peakRAM(TOM <- WGCNAplus::fastTOMsimilarity(adjacency, lowrank = 20))$Peak_RAM_Used_MiB

  nn <- c("WGCNA", "WGCNA+")
  mm <- paste0("Features = ",ncol(datExpr))
  barplot(c(st1, st2), names.arg = nn, cex.names = 2, ylab = "Time (s)", main = mm)
  barplot(c(MEM1, MEM2), names.arg = nn, cex.names = 2, ylab = "Peak RAM (MiB)", main = mm, col = "skyblue")
  
  message("Size: ", ncol(datExpr), ". Completed/n/n")
}
dev.off()
