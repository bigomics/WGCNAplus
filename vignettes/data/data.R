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
## multi-omics data
##--------------------------------------------------------------------

pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/mox-brca.pgx")
write.csv(pgx$X, file="brca/expression.csv")
write.csv(pgx$samples, file="brca/samples.csv")
write.csv(pgx$contrasts, file="brca/contrasts.csv")
save(pgx$GMT, file="brca/gmt.rda")


pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/mox-geiger.pgx")
write.csv(pgx$X, file="geiger/expression.csv")
write.csv(pgx$samples, file="geiger/samples.csv")
write.csv(pgx$contrasts, file="geiger/contrasts.csv")
dim(pgx$GMT)
object.size(pgx$GMT)
sel <- grep("^GO_|HALLMARK",colnames(pgx$GMT))
length(sel)
GMT <- pgx$GMT[,sel]
object.size(GMT)
save(GMT, file="geiger/gmt.rda")
