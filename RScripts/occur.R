
library(ccrepe)


otu <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
metadat[stoolsamp,]
spptab <- colnames(otu) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu.gut <- otu[-which(rowSums(otu[,spptab]) == 0),spptab]
dim(otu.gut)
ogut <- t(apply(otu.gut, 2, function(x) x/sum(x)))

test <- ccrepe(ogut)
