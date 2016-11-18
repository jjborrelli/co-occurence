library(ccrepe)


otu <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
metadat[stoolsamp,]
spptab <- colnames(otu) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu.gut <- otu[-which(rowSums(otu[,spptab]) == 0),spptab]
dim(otu.gut)
ogut <- t(apply(otu.gut, 2, function(x) x/sum(x)))

err <- c()
for(i in 1:100){
  m <- abs(matrix(rnorm(100*80), 100 , 80))
  m <- (apply(m, 2, function(x) x/sum(x)))
  test <- ccrepe(m)
  err[i] <- sum(test$p.values < 0.05, na.rm = T)/(80*79/2)
}

hist(err)
