### gpl341 validation by gds
# start with jpmEset.rdata: lapply(jpmEset, exprs): 15923 probe values per sample
load("~/GitHub/stevia/data/jpmEset.rdata")
rawEsets <- jpmEset[c("stromal", "r_hippo", "oculomotor", "spinal_cord")]
rm(jpmEset)
# apply impute::impute.knn k=10 nearest neihbors from max block size of 1500 genes, > 50% missing in block => column mean; $data is result
# apply PGSEA::aggregateExprs combining current probe SYMBOL annotations with "median":  10000 gene symbols per sample
# convert RNA expression levels to log2
# turn back to expression set with associated sample binary age variable: 0 = young, 1 = old

load("~/GitHub/stevia/data/symEsets.rdata")
gpl341Esets <- symEsets[c("stromal", "r_hippo", "oculomotor", "spinal_cord")]
rm(symEsets)
sapply(gpl341Esets, function(x) dim(exprs(x)))
#      stromal r_hippo oculomotor spinal_cord
# [1,]   10000   10000      10000       10000
# [2,]       3      78          9           9

library("ggplot2")
# genemeta analysis
library("GeneMeta")
load("~/GitHub/stevia/data/GMresults.rdata")
gpl341.fdr <- gpl.fdr$GPL341
gpl341.zscore <- gpl.zscore$GPL341
rm(gpl.fdr, gpl.zscore)
geneplotter::histStack(as.data.frame(gpl341.zscore[,1:4]), breaks = 100)
geneplotter::histStack(as.data.frame(gpl341.zscore[,2:4]), breaks = 100)
psych::multi.hist(gpl341.zscore[, 2:5], dcol = c("black", "red"), main = "differential z-score density\nnormal fit in red")

gpl341.z2 <- zScores(gpl341Esets, lapply(gpl341Esets, function(x) x$age))
geneplotter::histStack(as.data.frame(gpl341.z2[,1:4]), breaks = 100)
geneplotter::histStack(as.data.frame(gpl341.z2[,2:4]), breaks = 100)
psych::multi.hist(gpl341.z2[, 2:5], dcol = c("black", "red"), main = "differential z-score density\nnormal fit in red")

# r_hippo controls only
gpl341e2 <- gpl341Esets
gpl341e2$r_hippo <- gpl341e2$r_hippo[, c(1:10, 30:39)]
gpl341e2.zscore <- zScores(gpl341e2, lapply(gpl341e2, function(x) x$age))
geneplotter::histStack(as.data.frame(gpl341e2.zscore[,1:4]), breaks = 100)
geneplotter::histStack(as.data.frame(gpl341e2.zscore[,2:4]), breaks = 100)
psych::multi.hist(gpl341e2.zscore[, 2:5], dcol = c("black", "red"), main = "differential z-score density\nnormal fit in red")

# run without $stromal
gpl341e2.z2 <- zScores(gpl341e2[2:4], lapply(gpl341e2[2:4], function(x) x$age))
geneplotter::histStack(as.data.frame(gpl341e2.z2[,1:3]), breaks = 100)
psych::multi.hist(gpl341e2.z2[, 1:4], dcol = c("black", "red"), main = "differential z-score density\nnormal fit in red")
qq_plot(gpl341e2.z2[, "Qvals"], length(grep("zSco_Ex",colnames(gpl341e2.z2))), "GPL341b")

# make lattice plot comparing data sets and their combination
library(XDE)
GPL341e2.pos <- symbolsInteresting(rankingStatistic = pnorm(na.omit(gpl341e2.z2[, "zSco"])), percentile = .95)
png(file = "C:/Users/user/Documents/GitHub/XDE/GPL341e2_pairs.png", width = 1000, height = 1000)
pairs(gpl341e2.z2[GPL341e2.pos$order, 1:4], pch = GPL341e2.pos$pch, col = GPL341e2.pos$col, bg = GPL341e2.pos$bg, upper.panel = NULL, cex = GPL341e2.pos$cex, main = "GPL341b")
dev.off()

GPL341e2.neg <- symbolsInteresting(rankingStatistic = pnorm(1 - na.omit(gpl341e2.z2[, "zSco"])), percentile = .95)
png(file = "C:/Users/user/Documents/GitHub/XDE/GPL341e2_pairsN.png", width = 1000, height = 1000)
pairs(gpl341e2.z2[GPL341e2.neg$order, 1:4], pch = GPL341e2.neg$pch, col = GPL341e2.neg$col, bg = GPL341e2.neg$bg, upper.panel = NULL, cex = GPL341e2.neg$cex, main = "GPL341b")
dev.off()

# save symEsets with $r_hippo controls only
symEsets$r_hippo <- symEsets$r_hippo[, c(1:10, 30:39)]
save(symEsets, file = "~/GitHub/stevia/data/symEsets.rdata")
