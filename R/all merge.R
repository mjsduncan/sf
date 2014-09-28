### combine all data sets by gene homology

library("Biobase")
load("~/GitHub/stevia/data/probe2homo.rdata")
load("~/GitHub/stevia/data/symEsets.rdata")

probe2homo$kidney2 <- subset(probe2homo$kidney2, GEOgene != "--Control")

homoEsets <- list()
for(n in names(probe2homo)) {
  featureMap <- unique(na.omit(probe2homo[[n]][probe2homo[[n]]$SYMBOL %in% featureNames(symEsets[[n]]), c(3, 6)]))
  homoEsets[[n]] <- symEsets[[n]][featureMap$SYMBOL,]
  print(c(n, summary(duplicated(featureNames(homoEsets[[n]]))), summary(duplicated(featureMap$SYMBOL))))
#   featureNames(homoEsets[[n]]) <- featureMap$symbol
}

sapply(symEsets[7:26], dim)
#          muscle kidney1 kidney2 m_brain m_hippo liver m_heart  lung cochlea hemato_stem myo_progen r_hippo stromal
# Features     98    5106    3299   13016    8724  8724    8622 13016   13016       13016       8724   10000   10000
# Samples      10      10      10       6      23     7      12    15       6           8          4      78       3
#          spinal_cord oculomotor skeletal_ms extraoc_ms laryngeal_ms r_heart CA1_hipp2
# Features       10000      10000        4750       4750         4750    4750      4750
# Samples            9          9          12         12            9      11        29

sapply(homoEsets, dim)
#          muscle kidney1 kidney2 m_brain m_hippo liver m_heart  lung cochlea hemato_stem myo_progen r_hippo stromal
# Features    107    5437    4642   19924   10561 10561    9970 19924   19924       19924      10561   11298   11298
# Samples      10      10      10       6      23     7      12    15       6           8          4      78       3
#          spinal_cord oculomotor skeletal_ms extraoc_ms laryngeal_ms r_heart CA1_hipp2
# Features       11298      11298        6638       6638         6638    6638      6638
# Samples            9          9          12         12            9      11        29

for(n in names(homoEsets)) {
  HOMOesets[[n]] <- homoEsets[[n]]
}

for(n in names(homoEsets)) {
  featureMap[na.omit(match(featureNames(homoEsets[[n]]), featureMap$SYMBOL)), 2]
}

### rename & merge organism moses sets
allMoses <- lapply(orgMoses[2:3], function(x) x[, colnames(x) %in% featureMap$SYMBOL])
for(n in names(allMoses)) colnames(allMoses[[n]]) <- featureMap[match(colnames(allMoses[[n]]), featureMap$SYMBOL), 2]
allMoses <- mapply(function(x, y) cbind(x[, 1], y), orgMoses[2:3], allMoses)
allMoses[['Homo sapiens']] <- orgMoses[[1]]
for(i in 1:3) colnames(allMoses[[i]])[1] <- "age"
allMoses <- lapply(allMoses, function(x) x[, Reduce(intersect, lapply(allMoses, colnames))])
allMoses <- Reduce(rbind, allMoses)
