### XDE this uses jpmEset from "virtualArray.R"
library("XDE")

# make pheno dataframes.  fix lung age catagories!
age.cat$lung <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1 ,1)
age.pdata <- lapply(age.cat, function(x) data.frame(age = x, row.names = colnames(jpmEset[[names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]]]])))

# try imputing then combiningin
load("~/GitHub/stevia/data/jpmEset.rdata")
aggEsets <- lapply(jpmEset, exprs)

# check for missing expression values:  kidney1 & 2, CA1_hipp2 (1/3), m_heart (1/7) myo_pyogen (1/4) are missing a significant portion
# 13 of 26 data sets have 4 to 1900 rows missin > 50% with imputation by means rather than 10 nearest neihbors
sapply(aggEsets, function(x) summary(x))

# number of probes/rows
sapply(aggEsets, function(x) dim(x)[1])
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle      kidney1      kidney2      m_brain 
#        12625        22283        22645        22283        22645        12626        22690         6584         6595        45101 
#      m_hippo        liver      m_heart         lung      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord 
#        12488        12488        12654        45101        45101        45101        12488        15923        15923        15923 
#   oculomotor  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#        15923         8799         8799         8799         8799         8799 

# impute NA expression values
for(n in names(aggEsets)) aggEsets[[n]] <- impute::impute.knn(aggEsets[[n]])
sapply(aggEsets, function(x) summary(is.na(x$data)))
# no NAs! and row names are probes
sapply(aggEsets, function(x) dim(x)[1])
sapply(aggEsets, function(x) head(rownames(x)))

# get rid of random seed info supplied by impute package
aggEsets <- lapply(aggEsets, function(x) x$data)

#now aggregate imputed results
aggEsets <- mapply(function(x, y) PGSEA::aggregateExprs(x, package = y, using = "SYMBOL", FUN = median), aggEsets, dbpkg)

# much smaller & still no NAs, it worked!
sapply(aggEsets, function(x) dim(x)[1])
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle      kidney1      kidney2      m_brain 
#         8632        12494         9648        12494         9648         8616           98         5106         3299        13016 
#      m_hippo        liver      m_heart         lung      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord 
#         8724         8724         8622        13016        13016        13016         8724        10000        10000        10000 
#   oculomotor  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#        10000         4750         4750         4750         4750         4750 
# they are the same!
aggEsets <- lapply(aggEsets, log2)
save(aggEsets, file = "data/aggEsets.rdata")

# turn back into Esets
symEsets <- mapply(function(x, y) new("ExpressionSet", exprs = x, phenoData = y), aggEsets, lapply(age.pdata, function(x) as(x, "AnnotatedDataFrame")))
save(symEsets, file = "data/symEsets.rdata")

# make ExpressionSetLists by chip
chiplist <- lapply(chips, get)
names(chiplist) <- chips
chipEsets <- lapply(chiplist, function(x) symEsets[x])

# make expression set list objects
chipEsets <- lapply(chipEsets, function(x) new("ExpressionSetList", x))
sapply(chipEsets, validObject)
# GPL1261  GPL341   GPL81   GPL85   GPL96   GPL97 
#    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 

# final validity test from xde vignette to confirm all pData matches.  it works!
sapply(chipEsets, function(x) sapply(x, varLabels))

# make list of XdeParameter objects with empirical starting values
burnParams <- list()
for(n in names(chipEsets)) {
  burnParams[[n]] <- new("XdeParameter", esetList = chipEsets[[n]], phenotypeLabel = "age", iterations = 10,   directory = paste("data/", n, "burnLogs", sep = "")
, one.delta=FALSE, verbose = TRUE)
}

chipBurn <- mapply(xde, burnParams, chipEsets)

# graphs of c2 for each chipBurn set
lapply(chipBurn, function(x) plot.ts(x$c2, ylab = "c2", xlab = "iterations", plot.type = "single", main = names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]]))

for(n in names(burnParams)) {
  directory(burnParams[[n]]) <- paste("~/Github/XDE/", n, "chipLogs", sep = "")
  directory(chipBurn[[n]]) <- paste("~/Github/XDE/", n, "chipLogs", sep = "")
  
}

chipBurnBES <- lapply(chipBurn, calculateBayesianEffectSize)
chipBurnPosAv <- lapply(chipBurn, calculatePosteriorAvg, burnin = 8)

# z scores from GeneMeta and histogram of z values for arrays
chipBurnZ <- lapply(chipEsets, function(x) ssStatistic("z", "age", x))
lapply(chipBurnZ, function(x) hist(x[, "zSco"], main = names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], xlab = "z score"))
chip.pos <- lapply(chipBurnZ, function(x) symbolsInteresting(rankingStatistic = pnorm(x[, "zSco"]), percentile = .95))
 op.conc <- symbolsInteresting(rankingStatistic=postAvg[, "concordant"])

# GPL341 breaks down because stromal has all 0 for gene z scores !?!?
chipNsamples <- lapply(chipEsets, nSamples)
chipNsamples <- lapply(chipNsamples, sum)
chipBurnZx <- list()
names(chipBurnZx) <- names(chipEsets)
for(n in names(chipBurnZ)) {
  chipBurnZx[[n]] <- try(xsScores(chipBurnZ[[n]][, 1:length(chipNsamples[[n]])], chipNsamples[[n]]))
}


# this takes > 18 hours on celeron, waiting for bed stuy results
chipParams <- burnParams
for(n in names(chipEsets)) {
  iterations(chipParams[[n]]) <- 1000
  thin(chipParams[[n]]) <- 2
  directory(chipParams[[n]]) <- paste("~/Github/XDE/", n, "chipLogs", sep = "")
  firstMcmc(chipParams[[n]]) <- lastMcmc(chipBurn[[n]])
}

date()
chipRun <- mapply(xde, chipParams, chipEsets)
date()