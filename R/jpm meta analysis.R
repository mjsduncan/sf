### jpm meta analysis
# uses jpmPos, jpmNeg from 'jpm analysis.R' and jpmPos.df, jpmNeg.df from 'annotation mapping.R'

# make vector of all genes linearly associated with age from all data sets
HsPos <- unique(unlist(lapply(jpmPos[1:6], function(x) x$jpmGene)))
nonHsPos <- unique(unlist(lapply(jpmPos.df, function(x) x$homSym)))
Pos <- unique(c(HsPos, nonHsPos))

HsNeg <- unique(unlist(lapply(jpmNeg[1:6], function(x) x$jpmGene)))
nonHsNeg <- unique(unlist(lapply(jpmNeg.df, function(x) x$homSym)))
Neg <- unique(c(HsNeg, nonHsNeg))

# average fraction of individual data sets differentially expressed 2 sided p value < .05
# average of fraction of significant probes are close to jpm result except more negatively correlated than positive...
summary(unlist(lapply(jpmExp.slope[-c(19, 20)], function (x) (sum((x[, 3] <= .025 | x[, 3] >= .975) & x[, 1] > 0, na.rm = TRUE) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02203 0.02927 0.03548 0.03862 0.04633 0.08201 

summary(unlist(lapply(jpmExp.slope[-c(19, 20)], function (x) (sum((x[, 3] <= .025 | x[, 3] >= .975) & x[, 1] < 0, na.rm = TRUE) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01890 0.02451 0.03842 0.03984 0.04987 0.08399 
