### jpm meta analysis
# uses jpmPos, jpmNeg from 'jpm analysis.R' and jpmPos.df, jpmNeg.df from 'annotation mapping.R'

# make vector of unique gene symbols linearly associated with age from all data sets
HsPos <- unique(unlist(lapply(jpmPos[1:6], function(x) x$jpmGene)))
nonHsPos <- unique(unlist(lapply(jpmPos.df, function(x) x$homSym)))
Pos <- unique(c(HsPos, nonHsPos))
sapply(list(HsPos, nonHsPos, Pos), length)
# [1] 3623 3503 6342 -> 3623 + 3503 - 6342 = 784 duplicates

HsNeg <- unique(unlist(lapply(jpmNeg[1:6], function(x) x$jpmGene)))
nonHsNeg <- unique(unlist(lapply(jpmNeg.df, function(x) x$homSym)))
Neg <- unique(c(HsNeg, nonHsNeg))
sapply(list(HsNeg, nonHsNeg, Neg), length)
# [1] 4062 3516 6711 -> 4062 + 3516 - 6711 = 867 duplicates

# average fraction of individual data sets differentially expressed 2 sided p value < .05
# average of fraction of significant probes are close to jpm result except more are negatively correlated than positive...
summary(unlist(lapply(jpmExp.slope, function (x) (sum((x[, 3] <= .025 | x[, 3] >= .975) & x[,1] > 0, na.rm = TRUE) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02203 0.02742 0.03529 0.03773 0.04516 0.08202 

summary(unlist(lapply(jpmPos, function (x) sum(x[, 4] <= .025 | x[, 4] >= .975))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   187.0   309.8   502.0   646.8   946.0  1522.0 

summary(unlist(lapply(jpmExp.slope, function (x) (sum((x[, 3] <= .025 | x[, 3] >= .975) & x[,1] < 0, na.rm = TRUE) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01890 0.02512 0.03843 0.03940 0.04906 0.08400 

summary(unlist(lapply(jpmNeg, function (x) sum(x[, 4] <= .025 | x[, 4] >= .975))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   169.0   349.5   562.5   652.7  1016.0  1499.0 

# total expression levels in all experiments
sum(unlist(lapply(jpmExp.slope, function(x) dim(x)[1])))
# [1] 474305

# total expression levels linearly related to age
sum(unlist(lapply(jpmSig, function(x) dim(x)[1])))
# [1] 33790

# get k, the number of experiments in which a gene is over/under expressed.  
kpos <- numeric(length(Pos))
names(kpos) <- Pos
kpos.hs <- numeric(length(Pos))
names(kpos.hs) <- Pos
kpos.nhs <- numeric(length(Pos))
names(kpos.nhs) <- Pos
for(g in Pos) {
  kpos.hs[g] <- sum(sapply(jpmPos[1:6], function(x) g %in% x$jpmGene), na.rm = TRUE)
  kpos.nhs[g] <- sum(sapply(jpmPos.df, function(x) g %in% x$homSym), na.rm = TRUE)
  kpos[g] <- kpos.hs[g] + kpos.nhs[g]
}

# get n, the number of experiments in which a gene is measured
npos <- numeric(length(Pos))
names(npos) <- Pos
npos.hs <- numeric(length(Pos))
names(npos.hs) <- Pos
npos.nhs <- numeric(length(Pos))
names(npos.nhs) <- Pos
for(g in Pos) {
  npos.hs[g] <- sum(sapply(probe2gene[1:6], function(x) g %in% x$IDENTIFIER), na.rm = TRUE)
  npos.nhs[g] <- sum(mapply(function(x, y) x$SYMBOL[match(g, x$homSym)] %in% y$IDENTIFIER, jpmPos.df, probe2gene[c(7:26)]))
  npos[g] <- npos.hs[g] + npos.nhs[g]
}

npos.hom <- numeric(length(Pos))
names(npos.hom) <- Pos
for(g in Pos) {
  npos.hom[g] <- sum(mapply(function(x, y) x$PROBEID[match(g, x$homSym)] %in% y$ID_REF, jpmPos.df, probe2gene[c(7:26)]))
}

# p value for null hypothesis binomial dist
### functions
# calculate cummulative binomial distribution
cbd <- function(k, n, p) {
  pcb <- numeric(1)
  for(i in k:n){
    `<-`(pcb, pcb + choose(n, i)*p^i*(1-p)^(n-i))
  }
  return(pcb)
}

sigPos <- mapply(cbd, kpos, npos, 0.03773)
summary(sigPos)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.02634 0.10900 0.10490 0.17490 1.00000 
length(sigPos)
# [1] 6342

sigPos.min <- sigPos[sigPos < .05]
summary(sigPos.min)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001677 0.013190 0.012400 0.019290 0.037730 
length(sigPos.min)
# [1] 1861


# do the same for negatively associated genes
# get k, the number of experiments in which a gene is over/under expressed.  
kneg <- numeric(length(Neg))
names(kneg) <- Neg
kneg.hs <- numeric(length(Neg))
names(kneg.hs) <- Neg
kneg.nhs <- numeric(length(Neg))
names(kneg.nhs) <- Neg
for(g in Neg) {
  kneg.hs[g] <- sum(sapply(jpmNeg[1:6], function(x) g %in% x$jpmGene), na.rm = TRUE)
  kneg.nhs[g] <- sum(sapply(jpmNeg.df, function(x) g %in% x$homSym), na.rm = TRUE)
  kneg[g] <- kneg.hs[g] + kneg.nhs[g]
}

# get n, the number of experiments in which a gene is measured
nneg <- numeric(length(Neg))
names(nneg) <- Neg
nneg.hs <- numeric(length(Neg))
names(nneg.hs) <- Neg
nneg.nhs <- numeric(length(Neg))
names(nneg.nhs) <- Neg
for(g in Neg) {
  nneg.hs[g] <- sum(sapply(probe2gene[1:6], function(x) g %in% x$IDENTIFIER), na.rm = TRUE)
  nneg.nhs[g] <- sum(mapply(function(x, y) x$SYMBOL[match(g, x$homSym)] %in% y$IDENTIFIER, jpmNeg.df, probe2gene[c(7:26)]))
  nneg[g] <- nneg.hs[g] + nneg.nhs[g]
}

nneg.hom <- numeric(length(Neg))
names(nneg.hom) <- Neg
for(g in Neg) {
  nneg.hom[g] <- sum(mapply(function(x, y) x$PROBEID[match(g, x$homSym)] %in% y$ID_REF, jpmNeg.df, probe2gene[c(7:26)]))
}

sigNeg <- mapply(cbd, kneg, nneg, 0.03940)
summary(sigNeg)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.02856 0.11360 0.11010 0.18210 1.00000 
length(sigNeg)
# [1] 6711

sigNeg.min <- sigNeg[sigNeg < .05]
summary(sigNeg.min)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001899 0.014340 0.013300 0.020950 0.039400 
length(sigNeg.min)
# [1] 1861

# make data frame to compare with original results
agePos <- merge(jpmEz, data.frame(SYMBOL = names(sigPos.min), metaP = sigPos.min, stringsAsFactors = FALSE))
ageNeg <- merge(jpmEz, data.frame(SYMBOL = names(sigNeg.min), metaP = sigNeg.min, stringsAsFactors = FALSE))
# write.csv(agePos, file = "agePos.csv")
# write.csv(ageNeg, file = "ageNeg.csv")

# get original tables
over.expressed <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/aging signatures/over expressed.csv", stringsAsFactors = FALSE)
under.expressed <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/aging signatures/under expressed.csv", stringsAsFactors = FALSE)

write.csv(agePos[agePos$ENTREZID %in% over.expressed$EntrezGeneID,], file = "agePosInt.csv")
write.csv(ageNeg[agePos$ENTREZID %in% over.expressed$EntrezGeneID,], file = "ageNegInt.csv")


# add q values
library("qvalue")
qPos <- qvalue(sigPos, pi0.method="bootstrap")$qvalues
summary(qPos)
> summary(qPos)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 8.000e-10 1.456e-04 2.904e-04 2.267e-04 2.904e-04 1.472e-03 
# too small, can't be right!

library("RankProd")

