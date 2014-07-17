### jpm meta analysis, now using homologene homologue map and current probe annotations
# uses jpmPos, jpmNeg from 'jpm analysis.R' and jpmPosEz.homo, jpmNegEz.homo from 'homologene.R'

## original analysis:  duplicate gene log2 expression values averaged
# how many current GEO annotated significant probes are duplicate genes (possibly still more current than in original paper)
sapply(jpmSig, function(x)summary(duplicated(x$jpmGene)))
#      h_brain   muscle1   muscle2   muscle3   muscle4   muscle5   muscle    kidney1   kidney2  
# FALSE "2011"    "2397"    "1901"    "1951"    "1500"    "696"     "3002"    "420"     "431"    
# TRUE  "260"     "324"     "134"     "167"     "50"      "24"      "480"     "14"      "23"     
# 
#      m_brain   m_hippo   liver     m_heart   lung      cochlea   hemato_stem myo_progen r_hippo  
# FALSE "2206"    "1148"    "1408"    "482"     "2104"    "1185"    "3164"      "1336"     "2717"   
# TRUE  "113"     "35"      "68"      "10"      "184"     "26"      "165"       "93"       "149"    
# 
#     stromal   spinal_cord oculomotor skeletal_ms extraoc_ms laryngeal_ms r_heart   CA1_hipp2
# FALSE "972"     "1336"      "1223"     "430"       "350"      "956"        "595"     "757"    
# TRUE  "13"      "34"        "37"       "60"        "28"       "68"         "28"      "67"     

# how many significantly different probes map to more than one gene?
probe2many <- lapply(probe2gene, function(x) x$probe[x$probe %in% x$probe[duplicated(x$probe)]])
sapply(probe2many, function(x) length(unique(x)))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#          600         1148          692         1148          692          600            4 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#          272          266          936          557          557          645          936 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#          936          936          557          758          758          758          758 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#          561          561          561          561          561 

sapply(probe2many, function(x) length(unique(x)) / length(x))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#    0.3225806    0.3365582    0.4032634    0.3365582    0.4032634    0.3225806    0.4444444 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#    0.2921590    0.3159145    0.3280757    0.3175599    0.3175599    0.3429027    0.3280757 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#    0.3280757    0.3280757    0.3175599    0.3937662    0.3937662    0.3937662    0.3937662 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#    0.3337299    0.3337299    0.3337299    0.3337299    0.3337299 

# average fraction of individual data sets differentially expressed 2 sided p value < .05
# average of fraction of significant probes are close to jpm result except more are negatively correlated than positive...
summary(unlist(lapply(jpmExp.slope, function (x) (sum(x[, 3] <= .05 & x[,1] > 0) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01692 0.03093 0.04279 0.04398 0.05117 0.09408 

summary(unlist(lapply(jpmExp.slope, function (x) (sum(x[, 3] <= .05 & x[,1] < 0) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00997 0.02487 0.04636 0.04863 0.06555 0.10460 

# how many genes are duplicated using current annotation?  a lot!
slopeGenes <- mapply(function(x, y) x$GEOgene <- y$GEOgene[match(row.names(x), y$probe)], jpmExp.slope, probe2gene)
sapply(slopeGenes, function(x) summary(duplicated(x)))
#       h_brain   muscle1   muscle2   muscle3   muscle4   muscle5   muscle    kidney1   kidney2  
# FALSE "9450"    "14093"   "16846"   "14093"   "16846"   "9466"    "13949"   "4264"    "3821"   
# TRUE  "3110"    "8122"    "5731"    "8122"    "5731"    "3067"    "8741"    "543"     "222"    

#       m_brain   m_hippo   liver     m_heart   lung      cochlea   hemato_stem myo_progen
# FALSE "26764"   "9614"    "9614"    "8523"    "26764"   "26762"   "25914"     "7688"    
# TRUE  "18337"   "2874"    "2874"    "1520"    "18337"   "18275"   "16955"     "1787"    

#       r_hippo   stromal   spinal_cord oculomotor skeletal_ms extraoc_ms laryngeal_ms r_heart  
# FALSE "13444"   "13439"   "13444"     "13444"    "6081"      "6081"     "6081"       "6081"   
# TRUE  "2479"    "2474"    "2478"      "2478"     "2718"      "2718"     "2718"       "2718"   

#       CA1_hipp2
# FALSE "4852"   
# TRUE  "1758"   

# with duplicates removed -- percentage is halved!?!?!
jpmExp.slope.nd <- mapply(function(x, y) x[!duplicated(y),], jpmExp.slope, slopeGenes, SIMPLIFY = FALSE)

summary(unlist(lapply(jpmExp.slope.nd, function (x) (sum(x[, 3] <= .05 & x[,1] > 0) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01496 0.02024 0.02101 0.02185 0.02486 0.02859 

summary(unlist(lapply(jpmExp.slope.nd, function (x) (sum(x[, 3] <= .05 & x[,1] < 0) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01574 0.01841 0.02100 0.02095 0.02365 0.02795 
 
# total expression levels in all experiments
sum(unlist(lapply(jpmExp.slope, function(x) dim(x)[1])))
# [1] 474305

# total expression levels linearly related to age
sum(unlist(lapply(jpmSig, function(x) dim(x)[1])))
# [1] 39332

# make vector of unique gene symbols linearly associated with age from all data sets
Pos <- unique(unlist(lapply(jpmPosEz.homo, function(x) x$symbol)))
length(Pos)
# [1] 7750

Neg <- unique(unlist(lapply(jpmNegEz.homo, function(x) x$symbol)))
length(Neg)
# [1] 8395

length(intersect(Pos, Neg))
# [1] 4114  !!

# make list of vectors of gene named counts
geneCountPk <- vector("list", length(probe2gene))
names(geneCountPk) <- names(probe2gene)
for(n in names(geneCountPk)) {
  geneCountPk[[n]] <- table(jpmPosEz.homo[[n]]$symbol)
}

kpos <- numeric(length(Pos))
names(kpos) <- Pos
for(g in Pos) kpos[g] <- sum(sapply(geneCountPk, function(x) x[g]), na.rm = TRUE)

# get n, the number of experiments in which a gene is measured
# add homologene symbols to probes
homo2ezPos <- vector("list", 20)
names(homo2ezPos) <- names(probe2gene[7:26])
for(n in names(homo2ez)) {
  homo2ezPos[[n]] <- unique(merge(probe2gene[[n]][, c(1, 4)], jpmPosEz.homo[[n]][,c(3, 5)], by = "ENTREZID"))
}

geneCountPn <- vector("list", length(probe2gene))
names(geneCountPn) <- names(probe2gene)
for(n in names(geneCountPn[1:6])) geneCountPn[[n]] <- table(probe2gene[[n]]$SYMBOL)
for(n in names(geneCountPn[7:26])) geneCountPn[[n]] <- table(homo2ezPos[[n]]$symbol)

npos <- numeric(length(Pos))
names(npos) <- Pos
for(g in Pos) npos[g] <- sum(sapply(geneCountPn, function(x) x[g]), na.rm = TRUE)

# p value for null hypothesis binomial dist
# calculate cummulative binomial distribution
cbd <- function(k, n, p) {
  if(is.na(k) | is.na(n)) return(NA)
  pcb <- numeric(1)
  for(i in k:n){
    `<-`(pcb, pcb + choose(n, i)*p^i*(1-p)^(n-i))
  }
  return(pcb)
}

sigPos <- mapply(cbd, kpos, npos, 0.04398)
summary(sigPos)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.03504 0.16470 0.16720 0.27010 1.00000 
length(sigPos)
# [1] 7751

sigPos.min <- sigPos[sigPos < .05]
summary(sigPos.min)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001942 0.010940 0.017050 0.026020 0.048660 
length(sigPos.min)
# [1] 2391


# do the same for negatively associated genes
# get k, the number of experiments in which a gene is over/under expressed.  
# make list of vectors of gene named counts
geneCountNk <- vector("list", length(probe2gene))
names(geneCountNk) <- names(probe2gene)
for(n in names(geneCountNk)) {
  geneCountNk[[n]] <- table(jpmNegEz.homo[[n]]$symbol)
}

kneg <- numeric(length(Neg))
names(kneg) <- Neg
for(g in Neg) kneg[g] <- sum(sapply(geneCountNk, function(x) x[g]), na.rm = TRUE)

# get n, the number of experiments in which a gene is measured
# add homologene symbols to probes
homo2ezNeg <- vector("list", 20)
names(homo2ezNeg) <- names(probe2gene[7:26])
for(n in names(homo2ez)) {
  homo2ezNeg[[n]] <- unique(merge(probe2gene[[n]][, c(1, 4)], jpmNegEz.homo[[n]][,c(3, 5)], by = "ENTREZID"))
}

geneCountNn <- vector("list", length(probe2gene))
names(geneCountNn) <- names(probe2gene)
for(n in names(geneCountNn[1:6])) geneCountNn[[n]] <- table(probe2gene[[n]]$SYMBOL)
for(n in names(geneCountNn[7:26])) geneCountNn[[n]] <- table(homo2ezNeg[[n]]$symbol)

nneg <- numeric(length(Neg))
names(nneg) <- Neg
for(g in Neg) nneg[g] <- sum(sapply(geneCountNn, function(x) x[g]), na.rm = TRUE)

sigNeg <- mapply(cbd, kneg, nneg, 0.03883)
summary(sigNeg)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.0278  0.1264  0.1467  0.2421  1.0000 
length(sigNeg)
# [1] 8396

sigNeg.min <- sigNeg[sigNeg < .05]
summary(sigNeg.min)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001821 0.008585 0.014260 0.022430 0.049900 
length(sigNeg.min)
# [1] 2725

# add p values & slopes to homologues
jpmPos <- jpmPosEz.homo
jpmPos <- mapply(function (x, y) cbind(x, y[match(x$PROBEID, row.names(y)),]), jpmPos, jpmExp.slope, SIMPLIFY = FALSE)

jpmNeg <- jpmNegEz.homo
jpmNeg <- mapply(function (x, y) cbind(x, y[match(x$PROBEID, row.names(y)),]), jpmNeg, jpmExp.slope, SIMPLIFY = FALSE)

# # TODO: verify below, compare to original paper
# # add q values -- q values determined from probe p values from all complete data sets 
# pFslope <- unlist(sapply(jpmExp.slope, function(x) x$p.F))
# names(pFslope) <- unlist(sapply(jpmExp.slope, rownames))
# qFslope <- qvalue(pFslope, robust = TRUE)$qvalues #, pi0.method="bootstrap"
# summary(qFslope)
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# #  0.0029  0.6956  0.7899  0.7448  0.8435  0.8795       3 
# probeSig <- cbind(pFslope, qFslope)
# 
# for(n in names(jpmPos)) {
#   jpmPos[[n]]$q.F <- probeSig[match(jpmPos[[n]]$PROBEID, row.names(probeSig)), 2]
#   jpmNeg[[n]]$q.F <- probeSig[match(jpmNeg[[n]]$PROBEID, row.names(probeSig)), 2]
# }
# 
# sapply(jpmPos, function (x) dim(x[x$q.F < .1,])[1])
# #      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
# #           77           54           21           21            7            9            1 
# #      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
# #           15           32            7           33            6            0            3 
# #      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
# #            2           11            6          435           21           46          110 
# #  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
# #           18            2            5            0            1 
# 
# sapply(jpmNeg, function (x) dim(x[x$q.F < .1,])[1])
# #      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
# #          212           75           29           26            8            9            0 
# #      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
# #           17           19            2            7            0            0            1 
# #      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
# #            0            1            3          126           32           21           20 
# #  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
# #            1            0            0            1            1 

save(jpmPos, jpmNeg, file = "~/GitHub/stevia/data/jpmSig.rdata")
# TODO: complete below
# check against fisher's inverse chi-square: calculate sum of logs of p values for each gene
# function to calculate sum of logs of p values of a gene
logsum <- function (gene, df) {
  sum(log(df$p.F[df$symbol == gene]), na.rm = TRUE)
} 
invX2pos <- data.frame(gene = names(kpos), sigCount = kpos, slog = numeric(length(kpos)), stringsAsFactors = FALSE)
for(g in invX2pos$gene) {
  invX2pos$slog[invX2pos$gene == g] <- sum(sapply(jpmPos.uniq, function (x) logsum(g, x)), na.rm = TRUE)
}

invX2neg <- data.frame(gene = names(kneg), slog = numeric(length(kneg)), stringsAsFactors = FALSE)
jpmNeg.uniq <- lapply(jpmNeg, function(x) unique(x[!duplicated(x$PROBEID), c(1, 5, 8)]))
