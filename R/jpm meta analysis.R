### jpm meta analysis, now using homologene homologue map and current probe annotations
# uses jpmPos, jpmNeg from 'jpm analysis.R' and jpmPosEz.homo, jpmNegEz.homo from 'homologene.R'

## original analysis:  duplicate gene log2 expression values averaged
# how many current GEO annotated significant probes are duplicate genes (possibly still more current than in original paper)
sapply(jpmSig, function(x)summary(duplicated(x$jpmGene)))
#       h_brain   muscle1   muscle2   muscle3   muscle4   muscle5   muscle    kidney1   kidney2  
# FALSE "1624"    "1971"    "1566"    "1553"    "1312"    "693"     "2325"    "372"     "366"    
# TRUE  "187"     "194"     "93"      "100"     "44"      "16"      "287"     "10"      "21"     

#       m_brain   m_hippo   liver     m_heart   lung      cochlea   hemato_stem myo_progen r_hippo  
# FALSE "2291"    "937"     "1065"    "461"     "2241"    "1776"    "2630"      "969"      "2289"   
# TRUE  "74"      "29"      "43"      "5"       "142"     "67"      "113"       "52"       "120"    

#       stromal   spinal_cord oculomotor skeletal_ms extraoc_ms laryngeal_ms r_heart   CA1_hipp2
# FALSE "835"     "1086"      "1000"     "433"       "396"      "703"        "483"     "605"    
# TRUE  "7"       "18"        "22"       "45"        "23"       "42"         "13"      "41"     

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
summary(unlist(lapply(jpmExp.slope, function (x) (sum((x[, 3] <= .025 | x[, 3] >= .975) & x[,1] > 0, na.rm = TRUE) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02203 0.02742 0.03529 0.03773 0.04516 0.08202 

summary(unlist(lapply(jpmExp.slope, function (x) (sum((x[, 3] <= .025 | x[, 3] >= .975) & x[,1] < 0, na.rm = TRUE) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01890 0.02512 0.03843 0.03940 0.04906 0.08400 

# with duplicates removed TODO: map rownames/probes to gene symbols to get actual values
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

jpmExp.slope.nd <- mapply(function(x, y) x[!duplicated(y),], jpmExp.slope, slopeGenes, SIMPLIFY = FALSE)

summary(unlist(lapply(jpmExp.slope.nd, function (x) (sum((x[, 3] <= .025 | x[, 3] >= .975) & x[,1] > 0, na.rm = TRUE) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02145 0.02800 0.03507 0.03771 0.04178 0.08375 

summary(unlist(lapply(jpmExp.slope.nd, function (x) (sum((x[, 3] <= .025 | x[, 3] >= .975) & x[,1] < 0, na.rm = TRUE) / dim(x)[1]))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01612 0.02508 0.03686 0.03883 0.04887 0.07672 
 
# total expression levels in all experiments
sum(unlist(lapply(jpmExp.slope, function(x) dim(x)[1])))
# [1] 474305

# total expression levels linearly related to age
sum(unlist(lapply(jpmSig, function(x) dim(x)[1])))
# [1] 33790

# make vector of unique gene symbols linearly associated with age from all data sets
Pos <- unique(unlist(lapply(jpmPosEz.homo, function(x) x$symbol)))
length(Pos)
# [1] 7292

Neg <- unique(unlist(lapply(jpmNegEz.homo, function(x) x$symbol)))
length(Neg)
# [1] 7631

length(intersect(Pos, Neg))
# [1] 3651  !!

# get k, the number of experiments in which a gene is over/under expressed.  
kpos <- numeric(length(Pos))
names(kpos) <- Pos
for(g in Pos) kpos[g] <- sum(sapply(jpmPosEz.homo, function(x) table(x$symbol)[g]), na.rm = TRUE)

# get n, the number of experiments in which a gene is measured
# add homologene symbols to probes
homo2ez <- vector("list", 20)
names(homo2ez) <- names(probe2gene[7:26])
for(n in names(homo2ez)) {
  homo2ez[[n]] <- unique(merge(probe2gene[[n]][, c(1, 4)], jpmPosEz.homo[[n]][,c(3, 5)], by = "ENTREZID"))
}

npos <- numeric(length(Pos))
names(npos) <- Pos
npos.hs <- numeric(length(Pos))
names(npos.hs) <- Pos
npos.nhs <- numeric(length(Pos))
names(npos.nhs) <- Pos
for(g in Pos) {
  npos.hs[g] <- sum(sapply(probe2gene[1:6], function(x) table(x$SYMBOL)[g]), na.rm = TRUE)
  npos.nhs[g] <- sum(sapply(homo2ez, function(x) table(x$symbol)[g]), na.rm = TRUE)
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
kpos <- numeric(length(Pos))
names(kpos) <- Pos
for(g in Pos) kpos[g] <- sum(sapply(jpmPosEz.homo, function(x) g %in% x$symbol), na.rm = TRUE)

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
