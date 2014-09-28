### presentation
# genes positively and negatively expressed in reproduction of jpm analysis
jpmPosEz <- unique(unlist(sapply(jpmPosEz.homo, function(x) x$egID)))
jpmNegEz <- unique(unlist(sapply(jpmNegEz.homo, function(x) x$egID)))
jpmEz <- select(org.Hs.eg.db, unique(c(jpmPosEz, jpmNegEz)), columns = c("ENTREZID", "SYMBOL"))

# make data frame to compare with original results
agePos <- merge(jpmEz, data.frame(SYMBOL = names(sigPos.min), metaP = sigPos.min, stringsAsFactors = FALSE))
ageNeg <- merge(jpmEz, data.frame(SYMBOL = names(sigNeg.min), metaP = sigNeg.min, stringsAsFactors = FALSE))

dim(agePos)
# [1] 2390    3
dim(ageNeg)
# [1] 2725    3
length(ageNeg$ENTREZID[ageNeg$metaP < .0001])
# [1] 263
length(agePos$ENTREZID[agePos$metaP < .0001])
# [1] 207
length(intersect(agePos$ENTREZID[agePos$metaP < .0001], ageNeg$ENTREZID[ageNeg$metaP < .0001]))
# [1] 5

write.csv(agePos, file = "agePos.csv")
write.csv(ageNeg, file = "ageNeg.csv")

# get original tables
over.expressed <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/aging signatures/over expressed.csv", stringsAsFactors = FALSE)
under.expressed <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/aging signatures/under expressed.csv", stringsAsFactors = FALSE)

dim(over.expressed)
# [1] 232   7
dim(under.expressed)
# [1] 146   5

write.csv(agePos[agePos$ENTREZID %in% over.expressed$EntrezGeneID,], file = "agePosInt.csv")
write.csv(ageNeg[ageNeg$ENTREZID %in% under.expressed$EntrezGeneID,], file = "ageNegInt.csv")

dim(agePos[agePos$ENTREZID %in% over.expressed$EntrezGeneID,])
# [1] 169   3
dim(ageNeg[ageNeg$ENTREZID %in% under.expressed$EntrezGeneID,])
# [1] 110   3

#import gene lists significant using fishers chi2 & compare with F test gene lists
chi_overexpressed <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/aging signatures/chi_overexpressed.csv")
chi_underexpressed <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/aging signatures/chi_underexpressed.csv")
length(intersect(over.expressed$EntrezGeneID, chi_overexpressed$EntrezGeneID))
# [1] 34
sapply(list(over.expressed$EntrezGeneID, chi_overexpressed$EntrezGeneID), length)
# [1] 232  71
length(intersect(under.expressed$EntrezGeneID, chi_underexpressed$EntrezGeneID))
# [1] 4
sapply(list(under.expressed$EntrezGeneID, chi_underexpressed$EntrezGeneID), length)
# [1] 146  41


# export moses chip gene lists from xxx_cum_hx5.best to csv files
mosesCombos <- function(mlist) {
  out <- cbind(mlist[[3]], mlist[[1]], mlist[[2]])
  names(out) <- c(names(out)[1:13], "combo", "out")
  write.csv(out[, c(13:15, 1:12)], file = paste(deparse(substitute(mlist)), ".csv", sep = ""), row.names = FALSE)
}

for(n in ls(pattern = "_cum_")) print(paste("mosesCombos(", n, ")", sep = ""))
mosesCombos(gpl1261_cum_hx5.best)
mosesCombos(gpl341_cum_hx5.best)
mosesCombos(gpl81_cum_hx5.best)
mosesCombos(gpl85_cum_hx5.best)
mosesCombos(gpl96_cum_hx5.bestN5)
mosesCombos(gpl97_cum_hx5.best)

# export moses gene lists from xxx_cum_hx5.best to csv files
mosesLists <- function(mlist, db) {
  require(db, character.only = TRUE)
  df <- mlist[[4]][,1:3]
  if(substring(df$feature[1], 1, 1) == "X") probes <- substring(df$feature, 2) else probes <- df$feature
  df$feature <- probes
  annot <- select(get(db), probes, columns = c("SYMBOL","GENENAME", "ENTREZID"))
  out <- merge(annot, df, by.x = "PROBEID", by.y = "feature", all = TRUE)
#   write.csv(unique(out), file = paste(deparse(substitute(mlist)), "_features.csv", sep = ""), row.names = FALSE)
  return(out)
}

chipMlist <- list()

for(n in ls(pattern = "_cum_")) print(paste("mosesLists(", n, ", ", ")", sep = ""))
chipMlist$gpl1261 <- mosesLists(gpl1261_cum_hx5.best, "mouse430a2.db")
chipMlist$gpl341 <- mosesLists(gpl341_cum_hx5.best, "rae230a.db")
chipMlist$gpl81 <- mosesLists(gpl81_cum_hx5.best, "mgu74av2.db")
chipMlist$gpl85 <- mosesLists(gpl85_cum_hx5.best, "rgu34a.db")
chipMlist$gpl96 <- mosesLists(gpl96_cum_hx5.bestN5, "hgu133a.db")
chipMlist$gpl97 <- mosesLists(gpl97_cum_hx5.best, "hgu133b.db")

chipMlist <- lapply(chipMlist, function(x) return(unique(x[, -5])))

# compare moses genes with genAge build 17 genes
# human
genage_human <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/genage_human.csv", stringsAsFactors = FALSE)
humeInt <- lapply(chipMlist[c("gpl96", "gpl97")], function(x) intersect(x$SYMBOL, genage_human$symbol_hugo))
humeInt
# $gpl96
# [1] "HSPA1A" "HSPA1B" "HSPD1"  "NCOR1"  "PTPN1"  "IGFBP2" "FOXO1" 
# $gpl97
# character(0)

# mouse
genage_models <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/genage_models.csv", stringsAsFactors = FALSE)
mouseInt <- lapply(chipMlist[c("gpl81", "gpl1261")], function(x) intersect(x$ENTREZID, subset(genage_models, organism == "Mus musculus")$Entrez.gene.id))
mouseIntS <- lapply(chipMlist[c("gpl81", "gpl1261")], function(x) intersect(x$SYMBOL, subset(genage_models, organism == "Mus musculus")$symbol))
unlist(c(mouseInt, mouseIntS))
#   gpl81   gpl81 
# "16491" "Kcna3" 
select(org.Mm.eg.db, as.character(mouseInt[[1]]), columns = c("ENTREZID", "SYMBOL"))
#   ENTREZID SYMBOL
# 1    16491  Kcna3


# export GeneMeta/XDE gene lists to csv
xdelists <- function(df, fname, alpha = .025) {
  df <- cbind(df, pvalue = pnorm(abs(df[, "zSco"]), lower.tail = FALSE), neg = as.numeric(abs(df[, "zSco"]) != df[, "zSco"]))
  out <- df[df[, "pvalue"] < alpha,]
  write.csv(out, file = paste(fname, "_zSco.csv", sep = ""), row.names = TRUE)
  return(out)
}

chipZsigs <- list()
for(n in names(chipBurnZ)) chipZsigs[[n]] <- xdelists(chipBurnZ[[n]], n)

# calculate overlap of z diff and moses genes
chipGeneInt <- mapply(function(x, y) y[x$SYMBOL[x$SYMBOL %in% rownames(y)],], chipMlist, chipZsigs)
chipGeneInt <- lapply(chipGeneInt, unique)
chipGeneInt
# $gpl1261
#       zSco_Ex_1  zSco_Ex_2  zSco_Ex_3 zSco_Ex_4      zSco     pvalue neg
# Hif1a -1.837097 0.04258764 -0.6544122 -1.890809 -2.218862 0.01324804   1
# 
# $gpl341
#             zSco_Ex_1  zSco_Ex_2 zSco_Ex_3   zSco_Ex_4      zSco       pvalue neg
# B2m      0.000000e+00  1.7713353  7.067281  1.91000237  6.660742 1.362246e-11   0
# Ctsb     0.000000e+00  1.0952746  3.269044 -0.05824359  3.209052 6.658671e-04   0
# Cd74     0.000000e+00  2.9157313 10.570492  2.10999577  2.350764 9.367444e-03   0
# Cd63     0.000000e+00  1.7305293  5.744309  1.64435436  6.096451 5.422451e-10   0
# Rgs3     0.000000e+00  0.6852579  4.693605  0.03632444  2.793456 2.607403e-03   0
# Anxa3    0.000000e+00  1.9592885  1.656759  2.37392963  1.990490 2.326849e-02   0
# C3       0.000000e+00  3.1890516 10.534281  3.22123168  3.246058 5.850743e-04   0
# Laptm5   0.000000e+00  2.1089724  3.629847 -0.02806565  3.323994 4.436902e-04   0
# Cd53     0.000000e+00  2.9757338  5.173866  2.98165831  3.099311 9.698570e-04   0
# C1qb     0.000000e+00  2.1178451  6.937819  1.47747631  5.895144 1.871771e-09   0
# Rps6ka1  0.000000e+00 -0.4253993  4.370815  1.81582080  2.101614 1.779354e-02   0
# RT1-Db1 -4.532467e-17  3.1651728  8.047888  1.74405645  3.608644 1.539006e-04   0
# Dnase2   0.000000e+00  1.3725720  5.153990  0.82112017  5.211933 9.344171e-08   0
# Slc13a1  0.000000e+00 -1.3614964 -1.889538 -0.49088975 -2.223646 1.308613e-02   1
# Myo16    0.000000e+00 -1.1578368 -2.892006  0.36999309 -2.801269 2.545105e-03   1
# Hdc      0.000000e+00 -1.0596182 -2.375408 -0.22603041 -2.494176 6.312489e-03   1
# RT1-Da   0.000000e+00  3.3117041 10.884236  2.04157352  2.240233 1.253791e-02   0
# Ctsz     0.000000e+00  2.2440285  3.672433  0.54465407  4.051127 2.548582e-05   0
# Ass1     0.000000e+00  3.6899755  5.197604  0.89871718  2.117864 1.709328e-02   0
# Fcgr2b   0.000000e+00  3.2752899  8.929528  2.34669152  3.828740 6.440060e-05   0
# Ndn      0.000000e+00 -1.0754592 -6.756041 -1.54147000 -3.697238 1.089792e-04   1
# Fcer1g   0.000000e+00  2.9995490  7.475708  3.29024833  4.002811 3.129715e-05   0
# Cpped1   0.000000e+00  1.1903014  3.217871  1.15614070  3.560662 1.849604e-04   0
# Mbnl1    9.064933e-17 -0.2531191  2.564260  0.06125173  2.220534 1.319129e-02   0
# Cdh11    3.625973e-16 -0.9194865 -4.318196 -0.12649404 -4.128401 1.826470e-05   1
# 
# $gpl81
#           zSco_Ex_1  zSco_Ex_2 zSco_Ex_3      zSco       pvalue neg
# Mmp13    1.12053317  0.7513435  1.628173  2.057831 0.0198031999   0
# Ralgds   0.03034259 -1.9099671 -2.446406 -2.813205 0.0024525157   1
# Ap3b1   -1.85385546 -2.2080149 -2.277745 -3.421607 0.0003112606   1
# Foxred1 -0.16236723 -0.9611893 -2.468730 -2.493491 0.0063246944   1
# Cpn1     0.89719241  1.4382231  1.392349  2.095100 0.0180810326   0
# S100a9   0.10182927  1.9105154  3.358068  3.485221 0.0002458651   0
# Rasl2-9  0.01717541 -2.0646969 -1.671677 -2.121238 0.0169508696   1
# Kcna3    0.39422187 -0.8075787 -2.565378 -2.108130 0.0175098624   1
# 
# $gpl85
#          zSco_Ex_1  zSco_Ex_2  zSco_Ex_3  zSco_Ex_4  zSco_Ex_5      zSco       pvalue neg
# Sort1   -0.2779817  1.6062551  2.0133788  2.5781743  1.4360215  3.508835 2.250367e-04   0
# Plekha1 -1.6900939 -2.0131665 -0.5185432 -1.6583339 -0.7475168 -2.901667 1.855912e-03   1
# Apc     -0.1346913  1.2511556  2.2604265  2.4526976  0.8721922  3.189828 7.117873e-04   0
# Col1a1  -3.2218848 -3.4422036 -1.7840012 -2.9639169 -2.3406854 -4.615408 1.961623e-06   1
# Fnta    -1.6852922 -0.8044149 -0.4760663 -0.9020978 -1.0996987 -2.085922 1.849283e-02   1
# Hbb      1.7559490 -0.1199852  1.1610377  3.7676566  0.6949579  2.463248 6.884238e-03   0
# Ndn     -1.6693150 -1.0017654  0.5597292 -3.6944514 -1.1664877 -2.214890 1.338379e-02   1
# Stk24   -1.7454395 -0.3296910 -0.2346506 -1.9758002 -0.3039118 -2.140751 1.614704e-02   1
# Hmgcs1  -1.5889545 -0.9003278  0.3777407 -1.8123779 -3.1310330 -2.188014 1.433429e-02   1
# 
# $gpl96
#         zSco_Ex_1 zSco_Ex_2      zSco       pvalue neg
# ZNF207 -0.9344716 -1.959091 -2.025889 0.0213880976   1
# RAP1B   2.2522484  2.349379  3.253498 0.0005699682   0
# MAP4    2.4442860  2.187677  3.272864 0.0005323188   0
# CTSB    2.3390760  1.676553  2.826485 0.0023530957   0
# PRDX6   2.8652818  1.183481  2.061881 0.0196095528   0
# TMEM66  2.1201435  2.013030  2.922242 0.0017376065   0
# NCOR1  -1.1642059 -1.873883 -2.138028 0.0162572264   1
# FOXO1   2.6248338  2.594738  3.690753 0.0001117955   0
# HPRT1   2.2836935  2.157452  3.139816 0.0008452709   0
# EPHA2  -1.6764308 -1.954922 -2.565773 0.0051473083   1
# BPGM   -2.3782130 -1.112898 -2.429205 0.0075659919   1
# MET    -1.2398044 -2.115292 -2.354686 0.0092691867   1
# 
# $gpl97
#          zSco_Ex_1 zSco_Ex_2      zSco       pvalue neg
# C11orf57 -1.494400 -1.980112 -2.451167 0.0071196952   1
# PLBD1    -1.134591 -2.576460 -2.247421 0.0123065721   1
# PNPO     -1.650361 -2.013953 -2.587640 0.0048318005   1
# ASH1L     1.968908  1.310187  2.308965 0.0104727682   0
# OIP5-AS1 -1.580311 -1.828504 -2.408952 0.0079992029   1
# SCAF11    2.938310  1.263307  2.121940 0.0169213911   0
# SPAG9     2.635876  2.355557  3.526174 0.0002108049   0
# CAPRIN1   2.129494  2.138534  3.017949 0.0012724576   0

mapply(function(x, y) unique(x[x$SYMBOL %in% rownames(y),]), chipMlist, chipGeneInt, SIMPLIFY = FALSE)

# write files for eddie
mouse_homo <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/mouse_NoHiLo.csv", stringsAsFactors=FALSE)
mouse_homo$Hsym <- homoMap$mouseSYM[match(mouse_homo$feature, homoMap$mouseSYM)]
write.csv(mouse_homo, file = "mouse_homologyAll.csv", row.names = FALSE)
write.csv(na.omit(mouse_homo), file = "mouse_homology.csv", row.names = FALSE)
rat_homo <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/rat_NoHiLo.csv", stringsAsFactors=FALSE)
rat_homo$Hsym <- homoMap$ratSYM[match(rat_homo$feature, homoMap$ratSYM)]
write.csv(rat_homo, file = "rat_homologyAll.csv", row.names = FALSE)
write.csv(na.omit(rat_homo), file = "rat_homology.csv", row.names = FALSE)

# function to strip first n characters from character vector strings
strip <- function(char, n = 1) substr(char,n + 1, max(sapply(char, nchar)))

gpl81_homo <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl81_NoHiLo.csv", stringsAsFactors=FALSE)
gpl85_homo <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl85_NoHiLo.csv", stringsAsFactors=FALSE)
gpl96_homo <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl96_NoHiLo.csv", stringsAsFactors=FALSE)
gpl97_homo <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl97_NoHiLo.csv", stringsAsFactors=FALSE)
gpl341b_homo <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl341b_NoHiLo.csv", stringsAsFactors=FALSE)
gpl1261_homo <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl1261_NoHiLo.csv", stringsAsFactors=FALSE)

chipList <- lapply(as.list(ls(pattern = "gpl")), get)
names(chipList) <- ls(pattern = "gpl")
for(i in c(1:3, 5:6)) chipList[[i]]$feature <- strip(chipList[[i]]$feature)
db <- c("mouse430a2.db", "rae230a.db", "mgu74av2.db", "rgu34a.db", "hgu133a.db", "hgu133b.db")
for(i in 1:6) {
  require(db[i], character.only = TRUE)
  annot <- select(get(db[i]), chipList[[i]]$feature, columns = "SYMBOL")
  chipList[[i]] <- merge(annot, chipList[[i]], by.x = "PROBEID", by.y = "feature", all = TRUE)
}
for(i in 1:4) {
  if(i %% 2 == 0) chipList[[i]]$Hsym <- ratMap$symbol[match(chipList[[i]]$SYMBOL, ratMap$SYMBOL)]
  if(i %% 2 != 0) chipList[[i]]$Hsym <- mouseMap$symbol[match(chipList[[i]]$SYMBOL, mouseMap$SYMBOL)]
} 
# noHiLo didn't work?  clean up lists
chipList <- lapply(chipList, function(x) x[order(x$Freq, decreasing = TRUE),])
chipList <- lapply(chipList, na.omit)
chipList <- lapply(chipList, unique)

# check for HiLos & remove homologs with equal counts, keep highest count if unequal
> lapply(chipList[1:4], function(x) x[duplicated(x$Hsym) | duplicated(x$Hsym, fromLast = TRUE),])
$gpl1261_homo
#       PROBEID  SYMBOL Freq level    Hsym
# 8  1416029_at   Klf10   13  down   KLF10
# 28 1416048_at    Phc2    6  down    PHC2
# 17 1416033_at Tmem109    5    up TMEM109
# 27 1416048_at    Phc2    2    up    PHC2
# 9  1416029_at   Klf10    1    up   KLF10
# 13 1416032_at Tmem109    1    up TMEM109
# 14 1416032_at Tmem109    1  down TMEM109
# 20 1416036_at  Fkbp1a    1    up  FKBP1A
# 21 1416036_at  Fkbp1a    1  down  FKBP1A
# remove KLF10 up, PHC2 up, TMEM109 down & increase Freq by 1, both FKBP1A
chipList$gpl1261_homo <- chipList$gpl1261_homo[!(row.names(chipList$gpl1261_homo) %in% c("27", "9", "13", "14", "20", "21")),] 
chipList$gpl1261_homo["17", 3] <- 6

# $gpl341b_homo
#         PROBEID       SYMBOL Freq level   Hsym
# 91 1370487_a_at        Kalrn   10    up  KALRN
# 43   1367850_at       Fcgr2a    5    up FCGR2A
# 44   1367850_at    LOC498276    5    up FCGR2A
# 46   1367850_at LOC100912061    5    up FCGR2A
# 47   1367850_at LOC100912098    5    up FCGR2A
# 37 1367830_a_at         Tp53    3    up   TP53
# 38   1367831_at         Tp53    3    up   TP53
# 27   1367806_at          Gls    2  down    GLS
# 2    1367562_at        Sparc    1    up  SPARC
# 3    1367563_at        Sparc    1    up  SPARC
# 26   1367805_at          Gls    1  down    GLS
# 36 1367830_a_at         Tp53    1  down   TP53
# 40   1367835_at       Pcsk1n    1  down PCSK1N
# 41   1367835_at LOC100911286    1  down PCSK1N
# 92 1370487_a_at        Kalrn    1  down  KALRN
# remove Kalrn down, rm dup FCGR2A, rm TP53 hilo probe, merge GLS, merge SPARC, rm dup PCSK1N 
chipList$gpl341b_homo <- chipList$gpl341b_homo[!(row.names(chipList$gpl341b_homo) %in% c("92", "44", "46", "47", "36", "37", "26", "3", "41")),] 
chipList$gpl341b_homo["27", 3] <- 3
chipList$gpl341b_homo["2", 3] <- 2

# $gpl81_homo
#        PROBEID     SYMBOL Freq level    Hsym
# 59   103878_at      Ap3b1   34  down   AP3B1
# 54 103874_r_at    Ankrd42   13  down ANKRD42
# 56 103874_r_at    Ccdc90b   13  down CCDC90B
# 62   103881_at       Ppa2   10    up    PPA2
# 28 100557_g_at      Eif4b    7    up   EIF4B
# 14   100536_at       Mobp    4  down    MOBP
# 63   103881_at       Ppa2    4  down    PPA2
# 8    100497_at       Stx3    2  down    STX3
# 19   100539_at      Acot7    2  down   ACOT7
# 29   100561_at     Iqgap1    2    up  IQGAP1
# 51 103873_i_at    Ankrd42    2    up ANKRD42
# 52 103873_i_at    Ccdc90b    2    up CCDC90B
# 7    100497_at       Stx3    1    up    STX3
# 15   100536_at       Mobp    1    up    MOBP
# 18   100539_at      Acot7    1    up   ACOT7
# 27   100556_at      Eif4b    1  down   EIF4B
# 30   100561_at     Iqgap1    1  down  IQGAP1
# 53 103874_r_at    Ankrd42    1    up ANKRD42
# 55 103874_r_at    Ccdc90b    1    up CCDC90B
# 57   103875_at       Ngrn    1    up    NGRN
# 58   103875_at       Ngrn    1  down    NGRN
# 60   103878_at      Ap3b1    1    up   AP3B1
# 73 161663_f_at    Dnajc19    1  down DNAJC19
# 74 161663_f_at Dnajc19-ps    1  down DNAJC19
# remove AP3B1 up, rm ANKRD42 upx2, CCDC90Bx2, PPA2 down, EIF4B down, MOBP up, STX3 up, ACOT7 up, IQGAP1 down, NGRN updown, merge DNAJC19
chipList$gpl81_homo <- chipList$gpl81_homo[!(row.names(chipList$gpl81_homo) %in% c("60", "51", "53", "52", "55", "63", "27", "15", "7", "18", "30", "57", "58", "74")),] 
chipList$gpl81_homo["73", 3] <- 2

# $gpl85_homo
#             PROBEID       SYMBOL Freq level   Hsym
# 9       AF032666_at        Exoc2    6    up  EXOC2
# 52   rc_AA875035_at       Slc4a1    4  down SLC4A1
# 67 rc_AA875099_s_at       Npap60    4    up  NUP50
# 10    AF032666_g_at        Exoc2    3  down  EXOC2
# 51   rc_AA875035_at       Slc4a1    3    up SLC4A1
# 58   rc_AA875054_at         Tcp1    2    up   TCP1
# 54   rc_AA875047_at        Cct6a    1  down  CCT6A
# 55   rc_AA875047_at        Cct6a    1    up  CCT6A
# 59   rc_AA875054_at         Tcp1    1  down   TCP1
# 68 rc_AA875099_s_at       Npap60    1  down  NUP50
# 87        U92289_at       Ptgdrl    1  down  PTGDR
# 88        U92289_at        Ptgdr    1  down  PTGDR
# 89        U93092_at        Hoxa1    1  down  HOXA1
# 90        U93092_at LOC100911406    1  down  HOXA1
# remove EXOC2 down, SLC4A1 up, TCP1 down, NUP50 down, CCT6A downup, dup HOXA1, dup PTGDR
chipList$gpl85_homo <- chipList$gpl85_homo[!(row.names(chipList$gpl85_homo) %in% c("10", "51", "59", "68", "54", "55", "87", 90)),] 


lapply(chipList, function(x) write.csv(na.omit(x), file = paste(names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "logy.csv", sep = ""), row.names = FALSE))

## useful function to remove all packages
detachAllPackages <- function() {
  basic.packages <- 
    c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

# function to strip first n characters from character vector strings
strip <- function(char, n = 1) substr(char,n + 1, max(sapply(char, nchar)))
