### PSEA test
library("PGSEA")
msigSym <- readGmt("C:/Users/user/Desktop/biomind/artificial biologist/pgsea/msigdb.v4.0.symbols.gmt")
tftSym <- readGmt("C:/Users/user/Desktop/biomind/artificial biologist/pgsea/c3.tft.v4.0.symbols.gmt")
names(tftSym) <- sub("^.*\\$(.*)_.+$", "\\1", names(tftSym))
names(tftSym) <- gsub("_.*$", "", names(tftSym))

### functions
# count samples less than p value
PGSEAp <- function(pout, p.cut = .05) {
  rowSums(pout$p.results < p.cut)
}

# apply to human mash
Phum <- PGSEA(orgExprs[[1]], tftSym, ref = which(orgAges[[1]] == 0), p.value = TRUE)
Phum05 <- PGSEA(orgExprs[[1]], tftSym, ref = which(orgAges[[1]] == 0), p.value = .05)

dim(Phum[[1]])
# [1] 615 102
sum(orgAges[[1]])
# [1] 60

dim(Phum$results)

summary(PGSEAp(Phum))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   0.000   2.000   7.752  10.000  57.000     107 

# apply to human data sets
PhumI <- lapply(symEsets[1:6], function(x) PGSEA(x, tftSym, ref = which(x$age == 0), p.value = TRUE))
sapply(PhumI, function(x) dim(x$results))
#      h_brain muscle1 muscle2 muscle3 muscle4 muscle5
# [1,]     615     615     615     615     615     615
# [2,]      30      15      15      15      15      12
sapply(symEsets[1:6], function(x) sum(x$age))
# h_brain muscle1 muscle2 muscle3 muscle4 muscle5 
#      22       8       8       8       8       6 

sapply(PhumI, function(x) summary(PGSEAp(x)))
#         h_brain muscle1 muscle2 muscle3 muscle4 muscle5
# Min.       0.00   0.000  0.0000   0.000  0.0000  0.0000
# 1st Qu.    1.00   0.000  0.0000   0.000  0.0000  0.0000
# Median     2.00   1.000  1.0000   1.000  1.0000  1.0000
# Mean       2.79   1.191  0.9793   1.025  0.7259  0.8038
# 3rd Qu.    4.00   2.000  2.0000   2.000  1.0000  1.0000
# Max.      19.00   8.000  6.0000   7.000  4.0000  6.0000
# NA's      34.00  45.000 35.0000  45.000 35.0000 34.0000
