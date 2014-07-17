### presentation
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

