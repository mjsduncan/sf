### annotate bed-stuy moses runs
#M# team icog version
ex1 <- gpl341_cum_hx5.best$features[1]
ex2<- ex1[[1]]

#M# data frames are both a 2 dimensional array-like object that can be indexed with single brackets df[, y] and a list of vectors that can be indexed by double brackets  df[[y]] so you just need to do this:
ex2 <- gpl341_cum_hx5.best$features[[1]]


# extract Probe ID from result

#M# if you use a function from a package load the package in your code!
library("stringr")

#M# it is good to include your testing/exploratory code but include the results so it is clear what you you are doing
pattern = unlist(str_split(ex2, "X", 2))[2]
[1] "1370482_at"

ex3 <- ex2
ex4 <- ex3
# removing "X" from the probeids
for(i in 1:99) ex3[i]  <- unlist(str_split(ex2[i], "X", 2))[2]

#M# you can also use "substr(ex2, 2, 99)"

#map probeID to gene name from rae230a.db

#M# don't forget to load package to make code reproducible!
library(rae230a.db)

for(i in 1:99)
    { probemap <- select(rae230a.db, ex3[i] , columns = c("ENTREZID", "SYMBOL", "GENENAME" ))
                   ex4[i] <- probemap$GENENAME }

#M# select works on vectors so you can do this:
probemap <- select(rae230a.db, ex3, columns = c("ENTREZID", "SYMBOL", "GENENAME"), keytype = "PROBEID")
#M# this way you save all the information you got out of "rae230a.db"
head(probemap)
#        PROBEID ENTREZID  SYMBOL                                             GENENAME
# 1   1370482_at    24767  Scnn1b    sodium channel, non-voltage-gated 1, beta subunit
# 2 1370485_a_at    24888  Bcl2l1                                          Bcl2-like 1
# 3   1373924_at   302890  Cpped1 calcineurin-like phosphoesterase domain containing 1
# 4   1368000_at    24232      C3                               complement component 3
# 5 1370383_s_at   294270 RT1-Db1                              RT1 class II, locus Db1
# 6   1373945_at   314647   Celf5                     CUGBP, Elav-like family member 5

ex4.alt <- probemap[[4]]

#M# your way missed some of the map
length(ex4.alt)
# [1] 73
length(unique(ex4.alt))
# [1] 68
length(ex4)
# [1] 99
length(unique(ex4))
# [1] 60

# To keep the original data we copied it to var test
test <- gpl341_cum_hx5.best
test$features[1] <- ex4

### here is how i did it
# uses "gplXXX_cum_hx5.best" from "mash runs.R" in Rmoses and "probe2gene" from "data cleaning.R"
load("~/GitHub/stevia/data/gplCumX6.rdata")

# get list of feature dataframes from moses "best of" lists
jpmMftures <- lapply(ls(pattern = "_cum_"), function(x) get(x)[[4]])

# strip "X" from probe names
jpmMftures <- lapply(jpmMftures, function(x) x <- cbind(substr(x$feature, 2, 99), x[, -1], stringsAsFactors = FALSE))

names(jpmMftures) <- substr(ls(pattern = "_cum_"), 1, regexpr("_", ls(pattern = "_cum_")) - 1)

# fix gpl85: first probe character is original!
jpmMftures$gpl85 <- gpl85_cum_hx5.best[[4]]

# clean up results: fix probe column name and name list elements by microarray
for(i in seq_along(jpmMftures)) {
  names(jpmMftures[[i]])[1] <- "probe"
}

# add symbol column from "probe2gene" list from 
jpmMftures <- mapply(function(x, y) x <- merge(x, y[, c(1, 3, 4)], by = "probe", all.x = TRUE), jpmMftures, probe2gene[c(10, 19, 17, 24, 4, 3)], SIMPLIFY = FALSE)

#reduce array lists to unique symbols & frequencies, by high and low expression
library("plyr")

jpmMhigh <- lapply(jpmMftures, function(x) arrange(count(subset(x, low == FALSE, select = c(probe, SYMBOL, ENTREZID))), desc(freq)))
jpmMlow <- lapply(jpmMftures, function(x) arrange(count(subset(x, low == TRUE, select = c(probe, SYMBOL, ENTREZID))), desc(freq)))

# add homologue annotation using modified rown2nrows & dflist2nrows functions from homologene.R
jpmMhighHomo <- list()
jpmMlowHomo <- list()
for(n in names(jpmMhigh)[1:4]) {
  jpmMhighHomo[[n]] <- cbind(jpmMhigh[[n]], t(sapply(jpmMhigh[[n]]$ENTREZID, ez2homo)))
  jpmMlowHomo[[n]] <- cbind(jpmMlow[[n]], t(sapply(jpmMlow[[n]]$ENTREZID, ez2homo)))
}

# fix NAx10000 problem by adding "n > 100" clause to row2nrows in homologene.R
row2nrows2 <- function(row) {
  out <- row
  n = length(out$egID[[1]])
  if(n == 0 | n > 100) {  
    out$egID <- ""
    out$symbol <- ""
    return(out)
  }
  if(n == 1) {  
    out$egID <- as.character(unlist(out$egID))
    out$symbol <- unlist(out$symbol)
    return(out)
  }
  out3 <- data.frame(probe = character(0), SYMBOL = character(0), ENTREZID = character(0), freq = numeric(0), egID = character(0), symbol = character(0))
  for(i in 1:n) {
    out2 <- out[, 1:4]
    out2$egID <- as.character(out$egID[[1]][i])
    out2$symbol <- out$symbol[[1]][i]
    out3 <- rbind(out2, out3)
  }
  return(out3)
}

# apply row2nrows2 to whole dataframe
dflist2nrows2 <- function(df) {
  out <- data.frame(probe = character(0), SYMBOL = character(0), ENTREZID = character(0), freq = numeric(0), egID = character(0), symbol = character(0))
  for(i in seq_len(dim(df)[1])) {
    out <- rbind(out, row2nrows2(df[i,]))
  }
  return(out)
}
# add rows for rodent genes with multiple human gene homologues
jpmMhighHomo <- lapply(jpmMhighHomo, dflist2nrows2)
jpmMlowHomo <- lapply(jpmMlowHomo, dflist2nrows2)
