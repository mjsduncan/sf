# I have attached moses output run on each datasets.. and I did pull out perfect genes before running moses. This is how I did it

eachExprs <- lapply(symEsets,  exprs)

eachMoses <- list()
  
for(n in names(eachExprs)[-(7:9)]) {
  eachMoses[[n]] <- med.normalize(eachExprs[[n]])
  eachMoses[[n]] <- t(rbind(symEsets[[n]]$age , eachMoses[[n]]))
  colnames(eachMoses[[n]])[1] <- "age"
}
## removing perfect score genes
for(n in names(eachMoses)){
    perScorGens <- vector("list")
      j = 1
        for(i in 2:dim(eachMoses[[n]])[2]){
          if(length(unique(eachMoses[[n]][ , 1] == eachMoses[[n]][ ,i])) == 1){
            perScorGens[[j]] <- colnames(eachMoses[[n]])[i]
            j <- j+1
          }
         }
 if(length(perScorGens) > 0) {eachMoses[[n]] <- eachMoses[[n]][ , !colnames(eachMoses[[n]]) %in% perScorGens]}
}

sapply(eachExprs, dim)
#      h_brain muscle1 muscle2 muscle3 muscle4 muscle5 muscle kidney1 kidney2 m_brain m_hippo liver
# [1,]    8632   12494    9648   12494    9648    8616     98    5106    3299   13016    8724  8724
# [2,]      30      15      15      15      15      12     10      10      10       6      23     7
#      m_heart  lung cochlea hemato_stem myo_progen r_hippo stromal spinal_cord oculomotor
# [1,]    8622 13016   13016       13016       8724   10000   10000       10000      10000
# [2,]      12    15       6           8          4      20       3           9          9
#      skeletal_ms extraoc_ms laryngeal_ms r_heart CA1_hipp2
# [1,]        4750       4750         4750    4750      4750
# [2,]          12         12            9      11        29
# > sapply(eachMoses, dim)
#      h_brain muscle1 muscle2 muscle3 muscle4 muscle5 m_brain m_hippo liver m_heart  lung cochlea
# [1,]      30      15      15      15      15      12       6      23     7      12    15       6
# [2,]    8633   12495    9648   12495    9649    8617   12992    8725  8687    8623 13008   12982
#      hemato_stem myo_progen r_hippo stromal spinal_cord oculomotor skeletal_ms extraoc_ms
# [1,]           8          4      20       3           9          9          12         12
# [2,]       12975       8505   10000    9750        9990       9997        4750       4751
#      laryngeal_ms r_heart CA1_hipp2
# [1,]            9      11        29
# [2,]         4747    4751      4751

eachMoses <- list()
  
for(n in names(eachExprs)[-(7:9)]) {
  eachMoses[[n]] <- med.normalize(eachExprs[[n]])
  eachMoses[[n]] <- t(rbind(symEsets[[n]]$age , eachMoses[[n]]))
  colnames(eachMoses[[n]])[1] <- "age"
}
## removing perfect score genes
for(n in names(eachMoses)){
    perScorGens <- vector("list")
      j = 1
        for(i in 2:dim(eachMoses[[n]])[2]){
          if(length(unique(eachMoses[[n]][ , 1] != eachMoses[[n]][ ,i])) == 1){
            perScorGens[[j]] <- colnames(eachMoses[[n]])[i]
            j <- j+1
          }
         }
 if(length(perScorGens) > 0) {eachMoses[[n]] <- eachMoses[[n]][ , !colnames(eachMoses[[n]]) %in% perScorGens]}
}

sapply(eachMoses, dim)
#      h_brain muscle1 muscle2 muscle3 muscle4 muscle5 m_brain m_hippo liver m_heart  lung cochlea
# [1,]      30      15      15      15      15      12       6      23     7      12    15       6
# [2,]    8633   12495    9648   12495    9649    8617   12992    8725  8687    8623 13008   12982
#      hemato_stem myo_progen r_hippo stromal spinal_cord oculomotor skeletal_ms extraoc_ms
# [1,]           8          4      20       3           9          9          12         12
# [2,]       12975       8505   10000    9750        9990       9997        4750       4751
#      laryngeal_ms r_heart CA1_hipp2
# [1,]            9      11        29
# [2,]         4747    4751      4751

# count perfect features (add one to account for addition of "age" in eachMoses)
mapply(function(x, y) dim(x)[1] - dim(y)[2], eachExprs[-(7:9)], eachMoses)
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5      m_brain      m_hippo        liver 
#           -1           -1            0           -1           -1           -1           24           -1           37 
#      m_heart         lung      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#           -1            8           34           41          219            0          250           10            3 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#            0           -1            3           -1           -1 
# redo complete eachMoses
for(n in names(eachExprs)[-(7:9)]) {
  eachMoses[[n]] <- med.normalize(eachExprs[[n]])
  eachMoses[[n]] <- t(rbind(symEsets[[n]]$age , eachMoses[[n]]))
  colnames(eachMoses[[n]])[1] <- "age"
}

# function to pick out perfect features
getPerfect <- function(MosesMat, case = "age") {
  up <- colnames(MosesMat[, apply(MosesMat, 2, function(x) identical(x, MosesMat[, case]))])[-1]
  down <- colnames(MosesMat[, apply(MosesMat, 2, function(x) all(x == as.integer(!MosesMat[, case])))])
  return(list(up = up, down = down))
}

eachPerf <- lapply(eachMoses, getPerfect)

# make dataframes of perfect features for eddie
perFeat <- eachPerf
for(n in names(perFeat)) {
  if(length(c(perFeat[[n]]$up, perFeat[[n]]$down)) > 0) {
    perFeat[[n]] <- data.frame(perfectFeature = c(perFeat[[n]]$up, perFeat[[n]]$down),
                               score = dim(eachMoses[[n]])[1],
                               level = c(rep("up", length(perFeat[[n]]$up)), rep("down", length(perFeat[[n]]$down))),
                               stringsAsFactors = FALSE)
    if(n %in% names(probe2homo)) perFeat[[n]] <- merge(perFeat[[n]], probe2homo[[n]][match(perFeat[[n]][[1]], probe2homo[[n]]$SYMBOL), c(3, 4, 6, 5)], by.x = "perfectFeature", by.y = "SYMBOL", all.x = TRUE)

  }
}

EperFeat <- lapply(perFeat, function(x) if(length(x) == 6) x[, c(5, 2, 3)] else x)
for(n in names(EperFeat)) if(names(EperFeat[[n]])[1] == "symbol") names(EperFeat[[n]])[1] <- "perfectFeature"

# only 3 feature-level pairs are duplicated among the 23 data sets.
na.omit(Reduce(rbind, EperFeat[sapply(perFeat, length) >= 3]))[duplicated(na.omit(Reduce(rbind, EperFeat[sapply(perFeat, length) >= 3]))[, c(1, 3)]),]
#     perfectFeature score level
# 385          APAF1     3    up
# 451           FMOD     3    up
# 554       C1orf174     3  down

# combine the individual data set gene lists
EperFeat <- na.omit(Reduce(rbind, EperFeat[sapply(perFeat, length) >= 3]))

# aggregate duplicates
EperFeat <- aggregate(. ~ perfectFeature + level, EperFeat, sum)

write.csv(EperFeat[order(EperFeat$perfectFeature),c(1, 3, 2)], file = "indy_updown.csv", row.names = FALSE)
