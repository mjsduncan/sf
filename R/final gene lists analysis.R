# compile gene lists, just symbols, no expression levels
HomGeneList <- lapply(orgList, function(x) x[[1]])
HomGeneList <- c(HomGeneList, lapply(EchipList, function(x) x[[1]]))

# add moses genes to HomgeneList
HomGeneList$humanIndy <- singleList[[1]][[1]]
HomGeneList$mouseIndy <- na.omit(singleList[[2]][[3]])
HomGeneList$ratIndy <- na.omit(singleList[[3]][[3]])
HomGeneList$alldata <- na.omit(alldata[[2]][[1]])
HomGeneList <- lapply(HomGeneList, unique)

# what do we got?
# union of all
length(unique(unlist(HomGeneList)))
[1] 890
# by munge level:
length(unique(unlist(HomGeneList[1:3]))) # oranism
# [1] 470
length(unique(unlist(HomGeneList[4:9]))) # chip
# [1] 317
length(unique(unlist(HomGeneList[10:12]))) # wierd individual
# [1] 122
# by organism
length(unique(unlist(HomGeneList[c(1, 4, 5, 10)]))) # human
# [1] 307
length(unique(unlist(HomGeneList[c(2, 6, 7, 11)]))) # mouse
# [1] 257
length(unique(unlist(HomGeneList[c(3, 8, 9, 12)]))) # rat
# [1] 342
# by organism w/o wierd current individual results
length(unique(unlist(HomGeneList[c(1, 4, 5)]))) # human
# [1] 250
length(unique(unlist(HomGeneList[c(2, 6, 7)]))) # mouse
# [1] 224
length(unique(unlist(HomGeneList[c(3, 8, 9)]))) # rat
# [1] 310

# by intersection. function from http://codereview.stackexchange.com/questions/17905/compute-intersections-of-all-combinations-of-vectors-in-a-list-of-vectors-in-r
overlap <- function(l) {
  results <- lapply(l, unique) # combinations of m elements of list l
  for (m in seq(along=l)[-1]) {
    # generate and iterate through combinations of length m
    for (indices in combn(seq(length(l)), m, simplify=FALSE)) {
      # make name by concatenating the names of the elements
      # of l that we're intersecting
      name_1 <- paste(names(l)[indices[-m]], collapse="_")
      name_2 <- names(l)[indices[m]]
      name <- paste(name_1, name_2, sep="_")
      results[[name]] <- intersect(results[[name_1]], results[[name_2]])
    }
  }
  return(results)
}
interList <- overlap(HomGeneList)
for(n in names(interList)) {
  if(interList[[n]] == "") interList[[n]] <- NULL
}

intersect()