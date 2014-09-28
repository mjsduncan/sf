### organism report
load("~/GitHub/stevia/data/orgCumX4.rdata")
human <- combo2fcount(homo_cum_hx5.best$combo)
mouse <- combo2fcount(mouse_cum_hx5.best$combo)
rat <- combo2fcount(rat_cum_hx5.best$combo)
rat$Hsym <- ratMap$symbol[match(rat$feature, ratMap$SYMBOL)]
mouse$Hsym <- mouseMap$symbol[match(mouse$feature, mouseMap$SYMBOL)]
orgList <- list(human = human, mouse = mouse, rat = rat)
lapply(orgList, function(x) write.csv(x, file = paste(names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], ".csv", sep = ""), row.names = FALSE))

# make files for eddie with Escore
lapply(orgList, function(x) write.csv(Escore(x), file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_add.csv", sep = ""), row.names = FALSE))
lapply(orgList, function(x) write.csv(Escore(x, add = FALSE), file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_sub.csv", sep = ""), row.names = FALSE))
# convert to new preffered format
orgList[2:3] <- lapply(orgList[2:3], function(x) na.omit(x[, c(4, 2, 3)]))
lapply(orgList, function(x) write.csv(x, file = paste("../eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_updown.csv", sep = ""), row.names = FALSE))

# re-make chipList after deleting "_HiLo.csv" files
gpl81 <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl81.csv", stringsAsFactors=FALSE)
gpl85 <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl85.csv", stringsAsFactors=FALSE)
gpl96 <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl96.csv", stringsAsFactors=FALSE)
gpl97 <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl97.csv", stringsAsFactors=FALSE)
gpl341b <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl341b.csv", stringsAsFactors=FALSE)
gpl1261 <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/gpl1261.csv", stringsAsFactors=FALSE)

chipList <- list(gpl96, gpl97, gpl81, gpl1261, gpl85, gpl341b)
names(chipList) <- c("gpl96", "gpl97", "gpl81", "gpl1261", "gpl85", "gpl341b")
for(i in c(1:4, 6)) chipList[[i]]$feature <- strip(chipList[[i]]$feature)
db <- c("hgu133a.db", "hgu133b.db", "mgu74av2.db", "mouse430a2.db", "rgu34a.db", "rae230a.db")
for(i in 1:6) {
  require(db[i], character.only = TRUE)
  annot <- select(get(db[i]), chipList[[i]]$feature, columns = "SYMBOL")
  chipList[[i]] <- merge(annot, chipList[[i]], by.x = "PROBEID", by.y = "feature", all = TRUE)
}
for(i in 3:6) {
  if(i > 4) chipList[[i]]$Hsym <- ratMap$symbol[match(chipList[[i]]$SYMBOL, ratMap$SYMBOL)]
  if(i < 5) chipList[[i]]$Hsym <- mouseMap$symbol[match(chipList[[i]]$SYMBOL, mouseMap$SYMBOL)]
} 
# clean up lists
chipList <- lapply(chipList, function(x) x[order(x$SYMBOL, x$Freq),])
EchipList <- lapply(chipList, function(x) ifelse(length(x) == 5, return(na.omit(x[, c(5, 3, 4)])),return(na.omit(x[, 2:4]))))

# eddie files  TODO:  reformat chiplist for Escore....
lapply(chipList, function(x) write.csv(Escore(x[, 2:length(x)]), file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_add.csv", sep = ""), row.names = FALSE))
lapply(chipList, function(x) write.csv(Escore(x, add = FALSE), file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_sub.csv", sep = ""), row.names = FALSE))

# compile gene lists
HomGeneList <- lapply(orgList, function(x) x[[1]])
HomGeneList <- c(HomGeneList, lapply(EchipList, function(x) x[[1]]))

# get icog individual data set lists
Homo_sapiensgenesfromEach <- read.csv("C:/Users/user/Desktop/biomind/addis/icog/eachData/moses/Homo_sapiensgenesfromEach.csv", stringsAsFactors=FALSE)
Mus_musculusgenesfromEach <- read.csv("C:/Users/user/Desktop/biomind/addis/icog/eachData/moses/Mus_musculusgenesfromEach.csv", stringsAsFactors=FALSE)
Rattus_norvegicusgenesfromEach <- read.csv("C:/Users/user/Desktop/biomind/addis/icog/eachData/moses/Rattus_norvegicusgenesfromEach.csv", stringsAsFactors=FALSE)

# apply Escore to old, incorrect combo counts
singleList <- list(Homo_sapiensgenesfromEach, Mus_musculusgenesfromEach, Rattus_norvegicusgenesfromEach)
names(singleList) <- c("humanIndy", "mouseIndy", "ratIndy")
singleList <- lapply(singleList, function (x) x[-4][, c(1, 3, 2)])
for(i in 1:3) names(singleList[[i]]) <- c("feature", "Freq", "level")
singleList <- lapply(singleList, Escore)
singleList$mouseIndy$Hsym <- mouseMap$symbol[match(singleList$mouseIndy$feature, mouseMap$SYMBOL)]
singleList$ratIndy$Hsym <- ratMap$symbol[match(singleList$ratIndy$feature, ratMap$SYMBOL)]

# save add version of Escores
lapply(singleList, function(x) write.csv(x, file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_add.csv", sep = ""), row.names = FALSE))

# add to HomgeneList
HomGeneList$humanIndy <- singleList[[1]][[1]]
HomGeneList$mouseIndy <- na.omit(singleList[[2]][[3]])
HomGeneList$ratIndy <- na.omit(singleList[[3]][[3]])
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

# get complete data set gene list from 2 available moses log files
alldata <- getMout(lines = 11, drop = 1)
alldata <- Mout2str(alldata)

intersect()