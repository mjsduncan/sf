### final report

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
EorgList[2:3] <- lapply(orgList[2:3], function(x) na.omit(x[, c(4, 2, 3)]))
lapply(EorgList, function(x) write.csv(x, file = paste("../eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_updown.csv", sep = ""), row.names = FALSE))

# make csv files for org workbooks
annList <- list("rae230a.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db")
for(i in 1:3) require(annList[[i]], character = TRUE)

# export moses org gene lists from xxx_cum_hx5.best to csv files
mosesCombos <- function(mlist, name = deparse(substitute(mlist))) {
  out <- cbind(mlist[[3]], mlist[[1]], mlist[[2]])
  names(out) <- c(names(out)[1:13], "combo", "out")
  write.csv(out[, c(13:15, 1:12)], file = paste(name, ".csv", sep = ""), row.names = FALSE)
}

for(n in ls(pattern = "_cum_")) print(paste("mosesCombos(", n, ")", sep = ""))
mosesCombos(gpl341b_cum_hx5.best)
mosesCombos(homo_cum_hx5.best)
mosesCombos(mouse_cum_hx5.best)
mosesCombos(rat_cum_hx5.best)

# annotate and export moses gene lists from feature data frames to csv files
mosesLists <- function(olist, db) {
  require(db, character.only = TRUE)
  df <- olist
  if(substring(df$feature[1], 1, 1) == "X") {
    probes <- substring(df$feature, 2)
    df$feature <- probes
    annot <- select(get(db), probes, columns = c("SYMBOL","GENENAME", "ENTREZID"))
    out <- merge(annot, df, by.x = "PROBEID", by.y = "feature", all = TRUE)
  } else {
    annot <- select(get(db), df$feature, columns = c("SYMBOL","GENENAME", "ENTREZID"), keytype = "SYMBOL")
    out <- unique(merge(df, annot, by.x = "feature", by.y = "SYMBOL", all = TRUE))
  }
  write.csv(unique(out), file = paste(deparse(substitute(olist)), "_features.csv", sep = ""), row.names = FALSE)
  return(out)
}


gene.csv <- mapply(mosesLists, list(gpl341b = combo2fcount(gpl341b_cum_hx5.best$combo), human = human, mouse = mouse, rat = rat), annList)
gene.csv$gpl341b$Hsym <- ratMap$symbol[match(gene.csv$gpl341b$SYMBOL, ratMap$SYMBOL)]
lapply(gene.csv, function(x) write.csv(x, file = paste(names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_features.csv", sep = ""), row.names = FALSE))

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

# eddie files
# TODO:  reformat chiplist for Escore....
for(i in 1:6) names(EchipList[[i]])[1] <- "feature"
lapply(EchipList, function(x) write.csv(aggregate(. ~ feature + level, data = x, sum), file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_updown.csv", sep = ""), row.names = FALSE))

# lapply(chipList, function(x) write.csv(Escore(x[, 2:length(x)]), file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_add.csv", sep = ""), row.names = FALSE))
# lapply(chipList, function(x) write.csv(Escore(x, add = FALSE), file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_sub.csv", sep = ""), row.names = FALSE))
# robert not interested in above scoring

# get icog individual data set lists
Homo_sapiensgenesfromEach <- read.csv("C:/Users/user/Desktop/biomind/addis/icog/eachData/moses/Homo_sapiensgenesfromEach.csv", stringsAsFactors=FALSE)
Mus_musculusgenesfromEach <- read.csv("C:/Users/user/Desktop/biomind/addis/icog/eachData/moses/Mus_musculusgenesfromEach.csv", stringsAsFactors=FALSE)
Rattus_norvegicusgenesfromEach <- read.csv("C:/Users/user/Desktop/biomind/addis/icog/eachData/moses/Rattus_norvegicusgenesfromEach.csv", stringsAsFactors=FALSE)

# apply Escore to old, incorrect combo counts
singleList <- list(Homo_sapiensgenesfromEach, Mus_musculusgenesfromEach, Rattus_norvegicusgenesfromEach)
names(singleList) <- c("humanIndy", "mouseIndy", "ratIndy")
singleList <- lapply(singleList, function (x) x[-4][, c(1, 3, 2)])
for(i in 1:3) names(singleList[[i]]) <- c("feature", "Freq", "level")
for(i in 1:3) singleList[[i]]$level <- ifelse(singleList[[i]]$level, "up", "down")
singleList <- lapply(singleList, function(x) x[, c(4, 2, 3)])
singleList$mouseIndy$Hsym <- mouseMap$symbol[match(singleList$mouseIndy$feature, mouseMap$SYMBOL)]
singleList$ratIndy$Hsym <- ratMap$symbol[match(singleList$ratIndy$feature, ratMap$SYMBOL)]
for(i in 2:3) {
  singleList[[i]] <- na.omit(singleList[[i]][, c(4, 2, 3)])
  names(singleList[[i]])[1] <- "feature"
}

# save old new version of Escores
lapply(singleList, function(x) write.csv(x, file = paste("eddie files/", names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]], "_updown.csv", sep = ""), row.names = FALSE))

# get complete data set gene list from 2 available moses log files
alldata <- getMout(lines = 11, drop = 1)
alldata <- Mout2str(alldata)
write.csv(alldata[[2]], file = paste("eddie files/alldata_updown.csv", sep = ""), row.names = FALSE)

# use 2 meta runs from kelly, 3rd finished 4 days later
load("~/GitHub/alldata/bestOfRun_ 1 .rdata")
alldata_cum <- allMoses_1_hx5.bestN5
load("~/GitHub/alldata/bestOfRun_ 2 .rdata")
allCombos <- list(all1 = alldata_cum$combo, all2 = allMoses_1_hx5.bestN5$combo)
allCombos <- lapply(allCombos, combo2fcount)
mosesCombos(alldata_cum, "all1")
mosesCombos(allMoses_1_hx5.bestN5, "all2")

load("~/GitHub/alldata/bestOfRun_ 3 .rdata")
allCombos$all3 <- combo2fcount(allMoses_1_hx5.bestN5$combo)
mosesCombos(allMoses_1_hx5.bestN5, "all3")

# in 300 cross validated runs using the same seed (-r1), there were 418 unique features with only 32 common to all 3 meta-runs
Reduce(intersect, lapply(allCombos, function(x) x$feature))
#  [1] "ACPP"      "ADCYAP1R1" "ANGPT2"    "ATP1A1"    "BMP6"      "GCLM"      "GM2A"      "IFRD1"     "MAZ"       "MMP8"     
# [11] "MTDH"      "MYCN"      "NDUFS1"    "NF1"       "PFKFB1"    "SHH"       "SNAP23"    "SSTR2"     "TDO2"      "TRDMT1"   
# [21] "TXNL1"     "UBE2B"     "UGCG"      "UTRN"      "VAV2"      "VPS26A"    "VWF"       "WDR61"     "XRCC5"     "YBX3"     
# [31] "YWHAH"     "ZNF148"   
sapply(allCombos, dim)
#      all1 all2 all3
# [1,]  276  251  191
# [2,]    3    3    3
length(Reduce(union, lapply(allCombos, function(x) x$feature)))
# [1] 418

# make array table
array.data <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/aging signatures/array data.csv")
array.data <- array.data[order(array.data$organism, array.data$platform, array.data$sample.number),]
write.csv(array.data, file = "final.csv", row.names = FALSE)

