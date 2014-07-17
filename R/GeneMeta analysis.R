### test GeneMeta individual differential expression measured by t-test
library("GeneMeta")

# make expression sets from gds expression matrices
jpmEset <- lapply(gdsData, GDS2eSet)
for(n in names(jpmEset)) jpmEset[[n]] <- jpmEset[[n]][, colnames(jpmExp[[n]])]

# make catagory vectors from ages
age.cat <- lapply(ages, function(x) as.numeric(med.normalize(as.matrix(x))))
age.cat$h_brain[c(5:9, 23:24)] <- 1
age.cat$liver[1:4] <- 0
age.cat$stromal[1:2] <- 0

# make vectors to combine datasets from the same arrays
GPL1261 <- row.names(arrays)[arrays$to_acc == "GPL1261"]
GPL341 <- row.names(arrays)[arrays$to_acc == "GPL341"]
GPL85 <- row.names(arrays)[arrays$to_acc == "GPL85"]
GPL81 <- row.names(arrays)[arrays$to_acc == "GPL81"]
GPL97 <- row.names(arrays)[arrays$to_acc == "GPL97"]
GPL96 <- row.names(arrays)[arrays$to_acc == "GPL96"]

# apply wrapper function to lists data sets with the same array
gpl.zscore <- vector("list", 6)
names(gpl.zscore) <- unique(arrays$to_acc)[c(1, 4, 7, 9, 11:12)]
for(n in names(gpl.zscore)) gpl.zscore[[n]] <- zScores(jpmEset[get(n)], age.cat[get(n)])

# idr plots with ugly hack to provide titles
lapply(gpl.zscore, function(x) IDRplot(na.omit(x), main = names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]]))

# try fdr method
gpl.fdr <- vector("list", 6)
names(gpl.fdr) <- unique(arrays$to_acc)[c(1, 4, 7, 9, 11:12)]
for(n in names(gpl.fdr)) gpl.fdr[[n]] <- zScoreFDR(jpmEset[get(n)], age.cat[get(n)])
(gpl)
# idr plots with ugly hack to provide titles
lapply(gpl.fdr, function(x) IDRplot(na.omit(x), main = names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]]))
