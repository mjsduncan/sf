### XDE this uses jpmEset from "virtualArray.R"
library("XDE")

# make pheno dataframes.  fix lung age catagories!
age.cat$lung <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1 ,1)
age.pdata <- lapply(age.cat, function(x) data.frame(age = factor(x, labels = c("young", "old")), row.names = colnames(jpmExp[[names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]]]])))

# check for missing expression values:  kidney1 & 2, CA1_hipp2 (1/3), m_heart (1/7) myo_pyogen (1/4) are missing a significant portion
aggEsets <- lapply(jpmEset, exprs)
sapply(aggEsets, function(x) summary(x))
# 13 of 26 data sets have 4 to 1900 rows missin > 50% with imputation by means rather than 10 nearest neihbors
for(n in names(aggEsets)) aggEsets[[n]] <- impute::impute.knn(aggEsets[[n]])
sapply(aggEsets, function(x) summary(is.na(x$data)))
aggEsets <- lapply(aggEsets, function(x) x$data)

# aggregate by symbol
dbpkg <- sapply(names(jpmEset), function(x) arrays[x, "bioc_package"])
aggEsets <- mapply(function(x, y) PGSEA::aggregateExprs(x, package = y, using = "SYMBOL", FUN = median), jpmEset, dbpkg)

# turn back into Esets
aggEsets <- mapply(function(x, y) {exprs(x) <- y; return(invisible(x))}, jpmEset, aggEsets)
# this didn't work
aggEsets <- mapply(function(x, y) {pData(x) <- y; return(x)}, aggEsets, age.pdata, SIMPLIFY = FALSE)
dbpkg <- sapply(dbpkg, function(x) substr(x, 1, nchar(x) - 3))
aggEsets <- mapply(function(x, y) {annotation(x) <- y; return(x)}, jpmEset, dbpkg, SIMPLIFY = FALSE)
save(aggEsets, file = "data/aggEsets.rdata")
load("~/GitHub/stevia/data/aggEsets.rdata")

# make ExpressionSetLists by chip
chiplist <- lapply(chips, get)
names(chiplist) <- chips
chipEsets <- lapply(chiplist, function(x) aggEsets[x])

# add simplified age pheno data
chip.pdata <- lapply(chiplist, function(x) age.pdata[x])
chipEsets <- mapply(function(x, y) mapply(function(a, b) {pData(a) <- b; return(a)}, x, y, SIMPLIFY = FALSE), chipEsets, chip.pdata, SIMPLIFY = FALSE)
chipEsets <- lapply(chipEsets, function(x) new("ExpressionSetList", x))
sapply(chipEsets, validObject)
# GPL1261  GPL341   GPL81   GPL85   GPL96   GPL97 
#    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 

# impute missing values

# make list of XdeParameter objects
chipParams <- lapply(chipEsets, function(x) new("XdeParameter", esetList= x, phenotypeLabel = "age"))
