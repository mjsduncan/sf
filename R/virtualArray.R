### icog virtualArray test

xx<- jpmEset $r_hippo
yy<- jpmEset $r_heart
annotation(xx) <- "rae230a"
annotation(yy) <- "rgu34a"


virtArray_mouse <-NULL
virtArray_mouse[["wBatchEffects"]] <- virtualArrayExpressionSets(all_expression_sets=c(xx, yy))

### apply virtual array by chip
chips <- ls(pattern = "GPL")

# make pheno dataframes.  fix lung age catagories!
age.cat$lung <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1 ,1)
age.pdata <- lapply(age.cat, function(x) data.frame(age = factor(x, labels = c("young", "old")), row.names = colnames(jpmExp[[names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]]]])))
# chip.pdata <- lapply(chiplist, function(x) age.pdata[x])

# addd simplified age pheno data
# chipEsets <- mapply(function(x, y) mapply(function(a, b) {pData(a) <- b; return(a)}, x, y, SIMPLIFY = FALSE), chipEsets, chip.pdata, SIMPLIFY = FALSE)
load("~/GitHub/stevia/data/jpmEset.rdata")
jpmEset <- mapply(function(x, y) {pData(x) <- y; return(x)}, jpmEset, age.pdata, SIMPLIFY = FALSE)

# annotate with current chip annotation packages
dbpkg <- sapply(names(jpmEset), function(x) arrays[x, "bioc_package"])
dbpkg <- sapply(dbpkg, function(x) substr(x, 1, nchar(x) - 3))
jpmEset <- mapply(function(x, y) {annotation(x) <- y; return(x)}, jpmEset, dbpkg, SIMPLIFY = FALSE)

# log2 transform expression values
jpmEset <- lapply(jpmEset, function(x) {exprs(x) <- log2(exprs(x)); return(invisible(x))})

# apply virtualArray
attach(jpmEset)
chipVAs <- list()
for(c in chips) {
  for(n in get(c)) {
    chipVAs[[n]] <- virtualArrayExpressionSets(all_expression_sets = n, identifier = "ENTREZID")
  }
}

#### functions

rownames2col <- function(df, colname) {
  output <- cbind(row.names(df), df)
  colnames(output)[1] <- colname
  return(output)
  }

col2rownames <- function(df, colname, removecol=FALSE){
  row.names(df) <- df[,colname]
  if(removecol){df[,colname] <- NULL}
  return(df)
  }

# this doesn't work as a function
# list.names <-  function(x) names(eval(sys.call(1)[[2]]))[substitute(x)[[3]]]

### MergeMaid
library("MergeMaid")

# map gene symbols to probe names
jpmExpSym <- mapply(function(x, y) {rownames(x) <- y[match(row.names(x), y[, 1]), 3]; return(invisible(x))}, jpmExp, probe2gene, SIMPLIFY = FALSE)

# lots of NAs!
sapply(jpmExpSym, function(x) summary(is.na(row.names(x))))
#       h_brain   muscle1   muscle2   muscle3   muscle4   muscle5   muscle    kidney1   kidney2   m_brain   m_hippo   liver    
# FALSE "12079"   "21053"   "15775"   "21053"   "15775"   "12088"   "113"     "4572"    "3360"    "22209"   "11944"   "11944"  
# TRUE  "481"     "1162"    "6802"    "1162"    "6802"    "445"     "22577"   "235"     "683"     "22892"   "544"     "544"    
# 
# m_heart   lung      cochlea   hemato_stem myo_progen r_hippo   stromal   spinal_cord oculomotor skeletal_ms extraoc_ms
# FALSE "9686"    "22209"   "22196"   "21305"     "9156"     "12925"   "12922"   "12925"     "12924"    "7701"      "7701"    
# TRUE  "357"     "22892"   "22841"   "21564"     "319"      "2998"    "2998"    "2998"      "2998"     "1098"      "1098"    
# 
# laryngeal_ms r_heart   CA1_hipp2
# FALSE "7701"       "7701"    "5853"   
# TRUE  "1098"       "1098"    "757"    

# remove rows with NA row names
jpmExpSym <- lapply(jpmExpSym, function(x) x[!is.na(row.names(x)),])

# aggregate duplicate named rows
jpmExpSym <- lapply(jpmExpSym, function(x) aggregate(x, by=list(row.names(x)), FUN=median))

# turn back into matrices
jpmExpSym <- lapply(jpmExpSym, function(x) as.matrix(col2rownames(x, "Group.1", TRUE)))

merged <-mergeExprs(jpmExpSym[[1]], jpmExpSym[[2]])

jpmExpSymNCN <- lapply(jpmExpSym, function(x){colnames(x) <- NULL; return(invisible(x))})

merged <-mergeExprs(jpmExpSymNCN[[1]], jpmExpSymNCN[[2]])
