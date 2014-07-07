### make homologene homologue map
homologene <- read.delim("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data", header=F, stringsAsFactors = FALSE)
names(homologene) <- c("HID", "taxID", "egID", "symbol", "prot_gi", "prot_acc")
# keep mouse, rat, human
# 10090  Mus musculus
# 10116	Rattus norvegicus
# 9606  Homo sapiens
homologene <- subset(homologene, taxID %in% c(10090, 10116, 9606))
write.csv(homologene, file = "homo.csv", col.names = TRUE, row.names = FALSE)
homo <- read.csv("~/GitHub/stevia/data/homo.csv", stringsAsFactors = FALSE)
homo.hs <- subset(homo, taxID == 9606, select = c(HID, egID, symbol))
homo.ro <- subset(homo, taxID != 9606, select = c(HID, egID, symbol))

# make mapping homologene dfs
# load databases
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")
library("hom.Mm.inp.db")
library("hom.Rn.inp.db")


# make list of dataframes mapping affy probe to entrez id
probe2ez <- function(probes, db) {
  select(db, probes, columns = c("SYMBOL","ENTREZID"))
}

jpmPosEz.df <- vector("list", length(jpmPos))
names(jpmPosEz.df) <- names(jpmPos)
jpmNegEz.df <- vector("list", length(jpmNeg))
names(jpmNegEz.df) <- names(jpmNeg)
for(n in names(jpmPos)) {
  if(!is.null(jpmPos[[n]]))  {
  require(arrays[n, "bioc_package"], character.only = TRUE)
  jpmPosEz.df[[n]] <- probe2Ez.df(row.names(jpmPos[[n]]), get(arrays[n, "bioc_package"]))
  jpmNegEz.df[[n]] <- probe2Ez.df(row.names(jpmNeg[[n]]), get(arrays[n, "bioc_package"]))
  }
}
jpmPosEz.df <- lapply(jpmPosEz.df, function(x) x[!(is.na(x[, 2]) & is.na(x[, 3])),])
jpmNegEz.df <- lapply(jpmNegEz.df, function(x) x[!(is.na(x[, 2]) & is.na(x[, 3])),])

# add columns with human homologue ed & symbol
Ez.df2homo <- function(Ez.df) {
  homo.hs[homo.hs$HID == homo.ro$HID[homo.ro$egID == Ez.df], 2:3]
}
for(n in names(jpmPos)) {
  if(!is.null(jpmPos[[n]]))  {
  jpmPosEz.df[[n]] <- cbind(jpmPosEz.df[[n]], t(sapply(jpmPosEz.df[[n]]$ENTREz.dfID, Ez.df2homo)))
  jpmNegEz.df[[n]] <- cbind(jpmNegEz.df[[n]], t(sapply(jpmNegEz.df[[n]]$ENTREz.dfID, Ez.df2homo)))
  }
}

# copy human symbols & entrez ids to homo columns for symetry
for(n in names(human)) {
  jpmPosEz.df[[n]]$egID <- jpmPosEz.df[[n]]$ENTREZID
  jpmPosEz.df[[n]]$symbol <- jpmPosEz.df[[n]]$SYMBOL
  jpmNegEz.df[[n]]$egID <- jpmNegEz.df[[n]]$ENTREZID
  jpmNegEz.df[[n]]$symbol <- jpmNegEz.df[[n]]$SYMBOL
}

# convert list columns to multiple rows
row2nrows <- function(row) {
  out <- row
  n = length(out$egID[[1]])
  if(n == 0) {  
    out$egID <- ""
    out$symbol <- ""
    return(out)
  }
  if(n == 1) {  
    out$egID <- as.character(unlist(out$egID))
    out$symbol <- unlist(out$symbol)
    return(out)
  }
out3 <- data.frame(PROBEID = character(0), SYMBOL = character(0), ENTREZID = character(0), egID = character(0), symbol = character(0))
  for(i in 1:n) {
  out2 <- out[, 1:3]
  out2$egID <- as.character(out$egID[[1]][i])
  out2$symbol <- out$symbol[[1]][i]
  out3 <- rbind(out2, out3)
  }
  return(out3)
}

# apply row2nrows to whole dataframe
dflist2nrows <- function(df) {
  out <- data.frame(PROBEID = character(0), SYMBOL = character(0), ENTREZID = character(0), egID = character(0), symbol = character(0))
  for(i in seq_len(dim(df)[1])) {
    out <- rbind(out, row2nrows(df[i,]))
  }
  return(out)
}

jpmPosEz.homo <- lapply(jpmPosEz.df, dflist2nrows)
jpmNegEz.homo <- lapply(jpmNegEz.df, dflist2nrows)



# compare inparanoid and homologene maps
# homologene map
jpmPosEz.hs <- unique(unlist(sapply(jpmPosEz[1:6], function(x) x$ENTREZID)))
length(jpmPosEz.hs)
# [1] 3373

jpmPosEz.nhs <- unique(unlist(sapply(jpmPosEz[7:26], function(x) x$egID)))
length(jpmPosEz.nhs)
# [1] 5067

jpmPosEz <- unique(c(jpmPosEz.hs, jpmPosEz.nhs))
length(jpmPosEz)
# [1] 7291

jpmNegEz.hs <- unique(unlist(sapply(jpmNegEz[1:6], function(x) x$ENTREZID)))
length(jpmNegEz.hs)
# [1] 3842

jpmNegEz.nhs <- unique(unlist(sapply(jpmNegEz[7:26], function(x) x$egID)))
length(jpmNegEz.nhs)
# [1] 5118

jpmNegEz <- unique(c(jpmNegEz.hs, jpmNegEz.nhs))
length(jpmNegEz)
# [1] 7630

length(intersect(jpmPosEz, jpmNegEz))
# [1] 3650

jpmEz <- select(org.Hs.eg.db, unique(c(jpmPosEz, jpmNegEz)), columns = c("ENTREZID", "SYMBOL"))
dim(jpmEz)
# [1] 11271     2

# inparaniod map size
inparPos <- character(0)
lapply(jpmPos.df, function(x) inparPos <<- c(inparPos, x$homSym))
inparNeg <- character(0)
lapply(jpmNeg.df, function(x) inparNeg <<- c(inparNeg, x$homSym))
inpar <- unique(c(inparPos, inparNeg))

length(inpar)
# [1] 5633

length(unique(c(jpmPosEz.nhs, jpmNegEz.nhs)))
# [1] 7985

### homo map wins!!

# count mouse & rat homologues
# rat map size
ratPos <- character(4)
lapply(jpmPosEz.df[names(rat)], function(x) ratPos <<- rbind(ratPos, x[, 2:5]))
ratNeg <- character(4)
lapply(jpmNegEz.df[names(rat)], function(x) ratNeg <<- rbind(ratNeg, x[, 2:5]))
homoRat <- unique(rbind(ratPos, ratNeg))

# non-uniqueness of map
summary(as.factor(sapply(homoRat$egID, length)))
#    0    1    2    3    4 
#  806 4382   18    2    1 

# mouse map size
mousePos <- character(4)
lapply(jpmPosEz.df[names(mouse)], function(x) mousePos <<- rbind(mousePos, x[, 2:5]))
mouseNeg <- character(4)
lapply(jpmNegEz.df[names(mouse)], function(x) mouseNeg <<- rbind(mouseNeg, x[, 2:5]))
homoMouse <- unique(rbind(mousePos, mouseNeg))

# non-uniqueness of map
summary(as.factor(sapply(homoMouse$egID, length)))
#    0    1    2    5    7 
#  835 5641   24    1    1 

homoMap <- rbind(homoRat, homoMouse)
