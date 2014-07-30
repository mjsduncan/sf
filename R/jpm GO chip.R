### jpm paper GO analysis
# get GO codes for each organism from significant genes
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")

# # function to map vector of entrez ids to GO codes
# ez2go <- function(genes, db) {
#   select(db, genes, columns = c("SYMBOL","GO"))
# }
# 
# # go annotations per gene, TODO:  map by chip group
# humGOpos <- lapply(jpmPos[names(human)], function(x) ez2go(x$ENTREZID, org.Hs.eg.db))
# humGOneg <- lapply(jpmNeg[names(human)], function(x) ez2go(x$ENTREZID, org.Hs.eg.db))

### try GOfunction
library("GOFunction")

# try on moses genes
GOposM <- list()
GOnegM <- list()
org <- c("org.Mm.eg.db", "org.Rn.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Hs.eg.db", "org.Hs.eg.db")
names(org) <- names(jpmMhigh)
refGenes <- list()
posGenes <- list()
negGenes <- list()

for(n in names(jpmMhigh)) {
refGenes[[n]] <- unique(probe2gene[[get(toupper(n))[1]]]$ENTREZID)
posGenes[[n]] <- na.omit(jpmMhigh[[n]]$ENTREZID)
negGenes[[n]] <- na.omit(jpmMlow[[n]]$ENTREZID)
GOposM[[n]] <- list(
GOFunction(posGenes[[n]], refGenes[[n]], organism=org[n], ontology="BP", filename=paste(n, "_sigPosTermBP", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5),
GOFunction(posGenes[[n]], refGenes[[n]], organism=org[n], ontology="MF", filename=paste(n, "_sigPosTermMF", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5),
GOFunction(posGenes[[n]], refGenes[[n]], organism=org[n], ontology="CC", filename=paste(n, "_sigPosTermCC", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5)
)
GOnegM[[n]] <- list(
GOFunction(negGenes[[n]], refGenes[[n]], organism=org[n], ontology="BP", filename=paste(n, "_sigNegTermBP", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5),
GOFunction(negGenes[[n]], refGenes[[n]], organism=org[n], ontology="MF", filename=paste(n, "_sigNegTermMF", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5),
GOFunction(negGenes[[n]], refGenes[[n]], organism=org[n], ontology="CC", filename=paste(n, "_sigNegTermCC", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5)
)
}

### try GOstats
library("GOstats")

# function to make GO hypergeometric parameter object from GOstats package
GOHGPobj <- function(intGenes, refGenes, annDB, ont = "BP", pvalue = .05, cond = TRUE, testDir = "over") {
  new("GOHyperGParams", geneIds = intGenes, universeGeneIds = refGenes, annotation = annDB, ontology = ont, pvalueCutoff = pvalue, conditional = cond, testDirection = testDir)
}

# function to make list of 3 GO hypergeometric parameter objects, one for each GO ontology
GOstat3obj <- function(intGenes, refGenes, annDB, pvalue = .05, cond = TRUE, testDir = "over") {
  out <- vector("list", 3)
  names(out) <- c("BP", "MF", "CC")
  for(n in names(out)) {
    out[[n]] <- GOHGPobj(intGenes, refGenes, annDB, ont = n, pvalue, cond, testDir)
  }
  return(out)
}

# make lists for GOstats objects for high and low results from each moses chip run using gene lists from GOFunction tests above
chipHiObj <- vector("list", 6)
names(chipHiObj) <- names(refGenes)
chipLoObj <- vector("list", 6)
names(chipLoObj) <- names(refGenes)

for(n in names(chipHiObj)) {
chipHiObj[[n]] <- GOstat3obj(posGenes[[n]], refGenes[[n]], org[n])
chipLoObj[[n]] <- GOstat3obj(negGenes[[n]], refGenes[[n]], org[n])
}

# make lists of GOstat results for moses runs by hi & low expression, chip, & GO ontology
chipHi <- lapply(chipHiObj, function(x) lapply(x, hyperGTest))
chipGOhi <- lapply(chipHi, function(x) lapply(x, summary))

chipLo <- lapply(chipLoObj, function(x) lapply(x, hyperGTest))
chipGOlo <- lapply(chipLo, function(x) lapply(x, summary))

save(GOposM, GOnegM, refGenes, posGenes, negGenes, chipHiObj, chipHi, chipGOhi, chipLoObj, chipLo, chipGOlo, file = "mosesGOchip.rdata")

# function to rbind list of dfs adding a list name column
namedList2df <- function(list) {
  if(is.null(names(list))) return("error:  list needs names!")
  out <- list()
  for(n in names(list)) out <- rbind(cbind(list[[n]], name = n), out)
  return(out)
}

# change "GOxxID" to "GOID" so the 3 ontology dfs can be rbinded
chipGOhi <- lapply(chipGOhi, function(x) lapply(x, function(y) {names(y) <- c("GOID", names(y)[-1]); return(y)}))
chipGOlo <- lapply(chipGOlo, function(x) lapply(x, function(y) {names(y) <- c("GOID", names(y)[-1]); return(y)}))

# save rbinded dfs of the 3 ontologies for each chip
for(n in names(chipGOhi)) {
  write.csv(namedList2df(chipGOhi[[n]]), file = paste(n, "_chipGOhi.csv", sep = ""), row.names = FALSE)
  write.csv(namedList2df(chipGOlo[[n]]), file = paste(n, "_chipGOlo.csv", sep = ""), row.names = FALSE)
}
