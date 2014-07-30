### jpm paper GO analysis  uses probe2gene from "data cleaning.R" and agePos & ageNeg from "jpm analysis.R"
# get GO codes for each organism from significant genes
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")

### try GOfunction
library("GOFunction")
# compare only positively expressed genes with all human genes tested
refGenesJPM <- unique(unlist(lapply(probe2gene[names(human)], function(x) x$ENTREZID)))
refGenesJPM <- refGenesJPM[!is.na(refGenesJPM)]

# throw out genes both significantly increasing and decreasing in association with age 
onlyPos <- setdiff(agePos$ENTREZID, ageNeg$ENTREZID)

refGenesJPMPos <- c(refGenesJPM, onlyPos[!(onlyPos %in% refGenesJPM)])
sigPosTermBP <- GOFunction(onlyPos, refGenesJPMPos, organism="org.Hs.eg.db", ontology="BP", filename="sigPosTermBP")
sigPosTermMF <- GOFunction(onlyPos, refGenesJPMPos, organism="org.Hs.eg.db", ontology="MF", filename="sigPosTermMF")
sigPosTermCC <- GOFunction(onlyPos, refGenesJPMPos, organism="org.Hs.eg.db", ontology="CC", filename="sigPosTermCC")

onlyNeg <- setdiff(ageNeg$ENTREZID, agePos$ENTREZID)
refGenesJPMNeg <- c(refGenesJPM, onlyNeg[!(onlyNeg %in% refGenesJPM)])
sigNegTermBP <- GOFunction(onlyNeg, refGenesJPMNeg, organism="org.Hs.eg.db", ontology="BP", filename="sigNegTermBP")
sigNegTermMF <- GOFunction(onlyNeg, refGenesJPMNeg, organism="org.Hs.eg.db", ontology="MF", filename="sigNegTermMF")
sigNegTermCC <- GOFunction(onlyNeg, refGenesJPMNeg, organism="org.Hs.eg.db", ontology="CC", filename="sigNegTermCC")

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

# make lists for GOstats objects for high and low results from gene lists from GOFunction tests above
jpmHiObj <- GOstat3obj(onlyPos, refGenesJPMPos, "org.Hs.eg.db")
jpmLoObj <- GOstat3obj(onlyNeg, refGenesJPMNeg, "org.Hs.eg.db")

# make lists of GOstat results for moses runs by hi & low expression, chip, & GO ontology
jpmHi <- lapply(jpmHiObj, hyperGTest)
jpmGOhi <- lapply(jpmHi, summary)

jpmLo <- lapply(jpmLoObj, hyperGTest)
jpmGOlo <- lapply(jpmLo, summary)

save(refGenesJPMPos, onlyPos, jpmHiObj, jpmHi, jpmGOhi, refGenesJPMNeg, onlyNeg, jpmLoObj, jpmLo, jpmGOlo, file = "mosesGOjpm.rdata")

# function to rbind list of dfs adding a list name column
namedList2df <- function(list) {
  if(is.null(names(list))) return("error:  list needs names!")
  out <- list()
  for(n in names(list)) out <- rbind(cbind(list[[n]], name = n), out)
  return(out)
}

# change "GOxxID" to "GOID" so the 3 ontology dfs can be rbinded
jpmGOhi <- lapply(jpmGOhi, function(y) {names(y) <- c("GOID", names(y)[-1]); return(y)})
jpmGOlo <- lapply(jpmGOlo, function(y) {names(y) <- c("GOID", names(y)[-1]); return(y)})

# save rbinded dfs of the 3 ontologies for each chip
write.csv(namedList2df(jpmGOhi), file = "jpmGOhi.csv", row.names = FALSE)
write.csv(namedList2df(jpmGOlo), file = "jpmGOlo.csv", row.names = FALSE)

