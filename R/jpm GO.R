### jpm paper GO analysis
# get GO codes for each organism from significant genes
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")

# function to map vector of entrez ids to GO codes by organism
ez2go <- function(genes, db) {
  select(db, genes, columns = c("SYMBOL","GO"))
}

humGOpos <- lapply(jpmPos[names(human)], function(x) ez2go(x$ENTREZID, org.Hs.eg.db))
humGOneg <- lapply(jpmNeg[names(human)], function(x) ez2go(x$ENTREZID, org.Hs.eg.db))

### try GOfunction
library("GOFunction")
# compare only positively expressed genes with all human genes tested
refGenes <- unique(unlist(lapply(probe2gene[names(human)], function(x) x$ENTREZID)))
refGenes <- refGenes[!is.na(refGenes)]
refGenes <- c(refGenes, onlyPos[!(onlyPos %in% refGenes)])
onlyPos <- setdiff(agePos$ENTREZID, ageNeg$ENTREZID)
sigPosTerm <- GOFunction(onlyPos, refGenes, organism="org.Hs.eg.db", ontology="BP", filename="sigPosTerm")
