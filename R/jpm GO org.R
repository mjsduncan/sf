### jpm paper GO analysis
# get GO codes for each organism from significant genes
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")

# # function to map vector of entrez ids to GO codes by organism  TODO: do for rat & mouse arrays
# ez2go <- function(genes, db) {
#   select(db, genes, columns = c("SYMBOL","GO"))
# }
# 
# humGOpos <- lapply(jpmPos[names(human)], function(x) ez2go(x$ENTREZID, org.Hs.eg.db))
# humGOneg <- lapply(jpmNeg[names(human)], function(x) ez2go(x$ENTREZID, org.Hs.eg.db))

### try GOfunction
library("GOFunction")

# try on moses genes  REDO: use moses runs on organism array sets, not moses results from chip sets grouped by organism!
# TODO:  use better code from "jpm GO chip.R"
GOposM <- list()
GOnegM <- list()
org <- c("org.Mm.eg.db", "org.Rn.eg.db", "org.Hs.eg.db")
names(org) <- c("mouse", "rat", "human")

refGenes <-  list(unique(c(probe2gene[[10]]$ENTREZID, probe2gene[[17]]$ENTREZID)),
          unique(c(probe2gene[[19]]$ENTREZID, probe2gene[[22]]$ENTREZID)),
          unique(c(probe2gene[[3]]$ENTREZID, probe2gene[[4]]$ENTREZID))
                  )
posGenes <- list(na.omit(unique(c(jpmMhigh[[1]]$ENTREZID, jpmMhigh[[3]]$ENTREZID))),
                 na.omit(unique(c(jpmMhigh[[2]]$ENTREZID, jpmMhigh[[4]]$ENTREZID))),
                 na.omit(unique(c(jpmMhigh[[5]]$ENTREZID, jpmMhigh[[6]]$ENTREZID)))
)
negGenes <- list(na.omit(unique(c(jpmMlow[[1]]$ENTREZID, jpmMlow[[3]]$ENTREZID))),
                 na.omit(unique(c(jpmMlow[[2]]$ENTREZID, jpmMlow[[4]]$ENTREZID))),
                 na.omit(unique(c(jpmMlow[[5]]$ENTREZID, jpmMlow[[6]]$ENTREZID)))
)
for(i in 1:3) {
  GOposM[[names(org)[i]]] <- list(
    GOFunction(posGenes[[i]], refGenes[[i]], organism=org[i], ontology="BP", filename=paste(names(org)[i], "_sigPosTermBP", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5),
    GOFunction(posGenes[[i]], refGenes[[i]], organism=org[i], ontology="MF", filename=paste(names(org)[i], "_sigPosTermMF", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5),
    GOFunction(posGenes[[i]], refGenes[[i]], organism=org[i], ontology="CC", filename=paste(names(org)[i], "_sigPosTermCC", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5)
  )
  GOnegM[[names(org)[i]]] <- list(
    GOFunction(negGenes[[i]], refGenes[[i]], organism=org[i], ontology="BP", filename=paste(names(org)[i], "_sigNegTermBP", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5),
    GOFunction(negGenes[[i]], refGenes[[i]], organism=org[i], ontology="MF", filename=paste(names(org)[i], "_sigNegTermMF", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5),
    GOFunction(negGenes[[i]], refGenes[[i]], organism=org[i], ontology="CC", filename=paste(names(org)[i], "_sigNegTermCC", sep = ""), fdrth = 0.5, ppth =   0.5, pcth = 0.5, poth = 0.5, peth = 0.5)
  )
}

### try GOstats
mouseHiBP <- new("GOHyperGParams",
              geneIds = posGenes[[1]],
              universeGeneIds = refGenes[[1]],
              annotation = org[1],
              ontology = "BP",
              pvalueCutoff = .05,
              conditional = TRUE,
              testDirection = "over")
mouseHiMF <- mouseHiBP
ontology(mouseHiMF) <- "MF"
mouseHiCC <- mouseHiBP
ontology(mouseHiCC) <- "CC"
mouseHi <- list(hyperGTest(mouseHiBP), hyperGTest(mouseHiMF), hyperGTest(mouseHiCC))
names(mouseHi) <- c("BP", "MF", "CC")
mouseGOhi <- lapply(mouseHi, summary)

mouseLoBP <- new("GOHyperGParams",
              geneIds = negGenes[[1]],
              universeGeneIds = refGenes[[1]],
              annotation = org[1],
              ontology = "BP",
              pvalueCutoff = .05,
              conditional = TRUE,
              testDirection = "over")
mouseLoMF <- mouseLoBP
ontology(mouseLoMF) <- "MF"
mouseLoCC <- mouseLoBP
ontology(mouseLoCC) <- "CC"
mouseLo <- list(hyperGTest(mouseLoBP), hyperGTest(mouseLoMF), hyperGTest(mouseLoCC))
names(mouseLo) <- c("BP", "MF", "CC")
mouseGOlo <- lapply(mouseLo, summary)

# rat
ratHiBP <- new("GOHyperGParams",
              geneIds = posGenes[[2]],
              universeGeneIds = refGenes[[2]],
              annotation = org[2],
              ontology = "BP",
              pvalueCutoff = .05,
              conditional = TRUE,
              testDirection = "over")
ratHiMF <- ratHiBP
ontology(ratHiMF) <- "MF"
ratHiCC <- ratHiBP
ontology(ratHiCC) <- "CC"
ratHi <- list(hyperGTest(ratHiBP), hyperGTest(ratHiMF), hyperGTest(ratHiCC))
names(ratHi) <- c("BP", "MF", "CC")
ratGOhi <- lapply(ratHi, summary)

ratLoBP <- new("GOHyperGParams",
              geneIds = negGenes[[2]],
              universeGeneIds = refGenes[[2]],
              annotation = org[2],
              ontology = "BP",
              pvalueCutoff = .05,
              conditional = TRUE,
              testDirection = "over")
ratLoMF <- ratLoBP
ontology(ratLoMF) <- "MF"
ratLoCC <- ratLoBP
ontology(ratLoCC) <- "CC"
ratLo <- list(hyperGTest(ratLoBP), hyperGTest(ratLoMF), hyperGTest(ratLoCC))
names(ratLo) <- c("BP", "MF", "CC")
ratGOlo <- lapply(ratLo, summary)

# human
humanHiBP <- new("GOHyperGParams",
              geneIds = posGenes[[3]],
              universeGeneIds = refGenes[[3]],
              annotation = org[3],
              ontology = "BP",
              pvalueCutoff = .05,
              conditional = TRUE,
              testDirection = "over")
humanHiMF <- humanHiBP
ontology(humanHiMF) <- "MF"
humanHiCC <- humanHiBP
ontology(humanHiCC) <- "CC"
humanHi <- list(hyperGTest(humanHiBP), hyperGTest(humanHiMF), hyperGTest(humanHiCC))
names(humanHi) <- c("BP", "MF", "CC")
humanGOhi <- lapply(humanHi, summary)

humanLoBP <- new("GOHyperGParams",
              geneIds = negGenes[[3]],
              universeGeneIds = refGenes[[3]],
              annotation = org[3],
              ontology = "BP",
              pvalueCutoff = .05,
              conditional = TRUE,
              testDirection = "over")
humanLoMF <- humanLoBP
ontology(humanLoMF) <- "MF"
humanLoCC <- humanLoBP
ontology(humanLoCC) <- "CC"
humanLo <- list(hyperGTest(humanLoBP), hyperGTest(humanLoMF), hyperGTest(humanLoCC))
names(humanLo) <- c("BP", "MF", "CC")
humanGOlo <- lapply(humanLo, summary)
