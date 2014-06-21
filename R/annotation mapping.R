### get current human homologues for mouse & rat expression data sets
# uses jpmPos, jpmNeg from 'jpm analysis.R'

# load databases
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")
library("hom.Mm.inp.db")
library("hom.Rn.inp.db")


# make list of dataframes mapping affy probe to uniprot
probe2uniprot <- function(probes, db) {
  select(db, probes, columns = c("SYMBOL","UNIPROT"))
}

jpmPosUni <- vector("list", length(jpmPos))
names(jpmPosUni) <- names(jpmPos)
jpmNegUni <- vector("list", length(jpmNeg))
names(jpmNegUni) <- names(jpmNeg)
for(n in names(jpmPos)) {
  if(!is.null(jpmPos[[n]]))  {
  require(arrays[n, "bioc_package"], character.only = TRUE)
  jpmPosUni[[n]] <- probe2uniprot(row.names(jpmPos[[n]]), get(arrays[n, "bioc_package"]))
  jpmNegUni[[n]] <- probe2uniprot(row.names(jpmNeg[[n]]), get(arrays[n, "bioc_package"]))
  }
}
jpmPosUni <- lapply(jpmPosUni, function(x) x[!(is.na(x[, 2]) & is.na(x[, 3])),])
jpmNegUni <- lapply(jpmNegUni, function(x) x[!(is.na(x[, 2]) & is.na(x[, 3])),])

# make vector of inparanoid species
inp.species <- character(20)
names(inp.species) <- names(jpmPos)[7:26]
for(n in names(inp.species)) {
  if(arrays[n, "organism"] == "Mus musculus") `<-`(inp.species[[n]], "MUSMU") else `<-`(inp.species[[n]], "RATNO")
}

#make list of lists mapping uniprot to human homologue uniprot
jpmPosHom <- vector("list", length(jpmPos) - 6)
names(jpmPosHom) <- names(jpmPos)[7:26]
jpmNegHom <- vector("list", length(jpmNeg) - 6)
names(jpmNegHom) <- names(jpmNeg)[7:26]

for(n in names(jpmPosHom)) {
  if(!is.null(jpmPos[[n]]))  {
  jpmPosHom[[n]] <- inpIDMapper(jpmPosUni[[n]]$UNIPROT, inp.species[n], "HOMSA", "UNIPROT", "UNIPROT", TRUE, TRUE, TRUE)
  jpmNegHom[[n]] <- inpIDMapper(jpmNegUni[[n]]$UNIPROT, inp.species[n], "HOMSA", "UNIPROT", "UNIPROT", TRUE, TRUE, TRUE)
  }
}

# make dataframe from named list of character vectors with as.data.frame.list.R then combine by data set
jpmPosHom.df <- lapply(jpmPosHom, function(x) lapply(x, as.data.frame))
jpmPosHom.df <- lapply(jpmPosHom.df, function(x) lapply(x, function(y) `names<-`(y, "homo")))
for(n in names(jpmPosHom.df)) {
  if(length(jpmPosHom.df[[n]]) > 0) {
    for(i in 1:length(jpmPosHom.df[[n]])) {
      prot <- rep(names(jpmPosHom.df[[n]])[i], dim(jpmPosHom.df[[n]][[i]])[1])
      jpmPosHom.df[[n]][[i]] <- cbind(prot, jpmPosHom.df[[n]][[i]])
    }
  }
}
jpmPosHom.df <- lapply(jpmPosHom.df, function(x) if(length(x) > 0) do.call(rbind, x))

jpmNegHom.df <- lapply(jpmNegHom, function(x) lapply(x, as.data.frame))
jpmNegHom.df <- lapply(jpmNegHom.df, function(x) lapply(x, function(y) `names<-`(y, "homo")))
for(n in names(jpmNegHom.df)) {
  if(length(jpmNegHom.df[[n]]) > 0) {
    for(i in 1:length(jpmNegHom.df[[n]])) {
      prot <- rep(names(jpmNegHom.df[[n]])[i], dim(jpmNegHom.df[[n]][[i]])[1])
      jpmNegHom.df[[n]][[i]] <- cbind(prot, jpmNegHom.df[[n]][[i]])
    }
  }
}
jpmNegHom.df <- lapply(jpmNegHom.df, function(x) if(length(x) > 0) do.call(rbind, x))

# merge probe dfs with homologue dfs
jpmPos.df <- mapply(merge, jpmPosUni[7:26][!sapply(jpmPosUni[7:26], is.null)], jpmPosHom.df[!sapply(jpmPosHom.df, is.null)], MoreArgs = list(by.x = "UNIPROT", by.y = "prot", all = TRUE), SIMPLIFY = FALSE)
jpmNeg.df <- mapply(merge, jpmNegUni[7:26][!sapply(jpmNegUni[7:26], is.null)], jpmNegHom.df[!sapply(jpmNegHom.df, is.null)], MoreArgs = list(by.x = "UNIPROT", by.y = "prot", all = TRUE), SIMPLIFY = FALSE)

jpmPos.df <- lapply(jpmPos.df, function(x) x[!is.na(x[, 4]),])
jpmNeg.df <- lapply(jpmNeg.df, function(x) x[!is.na(x[, 4]),])

# make list of symbols of homologues
jpmPos.sym <- vector("list", length(jpmPos.df))
names(jpmPos.sym) <- names(jpmPos.df)
for(n in names(jpmPos.df)) {
  jpmPos.sym[[n]] <- intraIDMapper(as.character(jpmPos.df[[n]]$homo), "HOMSA", destIDType = "SYMBOL", keepMultGeneMatches = TRUE)
}
jpmPos.sym <- lapply(jpmPos.sym, unlist)
jpmPos.sym <- lapply(jpmPos.sym, function(x) data.frame(homo = names(x), homSym = x, stringsAsFactors = FALSE))

jpmNeg.sym <- vector("list", length(jpmNeg.df))
names(jpmNeg.sym) <- names(jpmNeg.df)
for(n in names(jpmNeg.df)) {
  jpmNeg.sym[[n]] <- intraIDMapper(as.character(jpmNeg.df[[n]]$homo), "HOMSA", destIDType = "SYMBOL", keepMultGeneMatches = TRUE)
}
jpmNeg.sym <- lapply(jpmNeg.sym, unlist)  
jpmNeg.sym <- lapply(jpmNeg.sym, function(x) data.frame(homo = names(x), homSym = x, stringsAsFactors = FALSE))

# merge jpmXxx.df & jpmXxx.sym and get final map from probes to human homologue symbol
jpmPos.df <- mapply(merge, jpmPos.df, jpmPos.sym, MoreArgs = list(all = TRUE, stringsAsFactors = FALSE), SIMPLIFY = FALSE)
jpmNeg.df <- mapply(merge, jpmNeg.df, jpmNeg.sym, MoreArgs = list(all = TRUE, stringsAsFactors = FALSE), SIMPLIFY = FALSE)
jpmPos.df <- lapply(jpmPos.df, function(x) x[!is.na(x[, 5]),])
jpmNeg.df <- lapply(jpmNeg.df, function(x) x[!is.na(x[, 5]),])

# gained 5% over just capitalizing mouse & rat gene symbols
summary(unlist(lapply(jpmPos.df, function(x) sum(toupper(x[, 4]) == x[, 5])/dim(x)[1])))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.9008  0.9463  0.9565  0.9557  0.9729  1.0000 
summary(unlist(lapply(jpmNeg.df, function(x) sum(toupper(x[, 4]) == x[, 5])/dim(x)[1])))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.6667  0.9342  0.9514  0.9373  0.9651  1.0000 

# 