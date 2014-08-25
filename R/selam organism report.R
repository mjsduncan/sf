Note: its all the same for human, just change the name from mouse to human 

# For mouses
#### save edited names for later 
mouse_MN <- Pgs_mouseBack_normalized
dimnames(mouse_MN)[[2]] <- make.names(dimnames(mouse_MN)[[2]])
edited_Mname  <- dimnames(mouse_MN)[[2]][dimnames(mouse_MN)[[2]] != dimnames(Pgs_mouseBack_normalized)[[2]]]

### load 10 best result for moses
for(i in 1:10){
  load(paste("~/icog/moses_R/moses_R_result /mouse/meta_run",i,"/mouses/mouses_",i,"_best.rdata",sep = ""))
  assign(paste("mouse_", i, "_hx5.best", sep = ""), mouse_1_hx5.best)
}

### bind best of runs 
bestofmouse <- mouse_1_hx5.best
for(i in 2:10){ bestofmouse <- lbind(bestofmouse , get(ls(pattern = "_hx5.best")[i]))}

#step 1 . generate mouse moses genes 
## generating result matrix from moses result
bestofmouse_res <- matrix(nrow = length(unique(bestofmouse$features$feature)), ncol = 5)
dimnames(bestofmouse_res)[[2]]  <- c("SYMBOL", "GENENAME","ENTREZID", "Score", "Low")

bestofmouse_res[ , "SYMBOL"] <- unique(bestofmouse$features$feature)
for (i in 1:dim(bestofmouse_res)[1]){  bestofmouse_res[ , "Low"][i] <- bestofmouse$features$low[bestofmouse$features$feature  ==  bestofmouse_res [ , "SYMBOL"][i]]}
for (i in 1:dim(bestofmouse_res)[1]){  bestofmouse_res[ , "Score"][i] <- bestofmouse$features$score[bestofmouse$features$feature  ==  bestofmouse_res [ , "SYMBOL"][i]]}

##### the rest parametres need symbol name correction 
#1# removing X added in symbol name by make.names 
for (i in 1:length(bestofmouse_res[ , "SYMBOL"])) {
  if (sum(bestofmouse_res[ , "SYMBOL"][i] %in% edited_Mname) > 0 ){
    if (sum (grep("^X", bestofmouse_res[ , "SYMBOL"][i]) )== 1) 
      bestofmouse_res[ , "SYMBOL"][i] <- substring(bestofmouse_res[ , "SYMBOL"][i], 2, nchar( bestofmouse_res[ , "SYMBOL"][i])) 
    else 
      bestofmouse_res[ , "SYMBOL"][i] <-  gsub("[.]","-",bestofmouse_res[ , "SYMBOL"][i])
  }
}

###### adding enterzID and gene name to the matrix

for (i in 1:dim(bestofmouse_res)[1]){
  if(is.na(bestofmouse_res[ , "ENTREZID"][i] )){
    if (length(mget(bestofmouse_res[ , "SYMBOL"],org.Mm.egALIAS2EG)[i][[1]]) > 1) {
      bestofmouse_res[ , "ENTREZID"][i] <-  mget(bestofmouse_res[ , "SYMBOL"],org.Mm.egALIAS2EG)[i][[1]][1] 
      for (j in 2: length(mget(bestofmouse_res[ , "SYMBOL"],org.Mm.egALIAS2EG)[i][[1]])){
        bestofmouse_res <- insertRow( bestofmouse_res, dim(bestofmouse_res)[1]+1,c( bestofmouse_res[ , "SYMBOL"][i] , bestofmouse_res[ , "GENENAME"][i] ,
                                                      bestofmouse_res[ , "ENTREZID"][i],bestofmouse_res[ , "Score"][i] , bestofmouse_res[ , "Low"][i] ))
        bestofmouse_res[ , "ENTREZID"][dim(bestofmouse_res)[1]] <- mget(bestofmouse_res[ , "SYMBOL"],org.Mm.egALIAS2EG)[i][[1]][j]
      } 
    }
    else 
      bestofmouse_res[ , "ENTREZID"][i] <-  mget(bestofmouse_res[ , "SYMBOL"][i],org.Mm.egALIAS2EG)[[1]]
  }
}

for (i in 1:dim(bestofmouse_res)[1]){ bestofmouse_res[ , "GENENAME"][i] <- mget(bestofmouse_res[ , "ENTREZID"],org.Mm.egGENENAME)[i][[1]]}

########## save files 
  best_lowT  <- subset(bestofmouse_res, bestofmouse_res[ , "Low"] == TRUE)
  best_lowF  <- subset(bestofmouse_res, bestofmouse_res[ , "Low"] != TRUE)
  bst  <- rbind(as.matrix(sortByCol(as.data.frame(best_lowF),"SYMBOL"))  ,as.matrix(sortByCol(as.data.frame(best_lowT),"SYMBOL")))
  write.csv(bst, file = "best_all.csv", row.names =FALSE)
 

##### step 2 ... generate best combo excel file 
##mike## it's not necessary to format berfore saving .csv files, save the whole thing and format in the spreadsheet!

bestcombo <- matrix(nrow = length(bestofmouse$combo), ncol = 15)
dimnames(bestcombo)[[2]]  <- c("Balanced Accuracy","combo",  "out",	"Accuracy",	"AccuracyLower",	"AccuracyUpper","AccuracyNull",	"AccuracyPValue",
                                    "Sensitivity",	"Specificity",	"Pos Pred Value",	"Neg Pred Value",	"Prevalence","Detection Rate",	"Detection Prevalence")

bestcombo[ , "combo"]                <- bestofmouse$combo
bestcombo[ , "Balanced Accuracy"]    <- format(round(bestofmouse$score[["Balanced Accuracy"]], digits=2), nsmall = 2)
bestcombo[ , "out"]                  <- format(round(bestofmouse$rank, 2), nsmall = 2)
bestcombo[ , "Accuracy"]             <- format(round(bestofmouse$score[["Accuracy"]], digits=2) , nsmall = 2)
bestcombo[ , "AccuracyLower"]        <- format(round(bestofmouse$score[["AccuracyLower"]], digits=2), nsmall = 2)
bestcombo[ , "AccuracyUpper"]        <- format(round(bestofmouse$score[["AccuracyUpper"]], digits=2), nsmall = 2)
bestcombo[ , "AccuracyNull"]         <- format(round(bestofmouse$score[["AccuracyNull"]], digits=2), nsmall = 2)
bestcombo[ , "AccuracyPValue"]       <- format(round(bestofmouse$score[["AccuracyPValue"]], digits=2) , nsmall = 2)
bestcombo[ , "Sensitivity"]          <- format(round(bestofmouse$score[["Sensitivity"]], digits=2), nsmall = 2)
bestcombo[ , "Specificity"]          <- format(round(bestofmouse$score[["Specificity"]], digits=2), nsmall = 2)
bestcombo[ , "Pos Pred Value"]       <- format(round(bestofmouse$score[["Pos Pred Value"]], digits=2), nsmall = 2)
bestcombo[ , "Neg Pred Value"]       <- format(round(bestofmouse$score[["Neg Pred Value"]], digits=2), nsmall = 2)
bestcombo[ , "Prevalence"]           <- format(round(bestofmouse$score[["Prevalence"]], digits=2), nsmall = 2)
bestcombo[ , "Detection Rate"]       <- format(round(bestofmouse$score[["Detection Rate"]], digits=2), nsmall = 2)
bestcombo[ , "Detection Prevalence"] <- format(round(bestofmouse$score[["Detection Prevalence"]], digits=2), nsmall = 2)


bestcombo <- subset(bestcombo, duplicated(bestcombo[ , "combo"] ) != TRUE)
write.csv(bestcombo, file = "bestcombo.csv", row.names =FALSE)


##### step: 3Gene Onlology 

#Gene to GO BP CC MF Conditional test for over-representation 

selectedEntrezIds   <- unlist( bestofmouse_res[ , "ENTREZID"])
entrezUniverse      <- unlist(mget(colnames(Pgs_mouseBack_normalized)[-1], org.Mm.egALIAS2EG))
names(entrezUniverse) <- NULL

GOparamBP <- new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=unique(entrezUniverse), 
                annotation="org.Mm.eg.db",ontology="BP",
                pvalueCutoff=0.001, conditional=TRUE, testDirection="over")
GOparamCC <- new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=unique(entrezUniverse), 
                 annotation="org.Mm.eg.db",ontology="CC",
                 pvalueCutoff=0.001, conditional=TRUE, testDirection="over")
GOparamMF <- new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=unique(entrezUniverse), 
                 annotation="org.Mm.eg.db",ontology="MF",
                 pvalueCutoff=0.001, conditional=TRUE, testDirection="over")

hgOverBP <- hyperGTest(GOparamBP)
hgOverCC <- hyperGTest(GOparamCC)
hgOverMF <- hyperGTest(GOparamMF)

Go_mouse <- rbind(as.matrix(summary(hgOverBP)), as.matrix(summary(hgOverCC)))
Go_mouse <- rbind(Go_mouse, as.matrix(summary(hgOverMF)))
Go_mouse <- insertCol(Go_mouse,  dim(Go_mouse)[2]+1)
dimnames(Go_mouse)[[2]][ dim(Go_mouse)[2]] <- "Ontology"
dimnames(Go_mouse)[[2]][1] <- "GOID"
Go_mouse[ , "Ontology"] <-  Ontology(Go_mouse[ , "GOID"])
## sort by pvalue

write.csv(as.matrix(sortByCol(as.data.frame(Go_mouse),"Pvalue")), file = "GO_OVER.csv", row.names =FALSE)

#Gene to GO BP Conditional test for under-representation 

GOparamBP <- new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=unique(entrezUniverse), 
                 annotation="org.Mm.eg.db",ontology="BP",
                 pvalueCutoff=0.001, conditional=TRUE, testDirection="under")
GOparamCC <- new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=unique(entrezUniverse), 
                 annotation="org.Mm.eg.db",ontology="CC",
                 pvalueCutoff=0.001, conditional=TRUE, testDirection="under")
GOparamMF <- new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=unique(entrezUniverse), 
                 annotation="org.Mm.eg.db",ontology="MF",
                 pvalueCutoff=0.001, conditional=TRUE, testDirection="under")

hgUnderBP <- hyperGTest(GOparamBP)
hgUnderCC <- hyperGTest(GOparamCC)   # has no value
hgUnderMF <- hyperGTest(GOparamMF)   # has no value


GoU_mouse <- rbind(as.matrix(summary(hgUnderBP)), as.matrix(summary(hgUnderCC)))
GoU_mouse <- rbind(GoU_mouse, as.matrix(summary(hgUnderMF)))
GoU_mouse <- insertCol(GoU_mouse,  dim(GoU_mouse)[2]+1)
dimnames(GoU_mouse)[[2]][ dim(GoU_mouse)[2]] <- "Ontology"
dimnames(GoU_mouse)[[2]][1] <- "GOID"
GoU_mouse[ , "Ontology"] <-  Ontology(GoU_mouse[ , "GOID"])

## sort by pval

write.csv(as.matrix(sortByCol(as.data.frame(GoU_mouse),"Pvalue")), file = "GO_UNDER.csv", row.names =FALSE)
