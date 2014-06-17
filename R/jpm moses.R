
### moses
mosesExp <- lapply(gdsExp, impute.matrix)
mosesExp <- lapply(mosesExp, med.normalize)

# add control binary (cases are transplants resulting in primary graft disfunction)
controls <- c(0,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1)
hlavtx.moses <- t(rbind(controls, hlavtx.moses))

# remove control spots
gpl96.controls <- as.character(fData(hlavtx)$ID[fData(hlavtx)$Platform_SPOTID == "--Control"])
hlavtx.moses <- hlavtx.moses[,!(colnames(hlavtx.moses) %in% gpl96.controls)]

# export file
write.csv(hlavtx.moses, file = "results/transplant_samples/hlavtx_moses.csv")


