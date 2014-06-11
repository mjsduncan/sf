### import datasets

# import geo data sets
library("GEOquery")
library("GEOmetadb")
# from de Magalhães 2009 supplement
human <- c(
  h_brain = "GDS707", # Lu, et al., 2004
  muscle1 = "GDS472", muscle2 = "GDS473", # Welle, et al., 2004
  muscle3 = "GDS287", muscle4 = "GDS288", # Welle, et al.,2003
  muscle5 = "GDS156" # Welle, et al., 2002
)
mouse <- c(
  muscle = "GDS2612", # Edwards, et al., 2007
  kidney1 = "GDS355", kidney2 = "GDS356",
  m_brain = "GDS1311", # Godbout, et al., 2005
  m_hippo = "GDS2082", # Verbitsky, et al., 2004
  liver = "GDS2019", # Papaconstantinou, et al., 2005
  m_heart = "GDS40",
  lung = "GDS2929", # Misra, et al., 2007
# not in GEO ftp site #  eye1 = "GDS396", eye2 = "GDS397", #Ida, et al., 2003
  cochlea = "GDS2681", # Someya, et al., 2007
  hemato_stem = "GDS1803", # Rossi, et al., 2005
  myo_progen = "GDS1079" # Beggs, et al., 2004
)
rat <- c(
  r_hippo = "GDS2639", # Rowe, et al., 2007
  stromal = "GDS2231", # Akavia, et al., 2006
  spinal_cord = "GDS1280",
  oculomotor = "GDS1280",
  skeletal_ms = "GDS1279",
  extraoc_ms = "GDS1279",
  laryngeal_ms = "GDS1278", # McMullen and Andrade, 2006
  r_heart = "GDS399", # Dobson, et al., 2003
# not in GEO ftp site # CA1_hipp1 = "GDS2315", # Burger, et al., 2007
  CA1_hipp2 = "GDS520" # Blalock, et al., 2003
)
geodsets <- c(human, mouse, rat)

# get annotation packages for geo datasets.  specify location of metadata mysql file and 
getSQLiteFile(destdir = 'C:/Users/user/Documents')
gpls <- geoConvert(geodsets, out_type = "gpl", sqlite = 'C:/Users/user/Documents/GEOmetadb.sqlite')[[1]]
row.names(gpls) <- names(geodsets[match(gpls[,1], geodsets)])
con <- dbConnect(SQLite(),'C:/Users/user/Documents/GEOmetadb.sqlite')
arrays <- merge(gpls, getBiocPlatformMap(con, bioc='all'), by.x = "to_acc", by.y = "gpl", all.x = TRUE)
dbDisconnect(con)
# file.remove('C:/Users/user/Documents/GEOmetadb.sqlite')
BiocInstaller::biocLite(paste(unique(arrays[["bioc_package"]]), "db", sep = "."))

# dowlnoad geo datasets into list of GEOData objects
gdsData <- vector("list",length(geodsets))
names(gdsData) <- names(geodsets)
for(i in 1:length(geodsets)){gdsData[[i]] <- getGEO(GEO = geodsets[i], destdir = '~/github/stevia/data')}

# make list of gds expression data matrices
gdsExp <- vector("list",length(geodsets))
names(gdsExp) <- names(geodsets)
for(i in 1:length(geodsets)){gdsExp[[i]] <- Table(gdsData[[i]])}
gdsExp <- lapply(gdsExp, GEOfix)
# save(gdsExp, file = "~/github/stevia/Exp/gdsExp.rdata")

# extract relevant samples with class labels & construct numeric age
library("stringr")

# make list of gds data tables.  note most columns are factors
gdsDt <- vector("list",length(geodsets))
names(gdsDt) <- names(geodsets)
for(i in 1:length(geodsets)){gdsDt[[i]] <- dataTable(gdsData[[i]])@columns}
 
# list sample variables to isolate age & controls vs treated if any
lapply(gdsDt, function(df) lapply(df, summary))
# # human
# $h_brain$age
#     106 y 20 - 39 y 40 - 59 y 60 - 79 y 80 - 99 y 
#         1         8         7         6         8 
gdsDt$h_brain$description[1] <- "Value for GSM27015: 26 year old male; src: human frontal cortex"
gdsDt$h_brain$years <- as.numeric(word(gdsDt$h_brain$description, 4))
# $muscle1$age
# 20-29 years 65-71 years 
#           7           8 
gdsDt$muscle1$years <- as.numeric(word(gdsDt$muscle1$description, 9))
# $muscle2$age
# 20-29 years 65-71 years 
#           7           8 
gdsDt$muscle2$years <- as.numeric(word(gdsDt$muscle2$description, 9))
# $muscle3$age
# 21-27 years 67-75 years 
#           7           8 
gdsDt$muscle3$years <- as.numeric(word(gdsDt$muscle3$description, 9))
# $muscle4$age
# 21-27 years 67-75 years 
#           7           8 
gdsDt$muscle4$years <- as.numeric(word(gdsDt$muscle4$description, 9))
# $muscle5$age
# 21-31 year 62-77 year 
#          6          6 
# # mouse
# $muscle$age
#   old young 
#    10     5 
levels(gdsDt$muscle$age) <- c("25 months", "5 months")
# $muscle$protocol
# calorie-restricted diet             normal diet 
#                       5                      10 
gdsExp$muscle <- gdsExp$muscle[,as.character(gdsDt$muscle$sample[gdsDt$muscle$protocol != 'calorie-restricted diet'])]
# $kidney1$protocol
# calorie-restricted        control fed 
#                  5                 10 
# $kidney1$age
# 30 month  5 month 
#       10        5 
# $kidney2$protocol
# calorie-restricted        control fed 
#                  5                 10 
# $kidney2$age
# 30 month  5 month 
#       10        5 
# $m_brain$development.stage
# adult  aged 
#     6     6  
# $m_brain$agent
#    LPS saline 
#      6      6 
# $m_hippo$age
# 15 months  2 months 
#        14         9 
# $liver$`genotype/variation`
#     Snell wild type 
#         7         7  
# $liver$age
# 22 m  6 m 
#    6    8 
# $m_heart$development.stage
#          1 week        12 month         3 month          4 week         5 month 
#               3               6               3               3               3 
# embryo day 12.5  neonatal day 1 
#               3               3 
# $lung$strain
# C57BL/6J   DBA/2J 
#        9        6  
# $lung$age
# 18 mo  2 mo 26 mo 
#     6     6     3 
# $cochlea$age
# 15 mo  4 mo 
#     6     3  
# $cochlea$protocol
# calorie restricted diet             normal diet 
#                       3                       6 
# $hemato_stem$age
#   2 - 3 mo 22 - 23 mo 
#          3          5 
# $myo_progen$age
# 23 mo  8 mo 
#     2     2 
# $r_hippo$age
#   old young 
#    49    29 
## rat 
# $r_hippo$disease.state
#       impaired not applicable     unimpaired 
#             19             20             39  
# $r_hippo$protocol
# 21 days post-training       5 days training             untrained 
#                    28                    30                    20 
# $stromal$age
# 15 m  3 m 
#    2    3  
# $stromal$agent
# dexamethasone     untreated 
#             2             3 
# $spinal_cord$tissue
# oculomotor nucleus        spinal cord 
#                  9                  9  
# $spinal_cord$age
# 18 m 30 m  6 m 
#    6    6    6 
# $oculomotor$tissue
# oculomotor nucleus        spinal cord 
#                  9                  9  
# $oculomotor$age
# 18 m 30 m  6 m 
#    6    6    6 
# $skeletal_ms$tissue
# extensor digitorum longus        extraocular muscle 
#                        12                        12  
# $skeletal_ms$age
# 18 m 30 m  6 m 
#    8    8    8 
# $extraoc_ms$tissue
# extensor digitorum longus        extraocular muscle 
#                        12                        12  
# $extraoc_ms$age
# 18 m 30 m  6 m 
#    8    8    8 
# $laryngeal_ms$age
# 18 m 30 m  6 m 
#    3    3    3 
# $r_heart$age
# 20-22 month   3-4 month 
#           6           5 
# $CA1_hipp2$age
# 14 months 24 months  4 months 
#        10        10         9 




















#   Boxplot for selected GEO samples (note:  the .jpg & .svg files in ~/openbiomind2/results/transplant_samples were exported from RStudio environment)

# order samples by group
ex <- exprs(hlavtx)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("acute rejection","no acute rejection")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(hlavtx)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(hlavtx)))/2),4,2,1))
title <- paste ("GSE2018"," (transplant bronchiolar lavage samples)", " log2 transformed expression levels", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, inset = c(.1, .2), fill=palette())

################################################################
#   construct moses dataset

# get probe x sample log2 normalized expression level matrix from expression set
hlavtx.moses <- exprs(hlavtx)

#median normalize with med.normalize() from "data cleaning.R"
hlavtx.moses <- med.normalize(hlavtx.moses)

# add control binary (cases are transplants resulting in primary graft disfunction)
controls <- c(0,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1)
hlavtx.moses <- t(rbind(controls, hlavtx.moses))

# remove control spots
gpl96.controls <- as.character(fData(hlavtx)$ID[fData(hlavtx)$Platform_SPOTID == "--Control"])
hlavtx.moses <- hlavtx.moses[,!(colnames(hlavtx.moses) %in% gpl96.controls)]

# export file
write.csv(hlavtx.moses, file = "results/transplant_samples/hlavtx_moses.csv")
# apply median norm to df columns
med.normalize <- function(mat) {
  out <- mat
  for (i in seq(dim(mat)[2])) { 
    vect <- mat[,i]
    med <- median(vect, na.rm = TRUE)
    out[,i] <- as.numeric(vect >= med)
  }
  return(out)
}

### functions

# convert  class "GEOData" Table method output to matrix
GEOfix <- function(df) {
 out <- lapply(df, function (vec) if(class(vec) == "character") vec <- as.numeric(vec) else vec <- as.character(vec))
 out <- as.data.frame(out, row.names = as.character(out[[1]]))[, -1]
 as.matrix(out[, -1])
}
