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
arrays <- rbind(arrays, arrays[8,], arrays[17,])
row.names(arrays)[25:26] <- c("oculomotor", "extraoc_ms")
row.names(arrays)[1:24] <- names(geodsets)[match(arrays[1:24,]$from_acc, geodsets)]
arrays$bioc_package <- paste(arrays$bioc_package, ".db", sep = "")

# install miacroarray annotation packages
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

# constrouct sample ages
ages <- vector("list", length(geodsets))
names(ages) <- names(geodsets)

# # human
# $h_brain$age
#     106 y 20 - 39 y 40 - 59 y 60 - 79 y 80 - 99 y 
#         1         8         7         6         8 
gdsDt$h_brain$description[1] <- "Value for GSM27015: 26 year old male; src: human frontal cortex"
ages$h_brain <- as.numeric(word(gdsDt$h_brain$description, 4))
# $muscle1$age -- females on U133A chip
# 20-29 years 65-71 years 
#           7           8 
ages$muscle1 <- as.numeric(word(gdsDt$muscle1$description, 9))
# $muscle2$age -- females on U133B chip
# 20-29 years 65-71 years 
#           7           8 
ages$muscle2 <- as.numeric(word(gdsDt$muscle2$description, 9))
# $muscle3$age -- males on U133A chip
# 21-27 years 67-75 years 
#           7           8 
ages$muscle3 <- as.numeric(word(gdsDt$muscle3$description, 9))
# $muscle4$age -- males on U133B chip
# 21-27 years 67-75 years 
#           7           8 
ages$muscle4 <- as.numeric(word(gdsDt$muscle4$description, 9))
# $muscle5$age
# 21-31 year 62-77 year 
#          6          6 
ages$muscle5 <- c(rep(21 + (31 - 21)/2, 6), rep(62 + (77 - 62)/2, 6))

for(x in seq_along(ages[1:6])) names(ages[1:6][[x]]) <- as.character(dataTable(gdsData[1:6][[x]])@columns$sample)

# mouse
# $muscle$age
#   old young 
#    10     5 
# "old" & "young" are 25 months & 5 months
ages$muscle <- c(rep(5, 5), rep(25, 5))
# $muscle$protocol
# calorie-restricted diet             normal diet 
#                       5                      10 
gdsExp$muscle <- gdsExp$muscle[,as.character(gdsDt$muscle$sample[gdsDt$muscle$protocol != 'calorie-restricted diet'])]

# $kidney1$protocol
# calorie-restricted        control fed 
#                  5                 10 
gdsExp$kidney1 <- gdsExp$kidney1[,as.character(gdsDt$kidney1$sample[gdsDt$kidney1$protocol != 'calorie-restricted'])]
# $kidney1$age
# 30 month  5 month 
#       10        5 
ages$kidney1 <- c(rep(5, 5), rep(30, 5))
# $kidney2$protocol
# calorie-restricted        control fed 
#                  5                 10 
gdsExp$kidney2 <- gdsExp$kidney2[,as.character(gdsDt$kidney2$sample[gdsDt$kidney2$protocol != 'calorie-restricted'])]
# $kidney2$age
# 30 month  5 month 
#       10        5 
ages$kidney2 <- c(rep(5, 5), rep(30, 5))

for(x in seq_along(ages[7:9])) names(ages[7:9][[x]]) <- 
  as.character(subset(dataTable(gdsData[7:9][[x]])@columns, word(protocol, 1) != "calorie-restricted", select = sample)[, 1])

# $m_brain$development.stage - 3- to 6-month-old (young adult) and 20- to 24-month-old (aged)
# adult  aged 
#     6     6  
ages$m_brain <- c(rep(4.5, 3), rep(22, 3))
# $m_brain$agent
#    LPS saline 
#      6      6 
names(ages$m_brain) <- as.character(subset(dataTable(gdsData$m_brain)@columns, agent != "LPS", select = sample)[, 1])
gdsExp$m_brain <- gdsExp$m_brain[,as.character(gdsDt$m_brain$sample[gdsDt$m_brain$agent != 'LPS'])]

# $m_hippo$age
# 15 months  2 months 
#        14         9 
ages$m_hippo <- c(rep(2, 9), rep(15, 14))
names(ages$m_hippo) <- as.character(dataTable(gdsData$m_hippo)@columns$sample)

# $liver$`genotype/variation`
#     Snell wild type 
#         7         7  
gdsExp$liver <- gdsExp$liver[,as.character(gdsDt$liver$sample[gdsDt$liver$`genotype/variation` != 'Snell'])]
# $liver$age
# 22 m  6 m 
#    6    8 
ages$liver <- c(rep(6, 4), rep(22, 3))
names(ages$liver) <- as.character(subset(dataTable(gdsData$liver)@columns, `genotype/variation` != "Snell", select = sample)[, 1])

# $m_heart$development.stage
#          1 week        12 month         3 month          4 week         5 month 
#               3               6               3               3               3 
# embryo day 12.5  neonatal day 1 
#               3               3 
gdsExp$m_heart <- gdsExp$m_heart[, 13:24]
ages$m_heart <- c(rep(3, 3), rep(5, 3), rep(12, 6))
names(ages$m_heart) <- as.character(dataTable(gdsData$m_heart)@columns$sample)[13:24]

# $lung$strain - c57 strain lives longer than dba strain
# C57BL/6J   DBA/2J 
#        9        6  
# $lung$age
# 18 mo  2 mo 26 mo 
#     6     6     3 
ages$lung <- c(rep(2, 3), rep(18, 3), rep(26, 6), rep(2, 3), rep(18, 3))
names(ages$lung) <- as.character(dataTable(gdsData$lung)@columns$sample)

# $cochlea$age
# 15 mo  4 mo 
#     6     3  
ages$cochlea <- c(rep(4, 3), rep(15, 3))
# $cochlea$protocol
# calorie restricted diet             normal diet 
#                       3                       6 
names(ages$cochlea) <- as.character(subset(dataTable(gdsData$cochlea)@columns, protocol == "normal diet", select = sample)[, 1])
gdsExp$cochlea <- gdsExp$cochlea[,as.character(gdsDt$cochlea$sample[gdsDt$cochlea$protocol == "normal diet"])]
# $hemato_stem$age
#   2 - 3 mo 22 - 23 mo 
#          3          5 
ages$hemato_stem <- c(rep(2.5, 3), rep(22.5, 5))
names(ages$hemato_stem) <- as.character(dataTable(gdsData$hemato_stem)@columns$sample)

# $myo_progen$age
# 23 mo  8 mo 
#     2     2 
ages$myo_progen <- c(8, 8, 23, 23)
names(ages$myo_progen) <- as.character(dataTable(gdsData$myo_progen)@columns$sample)

## rat 
# $r_hippo$age - Young adult (4–6 months of age) and aged (24–26 months of age)
#   old young 
#    49    29 
ages$r_hippo <- c(rep(5, 29), rep(25, 49))
names(ages$r_hippo) <- as.character(dataTable(gdsData$r_hippo)@columns$sample)

# $r_hippo$disease.state
#       impaired not applicable     unimpaired 
#             19             20             39  
# $r_hippo$protocol
# 21 days post-training       5 days training             untrained 
#                    28                    30                    20 

# construct remaining rat sample ages
ages[19:26] <- lapply(gdsData[19:26], function (x) as.numeric(word(as.character(dataTable(x)@columns$age, 1))))
ages$r_heart <- c(rep(3.5, 5), rep(21, 6))
for(x in seq_along(ages[19:26])) names(ages[19:26][[x]]) <- as.character(dataTable(gdsData[19:26][[x]])@columns$sample)
# $stromal$age
# 15 m  3 m 
#    2    3  
# $stromal$agent
# dexamethasone     untreated 
#             2             3 
gdsExp$stromal <- gdsExp$stromal[, -c(3, 5)]
ages$stromal <- ages$stromal[-c(3, 5)]

# $spinal_cord$tissue
# oculomotor nucleus        spinal cord 
#                  9                  9  
gdsExp$spinal_cord <- gdsExp$spinal_cord[, 1:9]
ages$spinal_cord <- ages$spinal_cord[1:9]
# $spinal_cord$age
# 18 m 30 m  6 m 
#    6    6    6 
# $oculomotor$tissue
# oculomotor nucleus        spinal cord 
#                  9                  9  
gdsExp$oculomotor <- gdsExp$oculomotor[, 10:18]
ages$oculomotor <- ages$oculomotor[10:18]
# $oculomotor$age
# 18 m 30 m  6 m 
#    6    6    6 
# $skeletal_ms$tissue
# extensor digitorum longus        extraocular muscle 
#                        12                        12  
gdsExp$skeletal_ms <- gdsExp$skeletal_ms[, 1:12]
ages$skeletal_ms <- ages$skeletal_ms[1:12]
# $skeletal_ms$age
# 18 m 30 m  6 m 
#    8    8    8 

# $extraoc_ms$tissue
# extensor digitorum longus        extraocular muscle 
#                        12                        12  
gdsExp$extraoc_ms <- gdsExp$extraoc_ms[, 13:24]
ages$extraoc_ms <- ages$extraoc_ms[13:24]
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

# extract probe and gene names from data sets -- lots ofduplicates!
probe2gene <- vector("list", length(geodsets))
names(probe2gene) <- names(geodsets)
for(i in 1:length(geodsets)) {probe2gene[[i]] <- dataTable(gdsData[[i]])@table[, 1:2]}
# probe2gene <- lapply(probe2gene, function(x) lapply(x, as.character))
probe2gene <- lapply(probe2gene, function(x) data.frame(probe = x[[1]], gene = x[[2]], stringsAsFactors = FALSE))
# count & percent of duplicates
sapply(probe2gene, function(x) sum(duplicated(x[[2]])))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle      kidney1      kidney2      m_brain      m_hippo        liver 
#         3137         8189         5796         8189         5796         3141         8741          897          445        18337         2874         2874 
#      m_heart         lung      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor  skeletal_ms   extraoc_ms laryngeal_ms 
#         2343        18337        18337        18337         2874         2479         2479         2479         2479         2718         2718         2718 
#      r_heart    CA1_hipp2 
#         2718         2718 
sapply(probe2gene, function(x) round(sum(duplicated(x[[2]]))/length(x[[2]]), 2))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle      kidney1      kidney2      m_brain      m_hippo        liver 
#         0.25         0.37         0.26         0.37         0.26         0.25         0.39         0.14         0.07         0.41         0.23         0.23 
#      m_heart         lung      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor  skeletal_ms   extraoc_ms laryngeal_ms 
#         0.19         0.41         0.41         0.41         0.23         0.16         0.16         0.16         0.16         0.31         0.31         0.31 
#      r_heart    CA1_hipp2 
#         0.31         0.31 

save(gdsData, gdsDt, gdsExp, probe2gene, file = "gdsData.rdata")

# add current anotations to probe2gene
probe2gene <- lapply(probe2gene, function(x) data.frame(probe = x$ID_REF, GEOgene = x$IDENTIFIER, stringsAsFactors = FALSE))

probe2ez <- function(probes, db) {
  select(db, probes, columns = c("SYMBOL","ENTREZID"))
}

for(n in names(probe2gene)) {
  if(!is.null(probe2gene[[n]]))  {
  require(arrays[n, "bioc_package"], character.only = TRUE)
  probe2gene[[n]] <- merge(probe2gene[[n]], probe2ez(probe2gene[[n]]$probe, get(arrays[n, "bioc_package"])), by.x = "probe", by.y = "PROBEID", all = TRUE)
  }
}

### functions

# convert  class "GEOData" Table method output to matrix
GEOfix <- function(df) {
 out <- lapply(df, function (vec) if(class(vec) == "character") vec <- as.numeric(vec) else vec <- as.character(vec))
 out <- as.data.frame(out, row.names = as.character(out[[1]]))[, -1]
 as.matrix(out[, -1])
}


