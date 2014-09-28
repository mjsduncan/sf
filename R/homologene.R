### make homologene homologue map
# uses jpmPos, jpmNeg from 'jpm analysis.R'

# download and clean homologene database
homologene <- read.delim("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data", header=F, stringsAsFactors = FALSE)
names(homologene) <- c("HID", "taxID", "egID", "symbol", "prot_gi", "prot_acc")
# keep mouse, rat, human
# 10090  Mus musculus
# 10116	Rattus norvegicus
# 9606  Homo sapiens
homologene <- subset(homologene, taxID %in% c(10090, 10116, 9606))
write.csv(homologene, file = "homo.csv", col.names = TRUE, row.names = FALSE)

#load saved homologene file and split into human and rodent parts
homo <- read.csv("~/GitHub/stevia/data/homo.csv", colClasses = "character", stringsAsFactors = FALSE)
homo.hs <- subset(homo, taxID == "9606", select = c(HID, egID, symbol))
homo.ro <- subset(homo, taxID != "9606", select = c(HID, egID, symbol))

# make mapping homologene dfs
# load databases
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")

# make lists of dfs mapping affy probe to entrez id for probes positively and negatively age associated
# function to map vector of probes to symbol and entrez id given array annotation db
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
  jpmPosEz.df[[n]] <- probe2ez(row.names(jpmPos[[n]]), get(arrays[n, "bioc_package"]))
  jpmNegEz.df[[n]] <- probe2ez(row.names(jpmNeg[[n]]), get(arrays[n, "bioc_package"]))
  }
}
jpmPosEz.df <- lapply(jpmPosEz.df, function(x) x[!(is.na(x[, 2]) & is.na(x[, 3])),])
jpmNegEz.df <- lapply(jpmNegEz.df, function(x) x[!(is.na(x[, 2]) & is.na(x[, 3])),])

# add columns with human homologue ed & symbol
ez2homo <- function(ez) {
  if(is.na(ez)) return(data.frame(egID = NA_character_, symbol = NA_character_))
  homo.hs[homo.hs$HID == homo.ro$HID[homo.ro$egID == ez], 2:3]
}

for(n in names(jpmPos)) {
  if(!is.null(jpmPos[[n]]))  {
  jpmPosEz.df[[n]] <- cbind(jpmPosEz.df[[n]], t(sapply(jpmPosEz.df[[n]]$ENTREZID, ez2homo)))
  jpmNegEz.df[[n]] <- cbind(jpmNegEz.df[[n]], t(sapply(jpmNegEz.df[[n]]$ENTREZID, ez2homo)))
  }
}

# copy human symbols & entrez ids to homo columns for symetry
for(n in names(human)) {
  jpmPosEz.df[[n]]$egID <- jpmPosEz.df[[n]]$ENTREZID
  jpmPosEz.df[[n]]$symbol <- jpmPosEz.df[[n]]$SYMBOL
  jpmNegEz.df[[n]]$egID <- jpmNegEz.df[[n]]$ENTREZID
  jpmNegEz.df[[n]]$symbol <- jpmNegEz.df[[n]]$SYMBOL
}

# duplicate rows with probes mapping to multiple genes
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

# function to apply row2nrows to whole dataframe
dflist2nrows <- function(df) {
  out <- data.frame(PROBEID = character(0), SYMBOL = character(0), ENTREZID = character(0), egID = character(0), symbol = character(0))
  for(i in seq_len(dim(df)[1])) {
    out <- rbind(out, row2nrows(df[i,]))
  }
  return(out)
}

jpmPosEz.homo <- lapply(jpmPosEz.df, dflist2nrows)
jpmNegEz.homo <- lapply(jpmNegEz.df, dflist2nrows)

# difference in row count after accounting for multiple genes per probe
mapply(function(x, y) dim(y)[1] - dim(x)[1], jpmPosEz.df, jpmPosEz.homo)
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#            0            0            0            0            0            0            0 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#           11            2            4            6            4            1            6 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#            4            9            3           15            2            0            4 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#            2            0            1            3            7 
mapply(function(x, y) dim(y)[1] - dim(x)[1], jpmNegEz.df, jpmNegEz.homo)
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#            0            0            0            0            0            0            0 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#            1            7            1            1            1            0            6 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#            0            2            5            4            4            3            1 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#            0            0            4            0            1 


# compare homologue count to original count

# compare duplicate counts
sapply(probe2gene, function (x) summary(duplicated(x$probe)))
#       h_brain   muscle1   muscle2   muscle3   muscle4   muscle5   muscle    kidney1   kidney2  
# FALSE "12625"   "22283"   "22645"   "22283"   "22645"   "12626"   "22690"   "6584"    "6595"   
# TRUE  "1260"    "2263"    "1024"    "2263"    "1024"    "1260"    "5"       "659"     "576"    
# 
# m_brain   m_hippo   liver     m_heart   lung      cochlea   hemato_stem myo_progen r_hippo  
# FALSE "45101"   "12488"   "12488"   "12654"   "45101"   "45101"   "45101"     "12488"    "15923"  
# TRUE  "1917"    "1197"    "1197"    "1236"    "1917"    "1917"    "1917"      "1197"     "1167"   
# 
# stromal   spinal_cord oculomotor skeletal_ms extraoc_ms laryngeal_ms r_heart   CA1_hipp2
# FALSE "15923"   "15923"     "15923"    "8799"      "8799"     "8799"       "8799"    "8799"   
# TRUE  "1167"    "1167"      "1167"     "1120"      "1120"     "1120"       "1120"    "1120"   

sapply(jpmPosEz.homo, function (x) summary(duplicated(x$PROBEID)))
#       h_brain   muscle1   muscle2   muscle3   muscle4   muscle5   muscle    kidney1   kidney2  
# FALSE "939"     "1334"    "728"     "767"     "547"     "290"     "14"      "193"     "109"    
# TRUE  "64"      "161"     "59"      "67"      "33"      "39"      "1"       "73"      "31"     

#       m_brain   m_hippo   liver     m_heart   lung      cochlea   hemato_stem myo_progen r_hippo  
# FALSE "683"     "544"     "635"     "264"     "780"     "307"     "935"       "629"      "1312"   
# TRUE  "73"      "77"      "143"     "50"      "82"      "12"      "89"        "66"       "207"    

#       stromal   spinal_cord oculomotor skeletal_ms extraoc_ms laryngeal_ms r_heart   CA1_hipp2
# FALSE "422"     "417"       "511"      "305"       "207"      "246"        "351"     "330"    
# TRUE  "45"      "34"        "50"       "51"        "8"        "72"         "42"      "108"    

sapply(jpmNegEz.homo, function (x) summary(duplicated(x$PROBEID)))
#       h_brain   muscle1   muscle2   muscle3   muscle4   muscle5   muscle    kidney1   kidney2  
# FALSE "1294"    "1296"    "834"     "1255"    "558"     "415"     "15"      "232"     "214"    
# TRUE  "101"     "69"      "47"      "111"     "37"      "48"      "1"       "19"      "39"     

#       m_brain   m_hippo   liver     m_heart   lung      cochlea   hemato_stem myo_progen r_hippo  
# FALSE "466"     "611"     "796"     "216"     "471"     "223"     "764"       "781"      "1060"   
# TRUE  "24"      "44"      "69"      "24"      "42"      "12"      "86"        "98"       "67"     

#       stromal   spinal_cord oculomotor skeletal_ms extraoc_ms laryngeal_ms r_heart   CA1_hipp2
# FALSE "379"     "702"       "518"      "145"       "124"      "678"        "194"     "418"    
# TRUE  "46"      "60"        "65"       "16"        "34"       "94"         "33"      "43"     


# count mouse & rat homologues -- TODO: silence output to increase speed
# rat map size
ratPos <- character(4)
lapply(jpmPosEz.df[names(rat)], function(x) ratPos <<- rbind(ratPos, x[, 2:5]))
ratNeg <- character(4)
lapply(jpmNegEz.df[names(rat)], function(x) ratNeg <<- rbind(ratNeg, x[, 2:5]))
homoRat <- unique(rbind(ratPos, ratNeg))

# non-uniqueness of map
summary(as.factor(sapply(homoRat$egID, length)))
#    0    1    2    3    4 
#  934 4965   20    1    1 

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

# add homologues to probe2gene
# add homologue to probe2gene element row
homoRow <- function(dfrow) {
  homo <- ez2homo(dfrow$ENTREZID)
  nRows <- dim(homo)[1]
  if(nRows == 0) return(cbind(dfrow, egID = NA_character_, symbol = NA_character_))
  cbind(dfrow[rep(1, nRows),], homo)
}

# step through each row of probe2gene data frame
homoDf <- function(df) {
  out <- data.frame(probe = character(0), GEOgene = character(0), SYMBOL = character(0), ENTREZID = character(0), egID = character(0), symbol = character(0))
  for(i in 1:dim(df)[1]) {
    out <- rbind(out, homoRow(df[i,]))
  }
  return(out)
}

probe2homo <- lapply(probe2gene[7:26], homoDf)
save(probe2homo, file = "probe2homo.rdata")

### annotationTools vignette example
library(annotationTools)
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")
homologene <- read.delim("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data", header=F)
homoMap <- select(org.Hs.eg.db, keys(org.Hs.eg.db, keytype = "SYMBOL"), columns = c("GENENAME", "ENTREZID"), keytype = "SYMBOL")
names(homologene) <- c("HID", "taxID", "egID", "symbol", "prot_gi", "prot_acc")
# keep mouse, rat, human
# 10090  Mus musculus
# 10116  Rattus norvegicus
# 9606  Homo sapiens
homoMap$mouseEG <- NA
for(i in seq_along(homoMap$mouseEG)) homoMap$mouseEG[i] <- unlist(getHOMOLOG(homoMap$ENTREZID[i], 10090, homologene), rec = FALSE)
homoMap$mouseEG <- as.character(homoMap$mouseEG)
mouseSYM <- select(org.Mm.eg.db, homoMap$mouseEG, columns = "SYMBOL", keytype = "ENTREZID")
homoMap$mouseSYM <- mouseSYM[match(homoMap$mouseEG, mouseSYM$ENTREZID), 2]

homoMap$ratEG <- NA
for(i in seq_along(homoMap$ratEG)) homoMap$ratEG[i] <- try(unlist(getHOMOLOG(homoMap$ENTREZID[i], 10116, homologene), rec = FALSE))
# stopped at human symbol index 10500, restarted.
for(i in 10500:length(homoMap$ratEG)) homoMap$ratEG[i] <- try(unlist(getHOMOLOG(homoMap$ENTREZID[i], 10116, homologene), rec = FALSE))
homoMap$ratEG <- as.character(homoMap$ratEG)
# this genrates an error:
#  Error in sqliteExecStatement(con, statement, bind.data) : 
#   RS-DBI driver: (error in statement: near "errors": syntax error) 
# 11 sqliteExecStatement(con, statement, bind.data) 
# 10 sqliteQuickSQL(conn, statement, ...) 
# 9 dbGetQuery(conn, SQL) 
# 8 dbGetQuery(conn, SQL) 
# 7 dbQuery(dbConn(x), sql) 
# 6 .extractData(x, cols = cols, keytype = keytype, keys = keys) 
# 5 .legacySelect(x, keys, cols, keytype, jointype) 
# 4 .select(x, keys, columns, keytype, jointype = jointype) 
# 3 .selectWarnJT(x, keys, columns, keytype, jointype = jointype, 
#     kt = kt, ...) 
# 2 select(org.Rn.eg.db, homoMap$ratEG, columns = "SYMBOL", keytype = "ENTREZID") 
# 1 select(org.Rn.eg.db, homoMap$ratEG, columns = "SYMBOL", keytype = "ENTREZID") 
# error messages mistakenly appended to string sent to mysql server apparently generated error.  rebooting fixed problem

ratSYM <- select(org.Rn.eg.db, homoMap$ratEG, columns = "SYMBOL", keytype = "ENTREZID")
homoMap$ratSYM <- ratSYM[match(homoMap$ratEG, ratSYM$ENTREZID), 2]

# alternate transform vector of human EGids then match to homoMap
# partial attempt didn't match above!
# 107 in below not in equivalent part above, 2 vice versa.
# TODO: followup
# ratEGpart <- getHOMOLOG(homoMap$ENTREZID[10500:length(homoMap$ENTREZID)], 10116, homologene)
# ratEGpart2many <- ratEGpart[sapply(ratEGpart, length) > 1]
# <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/org data/mouse.csv", stringsAsFactors=FALSE)
# mouse$eg <- unlist(mget(mouse$feature, org.Mm.egSYMBOL2EG, ifnotfound= list(NA)), rec = FALSE)
# mouse$homoEG <- getHOMOLOG(mouse$eg,10090,homologene)

# salvage using probe2homo homology map for data sets
load("~/GitHub/stevia/data/probe2homo.rdata")
arrays <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/arrays.tab", stringsAsFactors=FALSE)
ratData <- arrays$row.names[arrays$organism == "Rattus norvegicus"]
ratMap <- unique(Reduce(rbind, probe2homo[ratData]))
# remove rows without homology
ratMap <- subset(ratMap, !is.na(egID))
homoMap$ratSYM <- ratMap[match(homoMap$ratEG, ratMap$ENTREZID), 3]

# make mouse map & use rodent maps to get homolog symbols for chip gene lists
mouseData <- arrays$row.names[arrays$organism == "Mus musculus"]
mouseMap <- unique(Reduce(rbind, probe2homo[mouseData]))
mouseMap <- subset(mouseMap, !is.na(egID))

save(homoMap, mouseMap, ratMap, file = "~/GitHub/stevia/data/homoMaps.rdata")

