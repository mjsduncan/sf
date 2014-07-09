### make homologene homologue map
homologene <- read.delim("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data", header=F, stringsAsFactors = FALSE)
names(homologene) <- c("HID", "taxID", "egID", "symbol", "prot_gi", "prot_acc")
# keep mouse, rat, human
# 10090  Mus musculus
# 10116	Rattus norvegicus
# 9606  Homo sapiens
homologene <- subset(homologene, taxID %in% c(10090, 10116, 9606))
write.csv(homologene, file = "homo.csv", col.names = TRUE, row.names = FALSE)
homo <- read.csv("~/GitHub/stevia/data/homo.csv", stringsAsFactors = FALSE)
homo.hs <- subset(homo, taxID == 9606, select = c(HID, egID, symbol))
homo.ro <- subset(homo, taxID != 9606, select = c(HID, egID, symbol))

# make mapping homologene dfs
# load databases
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Rn.eg.db")
library("hom.Mm.inp.db")
library("hom.Rn.inp.db")


# make list of dataframes mapping affy probe to entrez id
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
  jpmPosEz.df[[n]] <- probe2Ez.df(row.names(jpmPos[[n]]), get(arrays[n, "bioc_package"]))
  jpmNegEz.df[[n]] <- probe2Ez.df(row.names(jpmNeg[[n]]), get(arrays[n, "bioc_package"]))
  }
}
jpmPosEz.df <- lapply(jpmPosEz.df, function(x) x[!(is.na(x[, 2]) & is.na(x[, 3])),])
jpmNegEz.df <- lapply(jpmNegEz.df, function(x) x[!(is.na(x[, 2]) & is.na(x[, 3])),])

# add columns with human homologue ed & symbol
Ez.df2homo <- function(Ez.df) {
  homo.hs[homo.hs$HID == homo.ro$HID[homo.ro$egID == Ez.df], 2:3]
}
for(n in names(jpmPos)) {
  if(!is.null(jpmPos[[n]]))  {
  jpmPosEz.df[[n]] <- cbind(jpmPosEz.df[[n]], t(sapply(jpmPosEz.df[[n]]$ENTREz.dfID, Ez.df2homo)))
  jpmNegEz.df[[n]] <- cbind(jpmNegEz.df[[n]], t(sapply(jpmNegEz.df[[n]]$ENTREz.dfID, Ez.df2homo)))
  }
}

# copy human symbols & entrez ids to homo columns for symetry
for(n in names(human)) {
  jpmPosEz.df[[n]]$egID <- jpmPosEz.df[[n]]$ENTREZID
  jpmPosEz.df[[n]]$symbol <- jpmPosEz.df[[n]]$SYMBOL
  jpmNegEz.df[[n]]$egID <- jpmNegEz.df[[n]]$ENTREZID
  jpmNegEz.df[[n]]$symbol <- jpmNegEz.df[[n]]$SYMBOL
}

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

# apply row2nrows to whole dataframe
dflist2nrows <- function(df) {
  out <- data.frame(PROBEID = character(0), SYMBOL = character(0), ENTREZID = character(0), egID = character(0), symbol = character(0))
  for(i in seq_len(dim(df)[1])) {
    out <- rbind(out, row2nrows(df[i,]))
  }
  return(out)
}

jpmPosEz.homo <- lapply(jpmPosEz.df, dflist2nrows)
jpmNegEz.homo <- lapply(jpmNegEz.df, dflist2nrows)



# compare inparanoid and homologene maps
# homologene map
jpmPosEz.hs <- unique(unlist(sapply(jpmPosEz[1:6], function(x) x$ENTREZID)))
length(jpmPosEz.hs)
# [1] 3373

jpmPosEz.nhs <- unique(unlist(sapply(jpmPosEz[7:26], function(x) x$egID)))
length(jpmPosEz.nhs)
# [1] 5067

jpmPosEz <- unique(c(jpmPosEz.hs, jpmPosEz.nhs))
length(jpmPosEz)
# [1] 7291

jpmNegEz.hs <- unique(unlist(sapply(jpmNegEz[1:6], function(x) x$ENTREZID)))
length(jpmNegEz.hs)
# [1] 3842

jpmNegEz.nhs <- unique(unlist(sapply(jpmNegEz[7:26], function(x) x$egID)))
length(jpmNegEz.nhs)
# [1] 5118

jpmNegEz <- unique(c(jpmNegEz.hs, jpmNegEz.nhs))
length(jpmNegEz)
# [1] 7630

length(intersect(jpmPosEz, jpmNegEz))
# [1] 3650

jpmEz <- select(org.Hs.eg.db, unique(c(jpmPosEz, jpmNegEz)), columns = c("ENTREZID", "SYMBOL"))
dim(jpmEz)
# [1] 11271     2

# inparaniod map size
inparPos <- character(0)
lapply(jpmPos.df, function(x) inparPos <<- c(inparPos, x$homSym))
inparNeg <- character(0)
lapply(jpmNeg.df, function(x) inparNeg <<- c(inparNeg, x$homSym))
inpar <- unique(c(inparPos, inparNeg))

length(inpar)
# [1] 5633

length(unique(c(jpmPosEz.nhs, jpmNegEz.nhs)))
# [1] 7985

### homo map wins!!

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

sapply(jpmPosEz.homo, function (x) length(x$PROBEID))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#          817         1168          611          729          526          372            7 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#          245          136          737          525          620          287          792 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#          509          820          530         1277          376          366          469 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#          341          233          258          302          378 
sapply(jpmPosEz.homo, function (x) sum(duplicated(x$PROBEID)))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#           79          130           34           71           27           53            0 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#           57           26           71           60          110           43           74 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#           33           73           58          142           28           26           34 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#           54           24           44           36           91 

sapply(jpmNegEz.homo, function (x) length(x$PROBEID))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#         1111         1100          745          990          507          415           10 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#          205          210          573          498          603          232          561 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#          480          712          588          924          351          613          470 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#          178          171          523          194          325 
sapply(jpmNegEz.homo, function (x) sum(duplicated(x$PROBEID)))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#           73           60           53           88           30           38            0 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#           18           44           40           26           42           18           65 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#           33           66           57           58           33           48           57 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#           21           17           67           23           27 


# count mouse & rat homologues
# rat map size
ratPos <- character(4)
lapply(jpmPosEz.df[names(rat)], function(x) ratPos <<- rbind(ratPos, x[, 2:5]))
ratNeg <- character(4)
lapply(jpmNegEz.df[names(rat)], function(x) ratNeg <<- rbind(ratNeg, x[, 2:5]))
homoRat <- unique(rbind(ratPos, ratNeg))

# non-uniqueness of map
summary(as.factor(sapply(homoRat$egID, length)))
#    0    1    2    3    4 
#  806 4382   18    2    1 

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
