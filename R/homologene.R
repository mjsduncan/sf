### make homologene homologue map
homologene <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/homologene/homologene.data", header=F, stringsAsFactors = FALSE)
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
