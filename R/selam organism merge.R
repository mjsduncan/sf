# 2014_08_08
# PGSEA::aggregateExprs -- package 
# gdsData.Rdata 
# library("impute")
# library("PGSEA")

Rattus_norvegicus <-row.names(arrays)[arrays$organism == "Rattus norvegicus"]
Mus_musculus      <- row.names(arrays)[arrays$organism == "Mus musculus"]
Homo_sapiens      <- row.names(arrays)[arrays$organism == "Homo sapiens"]
# Homo_sapiens
# [1] "h_brain" "muscle5" "muscle3" "muscle1" "muscle2" "muscle4"

pgsea_data <- gdsExp
pgsea_data <- lapply(pgsea_data, log2)

#### Removing NAs from each datasets 
for (n in names(pgsea_data)){
pgsea_data[[n]] <- (impute::impute.knn(pgsea_data[[n]]))$data
   }
## removes duplicates row names using  "SYMBOL"
for (n in names(pgsea_data)){
pgsea_data[[n]] <-  PGSEA::aggregateExprs( pgsea_data[[n]], package = arrays$bioc_package [row.names(arrays) == n] , using = "SYMBOL", FUN= median )
  }

#### fixing ages 
ages_norm <- age.pdata
ages_norm$h_brain[ c(5:9, 23:24), ] <- 1 
ages_norm$lung[ , ] <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1 ,1)
ages_norm$m_brain[1:3, ] <- 0
ages_norm$liver[1:4, ] <- 0
ages_norm$hemato_stem[, 1] <- c(0, 0, 0, 1, 1, 1, 1, 1)
ages_norm$stromal[1:2, ] <- 0
ages_norm$r_heart[1:5, ] <- 0

#### binding age by organism # 
homo_age <- ages_norm[[Homo_sapiens[1]]]
for(i in 2: length(Homo_sapiens)){
   homo_age <- rbind (homo_age, ages_norm[[Homo_sapiens[i]]])}

mus_age <- ages_norm[[Mus_musculus[1]]]
for(i in 2: length(Mus_musculus)){
  mus_age <- rbind (mus_age, ages_norm[[Mus_musculus[i]]])}

rat_age <- ages_norm[[Rattus_norvegicus[1]]]
for(i in 2: length(Rattus_norvegicus)){
  rat_age <- rbind (rat_age, ages_norm[[Rattus_norvegicus[i]]])}
#### merging by organism 
 ##### human #########
Pgs_homo <- pgsea_data[[Homo_sapiens[1]]]
for(i in 2: length(Homo_sapiens)){
  Pgs_homo <- merge(Pgs_homo, pgsea_data[[Homo_sapiens[i]]], by = "row.names", all = TRUE)
  row.names(Pgs_homo) <- Pgs_homo$Row.names
  Pgs_homo$Row.names <- NULL
}
all_H <- Pgs_homo
Pgs_homo <- as.matrix(rbind(age = t(homo_age), Pgs_homo))
dimnames( Pgs_homo)[[1]][1] <- "age"
Pgs_homoBack <- impute.matrix(Pgs_homo)
Pgs_homoBack_normalized   <- med.normalize(t(Pgs_homoBack))

 ###### mouse########
Pgs_mouse <- pgsea_data[[Mus_musculus[1]]]
for(i in 2: length(Mus_musculus)){
  Pgs_mouse <- merge(Pgs_mouse, pgsea_data[[Mus_musculus[i]]], by = "row.names", all = TRUE)
  row.names(Pgs_mouse) <- Pgs_mouse$Row.names
  Pgs_mouse$Row.names <- NULL
}
all_M <- Pgs_mouse
Pgs_mouse <- as.matrix(rbind(age= t(mus_age), Pgs_mouse))
dimnames( Pgs_mouse)[[1]][1] <- "age"
Pgs_mouseBack <- impute.matrix(Pgs_mouse)
Pgs_mouseBack_normalized   <- med.normalize(t(Pgs_mouseBack))

 #### rat############
Pgs_rat <- pgsea_data[[Rattus_norvegicus[1]]]
for(i in 2: length(Rattus_norvegicus)){
  Pgs_rat <- merge(Pgs_rat, pgsea_data[[Rattus_norvegicus[i]]], by = "row.names", all = TRUE)
  row.names(Pgs_rat) <- Pgs_rat$Row.names
  Pgs_rat$Row.names <- NULL
}
all_R <- Pgs_rat 
Pgs_rat <- as.matrix(rbind(age= t(rat_age), Pgs_rat))
dimnames( Pgs_rat)[[1]][1] <- "age"
Pgs_ratBack <- impute.matrix(Pgs_rat)
Pgs_ratBack_normalized   <- med.normalize(t(Pgs_ratBack))

#### all organisms together 

all_organism  <- merge(impute.matrix(all_H), impute.matrix(all_M) , by = "row.names", all = TRUE)
row.names(all_organism) <- all_organism$Row.names
all_organism$Row.names <- NULL
all_organism  <- merge(all_organism ,impute.matrix(all_R),  by = "row.names", all = TRUE)
row.names(all_organism) <- all_organism$Row.names
all_organism$Row.names <- NULL

all_age <- rbind(homo_age , mus_age, rat_age)
all_organism <- as.matrix(rbind(age= t(all_age), all_organism))
dimnames( all_organism)[[1]][1] <- "age"
all_organismBack <- impute.matrix(all_organism)
all_organismNorm <- med.normalize( t (all_organismBack))

######  tests

sapply(pgsea_data[Homo_sapiens], dim)
#      h_brain muscle5 muscle3 muscle1 muscle2 muscle4
# [1,]    8632    8616   12494   12494    9648    9648
# [2,]      30      12      15      15      15      15
sum(sapply(pgsea_data[Homo_sapiens], function(x) dim(x)[2]))
# [1] 102

sapply(pgsea_data[Mus_musculus], dim)
#      m_brain cochlea hemato_stem  lung m_heart muscle kidney1 kidney2 myo_progen liver m_hippo
# [1,]   13016   13016       13016 13016    8622     98    5106    3299       8724  8724    8724
# [2,]       6       6           8    15      12     10      10      10          4     7      23
sum(sapply(pgsea_data[Mus_musculus], function(x) dim(x)[2]))
# [1] 111

sapply(pgsea_data[Rattus_norvegicus], dim)
#      stromal spinal_cord r_hippo laryngeal_ms skeletal_ms r_heart CA1_hipp2 oculomotor extraoc_ms
# [1,]   10000       10000   10000         4750        4750    4750      4750      10000       4750
# [2,]       3           9      78            9          12      11        29          9         12
sum(sapply(pgsea_data[Rattus_norvegicus], function(x) dim(x)[2]))
# [1] 172

sapply(ls(pattern = "Pgs_"), function(x) dim(get(x)))
#      Pgs_homo Pgs_homoBack  ... Pgs_mouse Pgs_mouseBack ... Pgs_rat Pgs_ratBack ...
# [1,]    17985         8561  ...     13558          8421 ...   10263        4489 ...
# [2,]      102          102  ...       111           111 ...     172         172 ...




sapply(ls(pattern = "Back_normalized"), function(x) summary(colMeans(get(x))))
#         Pgs_homoBack_normalized Pgs_mouseBack_normalized Pgs_ratBack_normalized
# Min.                     0.5000                   0.5045                 0.5000
# 1st Qu.                  0.5000                   0.5045                 0.5000
# Median                   0.6275                   0.5495                 0.5000
# Mean                     0.6118                   0.5684                 0.5004
# 3rd Qu.                  0.6863                   0.6126                 0.5000
# Max.                     0.7843                   0.7838                 0.6395

sapply(ls(pattern = "all_org"), function(x) dim(get(x)))
#      all_organism all_organismBack all_organismNorm
# [1,]        18452             3007              385
# [2,]          385              385             3007
summary(colMeans(all_organismNorm))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.5013  0.5753  0.6286  0.6277  0.6805  0.7636 

# 