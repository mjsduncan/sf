##### recreate jpm paper analysis & compare with moses results
# uses gdsExp from 'data cleaning.R'
### jpm individual data set analysis recreation: linear model of individual gene expression with age

jpmExp <- lapply(gdsExp, impute.matrix)

# how many probe sets were eliminated?
unlist(lapply(gdsExp, function(x) dim(x)[1])) - unlist(lapply(jpmExp, function(x) dim(x)[1]))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle      kidney1 
#           65           68           68           68           68           93            0         1777 
#      kidney2      m_brain      m_hippo        liver      m_heart         lung      cochlea  hemato_stem 
#         2552            0            0            0         2611            0           64         2232 
#   myo_progen      r_hippo      stromal  spinal_cord   oculomotor  skeletal_ms   extraoc_ms laryngeal_ms 
#         3013            0            3            0            1            0            0            0 
#      r_heart    CA1_hipp2 
#            0         2189 

# check highest expression level values to insure log tranform not already aplied
unlist(lapply(jpmExp, max, na.rm = TRUE))
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle      kidney1 
#      10663.8      76224.0     150277.0      92814.0     171668.0      29071.1      43396.5     264948.0 
#      kidney2      m_brain      m_hippo        liver      m_heart         lung      cochlea  hemato_stem 
#     349472.0      22535.3      17548.9      48075.2      65959.1      30525.7     755370.0      10084.8 
#   myo_progen      r_hippo      stromal  spinal_cord   oculomotor  skeletal_ms   extraoc_ms laryngeal_ms 
#     160387.0      15736.0       6426.9      23463.4      20590.4      63871.6      46955.9      48487.6 
#      r_heart    CA1_hipp2 
#     307756.0      70316.6 

jpmExp <- lapply(jpmExp, log2)

# compute slopes and p values
jpmExp <- mapply(matrix.slope, ages, jpmExp)

#stromal & spinal_cord don't work
jpmExp.pb <- vector("list", 26)
names(jpmExp.pb) <- names(jpmExp)
for(i in 1:26) {
  try(jpmExp.pb[[i]] <- matrix.slope(ages[[i]], jpmExp[[i]]))
}

#get f statistics and p values
jpmExp.f <- vector("list", 26)
names(jpmExp.f) <- names(jpmExp)
for(i in 1:26) {
  try(jpmExp.f[[i]] <- matrix.poff(ages[[i]], jpmExp[[i]]))
}

save(list = c("jpmExp", "jpmExp.f", "jpmExp.pb"), file = "jpmData.rdata")

# fix matrix column names
for(i in 1:26) {
  if(!is.null(jpmExp.f[[i]])) {
    colnames(jpmExp.f[[i]])[last.ind(colnames(jpmExp.f[[i]]), 3)] <- "F"
    colnames(jpmExp.f[[i]])[last.ind(colnames(jpmExp.f[[i]]), 2)] <- "df1"
    colnames(jpmExp.f[[i]])[last.ind(colnames(jpmExp.f[[i]]))] <- "df2"
  }
}

jpmExp.fp <- lapply(jpmExp.f, function (x) try(apply(x, 1, function (y) pf(last.row(y, 3)[1], last.row(y, 3)[2], last.row(y, 3)[3]))))

# average of fraction of significant probes closely matches jpm result
 mean(unlist(lapply(jpmExp.fp[-c(19, 20)], function (x) (sum(x <= .025 | x >= .975, na.rm = TRUE) / length(x)))))
# [1] 0.07845494

# TODO use code from lm.fit and summary.lm to make efficient probe/row based model fitting
# see http://reliawiki.org/index.php/Simple_Linear_Regression_Analysis

# pick out genes with significant linear realtionship to age
jpmExp.pb[19:20] <- NA
jpmExp.b1p <- lapply(jpmExp.pb, function(x) if(!is.na(x)) x[, last.ind(colnames(x), 2):last.ind(colnames(x))])
jpmExp.slope <- mapply(cbind, jpmExp.b1p, jpmExp.fp)
for(i in 1:26) if(is.numeric(jpmExp.slope[[i]])) colnames(jpmExp.slope[[i]]) <- c("b1", "p of b1", "p of F")

# why don't F test signifcant probes have significant slopes? 
# $skeletal_ms
# age                    NA           NA        NA
# A01157cds_s_at 0.01849251  -0.06845535 0.9815075
#
# $m_heart                                         
# age              NA        NA         NA
# 100001_at 0.9712485  2.507187 0.02875146
#
# they are significant for 2 tailed tests!

jpmSig <- lapply(jpmExp.slope, function(x) if(is.numeric(x)) x[x[, 3] <= .025 | x[, 3] >= .975,])

# add original annotation
jpmSig <- lapply(jpmSig, function(x) x[-1,])
for(i in 1:26) {
  if(!is.null(jpmSig[[i]])) {
    jpmSig[[i]] <- cbind(jpmGene = probe2gene[[i]][match(rownames(jpmSig[[i]]), probe2gene[[i]][, 1]), 2], 
      as.data.frame(jpmSig[[i]]), stringsAsFactors = FALSE)
  }
}

# add homologue mapping
# add current annotation

# split into positively  and negatively correlated
jpmPos <- vector("list", 26)
names(jpmPos) <- names(jpmSig)
jpmNeg <- vector("list", 26)
names(jpmNeg) <- names(jpmSig)
for(n in names(jpmSig)) {
  if(!is.null(jpmSig[[n]]))  {
    jpmPos[[n]] <- jpmSig[[n]][jpmSig[[n]][, 2] > 0,]
    jpmNeg[[n]] <- jpmSig[[n]][jpmSig[[n]][, 2] < 0,]
  }
}

### functions

# if < 30% na replace with row average
impute.row <- function(row) {
  if(sum(is.na(row)) / length(row) > .3) return(NULL)
  row[is.na(row)] <- mean(row, na.rm = TRUE)
  return(row)
}

impute.matrix <- function(mat) {
  out <- apply(mat, 1, impute.row)
  out[sapply(out, is.null)] <- NULL
  if(class(out) != "list") {
    print(paste(deparse(substitute(mat)), "generated an error."))
    return(t(out))
  }
  do.call(rbind, out)
}

# get slope and p value from row regression
row.slope <- function (y, x) {
  summary(lm(y ~ x))$coefficients[2, c(1,4)]
}

# column naming doesn't work
matrix.slope <- function(y, x) {
  out <- t(apply(x, 1, function(r) row.slope(y, r)))
  out <- `colnames<-`(cbind(x, out), c(colnames(x), c("slope", "p")))
  `rownames<-`(rbind(c(y, NA, NA), out), c("age", rownames(out)))
}
  
# get p value of f statistic from row regression
row.poff <- function (y, x) {
  summary(lm(y ~ x))$fstatistic
}

# column naming doesn't work
matrix.poff <- function(y, x) {
  out <- t(apply(x, 1, function(r) row.poff(y, r)))
  out <- cbind(x, out)
  out <- rbind(c(y, NA, NA, NA), out)
  colnames(out) <- , c("age", rownames(out))
}
  
# fix failure to name columns in matrix.poff

last.row <- function (vec, n = 1) vec[(length(vec) - n + 1):length(vec)]
last.ind <- function (vec, n = 1) length(vec) - n + 1

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
  