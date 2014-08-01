##### recreate jpm paper analysis & compare with moses results
# uses gdsExp, probe2gene from 'data cleaning.R'
### jpm individual data set analysis recreation: linear model of individual gene expression with age

jpmExp <- lapply(gdsExp, impute.matrix)

# how many probe sets were eliminated because missing > 30%?
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

# NOTE: FIXED -- stromal & spinal_cord don't work - for some values model doesn't converge.  fixed with "tryCatch" in function "row.slope"
jpmExp.pb <- vector("list", 26)
names(jpmExp.pb) <- names(jpmExp)
for(i in 1:26) {
  try(jpmExp.pb[[i]] <- matrix.slope(ages[[i]], jpmExp[[i]]))
}

#get f statistics and p values
jpmExp.f <- vector("list", 26)
names(jpmExp.f) <- names(jpmExp)
for(i in 1:26) {
  try(jpmExp.f[[i]] <- poff(ages[[i]], jpmExp[[i]]))
}

jpmExp.fp <- lapply(jpmExp.f, function (x) try(apply(x, 1, function (y) pf(last.row(y, 3)[1], last.row(y, 3)[2], last.row(y, 3)[3], lower.tail = FALSE))))

# average of fraction of significant probes closely matches jpm result
 mean(unlist(lapply(jpmExp.fp, function (x) (sum(x <= .05, na.rm = TRUE) / length(na.omit(x))))))
# [1] 0.09261139

# TODO use code from lm.fit and summary.lm to make efficient probe/row based model fitting
# see http://reliawiki.org/index.php/Simple_Linear_Regression_Analysis

# pick out genes with significant linear realtionship to age
jpmExp.b1p <- lapply(jpmExp.pb, function(x) x[, last.ind(colnames(x), 2):last.ind(colnames(x))])

# fix stromal and spinal_cord
jpmExp.slope <- vector("list", length(jpmExp.b1p))
names(jpmExp.slope) <- names(jpmExp.b1p)
for(i in 1:26) {
  jpmExp.b1p[[i]] <- na.omit(jpmExp.b1p[[i]])
  jpmExp.slope[[i]] <- merge(as.data.frame(jpmExp.b1p[[i]]), as.data.frame(jpmExp.fp[[i]]), by = 0, all.x = TRUE)
  row.names(jpmExp.slope[[i]]) <- jpmExp.slope[[i]][, 1] 
  jpmExp.slope[[i]][, 1] <- NULL
  names(jpmExp.slope[[i]]) <- c("b1", "p.b1", "p.F")
}

save(list = c("jpmExp", "jpmExp.slope", "probe2gene"), file = "~/GitHub/stevia/data/jpmData.rdata")

jpmSig <- lapply(jpmExp.slope, function(x) x[x[, 3] <= .05,])

# add original annotation (untested after probe2gene modified to include current annotation)
for(i in 1:26) {
  jpmSig[[i]] <- cbind(jpmGene = probe2gene[[i]]$GEOgene[match(row.names(jpmSig[[i]]), probe2gene[[i]]$probe)], jpmSig[[i]], stringsAsFactors = FALSE)
}

# count probes measured
sapply(jpmExp.slope, function(x) dim(x)[1])
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#        12560        22215        22577        22215        22577        12533        22690 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#         4807         4043        45101        12488        12488        10043        45101 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#        45037        42869         9475        15923        15913        15922        15922 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#         8799         8799         8799         8799         6610 

# count probes significantly measured slopes
sapply(jpmSig, function(x) dim(x)[1])
#      h_brain      muscle1      muscle2      muscle3      muscle4      muscle5       muscle 
#         2271         2721         2035         2118         1550          720         3482 
#      kidney1      kidney2      m_brain      m_hippo        liver      m_heart         lung 
#          434          454         2319         1183         1476          492         2288 
#      cochlea  hemato_stem   myo_progen      r_hippo      stromal  spinal_cord   oculomotor 
#         1211         3329         1429         2866          985         1370         1260 
#  skeletal_ms   extraoc_ms laryngeal_ms      r_heart    CA1_hipp2 
#          490          378         1024          623          824 

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

# if < 30% na replace with row average, if > 30% missing delete row
impute.row <- function(row) {
  if(sum(is.na(row)) / length(row) > .3) return(NULL)
  row[is.na(row)] <- mean(row, na.rm = TRUE)
  return(row)
}

# apply impute.row to matrix
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
  tryCatch(summary(lm(y ~ x))$coefficients[2, c(1,4)], error = function (c) c(NA, NA))
}

# make nice matrix of values & results
matrix.slope <- function(y, x) {
  out <- t(apply(x, 1, function(r) row.slope(y, r)))
  out <- cbind(x, out)
  out <- rbind(c(y, Estimate = NA,'Pr(>|t|)' = NA), out)
  rownames(out) <- c("age", rownames(out)[-1])
  return(out)
}
  
# get p value of f statistic from row regression
row.poff <- function (y, x) {
  out <- try(summary(lm(y ~ x))$fstatistic)
  if(!is.numeric(out)) return(c(NA, NA, NA))
  return(out)
}

# column naming doesn't work
matrix.poff <- function(y, x) {
  out <- t(apply(x, 1, function(r) row.poff(y, r)))
  colnames(out) <- c("F", "df1", "df2")
  return(out)
}

  
# combine
poff <- function(y, x) {
  out <-matrix.poff(y, x)
  out <- cbind(x, out)
  out <- rbind(c(y, "F" = NA,"df1" = NA, "df2" = NA), out)
  rownames(out) <- c("age", rownames(out)[-1])
  return(out)
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
