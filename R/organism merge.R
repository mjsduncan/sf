### organism and complete data set mergeing
load("~/GitHub/stevia/data/symEsets.rdata")
# get r_hippo controls
symEsets$r_hippo <- symEsets$r_hippo[, c(1:10, 30:39)]

arrays <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/arrays.tab")
oEsets <- list()

# list expression sets by organism
for(n in levels(arrays$organism)) {
  oEsets[[n]] <- symEsets[as.character(arrays[arrays$organism == n, "row.names"])]
}
lapply(oEsets, function (x) sapply(x, dim))
# $`Homo sapiens`
#          h_brain muscle5 muscle3 muscle1 muscle2 muscle4
# Features    8632    8616   12494   12494    9648    9648
# Samples       30      12      15      15      15      15
# 
# $`Mus musculus`
#          m_brain cochlea hemato_stem  lung m_heart muscle kidney1 kidney2 myo_progen liver m_hippo
# Features   13016   13016       13016 13016    8622     98    5106    3299       8724  8724    8724
# Samples        6       6           8    15      12     10      10      10          4     7      23
# 
# $`Rattus norvegicus`
#          stromal spinal_cord r_hippo laryngeal_ms skeletal_ms r_heart CA1_hipp2 oculomotor extraoc_ms
# Features   10000       10000   10000         4750        4750    4750      4750      10000       4750
# Samples        3           9      20            9          12      11        29          9         12

# merge organism expression matrices by row
col2rownames <- function(df, colname = "Row.names", removecol = TRUE){
  row.names(df) <- df[,colname]
  if(removecol){df[,colname] <- NULL}
  return(df)
  }

mergeByRows <- function(x, y) col2rownames(merge(x, y, by = "row.names"))
orgExprs <- lapply(oEsets, function(x) Reduce(mergeByRows, lapply(x, exprs)))

sapply(orgExprs, dim)
#      Homo sapiens Mus musculus Rattus norvegicus
# [1,]         2463           18              4488
# [2,]          102          111               114

# throw out some mouse data sets to improve gene coverage
mouse2 <- Reduce(mergeByRows, lapply(oEsets[[2]][-6], exprs))
dim(mouse2)
# [1] 992 101

mouse3 <- Reduce(mergeByRows, lapply(oEsets[[2]][c(1:5, 7, 9:11)], exprs))
dim(mouse3)
# [1] 4733   91

mouse4 <- Reduce(mergeByRows, lapply(oEsets[[2]][c(1:5, 9:11)], exprs))
dim(mouse4)
# [1] 8319   81
# throwing out $muscle increases genes x50, $kidney2 x5, $kidney1 x2, with a sample drop of average 10% each time
orgExprs[[2]] <- mouse4

# add age catagory & make moses input files
orgAges <- lapply(oEsets, function(x) unlist(lapply(x, function(y) y$age)))
for(n in names(orgAges)) names(orgAges[[n]]) <- lapply(oEsets, function(x) unlist(lapply(x, sampleNames)))[[n]]
orgAges[[2]] <- orgAges[[2]][!(names(orgAges[[2]]) %in% unlist(lapply(oEsets[[2]][6:8], sampleNames)))]
identical(lapply(orgExprs, colnames), lapply(orgAges, names))
# [1] TRUE

# apply median norm to matrix columns
med.normalize <- function(mat) {
  out <- mat
  for (i in seq(dim(mat)[2])) { 
    vect <- mat[,i]
    med <- median(vect, na.rm = TRUE)
    out[,i] <- as.numeric(vect >= med)
  }
  return(out)
}

orgMoses <- list()
for(n in names(orgAges)) {
  orgMoses[[n]] <- t(med.normalize(orgExprs[[n]]))
  orgMoses[[n]] <- cbind(orgAges[[n]], orgMoses[[n]])
  colnames(orgMoses[[n]])[1] <- "age"
}

# TODO:  fix feature names

save(orgMoses, file = "~/GitHub/stevia/data/orgMoses.rdata")

