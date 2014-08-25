# program used to generate moses output  
### moses on merged dataset 
moses.flags <- "-j3 -W1 -u age --hc-crossover=1 -m 1000000 --balance=1 --enable-fs=1 -r1"
for(i in 1:10) {
  make.dir(getwd(),  paste("meta_run", i, sep = ""))
  dataset <- "mouses"
  cr <- mean( Pgs_mouseBack_normalized[, "age"])
  make.dir(getwd(), dataset)
  mouses_1.test <- makeMpartitions( Pgs_mouseBack_normalized, p = .5)
  mouses_1_hx5.train <- runMfolder(moses.flags)
  save(mouses_1.test, mouses_1_hx5.train, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  #   mouses_1_hx5 <- testClist(mouses_1_hx5.train, mouses_1.test, cr)
  #   mouses_1_hx5.best <- bestCombos(mouses_1_hx5)
  #   mouses_1_hx5.bestN5 <- bestCombos(mouses_1_hx5, N = 5)
  #   save(mouses_1_hx5.best, mouses_1_hx5.bestN5, file = paste("mouses_", i, "_best.rdata", sep = ""))
  setwd("..")
  
  dataset <- "rat"
  cr <- mean(Pgs_ratBack_normalized[, "age"])
  make.dir(getwd(), dataset)
  rat_1.test <- makeMpartitions(Pgs_ratBack_normalized, p = .5)
  rat_1_hx5.train <- runMfolder(moses.flags)
  save(rat_1.test, rat_1_hx5.train, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  #   rat_1_hx5 <- testClist(rat_1_hx5.train, rat_1.test, cr)
  #   rat_1_hx5.best <- bestCombos(rat_1_hx5)
  #   rat_1_hx5.bestN5 <- bestCombos(rat_1_hx5, N = 5)
  #   save(rat_1_hx5.best, rat_1_hx5.bestN5, file = paste("rat_", i, "_best.rdata", sep = ""))
  setwd("..")
  
  dataset <- " human"
  cr <- mean( Pgs_homoBack_normalized[, "age"])
  make.dir(getwd(), dataset)
  human_1.test <- makeMpartitions(Pgs_homoBack_normalized, p = .5)
  human_1_hx5.train <- runMfolder(moses.flags)
  save(human_1.test, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  
  #   human_1_hx5 <- testClist(human_1_hx5.train, human_1.test, cr)
  #   human_1_hx5.best <- bestCombos(human_1_hx5)
  #   human_1_hx5.bestN5 <- bestCombos(human_1_hx5, N = 5)
  #   save(human_1_hx5.best, human_1_hx5.bestN5, file = paste("human_", i, "_best.rdata", sep = ""))
  setwd("../..")
  
  #   save(mouses__1_hx5.best, mouses_1_hx5.bestN5, 
  #        rat_1_hx5.best, rat_1_hx5.bestN5, 
  #        human_1_hx5.bestN5, human_1_hx5.bestN5, 
  #        file = paste("bestOfRun_", i, ".rdata"))
}


###### edited make directory function 
make.dir <- function(mainDir =getwd(), subDir) {
  if (file.exists(paste(mainDir, subDir, "/", sep = "/", collapse = "/"))) {
    cat("subDir exists in mainDir and is a directory")
  } else if (file.exists(paste(mainDir, subDir, sep = "/", collapse = "/"))) {
    cat("subDir exists in mainDir but is a file")
    # you will probably want to handle this separately
  } else {
    cat("subDir does not exist in mainDir - creating")
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
  }  
  if (file.exists(paste(mainDir, subDir, "/", sep = "/", collapse = "/"))) {
    # By this point, the directory either existed or has been successfully created
    setwd(file.path(mainDir, subDir))
  } else {
    cat("subDir does not exist")
    # Handle this error as appropriate
  }
}

# run combos on testing data sets. 
# xxx_run1_moses.RData & xxx_run2_moses.RData contain identically named objects!
load("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run2/ human_run2_moses.Rdata")
human_2.test <- human_1.test
human_2_hx5.train <- human_1_hx5.train
load("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run2/mouses_run2_moses.Rdata")
mouses_2.test <- mouses_1.test
mouses_2_hx5.train <- mouses_1_hx5.train
load("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run2/rat_run2_moses.Rdata")
rat_2.test <- rat_1.test
rat_2_hx5.train <- rat_1_hx5.train
load("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run1/ human_run1_moses.Rdata")
load("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run1/mouses_run1_moses.Rdata")
load("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run1/rat_run1_moses.Rdata")

# load a training data set from each run for each organism to reconstruct complete sample age vector for each organism dataset
Pgs_homoBack_normalized_f11 <- read.csv("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run1/ human/Pgs_homoBack_normalized_f1.csv")
Pgs_homoBack_normalized_f12 <- read.csv("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run2/ human/Pgs_homoBack_normalized_f1.csv")

mouse_woutBack_normalized_f11 <- read.csv("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run1/mouses/mouse_woutBack_normalized_f1.csv")
mouse_woutBack_normalized_f12 <- read.csv("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run2/mouses/mouse_woutBack_normalized_f1.csv")


Pgs_ratBack_normalized_f11 <- read.csv("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run1/rat/Pgs_ratBack_normalized_f1.csv")
Pgs_ratBack_normalized_f12 <- read.csv("C:/Users/user/Desktop/biomind/addis/moses_R_result .tar/moses_R_result/meta_run2/rat/Pgs_ratBack_normalized_f1.csv")

# combine age vectors from testing and training sets by organism
ls(pattern = "_f11")
# [1] "mouse_woutBack_normalized_f11" "Pgs_homoBack_normalized_f11"   "Pgs_ratBack_normalized_f11"   
ls(pattern = "_1.test")[c(2, 1, 3)]
# [1] "mouses_1.test" "human_1.test"  "rat_1.test"   

age.train <- vector("list", 3)
names(age.train) <- c("mouse", "human", "rat")
for(i in 1:3) age.train[[i]] <- c(get(ls(pattern = "_f11")[i])$age, get(ls(pattern = "_1.test")[c(2, 1, 3)][i])[[1]][, "age"])

sapply(age.train, length)
# mouse human   rat 
#    81   102   172 

# get fraction of "1"s for balanced accuracy score
age.cr <- lapply(age.train, mean)
unlist(age.cr)
#     mouse     human       rat 
# 0.6049383 0.5882353 0.6395349 

# check with run2 training sets
age.train2 <- vector("list", 3)
names(age.train2) <- c("mouse", "human", "rat")
for(i in 1:3) age.train2[[i]] <- c(get(ls(pattern = "_f12")[i])$age, get(ls(pattern = "_2.test")[c(2, 1, 3)][i])[[1]][, "age"])

sapply(age.train2, length)
# mouse human   rat 
#    81   102   172 

# get fraction of "1"s for balanced accuracy score
age.cr2 <- lapply(age.train2, mean)
unlist(age.cr2)
#     mouse     human       rat 
# 0.6049383 0.5882353 0.6395349 

rm(list = ls(pattern = "_f1"))

# test combos
library("caret")

# 1810058I24Rik and 9 similiar in run1 f10, Krtap6-2 in f5, 1110057K04Rik in f7
mouses_1_hx5 <- testClist(mouses_1_hx5.train[c(2:5, 7, 9:10)], lapply(mouses_1.test[c(2:5, 7, 9:10)], as.data.frame), 0.6049383)
mouses_1_hx5.best <- bestCombos(mouses_1_hx5, metric = 13)
mouses_1_hx5.bestN5 <- bestCombos(mouses_1_hx5, N = 5, metric = 13)
save(mouses_1_hx5.best, mouses_1_hx5.bestN5, file = "mouses_1_best.rdata")

# Nkx2-3 in run2 f1, 4933415E08Rik in f6
mouses_2_hx5 <- testClist(mouses_2_hx5.train[c(1, 3:6, 8:10)], lapply(mouses_2.test[c(1, 3:6, 8:10)], as.data.frame), 0.6049383)
mouses_2_hx5.best <- bestCombos(mouses_2_hx5, metric = 13)
mouses_2_hx5.bestN5 <- bestCombos(mouses_2_hx5, N = 5, metric = 13)
save(mouses_2_hx5.best, mouses_2_hx5.bestN5, file = "mouses_2_best.rdata")

# Rpo1-3 in run1 f4
rat_1_hx5 <- testClist(rat_1_hx5.train[c(1:4, 6:10)], lapply(rat_1.test[c(1:4, 6:10)], as.data.frame), 0.6395349)
rat_1_hx5.best <- bestCombos(rat_1_hx5, metric = 13)
rat_1_hx5.bestN5 <- bestCombos(rat_1_hx5, N = 5, metric = 13)
save(rat_1_hx5.best, rat_1_hx5.bestN5, file = "rat_1_best.rdata")

rat_2_hx5 <- testClist(rat_2_hx5.train, lapply(rat_2.test, as.data.frame), 0.6395349)
rat_2_hx5.best <- bestCombos(rat_2_hx5, metric = 13)
rat_2_hx5.bestN5 <- bestCombos(rat_2_hx5, N = 5, metric = 13)
save(rat_2_hx5.best, rat_2_hx5.bestN5, file = "rat_2_best.rdata")

# WT1-AS in run1 f2, WT1-AS & NKX2-1 in run1 f7
human_1_hx5 <- testClist(human_1_hx5.train[c(1:2, 4:7, 9:10)], lapply(human_1.test[c(1:2, 4:7, 9:10)], as.data.frame), 0.5882353)
human_1_hx5.best <- bestCombos(human_1_hx5, metric = 13)
human_1_hx5.bestN5 <- bestCombos(human_1_hx5, N = 5, metric = 13)
save(human_1_hx5.best, human_1_hx5.bestN5, file = "human_1_best.rdata")

# NALCN-AS1 in run2 f10
human_2_hx5 <- testClist(human_2_hx5.train[c(2:10)], lapply(human_2.test[c(2:10)], as.data.frame), 0.5882353)
human_2_hx5.best <- bestCombos(human_2_hx5, metric = 13)
human_2_hx5.bestN5 <- bestCombos(human_2_hx5, N = 5, metric = 13)
save(human_2_hx5.best, human_2_hx5.bestN5, file = "human_2_best.rdata")

# bind together results function from "mash runs.R" in Rmoses:
lbind <- function(list1, list2) {
  if(sum(sapply(list1, class) != sapply(list2, class)) > 0) return("lists element classes do not match")
  out <- vector("list", length(list1))
  names(out) <- names(list1)
  for(i in seq_along(out)) {
    ifelse(is.atomic(list1[[i]]), out[[i]] <- c(list1[[i]], list2[[i]]), out[[i]] <- rbind(list1[[i]], list2[[i]]))
  }
    return(out)
}

best2human <- lbind(human_1_hx5.best, human_2_hx5.best)
best2mouse <- lbind(mouses_1_hx5.best, mouses_2_hx5.best)
best2rat <- lbind(rat_1_hx5.best, rat_2_hx5.best)


