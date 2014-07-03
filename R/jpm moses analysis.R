
# script variables
mainD <- "~/projects/sf_results/test_run1"
# moses.data <- list()
moses.flags <- "-j2 -W1 -u age --hc-crossover=1 -m 1000000 --balance=1 --enable-fs=1 -r1"

# main loop
for(i in 1:2) {
  make.dir(mainD, paste("meta_run", i, sep = ""))

  # # impute NAs and check for non-informative features
  # moses.data$gpl81 <- bin.impute.matrix(gpl81)
  # summary(colSums(moses.data$gpl81))
  # #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # #   16.00   18.00   18.00   18.07   18.00   22.00 
  # dim(moses.data$gpl81)
  # # [1]    34 12489
  
  ### dataset pgl81
  
  # gpl81 <- data.frame(moses.data$gpl81, check.names = TRUE)
  dataset <- "gpl81"
  cr <- mean(gpl81[, "age"])
  make.dir(getwd(), dataset)
  gpl81_1.test <- makeMpartitions(gpl81, p = .5)
  gpl81_1_hx5.train <- runMfolder(moses.flags)
  save(gpl81_1.test, gpl81_1_hx5.train, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  gpl81_1_hx5 <- testClist(gpl81_1_hx5.train, gpl81_1.test, cr)
  gpl81_1_hx5.best <- bestCombos(gpl81_1_hx5)
  gpl81_1_hx5.bestN5 <- bestCombos(gpl81_1_hx5, N = 5)
  save(gpl81_1_hx5.best, gpl81_1_hx5.bestN5, file = paste("gpl81_", i, "_best.rdata", sep = ""))
  setwd("..")
  
  # # impute NAs and check for non-informative features
  # moses.data$gpl85 <- bin.impute.matrix(gpl85)
  # summary(colSums(moses.data$gpl85))
  # #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # #   38.00   38.00   38.00   38.17   38.00   48.00 
  # dim(moses.data$gpl85)
  # # [1]   73 6611
  
  ### dataset gpl85
  
  # gpl85 <- data.frame(moses.data$gpl85, check.names = TRUE)
  dataset <- "gpl85"
  cr <- mean(gpl85[, "age"])
  make.dir(getwd(), dataset)
  gpl85_1.test <- makeMpartitions(gpl85, p = .5)
  gpl85_1_hx5.train <- runMfolder(moses.flags)
  save(gpl85_1.test, gpl85_1_hx5.train, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  gpl85_1_hx5 <- testClist(gpl85_1_hx5.train, gpl85_1.test, cr)
  gpl85_1_hx5.best <- bestCombos(gpl85_1_hx5)
  gpl85_1_hx5.bestN5 <- bestCombos(gpl85_1_hx5, N = 5)
  save(gpl85_1_hx5.best, gpl85_1_hx5.bestN5, file = paste("gpl85_", i, "_best.rdata", sep = ""))
  setwd("..")
  
  # # impute NAs and check for non-informative features
  # moses.data$gpl96 <- bin.impute.matrix(gpl96)
  # summary(colSums(moses.data$gpl96))
  # #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # #   16.00   16.00   16.00   16.09   16.00   20.00 
  # dim(moses.data$gpl96)
  # # [1]    30 22216
  
  ### dataset gpl96
  
  # gpl96 <- data.frame(moses.data$gpl96, check.names = TRUE)
  dataset <- "gpl96"
  cr <- mean(gpl96[, "age"])
  make.dir(getwd(), dataset)
  gpl96_1.test <- makeMpartitions(gpl96, p = .5)
  save(gpl96_1.test, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  gpl96_1_hx5.train <- runMfolder(moses.flags)
  gpl96_1_hx5 <- testClist(gpl96_1_hx5.train, gpl96_1.test, cr)
  gpl96_1_hx5.bestN5 <- bestCombos(gpl96_1_hx5, N = 5)
  save(gpl96_1_hx5.bestN5, gpl96_1_hx5.bestN5, file = paste("gpl96_", i, "_best.rdata", sep = ""))
  setwd("..")
  
  # # impute NAs and check for non-informative features
  # moses.data$gpl97 <- bin.impute.matrix(gpl97)
  # summary(colSums(moses.data$gpl97))
  # #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # #   16.00   18.00   18.00   18.07   18.00   22.00 
  # dim(moses.data$gpl97)
  # # [1]    34 12489
  
  ### dataset gpl97
  
  # gpl97 <- data.frame(moses.data$gpl97, check.names = TRUE)
  dataset <- "gpl97"
  cr <- mean(gpl97[, "age"])
  make.dir(getwd(), dataset)
  gpl97_1.test <- makeMpartitions(gpl97, p = .5)
  save(gpl97_1.test, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  gpl97_1_hx5.train <- runMfolder(moses.flags)
  gpl97_1_hx5 <- testClist(gpl97_1_hx5.train, gpl97_1.test, cr)
  gpl97_1_hx5.best <- bestCombos(gpl97_1_hx5)
  gpl97_1_hx5.bestN5 <- bestCombos(gpl97_1_hx5, N = 5)
  save(gpl97_1_hx5.best, gpl97_1_hx5.bestN5, file = paste("gpl97_", i, "_best.rdata", sep = ""))
  setwd("..")
  
  # # impute NAs and check for non-informative features
  # moses.data$gpl341 <- bin.impute.matrix(gpl341)
  # summary(colSums(moses.data$gpl341))
  # #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # #   50.00   51.00   51.00   51.13   51.00   99.00 
  # dim(moses.data$gpl341)
  # # [1]    99 15924
  
  # one sample with all positive features - leave in and TODO:  see if it shows up in results...
  # names(gpl341)[colSums(gpl341) == 99]
  # [1] "X1370846_at"
  # gpl341[, "X1370846_at"]
  #  [1] 1 1 1 1 1 ...
  
  ### dataset gpl341
  
  # gpl341 <- data.frame(moses.data$gpl341, check.names = TRUE)
  dataset <- "gpl341"
  cr <- mean(gpl341[, "age"])
  make.dir(getwd(), dataset)
  gpl341_1.test <- makeMpartitions(gpl341, p = .5)
  save(gpl341_1.test, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  gpl341_1_hx5.train <- runMfolder(moses.flags)
  gpl341_1_hx5 <- testClist(gpl341_1_hx5.train, gpl341_1.test, cr)
  gpl341_1_hx5.best <- bestCombos(gpl341_1_hx5)
  gpl341_1_hx5.bestN5 <- bestCombos(gpl341_1_hx5, N = 5)
  save(gpl341_1_hx5.best, gpl341_1_hx5.bestN5, file = paste("gpl341_", i, "_best.rdata", sep = ""))
  setwd("..")
  
  # # impute NAs and check for non-informative features
  # moses.data$gpl1261 <- bin.impute.matrix(gpl1261)
  # summary(colSums(moses.data$gpl1261))
  # #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # #   16.00   18.00   18.00   18.07   18.00   22.00 
  # dim(moses.data$gpl1261)
  # # [1]    34 12489
  
  ### dataset gpl1261
  
  # gpl1261 <- data.frame(moses.data$gpl1261, check.names = TRUE)
  dataset <- "gpl1261"
  cr <- mean(gpl1261[, "age"])
  make.dir(getwd(), dataset)
  gpl1261_1.test <- makeMpartitions(gpl1261, p = .5)
  save(gpl1261_1.test, file = paste("../", dataset, "_run", i, "_moses.Rdata", sep = ""))
  gpl1261_1_hx5.train <- runMfolder(moses.flags)
  gpl1261_1_hx5 <- testClist(gpl1261_1_hx5.train, gpl1261_1.test, cr)
  gpl1261_1_hx5.best <- bestCombos(gpl1261_1_hx5)
  gpl1261_1_hx5.bestN5 <- bestCombos(gpl1261_1_hx5, N = 5)
  save(gpl1261_1_hx5.best, gpl1261_1_hx5.bestN5, file = paste("gpl1261_", i, "_best.rdata", sep = ""))
  setwd("../..")
  save(gpl81_1_hx5.best, gpl81_1_hx5.bestN5, 
    gpl85_1_hx5.best, gpl85_1_hx5.bestN5, 
    gpl96_1_hx5.bestN5, gpl96_1_hx5.bestN5, 
    gpl97_1_hx5.best, gpl97_1_hx5.bestN5, 
    gpl341_1_hx5.best, gpl341_1_hx5.bestN5, 
    gpl1261_1_hx5.best, gpl1261_1_hx5.bestN5,
    file = paste("bestOfRun_", i, ".rdata"))
}


### functions

# if < 30% na replace with row average
bin.impute.row <- function(row) {
  if(sum(is.na(row)) / length(row) > .3) return(NULL)
  row[is.na(row)] <- rbinom(sum(is.na(row)), 1, mean(row, na.rm = TRUE))
  return(row)
}

bin.impute.matrix <- function(mat) {
  out <- apply(mat, 2, bin.impute.row)
  out[sapply(out, is.null)] <- NULL
  print(paste(deparse(substitute(mat)), "imputed!"))
  if(class(out) != "list") {
    print(paste(deparse(substitute(mat)), "generated an error."))
    return(out)
  }
  do.call(cbind, out)
}

make.dir <- function(mainDir = "~/projects/sf_results", subDir) {
  if (file.exists(paste(mainDir, subDir, "/", sep = "/", collapse = "/"))) {
    cat("subDir exists in mainDir and is a directory")
  } else if (file.exists(paste(mainDir, subDir, sep = "/", collapse = "/"))) {
    cat("subDir exists in mainDir but is a file")
    # you will probably want to handle this separately
  } else {
    cat("subDir does not exist in mainDir - creating")
    dir.create(file.path(mainDir, subDir))
  }  
  if (file.exists(paste(mainDir, subDir, "/", sep = "/", collapse = "/"))) {
    # By this point, the directory either existed or has been successfully created
    setwd(file.path(mainDir, subDir))
  } else {
    cat("subDir does not exist")
    # Handle this error as appropriate
  }
}
