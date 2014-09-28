### get test & train files, extract gene lists
library("stringr")
setwd("C:/Users/user/Desktop/biomind/addis/icog/eachData/moses/mose_each")

# get file list
getRdata <- function(dir = getwd()) {
  require(stringr)
  list.files(path = dir)[str_detect(list.files(path = dir), ignore.case("rdata$"))]
}
  
Mtest <- function(train, test, n = 10, trim = 5){
  out <- vector("list", 3)
  names(out) <- c("combo", "features", "scores")
  for(i in 1:n) {
    fold <- Mout2str(train[i], trimS = trim)
    out <- lbind2(out, testClist(fold[[1]], test[i]))
  }
  return(out)
}

mget(c("data_1.hx5.train", "data_1.test"), pos = paste("file:", files[i], sep = ""))



# bind together results from multiple moses runs using combo counting arguement
# this function will attempt to concatinate vectors and row bind higher dim objects in 2 lists if objects match by class in each position.  row bound objects are sum-aggregated by variables "feature" & "upregulated"
# TODO: generalize to arbitrary matching lists
lbindWfeatureDf <- function(list1, list2) {
  if(sum(sapply(list1, class) != sapply(list2, class)) > 0) return("lists element classes do not match")
  out <- vector("list", length(list1))
  names(out) <- names(list1)
  for(i in seq_along(out)) {
    if(is.atomic(list1[[i]])) {out[[i]] <- c(list1[[i]], list2[[i]])}
    if(!is.atomic(list1[[i]]))  {
      fdf <- rbind(list1[[i]], list2[[i]])
      out[[i]] <- aggregate(. ~ feature + upregulated, data = fdf, sum)
      out[[i]] <- out[[i]][order(out[[i]][, 1]),]
    }
  }
  return(out)
}
