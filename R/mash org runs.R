### load bed stuy moses results  TODO:  R.utils::use attachLocally
library("stringr")
setwd("~/GitHub/org_run3")

# attach output .rdata file from each run
for(i in 1:10) attach(paste("bestOfRun_", i, ".rdata"))

# get highest scoring combos
bests <- c("gpl341b_1_hx5.best", "homo_1_hx5.best", "mouse_1_hx5.best", "rat_1_hx5.best") 
bests.cum <- str_replace(bests, fixed("_1_"), "_cum_")


# get the best of the runs for each chip
getBlist <- function(obj, run) get(obj, pos = paste("file:bestOfRun_", run, ".rdata"))
for(i in 1:10){
  for(j in 1:4){
    rlist <- getBlist(bests[j], i)
    assign(str_replace(bests[j], fixed("_1_"), paste("_", i, "_", sep = "")), rlist, inherits = TRUE)
  }
}

for(n in grep("file", search(), value = TRUE)) detach(n, char = TRUE)

# bind together results by data set
lbind <- function(list1, list2) {
  if(sum(sapply(list1, class) != sapply(list2, class)) > 0) return("lists element classes do not match")
  out <- vector("list", length(list1))
  names(out) <- names(list1)
  for(i in seq_along(out)) {
    ifelse(is.atomic(list1[[i]]), out[[i]] <- c(list1[[i]], list2[[i]]), out[[i]] <- rbind(list1[[i]], list2[[i]]))
  }
    return(out)
}

# go through list of lists of results and mash them together
for(i in 1:4) assign(bests.cum[i], Reduce(lbind, lapply(ls(pos = 1, pattern = unlist(str_split(bests.cum[i], "_", 2))[1]), get)))

save(list = ls(pattern = "gpl341b"), file = "gpl341b")
save(list = ls(pattern = "homo"), file = "homo.rdata")
save(list = ls(pattern = "mouse"), file = "mouse.rdata")
save(list = ls(pattern = "rat"), file = "rat.rdata")
rm(list = ls(pattern = "[0123456789]_hx5"))

# there are 21 gpl341b, 15 human, 12 mouse , 17 gpl96 genes both hi & lo
sapply(lapply(ls(pattern = "_cum_"), get), function(x) c(length(x$combo), dim(x$features)[1]))
#      [,1] [,2] [,3] [,4]
# [1,]   56   29   84   47
# [2,]  131  193  196  231

sapply(lapply(ls(pattern = "_cum_"), get), function(x) c(length(unique(x$combo)), dim(unique(x$features))[1]))
#      [,1] [,2] [,3] [,4]
# [1,]   56   29   84   47
# [2,]  110  178  182  214

# look at ranks & accuracies of best combos
unlist(sapply(lapply(ls(pattern = "_cum_"), get), function(x) summary(as.factor(x$rank))))
#            gpl341b   human       mouse         rat
# rank        -1   0  -1   0  -2  -1   0  -2  -1   0 
# features    10  46   9  20   6  42  36   1  72  43

unlist(sapply(lapply(ls(pattern = "_cum_"), get), function(x) summary(as.factor(signif(x$score$Accuracy, 3)))))
#       gpl341b               human                   mouse                             rat
#   0.9  0.95 1   0.843 0.863 0.882   0.8 0.825  0.85 0.875   0.839 0.857 0.875 0.893 0.929 
#    19    32 5       9     9    11    12     7    43    22      56    27    23     9     1 

unlist(sapply(lapply(ls(pattern = "_cum_"), get), function(x) summary(as.factor(signif(x$score$`Balanced Accuracy`, 3)))))
#                   gpl341b                                                         human 
# 0.889 0.899 0.944 0.955 1   0.831  0.84 0.845 0.848 0.862 0.864 0.869 0.871 0.879 0.886 
#     6    13    18    14 5       4     1     5     5     2     1     1     1     2     7 
#                                                 mouse 
# 0.793 0.803 0.816 0.838 0.848 0.854 0.859 0.866 0.871 
#     2    10     7     6    10    19     8    15     7    
#                                                               rat
# 0.776 0.786   0.8 0.819 0.843 0.852 0.867 0.871 0.881  0.89 0.914
#    10     9     1     1     1    10     1     3     7     3     1 
    
save(list = ls(pattern = "_cum_"), file = "~/GitHub/stevia/data/orgCumX4.rdata")

# replace old feature lists with "Freq" (combo count) instead of "rank"
load("~/GitHub/stevia/data/gplCumX6.rdata")
rm(gpl341_cum_hx5.best)
# fix gpl96:  bestN5 -> best
gpl96_cum_hx5.best <- lapply(gpl96_cum_hx5.bestN5[1:2], function(x) x[gpl96_cum_hx5.bestN5$score[[13]] > .85])
gpl96_cum_hx5.best$score <- gpl96_cum_hx5.bestN5$score[gpl96_cum_hx5.bestN5$score[[1]] > .85,]
# uses combo2flist from "moses.R"
gpl96_cum_hx5.best$features <- combo2flist(gpl96_cum_hx5.best$combo, gpl96_cum_hx5.best$rank)
rm(gpl96_cum_hx5.bestN5)

# replace features data frame with combo2fcount
for(n in ls(pattern = "_cum_")) {
  old <- get(n)
  if(length(old$features[[1]]) == length(grep("^X", old$features[[1]]))) strip <- TRUE else strip <- FALSE
  old$features <- combo2fcount(old$combo, stripX = strip)
  assign(gsub("cum_hx5", "", n, fixed = TRUE), old, pos = 1)
}
rm(list = ls(pattern = "_cum_"))

# save .csv's with combo2Fcsv <- function(combo, name = deparse(substitute(combo)), dir = ".", strip = FALSE, ret = FALSE, noHiLo = FALSE)
for(n in ls(pattern = "best$")) {
combo2Fcsv(get(n)$combo, name = gsub(".best$", "", n))
combo2Fcsv(get(n)$combo, name = gsub(".best$", "_NoHiLo", n), noHiLo = TRUE)
}

