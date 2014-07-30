### icog virtualArray test

xx<- jpmEset $r_hippo
yy<- jpmEset $r_heart
annotation(xx) <- "rae230a"
annotation(yy) <- "rgu34a"


virtArray_mouse <-NULL
virtArray_mouse[["wBatchEffects"]] <- virtualArrayExpressionSets(all_expression_sets=c(xx, yy))

#### functions

rownames2col <- function(df, colname) {
  output <- cbind(row.names(df), df)
  colnames(output)[1] <- colname
  return(output)
  }

col2rownames <- function(df, colname, removecol=FALSE){
  row.names(df) <- df[,colname]
  if(removecol){df[,colname] <- NULL}
  return(df)
  }
