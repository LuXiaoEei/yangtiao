div <- function(a,b){
  x <- a/b
  x[is.nan(a/b)|is.infinite(a/b)]=0
  return(x)
}
