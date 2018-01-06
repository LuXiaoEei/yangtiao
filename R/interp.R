# 最近邻 因变量y是当x给定时的函数值。当自变量换为x0时，此函数返回因变量的值。
interp<- function(y,x,x0){
  n<- length(x0)
  yhat<- matrix(rep(0),n,1)
  n1<- length(x)
  for (i in 1:n) {
    a<- abs(x-x0[i])
    yhat[i]<- y[order(a)[1]]
  }
  return(yhat)
}
