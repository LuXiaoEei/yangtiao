# BIC准则选择节点和次数函数
BIC_knot_fit<-function(low,up,yj,uj,updegree=3,jumps,scale=1){
  n <- length(yj)
  lowdegree <- 1
  BIC<- matrix(0,length(lowdegree:updegree),up-low+1)
  for(k in low:up){
    for(degree in lowdegree:updegree){
      u0<- seq(0,1,length.out = k)*scale+0.5
      if(length(jumps)>0){
        # for (x in jumps){
        #   u0 <- u0[!((x-0.05*scale)<u0 & u0<(0.05*scale+x))]
        #   u0 <- sort(c(u0,rep(x,degree+1)))
        # }
        u0 <- setdiff(u0,jumps)
        u0 <- sort(c(u0,rep(jumps,degree+1)))
      }
      fit1 <- lm(yj ~ BaSplite1(uj,degree,u0)-1)
      y1<- fitted(fit1)
      RSS0<-sum(fit1$residuals^2)    # cat('RSS0',RSS0,'\n')
      BIC[degree-lowdegree+1,k-low+1] <- log(RSS0/n)+(k+degree-1)/n*log(n)
    }
  }
  tmp <- which(BIC==min(BIC),arr.ind = TRUE)[1,]
  degree <- tmp[1]-1+lowdegree
  k0 <- tmp[2]-1+low
  u0<- seq(0,1,length.out = k)*scale+0.5
  # cat(length(u0),degree,'\n')
  if(length(jumps)>0){
    # for (x in jumps){
    #   u0 <- u0[!((x-0.05*scale)<u0 & u0<(0.05*scale+x))]
    #   u0 <- sort(c(u0,rep(x,degree+1)))
    # }
    u0 <- setdiff(u0,jumps)
    u0 <- sort(c(u0,rep(jumps,degree+1)))
  }
  BIC_list<-list(u0,degree)
  return(BIC_list)
}
