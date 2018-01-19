# BIC准则选择节点函数
BIC_knot<-function(low,up,yj,uj,n,alpha,scale=1){
  lowdegree <- 1
  updegree <- 3
  BIC<- matrix(0,length(lowdegree:updegree),up-low+1)
  for(k in low:up){
    for(degree in lowdegree:updegree){
      u0<- seq(0,1.01,length.out = k)*scale
      # fit1<- lm(yj ~ bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
      # fit1<- lm(yj ~ sapply(0:(length(u0)-2+degree),BaSplite,x=uj,degree=degree,u=u0,n=length(u0))-1)
      fit1<- lm(yj ~ BaSplite1(uj,degree,u0)-1)
      y1<- fitted(fit1)
      RSS0<-sum(fit1$residuals^2)    # cat('RSS0',RSS0,'\n')
      BIC[degree-lowdegree+1,k-low+1]<- log(RSS0/n)+(k+degree-1)/n*log(n)
    }
  }
  tmp <- which(BIC==min(BIC),arr.ind = TRUE)[1,]
  degree <- tmp[1]-1+lowdegree
  k0 <- tmp[2]-1+low
  # q2 <- qchisq(1-alpha,degree+1)
  err <- 0.8/(k0-1)
  #cat('k0=',k0,'degree=',degree,'\n')
  u0<- seq(0,1.01,length.out = k0)*scale
  nu0<- length(u0)
  BIC_list<-list(u0,nu0,degree,k0,err)
  return(BIC_list)
}
