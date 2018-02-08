# BIC准则选择节点函数
BIC_knot<-function(low,up,yj,uj,updegree=3,scale=1){
  n <- length(yj)
  lowdegree <- 1
  BIC<- matrix(0,length(lowdegree:updegree),up-low+1)
  for(k in low:up){
    for(degree in lowdegree:updegree){
      # u0<- seq(0,1+0.5/scale,length.out = k)*scale
      u0<- seq(0,1,length.out = k)*scale+0.5
      # u0<- c(0,seq(1/(k-1),1,length.out = k-2),1+0.5/scale)*scale
      # fit1<- lm(yj ~ bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
      # fit1<- lm(yj ~ sapply(0:(length(u0)-2+degree),BaSplite,x=uj,degree=degree,u=u0,n=length(u0))-1)
      # k <- k-2
      # u0<- seq(1/(k+1),k/(k+1),length.out = k)*scale
      # B<-splineDesign(knots=c(rep(0,degree+1),u0,rep(1,degree+1)),ord=degree+1,x=uj,outer.ok = TRUE)

      fit1<- lm(yj ~ BaSplite1(uj,degree,u0)-1)
      # fit1<- lm(yj ~ B-1)
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
  # u0<- seq(0,1.01,length.out = k0)*scale
  u0<- seq(0,1,length.out = k)*scale+0.5
  BIC_list<-list(u0,degree,k0,err)
  return(BIC_list)
}
