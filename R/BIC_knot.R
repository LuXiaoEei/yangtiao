# BIC准则选择节点函数
BIC_knot<-function(low,up,yj,uj,n,alpha,scale=1)
{
  BIC<- matrix(0,3,up-low+1)
  for(k in low:up){
    for(degree in 1:3){
      u0<- seq(0,1,length.out = k)*scale
      # fit1<- lm(yj ~ bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
      # fit1<- lm(yj ~ sapply(0:(length(u0)-2+degree),BaSplite,x=uj,degree=degree,u=u0,n=length(u0))-1)
      fit1<- lm(yj ~ BaSplite1(uj,degree,u0)-1)
      y1<- fitted(fit1)
      RSS0<-sum(fit1$residuals^2)    # cat('RSS0',RSS0,'\n')
      BIC[degree,k-low+1]<- log(RSS0/n)+(k+degree+1)/n*log(n)
    }
  }
  k0<- ceiling(which.min(BIC)/3)+low-1
  degree<- (which.min(BIC)-1)%%3+1
  q2<- qchisq(1-alpha,degree+1)
  err<-  0.8/(k0-1)
  #cat('k0=',k0,'degree=',degree,'\n')
  u0<- seq(0,1,length.out = k0)*scale
  nu0<- length(u0)
  BIC_list<-list(u0,nu0,degree,k0,err)
  return(BIC_list)
}
