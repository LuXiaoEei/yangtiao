# GCV准则选择节点函数
# B<-bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)
GCV_knot<-function(low,up,yj,uj,updegree=3,scale=1){
  n <- length(yj)
  lowdegree <- 1
  C_lambda<-matrix(0,nrow=updegree-lowdegree+1,ncol=up-low+1)
  Lambda <- seq(0.00001,0.001,0.00001)
  GCV<-matrix(0,nrow = length(Lambda),ncol = 4)
  for (index in 1:length(Lambda)){
    lambda <- Lambda[index]
    for(k in low:up){
      for(degree in lowdegree:updegree){
        u0<- seq(0,1+0.5/scale,length.out = k)*scale
        B <- BaSplite1(uj,degree,u0)
        I <- diag(1,nrow = nrow(B))
        I1 <- lambda*nrow(B)*diag(rep(1, ncol(B)))
        C_lambda[degree-lowdegree+1,k-low+1]<- (t(yj)%*%(I-(2-nrow(B)*lambda)*B%*%solve(t(B)%*%B+I1)%*%t(B)+B%*%solve(t(B)%*%B+I1)%*%t(B)%*%B%*%solve(t(B)%*%B+I1)%*%t(B))%*%yj)/n
      }
    }

    tmp <- which(C_lambda==min(C_lambda),arr.ind = TRUE)[1,]
    degree <- tmp[1]-1+lowdegree
    k0 <- tmp[2]-1+low
    # cat('degree',degree,'k0',k0,'\n')
    err <- 0.8/(k0-1)
    #cat('k0=',k0,'degree=',degree,'\n')
    u0<- seq(0,1+0.5/scale,length.out = k0)*scale
    B <- BaSplite1(uj,degree,u0)
    I <- diag(1,nrow = nrow(B))
    t1 <- nrow(B)*(t(yj)%*%(I-B%*%solve(t(B)%*%B)%*%t(B))%*%yj)
    t2 <- sum(diag(I-B%*%solve(t(B)%*%B)%*%t(B)))^2

    # cat(index,'\n')
    GCV[index,]<- c(degree,k0,lambda,t1/t2)
  }
  lambda_op <- which.min(GCV[,4])

  degree <- GCV[lambda_op,1]
  k0 <- GCV[lambda_op,2]
  lambda <- GCV[lambda_op,3]
  # q2 <- qchisq(1-alpha,degree+1)
  err <- 0.8/(k0-1)
  #cat('k0=',k0,'degree=',degree,'\n')
  u0<- seq(0,1+0.5/scale,length.out = k0)*scale

  GCV_list<-list(u0,degree,k0,err,lambda)
  return(GCV_list)
}
