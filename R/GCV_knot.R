# GCV准则选择节点函数
# B<-bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)

GCV_knot<-function(low,up,yj,uj,updegree=3,scale=1){
  n <- length(yj)
  lowdegree <- 1

  Lambda <- c(0.0001,0.001,0.01,0.1,1,10,100,1000)
  GCV<-matrix(0,nrow = length(Lambda),ncol = 4)
  Args <- matrix(c(rep(low:up,each=length(lowdegree:updegree)),rep(lowdegree:updegree,length(low:up))),ncol = 2)
  C_lambda <- matrix(0,nrow=nrow(Args),ncol=length(Lambda))
  # for (index in 1:length(Lambda)){
  #   lambda <- Lambda[index]


  for(Index in 1:nrow(Args)){
    # print(Index)
    k <- Args[Index,1]
    degree <- Args[Index,2]
    u0<- seq(0,1+0.5/scale,length.out = k)*scale
    B <- BaSplite1(uj,degree,u0)
    # B<-bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,scale+1),intercept=T)
    # I <- diag(nrow = nrow(B))
    # II <- diag(nrow(B),nrow = ncol(B))
    C_lambda[Index,] <- sapply(Lambda,function(lambda){
      # a <- 2-n*lambda
      # I1 <- diag(lambda*n,nrow = ncol(B))
      # I1 <- diag(lambda,ncol(B))
      beta <- solve(t(B)%*%B+diag(lambda,nrow = ncol(B)))%*%t(B)%*%yj
      # P <- B%*%solve(t(B)%*%B+I1)%*%t(B)
      # return((t(yj)%*%(I-(2-nrow(B)*lambda)*B%*%solve(t(B)%*%B+I1)%*%t(B)+B%*%solve(t(B)%*%B+I1)%*%t(B)%*%B%*%solve(t(B)%*%B+I1)%*%t(B))%*%yj)/n)
      # return((sum(((B%*%solve(t(B)%*%B+diag(lambda*n,nrow = ncol(B)))%*%t(B)-diag(a/2,nrow = n))%*%yj)^2)+sum(yj^2)*(1-a^2/4))/n)
      return((sum((yj-B%*%beta)^2)+lambda*sum(beta^2))/n)
      })
  }

  for (index in 1:length(Lambda)){
    tmp <- which.min(C_lambda[,index])
    degree <- Args[tmp,2]
    k0 <- Args[tmp,1]
    lambda <- Lambda[index]
    # err <- 0.8/(k0-1)
    u0<- seq(0,1+0.5/scale,length.out = k0)*scale
    B <- BaSplite1(uj,degree,u0)
    # I <- diag(nrow = nrow(B))
    beta <- solve(t(B)%*%B)%*%t(B)%*%yj
    # t1 <- nrow(B)*(t(yj)%*%(I-P)%*%yj)
    t1 <- n*sum((yj-B%*%beta)^2)
    # t2 <- sum(diag(I-P))^2
    t2 <- sum(n-diag(B%*%solve(t(B)%*%B)%*%t(B)))^2
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
