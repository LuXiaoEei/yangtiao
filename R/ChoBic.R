ChoBic <- function(lower,upper,Post,Y,AllPost,AllY){

  if(missing(AllPost)) AllPost <- Post
  if(missing(AllY)) AllY <- Y

  updegree <- 3
  lowdegree <- 1
  lendegree <- length(lowdegree:updegree)

  Bic <- function(Xindex,Yindex){
    Bx <- Xlist[[Xindex%/%lendegree+1]][[Xindex%%lendegree+1]]
    By <- Ylist[[Yindex%/%lendegree+1]][[Yindex%%lendegree+1]]

    AllBx <- AllXlist[[Xindex%/%lendegree+1]][[Xindex%%lendegree+1]]
    AllBy <- AllYlist[[Yindex%/%lendegree+1]][[Yindex%%lendegree+1]]

    X <- Bx[,rep(1:ncol(Bx),ncol(By))]*By[,rep(1:ncol(By),each=ncol(Bx))]

    AllX <- AllBx[,rep(1:ncol(AllBx),ncol(AllBy))]*AllBy[,rep(1:ncol(AllBy),each=ncol(AllBx))]

    residual <- AllY-AllX %*%(MASS::ginv(t(X)%*%X)%*%t(X)%*%Y)
    return(log(mean(residual^2))+ncol(X)*log(nrow(X))/nrow(X))
  }

  Xlist <- list()
  AllXlist <- list()
  for (KnotNumX in lower:upper){
    KnotX <- seq(range(Post[,1])[1]-1,range(Post[,1])[2]+1,length.out = KnotNumX)
    Xlist[[KnotNumX-lower+1]] <- sapply(lowdegree:updegree, BaSplite1,x=Post[,1],u=KnotX)
    AllXlist[[KnotNumX-lower+1]] <- sapply(lowdegree:updegree,BaSplite1,x=AllPost[,1],u=KnotX)
  }

  Ylist <- list()
  AllYlist <- list()
  for (KnotNumY in lower:upper){
    KnotY <- seq(range(Post[,2])[1]-1,range(Post[,2])[2]+1,length.out = KnotNumY)
    Ylist[[KnotNumY-lower+1]] <- sapply(lowdegree:updegree, BaSplite1,x=Post[,2],u=KnotY)
    AllYlist[[KnotNumY-lower+1]] <- sapply(lowdegree:updegree, BaSplite1,x=AllPost[,2],u=KnotY)
  }

  len <- upper-lower+1
  # Bic <- matrix(NA,nrow = 4*len, 4*len)
  Arg <- matrix(c(rep(1:(lendegree*len),lendegree*len),rep(1:(lendegree*len),each=lendegree*len)),ncol = 2)

  res <- mapply(Bic,Arg[,1]-1,Arg[,2]-1)

  Xindex <- Arg[which.min(res)[1],][1]-1
  Yindex <- Arg[which.min(res)[1],][2]-1

  KnotNumX <- Xindex%/%lendegree+lower
  KnotNumY <- Yindex%/%lendegree+lower
  DegreeX <- Xindex%%lendegree+lowdegree
  DegreeY <- Yindex%%lendegree+lowdegree

  return(list(KnotNumX=KnotNumX,KnotNumY=KnotNumY,DegreeX=DegreeX,DegreeY=DegreeY))
}
