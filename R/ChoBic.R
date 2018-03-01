ChoBic <- function(lowerKnot,upperKnot,Post,Y,AllPost=NULL,AllY=NULL,updegree,lowdegree){

  if(missing(AllPost)) AllPost <- Post
  if(missing(AllY)) AllY <- Y

  # updegree <- 1
  # lowdegree <- 3
  lendegree <- length(lowdegree:updegree)

  Bic <- function(Xindex,Yindex){
    Bx <- Xlist[[Xindex%/%lendegree+1]][[Xindex%%lendegree+1]]
    By <- Ylist[[Yindex%/%lendegree+1]][[Yindex%%lendegree+1]]

    # AllBx <- AllXlist[[Xindex%/%lendegree+1]][[Xindex%%lendegree+1]]
    # AllBy <- AllYlist[[Yindex%/%lendegree+1]][[Yindex%%lendegree+1]]

    X <- Bx[,rep(1:ncol(Bx),ncol(By))]*By[,rep(1:ncol(By),each=ncol(Bx))]

    # AllX <- AllBx[,rep(1:ncol(AllBx),ncol(AllBy))]*AllBy[,rep(1:ncol(AllBy),each=ncol(AllBx))]

    residual <- Y-X %*%(MASS::ginv(t(X)%*%X)%*%t(X)%*%Y)
    return(log(mean(residual^2))+ncol(X)*log(nrow(X))/nrow(X))
    # return(mean(residual^2))
  }

  Xlist <- list()
  # AllXlist <- list()
  for (KnotNumX in lowerKnot:upperKnot){
    KnotX <- seq(range(Post[,1])[1]-1,range(Post[,1])[2]+1,length.out = KnotNumX)
    Xlist[[KnotNumX-lowerKnot+1]] <- sapply(lowdegree:updegree, BaSplite1,x=Post[,1],u=KnotX,simplify = FALSE)
    # AllXlist[[KnotNumX-lowerKnot+1]] <- sapply(lowdegree:updegree,BaSplite1,x=AllPost[,1],u=KnotX,simplify = FALSE)
  }

  Ylist <- list()
  # AllYlist <- list()
  for (KnotNumY in lowerKnot:upperKnot){
    KnotY <- seq(range(Post[,2])[1]-1,range(Post[,2])[2]+1,length.out = KnotNumY)
    Ylist[[KnotNumY-lowerKnot+1]] <- sapply(lowdegree:updegree, BaSplite1,x=Post[,2],u=KnotY,simplify = FALSE)
    # AllYlist[[KnotNumY-lowerKnot+1]] <- sapply(lowdegree:updegree, BaSplite1,x=AllPost[,2],u=KnotY,simplify = FALSE)
  }

  len <- upperKnot-lowerKnot+1
  # Bic <- matrix(NA,nrow = 4*len, 4*len)
  Arg <- matrix(c(rep(1:(lendegree*len),lendegree*len),rep(1:(lendegree*len),each=lendegree*len)),ncol = 2)

  res <- mapply(Bic,Arg[,1]-1,Arg[,2]-1)

  Xindex <- Arg[which.min(res)[1],][1]-1
  Yindex <- Arg[which.min(res)[1],][2]-1

  KnotNumX <- Xindex%/%lendegree+lowerKnot
  KnotNumY <- Yindex%/%lendegree+lowerKnot
  DegreeX <- Xindex%%lendegree+lowdegree
  DegreeY <- Yindex%%lendegree+lowdegree

  return(list(KnotNumX=KnotNumX,KnotNumY=KnotNumY,DegreeX=DegreeX,DegreeY=DegreeY))
}
