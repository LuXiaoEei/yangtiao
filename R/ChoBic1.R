ChoBic1 <- function(lowerKnot,upperKnot,Post,Y,AllPost=NULL,AllY=NULL,updegree=3,lowdegree=1){

  # if(missing(AllPost)) AllPost <- Post
  # if(missing(AllY)) AllY <- Y

  # lambda <- 0.0001

  # updegree <- 3
  # lowdegree <- 1
  lenDeg <- length(lowdegree:updegree)
  lenKnot <- length(lowerKnot:upperKnot)

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


  #分成两个方向的样条

  Arg <- matrix(c(
  rep(c(rep(lowdegree:updegree,lenKnot),rep(lowerKnot:upperKnot,each=lenDeg)),each=lenDeg*lenKnot),
  c(rep(rep(lowdegree:updegree,lenKnot),lenDeg*lenKnot),rep(rep(lowerKnot:upperKnot,each=lenDeg),lenDeg*lenKnot))
  ),ncol=4)
  Argx <- matrix(c(rep(lowdegree:updegree,lenKnot),rep(lowerKnot:upperKnot,each=lenDeg)),ncol = 2)


  Bic <- function(DegreeY,KnotNumY){
    KnotY <- seq(range(Post[,2])[1]-1,range(Post[,2])[2]+1,length.out = KnotNumY)
    G <- matrix(0,nrow = length(XPost),ncol = KnotNumY+DegreeY-1)
    for (ii in 1:length(XPost)){
      index <- XPost[ii]
      X <- Ylist[[KnotNumY-lowerKnot+1]][[DegreeY-lowdegree+1]][Post[,1]==index,,drop=FALSE]
      # X <- BaSplite1(x = Post[Post[,1]==index,2],degree = DegreeY,u = KnotY)
      if(sum(Post[,1]==index)>KnotNumY+DegreeY+5){
        # AllX <- BaSplite1(x = AllPost[AllPost[,1]==index,2],degree = DegreeY,u = KnotY)
        G[ii,] <- MASS::ginv(t(X)%*%X)%*%t(X)%*%Y[Post[,1]==index]
      }else{
        G[ii,1] <- Inf
      }
    }
    Gu <- XPost[!is.infinite(G[,1])]
    G <- G[!is.infinite(G[,1]),,drop=FALSE]


    BicX <- function(DegreeX,KnotNumX){
      # DegreeY <- Argx[index,1]
      # KnotNumY <- Argx[index,2]
      KnotX <- seq(range(Post[,1])[1]-1,range(Post[,1])[2]+1,length.out = KnotNumX)
      V <- matrix(0,nrow = KnotNumX+DegreeX-1,ncol = ncol(G))
      for (ii in 1:ncol(G)){
        X <- BaSplite1(x = Gu,degree = DegreeX,u = KnotX)
        V[,ii] <- MASS::ginv(t(X)%*%X)%*%t(X)%*%G[,ii]
      }

      # AllBx <- BaSplite1(x=AllPost[,1],degree = DegreeX,u = KnotX)
      # AllBy <- BaSplite1(x = AllPost[,2],degree = DegreeY,u = KnotY)
      # AllX <- AllBx[,rep(1:ncol(AllBx),ncol(AllBy))]*AllBy[,rep(1:ncol(AllBy),each=ncol(AllBx))]

      Bx <- BaSplite1(x=Post[,1],degree = DegreeX,u = KnotX)
      By <- BaSplite1(x = Post[,2],degree = DegreeY,u = KnotY)
      X <- Bx[,rep(1:ncol(Bx),ncol(By))]*By[,rep(1:ncol(By),each=ncol(Bx))]

      V <- matrix(V,ncol = 1)
      residual <- Y-X%*%V
      # (t(Y)%*%(I-(2-nrow(X)*lambda)*X%*%solve(t(X)%*%X+I1)%*%t(X)+X%*%solve(t(X)%*%X+I1)%*%t(X)%*%X%*%solve(t(X)%*%X+I1)%*%t(X))%*%Y)/length(Y)
      # return(mean(residual^2)+lambda*sum(V^2))
      return(log(mean(residual^2))+ncol(X)*log(nrow(X))/nrow(X))
      # return(mean(residual^2))
    }
    return(mapply(BicX,Argx[,1],Argx[,2]))
  }

  XPost <- unique(Post[,1])

  Res <- matrix(mapply(Bic,Argx[,1],Argx[,2]),ncol = 1)

  fin <- Arg[which.min(Res),]

  KnotNumX <- fin[4]
  KnotNumY <- fin[2]
  DegreeX <- fin[3]
  DegreeY <- fin[1]

  return(list(KnotNumX=KnotNumX,KnotNumY=KnotNumY,DegreeX=DegreeX,DegreeY=DegreeY))
}


