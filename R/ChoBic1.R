ChoBic1 <- function(lowerKnot,upperKnot,Post,Y,AllPost,AllY,updegree=3,lowdegree=1){

  if(missing(AllPost)) AllPost <- Post
  if(missing(AllY)) AllY <- Y

  # updegree <- 3
  # lowdegree <- 1
  lenDeg <- length(lowdegree:updegree)
  lenKnot <- length(lowerKnot:upperKnot)

  Xlist <- list()
  AllXlist <- list()
  for (KnotNumX in lowerKnot:upperKnot){
    KnotX <- seq(range(Post[,1])[1]-1,range(Post[,1])[2]+1,length.out = KnotNumX)
    Xlist[[KnotNumX-lowerKnot+1]] <- sapply(lowdegree:updegree, BaSplite1,x=Post[,1],u=KnotX,simplify = FALSE)
    AllXlist[[KnotNumX-lowerKnot+1]] <- sapply(lowdegree:updegree,BaSplite1,x=AllPost[,1],u=KnotX,simplify = FALSE)
  }

  Ylist <- list()
  AllYlist <- list()
  for (KnotNumY in lowerKnot:upperKnot){
    KnotY <- seq(range(Post[,2])[1]-1,range(Post[,2])[2]+1,length.out = KnotNumY)
    Ylist[[KnotNumY-lowerKnot+1]] <- sapply(lowdegree:updegree, BaSplite1,x=Post[,2],u=KnotY,simplify = FALSE)
    AllYlist[[KnotNumY-lowerKnot+1]] <- sapply(lowdegree:updegree, BaSplite1,x=AllPost[,2],u=KnotY,simplify = FALSE)
  }


  #分成两个方向的样条

#
#   DegreeY <- 1
#   KnotNumY <- 2
#   DegreeX <- 1
#   KnotNumX <- 2

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

      AllBx <- BaSplite1(x=AllPost[,1],degree = DegreeX,u = KnotX)
      AllBy <- BaSplite1(x = AllPost[,2],degree = DegreeY,u = KnotY)
      AllX <- AllBx[,rep(1:ncol(AllBx),ncol(AllBy))]*AllBy[,rep(1:ncol(AllBy),each=ncol(AllBx))]
      residual <- AllY-AllX%*%matrix(V,ncol = 1)
      return(log(mean(residual^2))+ncol(AllX)*log(nrow(AllX))/nrow(AllX))
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


