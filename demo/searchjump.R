# 带跳曲面拟合
require(yangtiao)
require(EBImage)

if(!exists('K')) K <- c(0) #K表示旋转的程度，旋转度数的正切值K

AllJumps3 <- matrix(NA,ncol = 2)
AllJumps3 <- AllJumps3[-1,]

message('######## ','Detecting...',' ########','\n')
for (k in K){
  message('######## k = ',k,' ########','\n')
  if(k==0){
    for (Index in c(FALSE,TRUE)){
      cat(Index,'\n')
      # 设置初始值 默认处理一列
      if (Index) {a <- COL;COL <- ROW;ROW <- a;Image_noise <- t(Image_noise)}
      n <- ROW   ##设置样本量
      # 设置节点选择时的参数
      alpha <- 0.05
      low <- 4
      up <- 14
      uj <- c(1:ROW)
      u1 <- c(seq(1,5,length.out = 15),seq(6,COL-5,1),seq(COL-4,COL,length.out = 15))#寻找跳点时的候选节点,加大边界寻找密度

      # 开始拟合 --------------------------------------------------------------------
      # Image_hat <- matrix(0,ROW,COL)
      AllJumps <- matrix(NA,ncol = 2)
      AllJumps <- AllJumps[-1,]

      for (index in 1:COL){
        # 选择sample
        yj <- Image_noise[,index]
        sigma2hat<- 1/2*1/(n-1)*sum((yj[-1]-yj[-n])^2)

        BIC_list<-BIC_knot(low,up,yj,uj,n,alpha,scale=ROW) # BIC准则选择节点
        u0<-BIC_list[[1]]         ##初始节点
        nu0<-BIC_list[[2]]        ##初始节点数
        degree<-BIC_list[[3]]     ##初始次数
        k0<-BIC_list[[4]]      ##初始节点数
        err<-BIC_list[[5]]

        jumps_find_list<-jumps_find(u0,u1,degree,yj,uj,err,sigma2hat,alpha,gamma = 0)
        u0<-jumps_find_list[[1]]
        jumps2<-jumps_find_list[[2]]

        if (length(jumps2)>0){
          AllJumps <- rbind(AllJumps,matrix(c(jumps2,rep(index,length(jumps2))),nrow = length(jumps2)))
        }
      }
      if(Index){
        AllJumpstmp2 <- AllJumps[,2:1]
        Image_noise <- t(Image_noise)
        a <- COL
        COL <- ROW
        ROW <- a
      }else{
        AllJumpstmp <- AllJumps
      }
    }
    # eval(parse(text = paste0('AllJumps',k,' <- round(rbind(AllJumpstmp2,AllJumpstmp))')))
    AllJumps3 <- rbind(AllJumps3,AllJumpstmp2,AllJumpstmp)
  }else{
    for (Index in c(FALSE,TRUE)){
      cat(Index,'\n')
      Num <- 1
      index <- 0
      AllJumps <- matrix(NA,ncol = 2)
      AllJumps <- AllJumps[-1,]
      if (Index) {a <- COL;COL <- ROW;ROW <- a;Image_noise <- t(Image_noise)}
      while (Num > 0){
        index <- index+1
        if(Index){
          allpoints <- matrix(c(1:ROW,1/k*c(1:ROW)-ROW/k+index),ncol = 2)
        }else{
          allpoints <- matrix(c(1:ROW,-1/k*c(1:ROW)+1/k+index),ncol = 2)
        }
        allpoints <- round(allpoints)
        IIndex <- allpoints[,2]<=COL&allpoints[,2]>=1&as.integer(allpoints[,2])==allpoints[,2]
        yj <- Image_noise[allpoints[IIndex,,drop=FALSE]]
        Num <- length(yj)
        if(Num>30){
          n <- Num   ##设置样本量
          # 设置节点选择时的参数
          alpha <- 0.05
          low <- 4
          up <- 14
          uj <- c(1:n)
          u1 <- seq(1,n,1)#寻找跳点时的候选节点

          sigma2hat<- 1/2*1/(n-1)*sum((yj[-1]-yj[-n])^2)

          BIC_list<-BIC_knot(low,up,yj,uj,n,alpha,scale=Num) # BIC准则选择节点
          u0<-BIC_list[[1]]         ##初始节点
          nu0<-BIC_list[[2]]        ##初始节点数
          degree<-BIC_list[[3]]     ##初始次数
          k0<-BIC_list[[4]]      ##初始节点数
          err<-BIC_list[[5]]

          jumps_find_list<-jumps_find(u0,u1,degree,yj,uj,err,sigma2hat,alpha,gamma = 0)
          u0 <- jumps_find_list[[1]]
          jumps2 <- jumps_find_list[[2]]

          if (length(jumps2)>0){
            Points <- allpoints[IIndex,,drop=FALSE]
            AllJumps <- rbind(AllJumps,Points[jumps2,,drop=FALSE])
          }
        }
      }
      if(Index){
        AllJumpstmp2 <- AllJumps[,2:1]
        Image_noise <- t(Image_noise)
        a <- COL
        COL <- ROW
        ROW <- a
      }else{
        AllJumpstmp <- AllJumps
      }
    }
  }
  # eval(parse(text = paste0('AllJumps',k,' <- round(rbind(AllJumpstmp2,AllJumpstmp))')))
  AllJumps3 <- rbind(AllJumps3,AllJumpstmp2,AllJumpstmp)
}

# 结果可视化

Jumps <- matrix(0,ROW,COL)
AllJumps <- unique.matrix(round(AllJumps3))

Jumps[AllJumps] <- 1

display(Jumps,method = 'raster')

AllJumps_clear <- AllJumps

# 去掉杂点
message('######## ','Processing...',' ########','\n')

ClearPoint <- function(hmin=2,count=1){
  clearpoint <- rep(0,nrow(Jumpstmp))
  for (index in 1:nrow(Jumpstmp)){
    x <- Jumpstmp[index,1]
    y <- Jumpstmp[index,2]
    nearpoint <- matrix(c(rep(c((x-hmin):(x+hmin)),2*hmin+1),rep(c((y-hmin):(y+hmin)),each=2*hmin+1)),ncol = 2)
    clearpoint[index] <- sum(apply(nearpoint,1,
              function(x){sign(Jumpstmp[,1]==x[1]&Jumpstmp[,2]==x[2])}))
    }

  Jumpstmp <- Jumpstmp[clearpoint>count,]
  return(Jumpstmp)
}
Jumpstmp <- AllJumps
Jumpstmp <- ClearPoint(hmin = 1,count = 1)
Jumpstmp <- ClearPoint(hmin = 2,count = 2)
Jumpstmp <- ClearPoint(hmin = 3,count = 3)
Jumpstmp <- ClearPoint(hmin = 4,count = 4)
Jumpstmp <- ClearPoint(hmin = 5,count = 5)

AllJumps_clear <- Jumpstmp


# 插值
Jumps_clear <- matrix(0,ROW,COL)
Jumps_clear[AllJumps_clear] <- 1
Jumps <- Jumps_clear

Jumps_edge <- Jumps_clear
Jumps_edge[c(1,ROW),] <- 1
Jumps_edge[,c(1,COL)] <- 1

display(Jumps_clear,method = 'r')
display(Jumps_edge,method = 'r')

Judg <- function(index,h=1){
  near <- GetArea(AllJumps_clear[index,],ROW,COL,h)
  points <- near[Jumps_edge[near]==1,,drop=FALSE]
  if (nrow(points)==1) return(Judg(index,h=h+1))
  if (max(dist(points))>2*h){
    return(h)
  }else{
    return(Judg(index,h=h+1))
  }
}

H <- matrix(rep(1,2*nrow(AllJumps_clear)),ncol = 2)

if (!exists('Count')) Count <- 1

for (index in 1:nrow(AllJumps_clear)){
  H[index,] <- c(index,Judg(index,h=1))
}
H <- H[order(H[,2],decreasing = TRUE),]
H <- H[H[,2]>Count,,drop=FALSE]

if (nrow(H)>0){
  for (iii in 1:nrow(H)){
  index <- H[iii,1]
  h <- Judg(index,h=1)
  if (h>Count){
    near <- GetArea(AllJumps_clear[index,],ROW,COL,h)
    points <- near[Jumps_edge[near]==1,,drop=FALSE]
    Dist <- as.matrix(dist(points))
    vertex <- points[which(Dist==max(Dist),arr.ind = TRUE)[,1],,drop=FALSE]
    x2 <- AllJumps_clear[index,]
    for (ii in 1:nrow(vertex)){
      x1 <- vertex[ii,]
      if(abs(x1[1]-x2[1])>=abs(x1[2]-x2[2])){
        res <- do.call(cbind,approx(x=c(x1[1],x2[1]),y=c(x1[2],x2[2]),xout = x1[1]:x2[1]))
      }else{
        res <- do.call(cbind,approx(x=c(x1[2],x2[2]),y=c(x1[1],x2[1]),xout = x1[2]:x2[2]))[,2:1]
      }
      Jumps_edge[round(res)] <- 1
      Jumps[round(res)] <- 1
    }
  }
  }
}


AllJumps_clear <- which(Jumps==1,arr.ind = TRUE)
display(Jumps,method = 'raster')

save(AllJumps,Jumps,Image_raw,Image_noise,sigma,AllJumps_clear,file = './/tmp//jumps.Rdata')



# if (clear){
#   # 补充点
#   Jumps_clear <- matrix(0,ROW,COL)
#   Jumps_clear[AllJumps_clear] <- 1
#   display(Jumps_clear,method = 'raster')
#
#   Jumps_clear <- matrix(0,ROW,COL)
#   for(index in 1:nrow(AllJumps_clear)){
#     hmax <- 1
#     nearpoint <- matrix(c(-hmax:hmax,rep(0,4*hmax+2),-hmax:hmax)+rep(AllJumps_clear[index,],each=4*hmax+2),ncol = 2)
#     # nearpoint <- matrix(c(rep(c((x-hmax):(x+hmax)),2*hmax+1),rep(c((y-hmax):(y+hmax)),each=2*hmax+1)),ncol = 2)
#     nearpoint <- nearpoint[nearpoint[,1]>=1&nearpoint[,1]<=ROW&nearpoint[,2]<=COL&nearpoint[,2]>=1,,drop=FALSE]
#     Jumps_clear[nearpoint] <- 1
#     }
#
#   # Jumps_clear <- matrix(0,ROW,COL)
#   # Jumps_clear[round(AllJumps[clearpoint>1,])] <- 1
#
#   display(Jumps_clear,method = 'raster')
#   Jumps <- Jumps_clear
# }

