# 带跳曲面拟合
require(yangtiao)
require(EBImage)

if(!exists('K')) K <- c(0) #K表示旋转的程度，旋转度数的正切值K
if(!exists('updegree')) updegree <- 3

AllJumps3 <- matrix(NA,ncol = 2)[-1,]
Sigma2hat <- c()
para <- matrix(NA,0,3) #记录degre和节点个数
colnames(para) <- c('degree','knotnum','lambda')
# 跳点寻找
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
      low <- 2
      uj <- c(1:ROW)
      u1 <- c(seq(1,5,length.out = 15),seq(6,COL-5,1),seq(COL-4,COL,length.out = 15))#寻找跳点时的候选节点,加大边界寻找密度

      # 开始拟合 --------------------------------------------------------------------
      # Image_hat <- matrix(0,ROW,COL)
      AllJumps <- matrix(NA,ncol = 2)
      AllJumps <- AllJumps[-1,]

      for (index in 1:COL){
        # 选择sample
        yj <- Image_noise[,index]
        if(exists('blur')){
          Num <- length(yj)
          yj <- blur[1]*yj+blur[2]*c(yj[1],yj[-Num])+blur[3]*c(yj[-1],yj[Num])
        }

        sigma2hat<- 1/2*1/(n-1)*sum((yj[-1]-yj[-n])^2)
        Sigma2hat <- c(Sigma2hat,sigma2hat)

        BIC_list<-BIC_knot(low,up,yj,uj,n,alpha,updegree,scale=ROW) # BIC准则选择节点
        BIC_list<-GCV_knot(low,up,yj,uj,updegree,scale=ROW) # BIC准则选择节点

        u0<-BIC_list[[1]]         ##初始节点
        degree<-BIC_list[[2]]     ##初始次数
        k0<-BIC_list[[3]]      #初始节点数
        err<-BIC_list[[4]]
        lambda <-BIC_list[[5]]
        para <- rbind(para,c(degree,k0,lambda))
        print(c(degree,k0,lambda))

        jumps_find_list <- jumps_find(u0,u1,degree,yj,uj,err,sigma2hat,alpha,gamma = 0)
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
          if(exists('blur')){
            yj <- blur[1]*yj+blur[2]*c(yj[1],yj[-Num])+blur[3]*c(yj[-1],yj[Num])
          }
          n <- Num   ##设置样本量
          # 设置节点选择时的参数
          alpha <- 0.05
          low <- 4
          up <- 14
          uj <- c(1:n)
          u1 <- seq(1,n,1)#寻找跳点时的候选节点

          sigma2hat<- 1/2*1/(n-1)*sum((yj[-1]-yj[-n])^2)
          Sigma2hat <- c(Sigma2hat,sigma2hat)

          BIC_list<-BIC_knot(low,up,yj,uj,n,alpha,updegree,scale=Num) # BIC准则选择节点
          u0<-BIC_list[[1]]         ##初始节点
          degree<-BIC_list[[2]]     ##初始次数
          k0<-BIC_list[[3]]      ##初始节点数
          err<-BIC_list[[4]]

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
#

index <- as.numeric(names(sort(table(para[,1]),decreasing = TRUE)[1]))
lowdegree <- max(index-1,1)
updegree <- max(index+1,3)

index <- as.numeric(names(sort(table(para[,2]),decreasing = TRUE)[1:2]))
lowerKnot <- max(min(index-1),2)
upperKnot <- max(index+1,5)

cat('lowdegree:',lowdegree,'updegree:',updegree,'lowerKnot:',lowerKnot,'upperKnot:',upperKnot,'\n')


# 结果可视化

Jumps <- matrix(0,ROW,COL)
AllJumps <- unique.matrix(round(AllJumps3))

Jumps[AllJumps] <- 1

display(Jumps,method = 'raster')

AllJumps_clear <- AllJumps
#
# # 去掉杂点
message('######## ','Processing...',' ########','\n')
#
# 在hmin大小的邻域内,跳点数量小于等于count的点被认为是杂点
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
# Jumpstmp <- AllJumps
#
# # 在hmin大小的邻域内,跳点数量小于等于count的点被认为是杂点
# for (ii in 1:(2*length(K)+1)){
#   Jumpstmp <- ClearPoint(hmin = ii,count = ii)
# }
#
# Jumps_clear <- matrix(0,ROW,COL)
# Jumps_clear[Jumpstmp] <- 1
# display(Jumps_clear,method = 'r')
#
# # 计算跳点的h邻域内跳点的个数，跳点个数小于均值减去2倍标准差的跳点被认为是杂点
# Freq <- sapply(1:nrow(Jumpstmp),function(x){
#   return(sum(Jumps_clear[GetArea(Jumpstmp[x,],Row = ROW,Col = COL,h = 4*length(K))]))
#   })
# Jumpstmp <- Jumpstmp[Freq>=mean(Freq)-2*sd(Freq),,drop=FALSE]
#
# AllJumps_clear <- Jumpstmp

Image_smooth <- matrix(0,ROW,COL)
for(i in 1:ROW){
  for(j in 1:COL){
    Image_smooth[i,j] <- mean(Image_noise[GetArea(c(i,j),ROW,COL,Smooth)])
  }
}
display(Image_smooth,method = 'r')


sigma2hat <- mean(Sigma2hat)
clear <- rep(0,nrow(AllJumps))
for (index in 1:nrow(AllJumps)){
  orig <- AllJumps[index,]
  near <- GetArea(AllJumps[index,],ROW,COL,Smooth)
  post <- t((t(near)-orig))
  # print(post)
  W <- diag(BivKernal(x = post[,1]/ROW,y=post[,2]/COL,h=0.035),nrow = nrow(near))
  Y <- Image_smooth[near]
  X <- cbind(1,post)
  beta <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%Y
  Area1 <- near[beta[2]*post[,1]+beta[3]*post[,2]>0,,drop=FALSE]
  Area2 <- near[beta[2]*post[,1]+beta[3]*post[,2]<0,,drop=FALSE]
  clear[index] <- abs(mean(Image_smooth[Area1])-mean(Image_smooth[Area2]))/
    sqrt(sigma2hat/(2*Smooth+1)^2*(1/nrow(Area1)+1/nrow(Area2)))
  # print(sqrt(sigma2hat/(2*Smooth+1)^2*(1/nrow(Area1)+1/nrow(Area2))))
}
# tmp <- AllJumps[clear>5,,drop=FALSE]
Class <- kmeans(clear,2) #聚类
Jumpstmp <- AllJumps[Class$cluster==which.max(Class$centers),,drop=FALSE]
Jumps_clear <- matrix(0,ROW,COL)
Jumps_clear[Jumpstmp] <- 1
display(Jumps_clear,method = 'r')

Jumpstmp <- ClearPoint(hmin = 10,count = 2)
AllJumps_clear <- Jumpstmp

# 插值
Jumps_clear <- matrix(0,ROW,COL)
Jumps_clear[AllJumps_clear] <- 1
Jumps <- Jumps_clear

Jumps_edge <- Jumps_clear

# for (i in 1:nrow(AllJumps_clear)){
#   tmp <- GetArea(point = AllJumps_clear[i,],Row = ROW,Col = COL,h = round(min(ROW,COL)*0.02))
#   tmp <- tmp[tmp[,1]==1|tmp[,1]==ROW|tmp[,2]==1|tmp[,2]==COL,,drop=FALSE]
#   Jumps_edge[tmp] <- 1
# }
#

display(Jumps_clear,method = 'r')
display(Jumps_edge,method = 'r')

# 寻找最小邻域h，在这邻域内相距最远的两个跳点的距离大于2h,当h超过hmax时候，h当做无穷处理
Judg <- function(index,h=1,hmax){
  if(h>hmax){return(Inf)}
  near <- GetArea(AllJumps_clear[index,],ROW,COL,h)
  points <- near[Jumps_edge[near]==1,,drop=FALSE]
  if (nrow(points)==1) return(Judg(index,h=h+1,hmax=hmax))
  if (max(dist(points))>2*h){
    return(h)
  }else{
    return(Judg(index,h=h+1,hmax=hmax))
  }
}

H <- matrix(rep(1,2*nrow(AllJumps_clear)),ncol = 2)

if (!exists('Count')) Count <- 1
if (!exists('gamma')) gamma <- 0.25

for (index in 1:nrow(AllJumps_clear)){
  H[index,] <- c(index,Judg(index,h=1,hmax=gamma*min(COL,ROW)))
}
H <- H[order(H[,2],decreasing = TRUE),]

# 将邻域集合H分为无穷和非无穷，并分开进行插值
Hinf <- H[is.infinite(H[,2]),,drop=FALSE]
H <- H[H[,2]>Count&!is.infinite(H[,2]),,drop=FALSE]

if (nrow(H)>0){
  for (iii in 1:nrow(H)){
  index <- H[iii,1]
  h <- Judg(index,h=1,hmax=gamma*min(COL,ROW))
  if (h>Count){
    near <- GetArea(AllJumps_clear[index,],ROW,COL,h)
    points <- near[Jumps_edge[near]==1,,drop=FALSE]
    Dist <- as.matrix(dist(points))
    vertex <- points[which(Dist==max(Dist),arr.ind = TRUE)[,1],,drop=FALSE] #邻域内相距最远的两个点
    x2 <- AllJumps_clear[index,]
    for (ii in 1:nrow(vertex)){ #中心点分别于上述两个点之间进行插值
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

display(Jumps,method = 'raster')

for (i in 1:nrow(AllJumps_clear)){
  tmp <- GetArea(point = AllJumps_clear[i,],Row = ROW,Col = COL,h = round(min(ROW,COL)*0.05))
  tmp <- tmp[tmp[,1]==1|tmp[,1]==ROW|tmp[,2]==1|tmp[,2]==COL,,drop=FALSE]
  Jumps_edge[tmp] <- 1
}

# Jumps_edge[c(1,ROW),] <- 1
# Jumps_edge[,c(1,COL)] <- 1
if(nrow(Hinf)>0){
  for (iii in 1:nrow(Hinf)){
    index <- Hinf[iii,1]
    h <- Judg(index,h=1,hmax = gamma*min(COL,ROW))
    if (h>Count&!is.infinite(h)){
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

display(Jumps,method = 'raster')
AllJumps_clear <- which(Jumps==1,arr.ind = TRUE)


save(AllJumps,Jumps,Image_raw,Image_noise,sigma,AllJumps_clear,lowdegree,updegree,lowerKnot,upperKnot,file = './/tmp//jumps.Rdata')



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

