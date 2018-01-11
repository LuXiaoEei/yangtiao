# 带跳曲面拟合
require(yangtiao)
require(EBImage)

rm(list = ls())
gc()

f2 <- function(x,y){
  z <- 0.25*(1-x)*y+(1+0.2*sin(2*pi*x))*sign(y>0.6*sin(pi*x)+0.2)
  return(z)
}
# 定义图像尺寸
ROW <- 200
COL <- 200
sigma=0.1 #噪声大小

# set.seed(111)
Image_raw <- matrix(f2(x=rep(c(1:ROW),COL)/ROW,y=rep(c(1:COL),each=ROW)/COL),nrow = ROW,ncol = COL)
Image_raw <- Image_raw/max(Image_raw)
Image_noise <- Image_raw+matrix(rnorm(ROW*COL,0,sigma),nrow = ROW,ncol = COL) #加噪

display(Image_raw,method = 'raster')
display(Image_noise,method = 'raster')

system.time(
  for (Index in c(FALSE,TRUE)){
    cat(Index,'\n')
  # 设置初始值 默认处理一列
    if (Index) {a <- COL;COL <- ROW;ROW <- a;Image_noise <- t(Image_noise)}
    n <- ROW   ##设置样本量
    # 设置节点选择时的参数
    alpha <- 0.05
    low <- 8
    up <- 16
    uj <- c(1:ROW)
    u1 <- c(seq(1,5,length.out = 25),seq(6,COL-5,1),seq(COL-4,COL,length.out = 25))#寻找跳点时的候选节点,加大边界寻找密度

    # 开始拟合 --------------------------------------------------------------------
    # Image_hat <- matrix(0,ROW,COL)
    AllJumps <- matrix(NA,ncol = 2)
    AllJumps <- AllJumps[-1,]

    for (index in 1:COL){
        # cat('Simulation ',index,'\n')
        # 选择sample
        yj <- Image_noise[,index]
        sigma2hat<- 1/2*1/(n-1)*sum((yj[-1]-yj[-n])^2)

        BIC_list<-BIC_knot(low,up,yj,uj,n,alpha,scale=ROW) # BIC准则选择节点
        u0<-BIC_list[[1]]         ##初始节点
        nu0<-BIC_list[[2]]        ##初始节点数
        degree<-BIC_list[[3]]     ##初始次数
        k0<-BIC_list[[4]]      ##初始内节点数
        err<-BIC_list[[5]]

        # cat('k0=',k0,'degree=',degree,'\n')

        jumps_find_list<-jumps_find(u0,u1,degree,yj,uj,err,sigma2hat,alpha,gamma = 0.05)
        u0<-jumps_find_list[[1]]
        jumps2<-jumps_find_list[[2]]

        # fit1<- lm(yj~ BaSplite1(uj,degree,u0)-1)
        # yjhat<- fitted(fit1)
        # Image_hat[,index] <- yjhat

        # cat('jumps',c(-1,jumps2),'\n')

        if (length(jumps2)>0){
          AllJumps <- rbind(AllJumps,matrix(c(rep(index,length(jumps2)),jumps2),nrow = length(jumps2)))
          }
        }
    if(!Index){
      AllJumpstmp <- AllJumps[,2:1]
    }else{
      AllJumpstmp2 <- AllJumps
    }
  }
)


# 结果可视化

Jumps <- matrix(0,ROW,COL)
AllJumps <- round(rbind(AllJumpstmp2,AllJumpstmp))

Jumps[AllJumps] <- 1

display(Jumps,method = 'raster')

save(AllJumps,Jumps,Image_raw,Image_noise,sigma,file = './/data//jumps.Rdata')

#
# # 去掉杂点
# clearpoint <- rep(0,nrow(AllJumps))
# for (index in 1:nrow(AllJumps)){
#   hmin <- 1
#   x <- AllJumps[index,1]
#   y <- AllJumps[index,2]
#   nearpoint <- matrix(c(rep(c(x-hmin,x,x+hmin),2*hmin+1),rep(c(y-hmin,y,y+hmin),each=2*hmin+1)),ncol = 2)
#   nearpoint <-  nearpoint[nearpoint[,1]>=1&nearpoint[,1]<=ROW&nearpoint[,2]<=COL&nearpoint[,2]>=1,,drop=FALSE]
#   clearpoint[index] <- sum(Jumps[nearpoint])
# }
#
#
# # 补充点
# AllJumps2 <- AllJumps[clearpoint>1,]
#
# Jumps_clear <- matrix(0,ROW,COL)
# Jumps_clear[AllJumps2] <- 1
# display(Jumps_clear,method = 'raster')
#
# Jumps_clear <- matrix(0,ROW,COL)
# for(index in 1:nrow(AllJumps2)){
#   hmax <- 1
#   x <- AllJumps2[index,1]
#   y <- AllJumps2[index,2]
#   nearpoint <- matrix(c(rep(c(x-hmax,x,x+hmax),2*hmax+1),rep(c(y-hmax,y,y+hmax),each=2*hmax+1)),ncol = 2)
#   nearpoint <- nearpoint[nearpoint[,1]>=1&nearpoint[,1]<=ROW&nearpoint[,2]<=COL&nearpoint[,2]>=1,,drop=FALSE]
#   Jumps_clear[nearpoint] <- 1
#   }
#
# # Jumps_clear <- matrix(0,ROW,COL)
# # Jumps_clear[round(AllJumps[clearpoint>1,])] <- 1
#
# display(Jumps_clear,method = 'raster')

