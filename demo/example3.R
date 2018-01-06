# 带跳曲面拟合
require(yangtiao)
require(EBImage)

f2 <- function(x,y){
  z <- 0.25*(1-x)*y+(1+0.2*sin(2*pi*x))*sign(y>0.6*sin(pi*x)+0.2)
  return(z)
}
# 定义图像尺寸
ROW <- 200
COL <- 200
sigma=0.1 #噪声大小

set.seed(111)
Image_raw <- matrix(f2(x=rep(c(1:ROW),COL)/ROW,y=rep(c(1:COL),each=ROW)/COL),nrow = ROW,ncol = COL)
Image_raw <- Image_raw/max(Image_raw)
Image_noise <- Image_raw+matrix(rnorm(ROW*COL,0,sigma),nrow = ROW,ncol = COL) #加噪

display(Image_raw,method = 'raster')
display(Image_noise,method = 'raster')


# 设置初始值 默认处理一列
if (FALSE) {a <- COL;COL <- ROW;ROW <- a;Image_noise <- t(Image_noise)}
n<- ROW   ##设置样本量
# 设置节点选择时的参数
alpha<- 0.05
low<- 8
up<- 16
uj=c(1:ROW)
u1<- seq(1,COL,1) #寻找跳点时的候选节点

# 开始拟合 --------------------------------------------------------------------
Image_hat <- matrix(0,ROW,COL)
AllJumps <- matrix(NA,ncol = 2)
AllJumps <- AllJumps[-1,]

system.time(
  for (index in 1:COL){
      cat('Simulation ',index,'\n')
      # 选择sample
      yj <- Image_noise[,index]

      sigma2hat<- 1/2*1/(n-1)*sum((yj[-1]-yj[-n])^2)

      BIC_list<-BIC_knot(low,up,yj,uj,n,alpha,scale=ROW) # BIC准则选择节点
      u0<-BIC_list[[1]]         ##初始节点
      nu0<-BIC_list[[2]]        ##初始节点数
      degree<-BIC_list[[3]]     ##初始次数
      k0<-BIC_list[[4]]      ##初始内节点数
      err<-BIC_list[[5]]

      cat('k0=',k0,'degree=',degree,'\n')

      ###############拟合插入跳点jumps_option之后的模型
      jumps_find_list<-jumps_find(u0,u1,degree,yj,uj,err,sigma2hat,alpha)
      u0<-jumps_find_list[[1]]

      jumps2<-jumps_find_list[[2]]
      # fit1_1<- lm(yj_1 ~ bs(uj, degree = degree_1, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
      # fit1_1<- lm(yj_1 ~ sapply(0:(length(u0)-2+degree_1),BaSplite,x=uj,degree=degree_1,u=u0,n=length(u0))-1)
      fit1<- lm(yj~ BaSplite1(uj,degree,u0)-1)
      yjhat<- fitted(fit1)
      Image_hat[,index] <- yjhat

      cat('jumps',c(-1,jumps2),'\n')

      if (length(jumps2)>0){
        AllJumps <- rbind(AllJumps,matrix(c(rep(index,length(jumps2)),jumps2),nrow = length(jumps2)))
        }
      }
  )

# 结果可视化
Jumps <- matrix(0,ROW,COL)
Jumps[round(AllJumps)] <- 1
display(t(Jumps),method = 'raster')
display(Image_hat,method = 'raster')
# display(t(Image_hat))
# display((Image_hat2+t(Image_hat))/2)
# Image_hat1_1 <- Image_hat
# AllJumps1_1 <- AllJumps

# save(AllJumps0_1,AllJumps0_2,AllJumps1_1,AllJumps1_2,
#       Image_hat0_1,Image_hat0_2,Image_hat1_1,Image_hat1_2,
#       file = 'result.Rdata')
