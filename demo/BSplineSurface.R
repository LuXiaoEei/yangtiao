if(FALSE) source('.//demo//BSplineSurface.R')

require(EBImage)
require(yangtiao)

# rm(list=ls())
# gc()


if (dir.exists('.//tmp')) unlink('.//tmp', recursive=TRUE)
dir.create('.//tmp')

# 构造数据 --------------------------------------------------------------------

f2 <- function(x,y){
  z <- 0.25*(1-x)*y+(1+0.2*sin(2*pi*x))*sign(y>0.6*sin(pi*x)+0.2)
  return(z)
}

# f2 <- function(x,y){
#   z <- 1/3*sign(y>0.3*sin(x*pi)+0.4)+1/4*sign(x>0.2*sin(y*pi)+0.5)+(1-x)*y*5/12
#   return(z)
# }
#
# f2 <- function(x,y){
#   z <- 1/2*sign(y>0.6*sin(pi*x)+0.2)+1/2*sign(x>0.6*sin(pi*y)+0.2)
#   return(z)
# }
# f2 <- function(x,y){
#   z <- -2*(x-0.5)^2-2*(y-0.5)^2+sign((x-0.5)^2+(y-0.5)^2<0.25^2)
#   return(z)
# }
#
# f2 <- function(x,y){
#   z <-cos(4*pi*(1-x-y))-2*cos(4*pi*(1-x-y))*sign(x+y-1>0)
#   return(z)
# }

# 定义图像尺寸
ROW <- 256
COL <- 256

sigma <- 0.5 #噪声大小

# set.seed(1234)

Image_raw <- matrix(f2(x=rep(c(1:ROW),COL)/ROW,y=rep(c(1:COL),each=ROW)/COL),nrow = ROW,ncol = COL)
# Image_raw <- Image_raw/max(Image_raw)
Image_noise <- Image_raw+matrix(rnorm(ROW*COL,0,sigma),nrow = ROW,ncol = COL) #加噪

display(Image_raw,method = 'raster')
display(Image_noise,method = 'raster')

# 跳点检测 --------------------------------------------------------------------
message('######## ','Detect Jumps',' ########','\n')
startime <- Sys.time()

if(exists('blur')) rm(blur)
# 寻找跳点中节点和次数的上限设置
up <- 10
updegree <- 3

Smooth <- 2
K <- c(0)
# if(length(K)>2) blur <- c(1/3,1/3,1/3)
Count <- 1
gamma <- 0.15
source('.//demo//searchjump.R')
print(Sys.time()-startime)


# 图形分块 --------------------------------------------------------------------
message('######## ','Block Image',' ########','\n')
# rm(list=ls())
# gc()
startime <- Sys.time()

hmax <- 5
source('.//demo//block.R')
print(Sys.time()-startime)


# 分块拟合 --------------------------------------------------------------------
message('######## ','Fit Image',' ########','\n')
# rm(list=ls())
# gc()
startime <- Sys.time()

# denoise <- FALSE

source('.//demo//fit.R')
print(Sys.time()-startime)
