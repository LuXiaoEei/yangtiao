if(FALSE) source('.//demo//BSplineSurface.R')

require(EBImage)
require(yangtiao)

rm(list=ls())
gc()

set.seed(1234)

if (dir.exists('.//tmp')) unlink('.//tmp', recursive=TRUE)
dir.create('.//tmp')

# 构造数据 --------------------------------------------------------------------

f2 <- function(x,y){
  z <- 0.25*(1-x)*y+(1+0.2*sin(2*pi*x))*sign(y>0.6*sin(pi*x)+0.2)
  return(z)
}


# 定义图像尺寸
ROW <- 256
COL <- 256

sigma <- 0.2 #噪声大小

Image_raw <- matrix(f2(x=rep(c(1:ROW),COL)/ROW,y=rep(c(1:COL),each=ROW)/COL),nrow = ROW,ncol = COL)
# Image_raw <- Image_raw/max(Image_raw)
Image_noise <- Image_raw+matrix(rnorm(ROW*COL,0,sigma),nrow = ROW,ncol = COL) #加噪

display(Image_raw,method = 'raster')
display(Image_noise,method = 'raster')

# 跳点检测 --------------------------------------------------------------------
message('######## ','Detect Jumps',' ########','\n')
startime <- Sys.time()

K <- c(0)
Count <- 1
source('.//demo//searchjump.R')
print(Sys.time()-startime)


# 图形分块 --------------------------------------------------------------------
message('######## ','Block Image',' ########','\n')
# rm(list=ls())
# gc()
startime <- Sys.time()

hmax <- 5
lambda <- 3
source('.//demo//block.R')
print(Sys.time()-startime)


# 分块拟合 --------------------------------------------------------------------
message('######## ','Fit Image',' ########','\n')
# rm(list=ls())
# gc()
startime <- Sys.time()

denoise <- FALSE
lowerKnot <- 2
upperKnot <- 8
source('.//demo//fit.R')
print(Sys.time()-startime)


