require(EBImage)
require(yangtiao)

# 图像分块 --------------------------------------------------------------------

#分块主要用到，GetArea,Inita,SearPoints三个函数

## 由于函数里引用了一定数量的全局变量，所有必须保证如下变量在全局环境中存在
## Jumps,ROW,COL,ii,maxi,hax,mark,searched,Image_noise

## 随着噪声的增加，为了相似度判断的精确，需要减少lambda的值，
## 这样会导致更多的点无法标记，这时候要适当的减少theta值，较小的theta的值会使得一些小的快难以识别出来.

# rm(list = ls())
# gc()

# load('.//tmp//jumps.Rdata')


def1 <- function(Jumps,hmax=5,theta=0.95,lambda=3){

  Jumps <<- Jumps
  ROW <<- dim(Jumps)[1]
  COL <<- dim(Jumps)[2]
  JumpsPro <<- sum(Jumps==0)/(ROW*COL)
  tag <<- 0
  h <<- 60 #Inita中最大窗宽
  hmax <<- hmax #Searpoints中最大窗宽

  N <<- 1000 #Inita中随机取点的个数
  n <<- 5 #Inita中h递减速度
  theta <<- theta  #随着噪声的增加，为了相似度判断的精确，需要减少lambda的值，这样会导致更多的点无法标记，这时候要适当的减少theta值
  lambda <<- lambda

  searched <<- matrix(0,ROW,COL) #搜索过的点记录矩阵
  mark <<- matrix(0,ROW,COL) #标记矩阵

  prob <<- 1-sum(mark==0)/(ROW*COL) #已标记比例
} #初始化参数

# Jumps <- readImage('trypict.png')@.Data[,,1]
# Jumps[Jumps>0] <- 1

# def2(Jumps = Jumps_clear,hmax = 5)

if (!exists('hmax')){hmax <- 5}
if (exists('lambda')){rm(lambda);lambda <- 1000000}

def1(Jumps = Jumps,hmax = hmax,theta = 0.95,lambda = lambda) #初始化参数

system.time(
  while (prob<JumpsPro*theta){
    tag <- tag+1
    cat('######### tag',tag,'##########','\n')
    inint <- Inita(h = h,mark = mark,N = N,n = n)

    origpoint <- inint$origpoint
    h <- inint$h

    Area <- GetArea(point = origpoint,Row = ROW,Col = COL,h = h)
    Area <- Area[mark[Area]==0,,drop=FALSE] #未标记
    Area <- Area[sample(1:nrow(Area),nrow(Area),replace = FALSE),]#打乱顺序

    searched[origpoint] <- 1
    mark[Area] <- tag

    newpoints <- ROW*COL
    chaPro <- 0
    while (newpoints!=0){
      Areatmp <- matrix(NA,0,2)[-1,]
      for (index in 1:nrow(Area)){
        post <- Area[index,,drop=FALSE]
        Areatmp <- rbind(Areatmp,SearPoints1(post,h = hmax,lambda = lambda))
      }
      Area2 <- Area
      Area <- Areatmp[!is.na(Areatmp[,1]),,drop=FALSE]
      oldpoints <- newpoints
      newpoints <- nrow(Area)
      chaPro <- newpoints/oldpoints
      # cat(':',newpoints,'\n')
      if(chaPro>=3){
        Area22 <- Area2
        Area23 <- Area
        mark[rbind(Area,Area2)] <- 0
        newpoints <- 0
      }
      # cat(ii,':',oldpoints,'\n')
    }
    prob <- 1-sum(mark==0)/(ROW*COL)
  }
)

display(mark/max(mark),method = 'raster')
display(searched,method = 'raster')
display(Jumps,method = 'raster')

save(mark,searched,tag,lowdegree,updegree,lowerKnot,upperKnot,file = './/tmp//mark.Rdata')

# 分块去噪 --------------------------------------------------------------------
cat('######### block denoise ##########','\n')

# load('.//tmp//mark.Rdata')
# load('.//tmp//jumps.Rdata')

ModeNum <- tag
ROW <- nrow(Image_noise)
COL <- ncol(Image_noise)
AllJumps <- which(mark==0,arr.ind = TRUE)

JumpValue <- matrix(Image_noise[AllJumps],ncol = ModeNum+1,nrow = nrow(AllJumps))
# colnames(JumpValue) <- c('orig',paste0('tag_',1:ModeNum))

for (index in 1:nrow(AllJumps)){
  h <- 5
  Num <- 0
  while(Num==0){
    Area <- GetArea(AllJumps[index,],Row = ROW,Col = COL,h = h)
    Num <- sum(mark[Area]!=0)
    h <- h+5
  }
  Values <- cbind(mark[Area],Image_noise[Area])
  JumpValue[index,2:(ModeNum+1)] <- sapply(1:ModeNum,function(x){
  mean(Values[Values[,1]==x,2,drop=FALSE]) #计算邻域内多个模式的点的平均值
  # JumpValue[index,2:(ModeNum+1)] <- sapply(1:ModeNum,function(x){
  #   nrow(Values[Values[,1]==x,2,drop=FALSE]) #计算邻域内多个点的个数
  }
  )
}


JumpValue[is.nan(JumpValue)] <- Inf
res <- apply(abs(JumpValue[,-1,drop=FALSE]-JumpValue[,1]),1,which.min) #通过相似判断
# res <- apply(JumpValue[,-1,drop=FALSE],1,which.max) #通过点位置判断，比较光滑

mark[AllJumps] <- res
display(mark/max(mark),method = 'raster')
save(AllJumps,mark,Image_noise,Image_raw,sigma,ModeNum,AllJumps_clear,lowdegree,updegree,lowerKnot,upperKnot,file = './/tmp//block.Rdata')
