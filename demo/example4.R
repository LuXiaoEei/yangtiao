require(EBImage)
require(yangtiao)
options(expressions = 5e5) #设置递归最大深度



# 开始计算 --------------------------------------------------------------------

#分块主要用到，GetArea,Inita,SearPoints三个函数

## 由于函数里引用了一定数量的全局变量，所有必须保证如下变量在全局环境中存在
## Jumps,ROW,COL,ii,maxi,hax,mark,searched,ROW,COL




def <- function(Jumps,hmax=5){
  # Jumps <- Jumps
  Jumps <<- Jumps
  # Jumps <- Jumps

  # 图像边界
  ROW <<- dim(Jumps)[1]
  COL <<- dim(Jumps)[2]
  JumpsPro <<- sum(Jumps==0)/(ROW*COL)
  tag <<- 0
  h <<- 60 #Inita中最大窗宽
  maxi <<- 20000 #Searpoints中递归最大次数
  hmax <<- hmax #Searpoints中最大窗宽
  N <<- 1000 #Inita中随机取点的个数
  n <<- 5 #Inita中h递减速度

  mark <<- matrix(0,ROW,COL) #标记矩阵
  searched <<- matrix(0,ROW,COL) #搜索过的点记录矩阵

  prob <<- 1-sum(mark==0)/(ROW*COL) #已标记比例
}

# Jumps <- readImage('trypict.png')@.Data[,,1]
# Jumps[Jumps>0] <- 1

# def(Jumps = Jumps_clear,hmax = 10) #参数初始化
def(Jumps = Jumps,hmax = 3) #参数初始化

system.time(
  while (prob<JumpsPro*0.95){
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

    record <- rep(1,100)
    for (index in 1:nrow(Area)){
      post <- Area[index,,drop=FALSE]
      ii <<- 1
      tag <- mark[post]
      SearPoints(post = post,h=hmax)
      # cat('############## tag',index,'end###################','\n')
      record[index%%100+1] <- sum(mark==0)
      if (var(record)==0){ #听着准则
        break()
      }
    }
    prob <- 1-sum(mark==0)/(ROW*COL)
    print(prob)
  }
)


display(mark/max(mark),method = 'raster')
display(searched,method = 'raster')
display(Jumps,method = 'raster')
