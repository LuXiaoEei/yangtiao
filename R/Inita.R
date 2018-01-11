# 初始化
# 寻找一个尽量大的的窗宽对应的起始点,从给定的最大h开始寻找，随机寻找N个点，
# 如果N个点的h邻域内都没有满足条件的点，则减少h继续寻找，直到满足条件
# 条件：邻域内没有跳点，
# 返回: 邻域的大小和初始点，
# h：最大的窗宽，N：每个窗宽内随机点的个数，Jumps，带跳的0-1矩阵（1代表跳）;n:每次窗宽减少的量
# mark：不作为初始点的矩阵，默认是0矩阵，对于分图任务中，mark表示已经标记了分类标签的点
Inita <- function(h,mark,N,n=5){
  if (missing(mark)){
    mark <- Jumps
    mark <- 0
  }
  origpoint <- matrix(c(sample(1:ROW,N,replace = TRUE),sample(1:COL,N,replace = TRUE)),ncol  = 2) #61,180
  origpoint <- origpoint[mark[origpoint]==0,,drop=FALSE]
  CountJumps <- apply(origpoint,1,function(x){
    Area <- GetArea(point = x,Row = ROW,Col = COL,h = h)
    return(sum(Jumps[Area]))
  })
  if(sum(CountJumps==0)>0){
    return(list(origpoint=origpoint[match(TRUE,CountJumps==0),,drop=FALSE],h=h))
  }else{
    return(Inita(h-n,mark,N))
  }
}
