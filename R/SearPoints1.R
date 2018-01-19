# 在不跨越跳点的时候寻找同一区域的点并标记的树递归算法
# h：搜索窗宽;H:最大窗宽，post；搜索中心，Jumps：跳0-1矩阵（1表示跳点）
# ii和maxi表示计数，防止递归太深；ROW，COL：图像边界;searched:记录搜索过的点
SearPoints1 <- function(post,h,lambda=2){
    if(searched[post]==0){ #确认post点未搜索过
      Area <- GetArea(point = post,Row = ROW,Col = COL,h = h)
      tmpmark <- mark[Area]
      if(sum(Jumps[Area])>0){#保证h范围内没有跳点，如果有就减少h
        return(SearPoints1(post,h-1,lambda = lambda))
      }else{
        tag <- mark[post]
        Area_marked <- Area[tmpmark==tag,,drop=FALSE] #取出已经标记的点
        Area <- Area[tmpmark==0,,drop=FALSE] #去掉已经标记的点
        if(nrow(Area)>0){#保证还有未标记的点
          if(nrow(Area_marked)>2){  #相似度判断是否越界
            orig <- Image_noise[Area_marked]
            Max <- mean(orig)+lambda*sd(orig)
            Min <- mean(orig)-lambda*sd(orig)
            new <- Image_noise[Area]
            Area <- Area[new<Max&new>Min,,drop=FALSE]
          }
          mark[Area] <<- tag #标记
          searched[post] <<- 1 #加入搜索过矩阵
          # Area <- Area[sample(1:nrow(Area),nrow(Area),replace = FALSE),]
          # print(Area)
          return(Area)
        }
      }
      searched[post] <<- 1
    }
  return(matrix(NA,0,2))
}

