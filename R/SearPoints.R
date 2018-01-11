# 在不跨越跳点的时候寻找同一区域的点并标记的树递归算法
# h：搜索窗宽;H:最大窗宽，post；搜索中心，Jumps：跳0-1矩阵（1表示跳点）
# ii和maxi表示计数，防止递归太深；ROW，COL：图像边界;searched:记录搜索过的点
SearPoints <- function(post,h){
  # cat('ii:',ii,'\n')
  if(ii < maxi){
    ii <<- ii+1
    # cat('mark:',sum(mark==tag),'\n')
    # cat('searched:',sum(searched==1),'\n')
    if(searched[post]==0){ #确认post点未搜索过
      Area <- GetArea(point = post,Row = ROW,Col = COL,h = h)
      tmpmark <- mark[Area]
      Area <- Area[tmpmark==0,,drop=FALSE] #去掉已经标记的点
      if(sum(Jumps[Area])>0){#保证h范围内没有跳点，如果有就减少h
        SearPoints(post,h-1)
      }else{
        if(nrow(Area)>0){#保证还有未标记的点
          tag <- mark[post]
          mark[Area] <<- tag #标记
          searched[post] <<- 1 #加入搜索过矩阵
          Area <- Area[sample(1:nrow(Area),nrow(Area),replace = FALSE),,drop=FALSE] #打乱顺序
          for (index in 1:nrow(Area)){ #对标记的每个点再次展开算法
            post <- Area[index,,drop=FALSE]
            SearPoints(post,h=hmax)
          }
        }
      }
      searched[post] <<- 1
    }
    # cat('searched:',sum(searched==1),'\n')
  }
}
