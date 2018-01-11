# 返回给定点正方形窗宽邻域点坐标; point:给定点，Row,Col:边界,h:窗宽
GetArea <- function(point,Row,Col,h){
  x <- point[1]
  y <- point[2]
  Area <- matrix(c(rep((x-h):(x+h),2*h+1),
                   rep((y-h):(y+h),each=2*h+1)),ncol = 2)
  Area <- Area[Area[,1]%in%c(1:Row)&Area[,2]%in%c(1:Col),,drop=FALSE]
  return(Area)
}
