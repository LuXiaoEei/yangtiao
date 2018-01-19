# 生成曲面的的样条基函数两两相乘的矩阵
# Knot:节点向量 Degree:次数  Value:求样条的点的位置值
BaValue <- function(KnotX,KnotY,DegreeX,DegreeY,ValueX,ValueY){
  Bx <- BaSplite1(x = ValueX,degree = DegreeX,u = KnotX)
  By <- BaSplite1(x = ValueY,degree = DegreeY,u = KnotY)
  for (index in 1:ncol(Bx)){
    if (index==1){
      X <- Bx[,index]*By
    }else{
      X <- cbind(X,Bx[,index]*By)
    }
  }
  return(X)
}
