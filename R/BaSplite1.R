## 计算B样条基的非递归版本
## x:value   degree:degree(次数,从0开始) u:knots
## 总节点个数(length(u)+2*degree)=样条基个数(max(i)+1)+阶数(degree+1)
## 样条基个数=内节点个数(length(u)-2)+阶数(degree+1)
BaSplite1 <- function(x,degree,u){
  u <- c(rep(u[1],degree),u,rep(u[length(u)],degree))
  Bas <- sapply(1:(length(u)-1),function(i) sign(x>=u[i])*sign(x<u[i+1])) #0次
  if (degree>0){
    UpdBas <- function(index){# 更新样条基
      div(x-u[index],u[index+p]-u[index])*Bas[,index]+
        div(u[index+p+1]-x,u[index+p+1]-u[index+1])*Bas[,index+1]
    }
    for (p in 1:degree){
      Bas <- sapply(1:(length(u)-p-1),UpdBas)
      # Bas <- matrix(0,length(x),length(u)-i)
    }
  }
  return(Bas)
}

