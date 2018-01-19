## 通过递归的方式计算B样条基
## x:value i:ith basis  degree:degree(次数,从0开始) u:knots i:0-max(i)
## n: the length of origin knots
## !note length(u)-1=max(i)+1+n
## 总节点个数(length(u)+2*degree)=样条基个数(max(i)+1)+阶数(degree+1)
## max(i)=length(u)+degree-2
## 样条基个数(max(i)+1)=内节点个数(length(u)-2)+阶数(degree+1)
BaSplite <- memoise::memoise(function(x,i,degree,u,n){
  if (degree==0){
    return(sign(x>=u[i+1])*sign(x<u[i+2]))
    }else{
      if (length(u)==n){
        u <- c(rep(u[1],degree),u,rep(u[length(u)],degree))
      }
      return(div(BaSplite(x,i,degree-1,u,n),1/div(x-u[i+1],u[i+degree+1]-u[i+1]))+
               div(BaSplite(x,i+1,degree-1,u,n),1/div(u[i+degree+2]-x,u[i+degree+2]-u[i+2])))
    }
})

