## x:value i:ith basis  n:degree u:knots i:0-max(i)
## !note length(u)-1=max(i)+1+n
BaSplite <- function(x,i,n,u){
  if (n==0){
    return(sign(x>=u[i+1]&x<u[i+2]))
    }else{
      return(div(BaSplite(x,i,n-1,u),1/div(x-u[i+1],u[i+n+1]-u[i+1]))+
               div(BaSplite(x,i+1,n-1,u),1/div(u[i+n+2]-x,u[i+n+2]-u[i+2])))
    }
}


