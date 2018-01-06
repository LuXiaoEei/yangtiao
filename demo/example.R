require(plotly)
require(MASS)
require(ggplot2)
require(yangtiao)
require(rgl)

# 二维曲线——圆 -----------------------------------------------------------------
degree=3
N=7
t <- seq(0,2*pi,0.01)
x <- cos(t)
y <- sin(t)
# t <- c(1:length(t))
# BaSplite(t,15,0,u)
# i=15
# sign(t>=u[i+1]&t<u[i+2])

# u <- seq(-10,3*pi,length.out = N)
u <- c(seq(0,2*pi,length.out = N))
# u <- seq(-80,850,length.out = N)
X <- sapply(0:(length(u)+degree-2),BaSplite,x=t,degree=degree,u=u,n=length(u))

mu <- tryCatch({
  solve(t(X)%*%X)%*%t(X)},
  error=function(e){
    print(e)
    MASS::ginv(t(X)%*%X)%*%t(X) #广义逆
  }
)
betax <- mu%*%x
betay <- mu%*%y
xx <- runif(1000,0,2*pi) #产生新的点
new <- sapply(0:(length(u)+degree-2),BaSplite,x=xx,degree=degree,u=u,n=length(u))
# new <- sapply(0:(length(u)-2-n),BaSplite,x=t,n=n,u=u)
newx <- new%*%betax #产生新的点
newy <- new%*%betay #产生新的点
p <- ggplot(data.frame(newx,newy),aes(newx,newy))+
  geom_point(size=0.5,color='red')
p
pp <- p+geom_path(data=data.frame(x=betax,y=betay),aes(x,y))+ # 控制点
  geom_text(data=data.frame(x=betax,y=betay),aes(x,y),label=c(0:(length(u)+degree-2)),colour="blue",size=8)
pp
pp+geom_point(data = data.frame(x,y),aes(x,y),size=0.5)


# res <- data.frame(newx,newy)
# plot(xx,newx)

# 二维曲线——sin(t) ------------------------------------------------------------
  m=0
  N=8
  u <- seq(0,2*pi,length.out = N)
  t <- seq(0,2*pi,0.01)
  degree=3
  X <- sapply(0:(length(u)+degree-2-m),BaSplite,x=t,degree,u,length(u))
  y=sin(t)
  betay <- solve(t(X)%*%X)%*%t(X)%*%y
  betax <- solve(t(X)%*%X)%*%t(X)%*%t

  new <- X%*%betay

  ggplot(data = data.frame(x=t,y=new),aes(x,y))+geom_line()+
    geom_path(data=data.frame(t,y),aes(t,y),color='red')+
    geom_point(data=data.frame(betax,betay),aes(betax,betay))+
    geom_path(data=data.frame(betax,betay),aes(betax,betay))+
    geom_text(data=data.frame(x=betax,y=betay),aes(x,y),label=c(0:(length(u)+degree-2-m)),colour="blue",size=5)


# plot(t,new,col='red',cex=0.05)
# points(t,y,cex=0.05)
# points(betax,betay)
# sum((new-y)^2)


# 三维曲面 --------------------------------------------------------------------

# 生成对应的样条基函数
BaValue <- function(ux,uy,nx,ny,x,y){
  Bx <- sapply(0:(length(ux)+nx-2),BaSplite,x=x,degree=nx,u=ux,n=length(ux))
  By <- sapply(0:(length(uy)+ny-2),BaSplite,x=y,degree=ny,u=uy,n=length(uy))

  for (index in 1:ncol(Bx)){
    if (index==1){
      X <- Bx[,index]*By
    }else{
      X <- cbind(X,Bx[,index]*By)
    }
  }
  return(X)
}

data(geyser)
kd<-with(geyser,kde2d(duration,waiting,n=50))
dat <- data.frame(x=rep(kd$x,each=50),y=rep(kd$y,50),z=matrix(kd$z,2500,1))
Nx=20
Ny=20
nx=3
ny=3
ux=seq(-2,8,length.out = Nx)
uy=seq(30,130,length.out = Ny)


X <- BaValue(ux,uy,nx,ny,dat$x,dat$y)

mu <- tryCatch({
  solve(t(X)%*%X)%*%t(X)},
  error=function(e){
    print(e)
    MASS::ginv(t(X)%*%X)%*%t(X) #广义逆
  }
)

# 控制点
betax <- mu%*%dat$x
betay <- mu%*%dat$y
beta <- mu%*%dat$z

newx <- runif(2500,min(dat$x),max(dat$x))
newy <- runif(2500,min(dat$y),max(dat$y))
newX <- BaValue(ux,uy,nx,ny,newx,newy)
newz <- newX%*%beta

plot_ly(data=data.frame(newx,newy,newz),x=~newx,y=~newy,z=~newz,sizes=c(3,6))%>%
        add_markers(size=~c(1),color=I('#4AC6B7'),name='fitted points')%>%#拟合
        add_mesh(x=dat$x,y=dat$y,z=dat$z,opacity=0.7,color=I('#C61951'))%>%
        add_markers(data = data.frame(betax,betay,beta),
                    x=~betax,y=~betay,z=~beta,size=~c(2),
                    color=I('#1972A4'),name='control points') #加入控制点


# 三维曲面——球 -----------------------------------------------------------------

# 生成对应的样条基函数
BaValue <- function(ux,uy,nx,ny,x,y){
  Bx <- sapply(0:(length(ux)+nx-2),BaSplite,x=x,degree=nx,u=ux,n=length(ux))
  By <- sapply(0:(length(ux)+ny-2),BaSplite,x=y,degree=ny,u=uy,n=length(uy))

  for (index in 1:ncol(Bx)){
    if (index==1){
      X <- Bx[,index]*By
    }else{
      X <- cbind(X,Bx[,index]*By)
    }
  }
  return(X)
}


Fai <- seq(0,pi,length.out = 100)
Theta<- seq(0,2*pi,length.out = 100)
fai <- rep(Fai,each=length(Theta))
theta <- rep(Theta,length(Fai))
x <- sin(fai)*cos(theta)
y <- sin(fai)*sin(theta)
z <- cos(fai)
# require(rgl)
# plot3d(x,y,z,type = 'l')
Nfai=20
Ntheta=20
nfai=3
ntheta=3
ufai=seq(-2,2*pi,length.out = Nfai)
utheta=seq(-2,3*pi,length.out = Ntheta)
X <- BaValue(ufai,utheta,nfai,ntheta,fai,theta)
mu <- tryCatch({
  solve(t(X)%*%X)%*%t(X)},
  error=function(e){
    print(e)
    MASS::ginv(t(X)%*%X)%*%t(X) #广义逆
  }
)
betax <- mu%*%x
betay <- mu%*%y
betaz <- mu%*%z
# newx <- runif(2500,min(dat$x),max(dat$x))
# newy <- runif(2500,min(dat$y),max(dat$y))
# newX <- BaValue(ux,uy,nx,ny,newx,newy)
newfai <- runif(1000,0,pi)
newtheta <- runif(1000,0,2*pi)
newX <- BaValue(ufai,utheta,nfai,ntheta,newfai,newtheta)
newx <- newX%*%betax
newy <- newX%*%betay
newz <- newX%*%betaz
new <- data.frame(newx,newy,newz)
TrVal <- data.frame(x=sin(newfai)*cos(newtheta),y=sin(newfai)*sin(newtheta),z=cos(newfai))

plot_ly(data.frame(x,y,z),x=x,y=y,z=z,sizes=c(5,6))%>%
  add_markers(size=~c(1),color=I('red'),name='origin points')%>%
  add_markers(data = new,x=~newx,y=~newy,z=~newz,size=~c(1),color=I('green'),name='fitted points')%>%
  add_markers(data = data.frame(betax,betay,betaz),
              x=~betax,y=~betay,z=~betaz,size=~c(2),
              color=I('yellow'),name='control points')%>%
  add_markers(data=TrVal,x=~x,y=~y,z=~z,size=~c(1),color=I('black'),name='true points')

# require(rgl)
# plot3d(betax,betay,betaz)
# plot3d(x,y,z)
