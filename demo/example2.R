require(yangtiao)
## 带有限跳曲线拟合
# 函数列表
# 设置带跳函数的具体形式
f1<-function(u)
{ k<-length(u)
f1<-rep(0,k)
for (i in 1:k) {
  if (u[i]<0.3) {
    f1[i]<--3*u[i]+2
  }else {
    if(u[i]>=0.7) {
      f1[i]<-u[i]/2+1.55
    }else {
      f1[i]<--3*u[i]+3-sin((u[i]-0.3)*pi/0.2)
    }
  }
}
return(f1)
}

f2<-function(u)
{ k<-length(u)
f2<-rep(0,k)
for (i in 1:k) {
  if (u[i]<0.3) {
    f2[i]<-0
  }else {
    if(u[i]>=0.7) {
      f2[i]<-4*u[i]^2+1.24
    }else {
      f2[i]<-3*u[i]^2+0.93
    }
  }
}
return(f2)
}

f3<-function(u)
{ k<-length(u)
f3<-rep(0,k)
for (i in 1:k) {
  if (u[i]<0.5 & u[i]>=0) {
    f3[i]<-cos(8*pi*(0.5-u[i]))
  }
  if(u[i]>=0.5 & u[i]<=1) {
    f3[i]<- -cos(8*pi*(0.5-u[i]))
  }
}

return(f3)
}

# plot(ufix,f1(ufix),'l')
# plot(ufix,f2(ufix),'l')
# plot(ufix,f3(ufix),'l')

# 主程序 ---------------------------------------------------------------------
# 设置初值

n<- 200   ##设置样本量
sigmaeps_1<- 0.1
sigmaeps_2<-0.2
sigmaeps_4<-0.4##设置误差项的标准差
N<- 10   # 设置模拟次数
delta<-0.001
# 设置节点选择时的参数
alpha<- 0.05
low<- 8
up<- 16
ufix<- seq(0,1,length.out = 201)
yfix_1<-yfix_2<-yfix_4<- matrix(0,201,N)
ISE_1<- ISE3_1<- ISE7_1<-lambdaM_1_1<-lambdaM_2_1<- rep(0,N)
ISE_2<- ISE3_2<- ISE7_2<-lambdaM_1_2<-lambdaM_2_2<- rep(0,N)
ISE_4<- ISE3_4<- ISE7_4<-lambdaM_1_4<-lambdaM_2_4<- rep(0,N)
rt<- more<- less<- error<-crt<- 0
Nt<- 1
jumps.list_1<-jumps.list_2<- jumps.list_4<- list()
jumps_1<-jumps_2<-jumps_4<- c()

# 开始计算
system.time(
  while(Nt<N+1){

    cat('Simulation ', Nt,'\n')
    ####   设置取点步长
    u1<- seq(0,1,by=0.001)
    DRSStotal<- matrix(rep(0),length(u1),1)
    # 生成sample
    eps_1<- rnorm(n,0,sigmaeps_1)  ##Normal random error with sample size 'n' and parameters '0','sigmaeps'
    eps_2<- rnorm(n,0,sigmaeps_2)  ##Normal random error with sample size 'n' and parameters '0','sigmaeps'
    eps_4<- rnorm(n,0,sigmaeps_4)
    ##Normal random error with sample size 'n' and parameters '0','sigmaeps'
    uj<- sort(runif(n,0,1))  ##Uniform random design with sample size 'n' on internal [0,1]

    yj_1<- f1(uj)+eps_1      ## random response
    yj_2<- f1(uj)+eps_2      ## random response
    yj_4<- f1(uj)+eps_4      ## random response
    #估计sigma2

    # plot(uj,yj_1)
    # plot(uj,yj_2)
    # plot(uj,yj_4)

    sigma2hat_1<- 1/2*1/(n-1)*sum((yj_1[-1]-yj_1[-n])^2)
    sigma2hat_2<- 1/2*1/(n-1)*sum((yj_2[-1]-yj_2[-n])^2)
    sigma2hat_4<- 1/2*1/(n-1)*sum((yj_4[-1]-yj_4[-n])^2)
    cat('sigma2hat_1',sigma2hat_1,'sigma2hat_2',sigma2hat_2,'sigma2hat_4',sigma2hat_4,'\n')


    BIC_list_1<-BIC_knot(low,up,yj_1,uj,n,alpha) ####### sigma=0.1 BIC准则选择节点结果
    u0_1<-BIC_list_1[[1]]         ##初始节点
    nu0_1<-BIC_list_1[[2]]        ##初始节点数
    degree_1<-BIC_list_1[[3]]     ##初始自由度
    k0_1<-BIC_list_1[[4]]      ##初始内节点数
    err_1<-BIC_list_1[[5]]


    BIC_list_2<-BIC_knot(low,up,yj_2,uj,n,alpha) #######sigma=0.2 BIC准则选择节点结果
    u0_2<-BIC_list_2[[1]]
    nu0_2<-BIC_list_2[[2]]
    degree_2<-BIC_list_2[[3]]
    k0_2<-BIC_list_2[[4]]
    err_2<-BIC_list_2[[5]]


    BIC_list_4<-BIC_knot(low,up,yj_4,uj,n,alpha) ####### sigma=0.4 BIC准则选择节点结果
    u0_4<-BIC_list_4[[1]]
    nu0_4<-BIC_list_4[[2]]
    degree_4<-BIC_list_4[[3]]
    k0_4<-BIC_list_4[[4]]
    err_4<-BIC_list_4[[5]]


    cat('k0_1=',k0_1,'degree_1=',degree_1,'k0_2=',k0_2,'degree_2=',degree_2,  'k0_4=',k0_4,'degree_4=',degree_4,'\n')

    ###############拟合插入跳点jumps_option之后的模型
    jumps_find_list_1<-jumps_find(u0_1,u1,degree_1,yj_1,uj,err_1,sigma2hat_1,alpha,lambdaM_1_1,lambdaM_2_1)
    u0<-jumps_find_list_1[[1]]

    jumps2<-jumps_find_list_1[[2]]
    # fit1_1<- lm(yj_1 ~ bs(uj, degree = degree_1, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
    # fit1_1<- lm(yj_1 ~ sapply(0:(length(u0)-2+degree_1),BaSplite,x=uj,degree=degree_1,u=u0,n=length(u0))-1)
    fit1_1<- lm(yj_1 ~ BaSplite1(uj,degree_1,u0)-1)

    y1_1<- fitted(fit1_1)
    yfix_1[,Nt]<- interp(y1_1,uj,ufix)

    # jumps.list[[Nt]]<- jumps_option
    # jumps<-c(jumps,jumps_option)
    jumps.list_1[[Nt]]<-c(-1,jumps2)

    jumps_1<-c(jumps_1,jumps2)
    cat('jumps.list_1',jumps.list_1[[Nt]],'\n')

    #计算MISE
    ISE_1[Nt]<- sum((y1_1-f1(uj))^2)/n
    id1<- which((uj<0.35)&(uj>0.25))
    ISE3_1[Nt]<- sum((y1_1[id1]-f1(uj[id1]))^2)/n
    id2<- which((uj<0.75)&(uj>0.65))
    ISE7_1[Nt]<- sum((y1_1[id2]-f1(uj[id2]))^2)/n
    ###############拟合插入跳点jumps_option之后的模型
    jumps_find_list_2<-jumps_find(u0_2,u1,degree_2,yj_2,uj,err_2,sigma2hat_2,alpha,lambdaM_1_2,lambdaM_2_2)
    u0<-jumps_find_list_2[[1]]
    jumps2<-jumps_find_list_2[[2]]
    # fit1_2<- lm(yj_2 ~ bs(uj, degree = degree_2, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
    # fit1_2<- lm(yj_2 ~ sapply(0:(length(u0)-2+degree_2),BaSplite,x=uj,degree=degree_2,u=u0,n=length(u0))-1)
    fit1_2<- lm(yj_2 ~ BaSplite1(uj,degree_2,u0)-1)


    y1_2<- fitted(fit1_2)
    yfix_2[,Nt]<- interp(y1_2,uj,ufix)

    # jumps.list[[Nt]]<- jumps_option
    # jumps<-c(jumps,jumps_option)
    jumps.list_2[[Nt]]<-c(-1,jumps2)

    jumps_2<-c(jumps_2,jumps2)
    cat('jumps.list_2',jumps.list_2[[Nt]],'\n')

    #计算MISE
    ISE_2[Nt]<- sum((y1_2-f1(uj))^2)/n
    id1<- which((uj<0.35)&(uj>0.25))
    ISE3_2[Nt]<- sum((y1_2[id1]-f1(uj[id1]))^2)/n
    id2<- which((uj<0.75)&(uj>0.65))
    ISE7_2[Nt]<- sum((y1_2[id2]-f1(uj[id2]))^2)/n

    ###############拟合插入跳点jumps_option之后的模型
    jumps_find_list_4<-jumps_find(u0_4,u1,degree_4,yj_4,uj,err_4,sigma2hat_4,alpha,lambdaM_1_4,lambdaM_2_4)
    u0<-jumps_find_list_4[[1]]
    jumps2<-jumps_find_list_4[[2]]
    # fit1_4<- lm(yj_4 ~ bs(uj, degree = degree_4, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
    # fit1_4<- lm(yj_4 ~ sapply(0:(length(u0)-2+degree_4),BaSplite,x=uj,degree=degree_4,u=u0,n=length(u0))-1)
    fit1_4<- lm(yj_4 ~ BaSplite1(uj,degree_4,u0)-1)


    y1_4<- fitted(fit1_4)
    yfix_4[,Nt]<- interp(y1_4,uj,ufix)

    # jumps.list[[Nt]]<- jumps_option
    # jumps<-c(jumps,jumps_option)
    jumps.list_4[[Nt]]<-c(-1,jumps2)

    jumps_4<-c(jumps_4,jumps2)
    cat('jumps.list_4',jumps.list_4[[Nt]],'\n')

    #计算MISE
    ISE_4[Nt]<- sum((y1_4-f1(uj))^2)/n
    id1<- which((uj<0.35)&(uj>0.25))
    ISE3_4[Nt]<- sum((y1_4[id1]-f1(uj[id1]))^2)/n
    id2<- which((uj<0.75)&(uj>0.65))
    ISE7_4[Nt]<- sum((y1_4[id2]-f1(uj[id2]))^2)/n

    Nt<- Nt+1
  }
)


# 绘图 ----------------------------------------------------------------------

MISE_1<- mean(ISE_1)
MISE3_1<- mean(ISE3_1)
MISE7_1<- mean(ISE7_1)
cat('MISE_1=',MISE_1,'MISE0.3_1=',MISE3_1,'MISE0.7_1=',MISE7_1,'LMISE_1=',MISE3_1+MISE7_1,'\n')

MISE_2<- mean(ISE_2)
MISE3_2<- mean(ISE3_2)
MISE7_2<- mean(ISE7_2)
cat('MISE_2=',MISE_2,'MISE0.3_2=',MISE3_2,'MISE0.7_2=',MISE7_2,'LMISE_2=',MISE3_2+MISE7_2,'\n')

MISE_4<- mean(ISE_4)
MISE3_4<- mean(ISE3_4)
MISE7_4<- mean(ISE7_4)
# cat('Mean lambda=',mean(lambdaM))
cat('MISE_4=',MISE_4,'MISE0.3_4=',MISE3_4,'MISE0.7_4=',MISE7_4,'LMISE_4=',MISE3_4+MISE7_4,'\n')


yfinal_1<- apply(yfix_1,1,mean)  # apply函数只能用于处理矩阵类型的数据，也就是说所有的数据必须是同一类型。apply函数一般有三个参数，第一个参数代表矩阵对象，第二个参数代表要操作矩阵的维度，1表示对行进行处理，2表示对列进行处理。第三个参数就是处理数据的函数。apply会分别一行或一列处

#####  绘图  sigmaeps=0.1#########
mytitle<- paste('N=',N,'n=',n,'sigma=',sigmaeps_1)
plot(ufix,yfinal_1,type = 'l',lty=2,col='blue',xlim=c(0,1),ylim=c(0.5,2.5),
     xlab="x",ylab="y",main=mytitle)
lines(ufix,f1(ufix),col='red')
# 置信区间
mytitle<- paste('N=',N,'n=',n,'sigma=',sigmaeps_1)
plot(ufix,yfinal_1,type = 'l',lty=2,col='blue',xlim=c(0,1),ylim=c(0.5,2.5),
     xlab="x",ylab="y",main=mytitle)
lines(ufix,f1(ufix),col='red')
fL<-apply(yfix_1,1,quantile,probs = c(0.025,0.975))
lines(ufix,fL[1,],col='blue')
lines(ufix,fL[2,],col='blue')
####### 核密度估计
plot(density(jumps_1,bw = 0.015,kernel = 'gaussian'),xlim=c(0,1))

#####  绘图  sigmaeps=0.2#########
yfinal_2<- apply(yfix_2,1,mean)


mytitle<- paste('N=',N,'n=',n,'sigma=',sigmaeps_2)
plot(ufix,yfinal_2,type = 'l',lty=2,col='blue',xlim=c(0,1),ylim=c(0.5,2.5),
     xlab="x",ylab="y",main=mytitle)
lines(ufix,f1(ufix),col='red')
# 置信区间
mytitle<- paste('N=',N,'n=',n,'sigma=',sigmaeps_2)
plot(ufix,yfinal_2,type = 'l',lty=2,col='blue',xlim=c(0,1),ylim=c(0.5,2.5),
     xlab="x",ylab="y",main=mytitle)
lines(ufix,f1(ufix),col='red')
fL<-apply(yfix_2,1,quantile,probs = c(0.025,0.975))
lines(ufix,fL[1,],col='blue')
lines(ufix,fL[2,],col='blue')
####### 核密度估计
plot(density(jumps_2,bw = 0.015,kernel = 'gaussian'),xlim=c(0,1))

#####  绘图  sigmaeps=0.4#########
yfinal_4<- apply(yfix_4,1,mean)
mytitle<- paste('N=',N,'n=',n,'sigma=',sigmaeps_4)
plot(ufix,yfinal_4,type = 'l',lty=2,col='blue',xlim=c(0,1),ylim=c(0.5,2.5),
     xlab="x",ylab="y",main=mytitle)
lines(ufix,f1(ufix),col='red')
# 置信区间
mytitle<- paste('N=',N,'n=',n,'sigma=',sigmaeps_4)
plot(ufix,yfinal_4,type = 'l',lty=2,col='blue',xlim=c(0,1),ylim=c(0.5,2.5),
     xlab="x",ylab="y",main=mytitle)
lines(ufix,f1(ufix),col='red')
fL<-apply(yfix_4,1,quantile,probs = c(0.025,0.975))
lines(ufix,fL[1,],col='blue')
lines(ufix,fL[2,],col='blue')
####### 核密度估计
plot(density(jumps_4,bw = 0.015,kernel = 'gaussian'),xlim=c(0,1))



##########sigmaeps=0.1
J1<- which(jumps_1<0.35&jumps_1>0.25)
mean3<- mean(jumps_1[J1])
sd3<- sd(jumps_1[J1])
J2<- which(jumps_1<0.75&jumps_1>0.65)
mean7<- mean(jumps_1[J2])
sd7<- sd(jumps_1[J2])
cat('at 0.3, bias=',mean3-0.3,'sd=',sd3,'power=',length(jumps_1[J1])/N,'\n')
cat('at 0.7, bias=',mean7-0.7,'sd=',sd7,'power=',length(jumps_1[J2])/N,'\n')
Len<- rep(0,N)
for(i in 1:N){
  Len[i]<- length(jumps.list_1[[i]])-1
}
cat(sum(Len==0),sum(Len==1),sum(Len==2),sum(Len==3),sum(Len>3),'\n')
cat('Total jumps=',length(jumps_1),'\n')

##########sigmaeps=0.2
J1<- which(jumps_2<0.35&jumps_2>0.25)
mean3<- mean(jumps_2[J1])
sd3<- sd(jumps_2[J1])
J2<- which(jumps_2<0.75&jumps_2>0.65)
mean7<- mean(jumps_2[J2])
sd7<- sd(jumps_2[J2])
cat('at 0.3, bias=',mean3-0.3,'sd=',sd3,'power=',length(jumps_2[J1])/N,'\n')
cat('at 0.7, bias=',mean7-0.7,'sd=',sd7,'power=',length(jumps_2[J2])/N,'\n')
Len<- rep(0,N)
for(i in 1:N){
  Len[i]<- length(jumps.list_2[[i]])-1
}
cat(sum(Len==0),sum(Len==1),sum(Len==2),sum(Len==3),sum(Len>3),'\n')
cat('Total jumps=',length(jumps_2),'\n')
#save.image('RSSout0613RData')

##########sigmaeps=0.4
J1<- which(jumps_4<0.35&jumps_4>0.25)
mean3<- mean(jumps_4[J1])
sd3<- sd(jumps_4[J1])
J2<- which(jumps_4<0.75&jumps_4>0.65)
mean7<- mean(jumps_4[J2])
sd7<- sd(jumps_4[J2])
cat('at 0.3, bias=',mean3-0.3,'sd=',sd3,'power=',length(jumps_4[J1])/N,'\n')
cat('at 0.7, bias=',mean7-0.7,'sd=',sd7,'power=',length(jumps_4[J2])/N,'\n')
Len<- rep(0,N)
for(i in 1:N){
  Len[i]<- length(jumps.list_4[[i]])-1
}
cat(sum(Len==0),sum(Len==1),sum(Len==2),sum(Len==3),sum(Len>3),'\n')
cat('Total jumps=',length(jumps_4),'\n')
#save.image('RSSout0613RData')

