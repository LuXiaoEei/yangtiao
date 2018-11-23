#B样条检测跳点，2016年06月检测双跳函数，单次检测一个跳点

library(splines)
library(MASS)
library(yangtiao)
#清空内存变量
rm(list=ls())
gc()


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

# begin -------------------------------------------------------------------

# 最近邻 因变量y是当x给定时的函数值。当自变量换为x0时，此函数返回因变量的值。
interp<- function(y,x,x0){
  n<- length(x0)
  yhat<- matrix(rep(0),n,1)
  n1<- length(x)
  for (i in 1:n) {
    a<- abs(x-x0[i])
    yhat[i]<- y[order(a)[1]]
  }
  return(yhat)
}

# BIC准则选择节点函数
BIC_knot<-function(low,up,yj,uj,n,alpha)
{
  BIC<- matrix(0,3,up-low+1)
  for(k in low:up){
    for(degree in 1:3){
      u0<- seq(0,1,length.out = k)
      # fit1<- lm(yj ~ bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
      # fit1<- lm(yj ~ sapply(0:(length(u0)-2+degree),BaSplite,x=uj,degree=degree,u=u0,n=length(u0))-1)
      fit1<- lm(yj ~ BaSplite1(uj,degree,u0)-1)
      y1<- fitted(fit1)
      RSS0<-sum(fit1$residuals^2)    # cat('RSS0',RSS0,'\n')
      BIC[degree,k-low+1]<- log(RSS0/n)+(k+degree+1)/n*log(n)
    }
  }
  k0<- ceiling(which.min(BIC)/3)+low-1
  degree<- (which.min(BIC)-1)%%3+1
  q2<- qchisq(1-alpha,degree+1)
  err<-  0.8/(k0-1)
  #cat('k0=',k0,'degree=',degree,'\n')
  u0<- seq(0,1,length.out = k0)
  nu0<- length(u0)
  BIC_list<-list(u0,nu0,degree,k0,err)
  return(BIC_list)
}

# 利用惩罚选择跳点个数函数
penalty<-function(d,k0,degree,uj,yj,n,jumps2,err){
  penalty<-matrix(rep(0),1,length(d)+1)
  BIC_jump<-matrix(rep(0),1,length(d)+1)
  jumps_option<-c()

  u0<- seq(0,1,length.out = k0)
  # fit1<- lm(yj ~ bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
  # fit1<- lm(yj ~ sapply(0:(length(u0)-2+degree),BaSplite,x=uj,degree=degree,u=u0,n=length(u0))-1)
  fit1<- lm(yj ~ BaSplite1(uj,degree,u0)-1)

  y1<- fitted(fit1)
  RSS0<-sum(fit1$residuals^2)
  BIC_jump[1]<-n*log(RSS0/n)


  if(length(d)>0){
  d1<-d
  # cat('u0',u0,'\n')
  P_n_s=sqrt(n*((0.9-0.1)/(k0+1))^2*(log(n))^2)
  P_n_m=sqrt(n/k0*log(n))
  P_n_l=sqrt(n/k0)*log(n)

  for (i in 1:length(d))
  {


    du<- abs(u0-jumps2[which(d==max(d1))])
    m<- which(du==min(du))
    m<- m[1]
    if(du[m]<err){
      u0<- c(u0[-m],rep(jumps2[which(d==max(d1))],degree+1))
      u0<- sort(u0)
    }else{
      u0<- c(u0,rep(jumps2[which(d==max(d1))],degree+1))
      u0<- sort(u0)
    }
    # fit1<- lm(yj ~ bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
    # fit1<- lm(yj ~ sapply(0:(length(u0)-2+degree),BaSplite,x=uj,degree=degree,u=u0,n=length(u0))-1)
    fit1<- lm(yj ~ BaSplite1(uj,degree,u0)-1)

    y1<- fitted(fit1)
    RSS0<-sum(fit1$residuals^2)
    #cat('RSS0',RSS0,'\n')
    BIC_jump[i+1]<-n*log(RSS0/n)
    #cat('BIC_jump',BIC_jump,'\n')
    penalty[i+1]<-P_n_s/abs(d[which(d==max(d1))])
    # cat('penalize',penalize,'\n')
    jumps_option<-c(jumps_option,jumps2[which(d==max(d1))])


    d1<-d1[-which(d==max(d1))]
    # cat('u0',u0,'\n')

  }
  }
  BIC_jump_1<- BIC_jump+cumsum(penalty)
  cat('BIC_jump_1',BIC_jump_1,'\n')
  opt_k<-which(BIC_jump_1==min(BIC_jump_1))-1
  cat('opt_k',opt_k,'\n')
  if(opt_k>=1){
    jumps_option<-jumps_option[-(opt_k+1):-(length(d)+1)]}
  else{jumps_option<-c()}
  cat('jumps_option',jumps_option,'\n')

  return(jumps_option)
}



jumps_find<-function(u0,u1,degree,yj,uj,err,sigma2hat,alpha,lambdaM_1,lambdaM_2){
  jumps2 <- c()
  n <- length(yj)
  k0 <- length(u0)
  nu1 <- length(u1)

  # fit1<- lm(yj ~ bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,1),intercept=T)-1)
  # fit1<- lm(yj ~ sapply(0:(length(u0)-2+degree),BaSplite,x=uj,degree=degree,u=u0,n=length(u0))-1)
  fit1<- lm(yj ~ BaSplite1(uj,degree,u0)-1)
  yfit1<- fitted(fit1)


  DRSSu<-DRSSu2<-DRSSu3<-DRSSuT<- rep(0,nu1)
  for(i in 1:nu1){
    us<- u0
    m=which(u0==u1[i])
    if(length(m)>0){
      us<- c(u0[-m],rep(u1[i],degree+1))
      us<- sort(us)
      DRSStotal[i]<- NA
    }else{
      us<- c(u0,rep(u1[i],degree+1))
      us<- sort(us)
      # fit2<- lm(yj ~ bs(uj, degree = degree, knots = us,Boundary.knots = c(0,1),intercept=T)-1)
      # fit2<- lm(yj ~ sapply(0:(length(us)-2+degree),BaSplite,x=uj,degree=degree,u=us,n=length(us))-1)
      fit2<- lm(yj ~ BaSplite1(uj,degree,us)-1)

      yfit2<- fitted(fit2)
      DRSStotal[i]<- sum((yfit1-yfit2)^2)
      #确定x的区间
      it<- order(abs(u0-u1[i]))[1:2]
      it1<- min(it)
      it2<- max(it)
      idM2<- which(uj>=u0[it2]+err|uj<u0[it1])
      RSS2<- sum((yfit1[idM2]-yfit2[idM2])^2)*n/(length(idM2))
      idM3<- which(uj>=u0[it2]|uj<u0[it1]-err)
      RSS3<- sum((yfit1[idM3]-yfit2[idM3])^2)*n/(length(idM3))
      DRSSuT[i]<- max(RSS2,RSS3)

    }
  }
  lambda_1<- mean(DRSStotal,na.rm = TRUE)/sigma2hat-degree-1

  lambda_1<- max(lambda_1,0)
  lambdaM_1[Nt]<- lambda_1
  q2_1<- qchisq(1-alpha,degree+1,lambda_1)*sigma2hat


  lambda_2<- max(mean(DRSSuT)/sigma2hat-degree-1,0)
  lambdaM_2[Nt]<- lambda_2
  q2_2<- qchisq(1-alpha,degree+1,lambda_2)*sigma2hat


  #plot(u1,DRSStotal,main = paste('Simulation',Nt))
  #abline(h=q2,col='red')
  ### 找跳点
  id_1<- which(DRSStotal>q2_1)
  id_2<- which(DRSStotal>q2_2)
  #cat('id_1',id_1,'\n')
  #cat('id_2',id_2,'\n')
  idtemp_1<-id_1
  idtemp_2<-id_2
  DRSStotal_temp<-DRSStotal

  while (length(idtemp_2)>0)

    ####从id_1中把包含在nx的（2p+1）hn的邻域去掉

  {
    nx<- which(DRSStotal_temp==max(DRSStotal_temp,na.rm = TRUE))
    # cat('nx',nx,'\n')
    mean(nx)
    u.jump<- mean(u1[nx])
    jumps2<- c(jumps2,u.jump)
    # cat(u.jump,'\n')
    #将跳点插入
    if(length(u.jump)>0)
    { du<- abs(u0-u.jump)
    m<- which(du==min(du))
    m<- m[1]
    if(du[m]<err/2)
    {
      u0<- c(u0[-m],rep(u.jump,degree+1))
      u0<- sort(u0)
    }
    else{
      u0<- c(u0,rep(u.jump,degree+1))
      u0<- sort(u0)
    }
    }


    idinterval<-  c( max(1,ceiling((mean(nx)))-ceiling(1000*2*(degree+0.5)/(k0+1))):min(ceiling((mean(nx)))+ceiling(1000*2*(degree+0.5)/(k0+1)),length(u1)))
    #cat('idinterval',idinterval,'\n')
    if (length(idtemp_1)>0)
    {   id1_1<-idtemp_1
    for (j in 1:length(id1_1))
    {
      for (i in 1:length(idinterval))
      {
        if (id1_1[j]==idinterval[i]) {idtemp_1<-idtemp_1[-(which(idtemp_1==id1_1[j]))]}
      }
    }
    }
    id1_2<-idtemp_2
    for (j in 1:length(id1_2))
    {
      for (i in 1:length(idinterval))
      {
        if (id1_2[j]==idinterval[i]) {idtemp_2<-idtemp_2[-(which(idtemp_2==id1_2[j]))]}
      }
    }
    # cat('idtemp_1',idtemp_1,'\n')
    # cat('idtemp_2',idtemp_2,'\n')
    DRSStotal_temp[max(1,ceiling((mean(nx)))-ceiling(1000*2*(degree+0.5)/(k0+1))):min(ceiling((mean(nx)))+ceiling(1000*2*(degree+0.5)/(k0+1)),length(u1))]<-0
    # cat('DRSStotal_temp',DRSStotal_temp,'\n')
  }
  # cat('jumps2',jumps2,'\n')
  jumps_find_list<-list(u0,jumps2)
  return(jumps_find_list)
  }


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
    jumps_find_list_1<-jumps_find(u0_1,u1,degree_1,yj_1,uj,err_1,sigma2hat_1,alpha)
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
    jumps_find_list_2<-jumps_find(u0_2,u1,degree_2,yj_2,uj,err_2,sigma2hat_2,alpha)
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
    jumps_find_list_4<-jumps_find(u0_4,u1,degree_4,yj_4,uj,err_4,sigma2hat_4,alpha)
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

