jumps_find<-function(u0,u1,degree,yj,uj,err,sigma2hat,alpha){
  jumps2 <- c()
  n <- length(yj)
  k0 <- length(u0)
  nu1 <- length(u1)
  DRSStotal<- matrix(rep(0),length(u1),1)
  # fit1<- lm(yj ~ bs(uj, degree = degree, knots = u0,Boundary.knots = c(0,200),intercept=T)-1)
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
      # fit2<- lm(yj ~ bs(uj, degree = degree, knots = us,Boundary.knots = c(0,200),intercept=T)-1)
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
  # lambdaM_1[Nt]<- lambda_1
  q2_1<- qchisq(1-alpha,degree+1,lambda_1)*sigma2hat

  lambda_2<- max(mean(DRSSuT)/sigma2hat-degree-1,0)
  # lambdaM_2[Nt]<- lambda_2
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

  while (length(idtemp_2)>0){
    #从id_1中把包含在nx的（2p+1）hn的邻域去掉
    nx<- which(DRSStotal_temp==max(DRSStotal_temp,na.rm = TRUE))
    u.jump<- mean(u1[nx])
    jumps2<- c(jumps2,u.jump)
    #将跳点插入
    if(length(u.jump)>0){
      du<- abs(u0-u.jump)
      m<- which(du==min(du))
      m<- m[1]
      if(du[m]<err/2){
        u0<- c(u0[-m],rep(u.jump,degree+1))
        u0<- sort(u0)
      }else{
        u0<- c(u0,rep(u.jump,degree+1))
        u0<- sort(u0)
      }
      }
    idinterval<-  c(max(1,ceiling((mean(nx)))-ceiling(nu1*2*(degree+0.5)/(k0+1))):min(ceiling((mean(nx)))+ceiling(nu1*2*(degree+0.5)/(k0+1)),length(u1)))
    #cat('idinterval',idinterval,'\n')
    if (length(idtemp_1)>0){
      id1_1<-idtemp_1
      for (j in 1:length(id1_1)){
        for (i in 1:length(idinterval)){
          if (id1_1[j]==idinterval[i]){
            idtemp_1<-idtemp_1[-(which(idtemp_1==id1_1[j]))]
          }
        }
      }
      }
    id1_2<-idtemp_2
    for (j in 1:length(id1_2)){
      for (i in 1:length(idinterval)){
        if (id1_2[j]==idinterval[i]){
          idtemp_2<-idtemp_2[-(which(idtemp_2==id1_2[j]))]
        }
      }
      }
    # cat('idtemp_1',idtemp_1,'\n')
    # cat('idtemp_2',idtemp_2,'\n')
    DRSStotal_temp[max(1,ceiling((mean(nx)))-ceiling(nu1*2*(degree+0.5)/(k0+1))):min(ceiling((mean(nx)))+ceiling(nu1*2*(degree+0.5)/(k0+1)),length(u1))]<-0
    # cat('DRSStotal_temp',DRSStotal_temp,'\n')
    }
  # cat('jumps2',jumps2,'\n')
  jumps_find_list<-list(u0,jumps2)
  return(jumps_find_list)
  }
