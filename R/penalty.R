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

    for (i in 1:length(d)){
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
