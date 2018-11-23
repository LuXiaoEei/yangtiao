require(yangtiao)
require(ggplot2)

x <- sort(runif(200,0,1))
y <- sin(x*2*pi)+sign(x>0.5)

u1 <- c(0,0.3,0.6,0.9,1.0001)
X1 <- BaSplite1(x = x,degree = 2,u = u1)
fit1 <- X1%*%solve(t(X1)%*%X1)%*%t(X1)%*%y


u2 <- c(0,0.2,0,2,0.2,0.6,0.9,1.0001)
X2 <- BaSplite1(x = x,degree = 2,u = u2)
fit2 <- X2%*%MASS::ginv(t(X2)%*%X2)%*%t(X2)%*%y



u3 <- c(0,0.3,0,49,0,49,0.49,0.6,0.9,1.0001)
X3 <- BaSplite1(x = x,degree = 2,u = u3)
fit3 <- X3%*%MASS::ginv(t(X3)%*%X3)%*%t(X3)%*%y

pic <- data.frame(x=x,y=c(y,fit1,fit2,fit3),tag=rep(c('true values','no t_0','t_0=0.2','t_0=0.5'),each=200))


ggplot(pic,aes(x = x,y = y))+geom_line(aes(color=tag,linetype=tag),show.legend=F,size=0.8)
  # scale_colour_manual(values=c('red','blue','yellow','black'))+
  # scale_colour_discrete(name='tag',breaks=c('true values','no t_0','t_0=0.2','t_0=0.5'),
   # labels=c('true values','no duplicate',expression(paste(t[0],'=0.2')), expression(paste(t[0],'=0.5'))))
ggplot(pic[pic$tag=='true values',],aes(x = x,y = y))+geom_line()


plot(x,y,type = 'l')
lines(x,fit1,type = 'l')
lines(x,fit2,type = 'l')
lines(x,fit3,type = 'l')
