require(EBImage)
require(yangtiao)

# rm(list = ls())
# gc()

# load('.//tmp//block.Rdata')
# display(mark/max(mark))

ROW <- nrow(Image_noise)
COL <- ncol(Image_noise)

IgnPoints <- lapply(1:nrow(AllJumps),function(x) GetArea(point = AllJumps[x,],Row = ROW,Col = COL,h = 1))
IgnPoints <- unique.matrix(do.call(rbind,IgnPoints)) #这些点由于噪声或者跳的原因不用于拟合
mark_noise <- mark
mark_noise[IgnPoints] <- 0
display(mark_noise/max(mark_noise),method = 'r')

if(!exists('denoise')) denoise <- FALSE

Image_fit <- matrix(0,ROW,COL)

system.time(
  for(Mode in 1:ModeNum){
    cat('######## Mode ',Mode,' ########','\n')
    AllPost <- which(mark==Mode,arr.ind = TRUE)
    AllY <- Image_noise[AllPost]
    Post <- AllPost
    Ynoise <- AllY

    if(denoise){
      Post <- which(mark_noise==Mode,arr.ind = TRUE)
      Ynoise <- Image_noise[Post]
    }
    res <- ChoBic1(lowerKnot = lowerKnot,upperKnot = upperKnot,Post = Post,Y = Ynoise,AllPost=AllPost,AllY=AllY,updegree = updegree,lowdegree = lowdegree)
    # res <- ChoBic(lower = lowerKnot,upper = upperKnot,Post = Post,Y = Ynoise,AllPost=AllPost,AllY=AllY)
    print(res)

    # 节点
    KnotNumX <- res$KnotNumX
    KnotNumY <- res$KnotNumY

    # 次数
    DegreeX <- res$DegreeX
    DegreeY <- res$DegreeY

    # DegreeX <- 5
    # DegreeY <- 5

    KnotX <- seq(range(Post[,1])[1]-1,range(Post[,1])[2]+1,length.out = KnotNumX)
    KnotY <- seq(range(Post[,2])[1]-1,range(Post[,2])[2]+1,length.out = KnotNumY)

    X <- BaValue(KnotX,KnotY,DegreeX,DegreeY,ValueX=Post[,1],ValueY=Post[,2])
    AllX <- BaValue(KnotX,KnotY,DegreeX,DegreeY,ValueX=AllPost[,1],ValueY=AllPost[,2])


    Coef <- coef(lm(Ynoise~X-1))
    Coef[is.na(Coef)] <- 0
    Image_fit[AllPost] <- AllX%*%Coef
  }
)


if(exists('Image_raw')){
  Diff <- abs(Image_fit-Image_raw)
  MISE <- mean(Diff^2)
  print(MISE)

  AllNear <- lapply(1:nrow(AllJumps_clear),function(x) GetArea(point = AllJumps_clear[x,],Row = ROW,Col = COL,h = round(min(ROW,COL)*0.04)))
  AllNear <- unique.matrix(do.call(rbind,AllNear))

  MISEe <- mean((Image_fit[AllNear]-Image_raw[AllNear])^2)
  print(MISEe)

  MISEE <- matrix(0,ROW,COL)
  MISEE[AllNear] <- 1
  display(MISEE,method = 'raster')
  display(Image_raw,method = 'raster')
  display(Diff,method = 'raster')
  # Image_fit[200,67]
}

display(Image_fit,method = 'raster')
save.image(paste0('.//tmp//',MISE*10^4,'.Rdata'))
# save.image(paste0('.//res//',MISE*10^4,'.Rdata'))
mean((Image_fit[AllPost]-Image_raw[AllPost])^2)



# mu <- tryCatch({
#   solve(t(X)%*%X)%*%t(X)},
#   error=function(e){
#     print(e)
#     MASS::ginv(t(X)%*%X)%*%t(X) #广义逆
#   }
# )
#
# # 控制点
# betax <- mu%*%dat$x
# betay <- mu%*%dat$y
# beta <- mu%*%dat$z
#
# newx <- runif(2500,min(dat$x),max(dat$x))
# newy <- runif(2500,min(dat$y),max(dat$y))
# newX <- BaValue(ux,uy,nx,ny,newx,newy)
# newz <- newX%*%beta
#
# plot_ly(data=data.frame(newx,newy,newz),x=~newx,y=~newy,z=~newz,sizes=c(3,6))%>%
#   add_markers(size=~c(1),color=I('#4AC6B7'),name='fitted points')%>%#拟合
#   add_mesh(x=dat$x,y=dat$y,z=dat$z,opacity=0.7,color=I('#C61951'))%>%
#   add_markers(data = data.frame(betax,betay,beta),
#               x=~betax,y=~betay,z=~beta,size=~c(2),
#               color=I('#1972A4'),name='control points') #加入控制点
