# display(Jumps)
# AllJumps
Smooth <- 2
Image_smooth <- matrix(0,ROW,COL)
for(i in 1:ROW){
  for(j in 1:COL){
    Image_smooth[i,j] <- mean(Image_noise[GetArea(c(i,j),ROW,COL,Smooth)])
  }
}
display(Image_smooth,method = 'r')

# Image_smooth <- Image_noise
sigma2hat <- mean(Sigma2hat)
clear <- rep(0,nrow(AllJumps))
for (index in 1:nrow(AllJumps)){
  orig <- AllJumps[index,]
  near <- GetArea(AllJumps[index,],ROW,COL,Smooth)
  post <- t((t(near)-orig))
  # print(post)
  W <- diag(BivKernal(x = post[,1]/ROW,y=post[,2]/COL,h=0.035),nrow = nrow(near))
  Y <- Image_smooth[near]
  X <- cbind(1,post)
  beta <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%Y
  Area1 <- near[beta[2]*post[,1]+beta[3]*post[,2]>0,,drop=FALSE]
  Area2 <- near[beta[2]*post[,1]+beta[3]*post[,2]<0,,drop=FALSE]
  clear[index] <- abs(mean(Image_smooth[Area1])-mean(Image_smooth[Area2]))/
    sqrt(sigma2hat/(2*Smooth+1)^2*(1/nrow(Area1)+1/nrow(Area2)))
  # print(sqrt(sigma2hat/(2*Smooth+1)^2*(1/nrow(Area1)+1/nrow(Area2))))
}
plot(clear)
# tmp <- AllJumps[clear>1.96,,drop=FALSE]
lxlx <- kmeans(clear,2,nstart = 10)
lxlx$centers
tmp <- AllJumps[lxlx$cluster==which.max(lxlx$centers),,drop=FALSE]
lxl <- matrix(0,ROW,COL)
lxl[tmp] <- 1
lxl[near] <- 1
display(lxl,method = 'r')

display(Image_noise)

mydata <- clear
wss <- rep(0,14)
for (i in 2:15) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")



library(fpc)
pamk.best <- pamk(mydata)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")




require(vegan)
fit <- cascadeKM(scale(mydata, center = TRUE,  scale = TRUE), 1, 10, iter = 1000)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

d <- clear
library(mclust)
d_clust <- Mclust(as.matrix(d), G=1:20)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")   # 4 clusters
plot(d_clust)


library(apcluster)
d.apclus <- apcluster(negDistMat(r=2), matrix(d,ncol = 1))
cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")    # 4
heatmap(d.apclus)
plot(d.apclus, matrix(d,ncol = 1))

library(cluster)
clusGap(matrix(d,ncol = 1), kmeans, 10, B = 100, verbose = interactive())
