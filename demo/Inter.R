# save(AllJumps_clear,file = 'example.Rdata')
require(EBImage)
require(yangtiao)
load('example.Rdata')

ROW=200
COL=200


Jumps_clear <- matrix(0,ROW,COL)
Jumps_clear[AllJumps_clear] <- 1
Res <- Jumps_clear

Jumps_edge <- Jumps_clear
Jumps_edge[c(1,ROW),] <- 1
Jumps_edge[,c(1,COL)] <- 1
Jumpstmp <- Jumps_edge

display(Jumps_clear)
display(Jumps_edge)


Judg <- function(index,h=1){
  near <- GetArea(AllJumps_clear[index,],ROW,COL,h)
  points <- near[Jumps_edge[near]==1,,drop=FALSE]
  if (nrow(points)==1) return(Judg(index,h=h+1))
  if (max(dist(points))>=2*h){
    return(h)
  }else{
    return(Judg(index,h=h+1))
  }
}

H <- matrix(rep(1,2*nrow(AllJumps_clear)),ncol = 2)

for (index in 1:nrow(AllJumps_clear)){
  H[index,] <- c(index,Judg(index,h=1))
}
H <- H[order(H[,2],decreasing = TRUE),]
H <- H[H[,2]>1,]


for (iii in 1:nrow(H)){
  index <- H[iii,1]
  h <- Judg(index,h=1)
  if (h>1){
    near <- GetArea(AllJumps_clear[index,],ROW,COL,h)
    points <- near[Jumps_edge[near]==1,,drop=FALSE]
    # browser()
    # lxl <- matrix(0,ROW,COL)
    # lxl[points] <- 1
    # display(lxl)
    # display(lxl[1:10,41:50])
    Dist <- as.matrix(dist(points))
    vertex <- points[which(Dist==max(Dist),arr.ind = TRUE)[,1],,drop=FALSE]
    x2 <- AllJumps_clear[index,]
    for (ii in 1:nrow(vertex)){
      x1 <- vertex[ii,]
      if(abs(x1[1]-x2[1])>=abs(x1[2]-x2[2])){
        res <- do.call(cbind,approx(x=c(x1[1],x2[1]),y=c(x1[2],x2[2]),xout = x1[1]:x2[1]))
      }else{
        res <- do.call(cbind,approx(x=c(x1[2],x2[2]),y=c(x1[1],x2[1]),xout = x1[2]:x2[2]))[,2:1]
      }
      Jumps_edge[round(res)] <- 1
      Res[round(res)] <- 1
    }
  }
}

display(Res)
display(Jumps_edge)
display(Jumps_edge-Jumpstmp+Jumps_clear)
display(Jumps_edge-Jumpstmp)


# Jumps_clear <- Jumps_edge-Jumpstmp+Jumps_clear
AllJumps_clear <- which(Res==1,arr.ind = TRUE)
for (i in 1:nrow(AllJumps_clear)){
  if(Judg(i,1)>1) {print(AllJumps_clear[i,]);print(Judg(i,1))}
}
# row col
# 5  48
# row col
# 195  50
# row col
# 194  52
# row col
# 185  70 a
# row col
# 175  86 a
# row col
# 152 124 a
# row col
# 149 126
# row col
# 148 128
# AllJumps_clear[183,]
