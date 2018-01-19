require(parallel)
require(yangtiao)
setwd('//gpfsTMP//math//xllu//R//Rfile//parallel//sigma2')

if (dir.exists('.//res')) unlink('.//res', recursive=TRUE)
dir.create('.//res')

core <- detectCores()
core
cl <-  makeCluster(core)
cl
# 
AllMISE <- mclapply(1:1280,function(x){
  message('######## ',x,' ########')
  source('.//demo//BSplineSurface.R')
  return(c(MISE,MISEe))
  }, mc.cores = core)

print(AllMISE)
save(AllMISE,file = './/res//AllMISE.Rdata')
stopCluster(cl = cl)



  