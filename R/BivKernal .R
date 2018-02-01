### bivariate kernal二元 核函数支撑要求输入值在-1到1之间
BivKernal <- function(x,y,h=0.01){
  x <- x/h
  y <- y/h
  return(1/h^2*exp(-(x^2+y^2)/2)-exp(-0.5))/(2*pi-3*pi*exp(-0.5))
}
