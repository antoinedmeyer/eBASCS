normalize_psf <- function(bound){

  norm_fun <- function(pts){
    return(psf(pts,c(0,0)))
  }
  
  value <- adaptIntegrate(norm_fun,c(-bound,-bound),c(bound,bound))
  return(value)
}

# func_x0 <- function(y) {
#   pts <- c(0,y)
#   return(0.04470557-psf(pts,c(0,0)))
# }
# 
# func_y0 <- function(x) {
#   pts <- c(x,0)
#   return(0.04470557-psf(pts,c(0,0)))
# }
# 
# uniroot(func_x0,c(0,5))
# uniroot(func_y0,c(0,5))

