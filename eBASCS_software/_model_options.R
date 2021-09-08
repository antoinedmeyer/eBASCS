# Choose model options for spatial, spectral and temporal submodels

###############
# Spatial Model
###############
psf.form <- "parametric"

home = getwd()
setwd(paste(home,"/psf/",sep="")) 
source("psf.R") # 2D King profile density (King 1962)

if (psf.form == "parametric"){
  
  function_load_name <- ".R"
  slope <- 1.5  # power-law slope
  r0 <- 0.6  # core radius
  ellip <- 0.00573533  # ellipticity
  off.angle <- 0  # off-axis angle
  
  # Calculate normalizing constant
  source("normalize_psf.R")
  #source("/home/antoinem/Documents/astrostats/sim_study/psf/normalize_psf.R")
  psf.norm <- 1
  psf.norm <- 1/normalize_psf(100)$integral 
  
  # Load C++ version PSF (requires psf.norm calculated above)
  #sourceCpp('psf_cpp.cpp')
  #sourceCpp('/home/antoinem/Documents/astrostats/sim_study/psf/psf_cpp.cpp')
} else {
  
  # Alternatively load an image of the PSF 
  # User must supply image
  # List of length three:
  # psf_image[[1]] x vector
  # psf_image[[2]] y vector
  # psf_image[[3]] grid of psf values 
  function_load_name <- "_imagepsf.R"
  setwd(paste(home,"/psf/psf_image/",sep=""))
  load("psf_image.RData")
  source("psf_image_function.R") 
  source("adjust_grids.R")
  
  xpsf_grid <- psf_image[[1]] 
  ypsf_grid <- psf_image[[2]] 
  psf_image_values <- psf_image[[3]] 
  
  midpoint_alignment <- adjust_grids()
  xpsf_grid <- midpoint_alignment[[1]]
  ypsf_grid <- midpoint_alignment[[2]]
  psf_image_values <- midpoint_alignment[[3]]
  
}



################
# Temporal Model
################
time_spacing <- "none"
if (model == 'spacetime' | model == 'eBASCS') {
  # (0) Time Model: dirichlet
  time_model <- "dirichlet"
  
  # (1) Choose type of temporal segmentation. Options: "equal" (choose time breaks), "algo.UVCet"
  num_time_breaks <- 4
  time_spacing <- "algo.UVCet"
  if (time_spacing == 'equal') {
    num_time_breaks <- 4
  }
}



################
# Spectral Model
################

if (model == 'BASCS' | model == 'eBASCS') {
  # Choose spectral model. Options: "full", "extended", "extended_full"
  spectral_model <- "extended"
} else {
  spectral_model <- 'none'
}





setwd(home)