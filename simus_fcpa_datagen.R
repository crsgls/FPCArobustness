library(mvtnorm)
library(nlraa)
library(fdapace)

# data generating function
gene_data <- function(n, nmes, tmax){
  
  # initialization
  y <- matrix(NA, nrow = nmes, ncol = n)
  t <- matrix(0, nrow = nmes, ncol = n)
  tdef <- seq(0, tmax, length.out = nmes)
  
  # noise around planned time measurement
  for (i in seq(n)){
    t[,i] <- tdef + c(0, runif(nmes-1,-tmax/(2*(nmes-1)), tmax/(2*(nmes-1))))
  }
  
  for (i in seq(n)){
    # random effects
    bi <- rmvnorm(1, mean = rep(0,4), sigma = diag(c(3,5,2,2)))
    # outcome generation (scenario 2)
    y[, i] <- SSlogis5(t[,i], 10+bi[1], 35+bi[2], 10+bi[3], 5+bi[4], 1 + bi[4]/3) + rnorm(nmes, 0, 3) 
  }
  
  return(list(t, y))
  
}

# simulation hyperparameters
n <- 200
dt_preds <- list()
cpt <- 1

for (nmes in c(5, 7, 13)){

  # generating full data
  dt <- gene_data(n, nmes, 12)
  
  # preprocessing
  dt_Lt <- as.list(as.data.frame(dt[[1]]))
  dt_Ly <- as.list(as.data.frame(dt[[2]]))
  
  # estimation
  fpca_est_K2 <- FPCA(dt_Ly, dt_Lt, list(dataType = "Sparse", methodSelectK = 2, usergrid = TRUE))
  fpca_est_K2_grid <- FPCA(dt_Ly, dt_Lt, list(dataType = "Sparse", methodSelectK = 2, usergrid = FALSE))
  
  # extraction of FPC & scores
  mu <- fpca_est_K2$mu
  fpc <- fpca_est_K2$phi
  scores <- fpca_est_K2$xiEst
  times <- fpca_est_K2$obsGrid
  
  mu_grid <- fpca_est_K2_grid$mu
  fpc_grid <- fpca_est_K2_grid$phi
  times_grid <- fpca_est_K2_grid$workGrid
  
  # compute prediction for each subject
  # dt_pred <- matrix(NA, nrow = nmes, ncol = n)
  # for (i in seq(n)){
  #   ti <- dt_Lt[[i]]
  #   idxti <- which(round(times,8) %in% round(ti,8))
  #   dt_pred[,i] <- mu[idxti] + scores[i,] %*% t(fpc[idxti,])
  # }
  
  dt_preds[[cpt]] <- list(mu, fpc, scores, times, dt_Lt, mu_grid, fpc_grid, times_grid)
  print(nmes)
  
  cpt <- cpt + 1
}

save(dt_preds, file = "dt_preds_n200.Rdata")
