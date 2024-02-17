# load packages ================================================================
library(mvtnorm)
library(nlraa)
library(fdapace)

# function for simulations =====================================================

# simulation function
simu_fcpa <- function(n, nmes, txdo, missing){

  # create noisy data with dropout from FCPA prediction ========================
  
  # generate dropout
  if (missing == "fixed"){
    dtpast <- rbind(rep(10, n), dt_noisy[-nmes,])
    dropout <- matrix(NA, ncol = n, nrow = nmes)
    dropout <- ifelse(dtpast > 24 + 4*(txdo==0.3) - 3.5*(nmes == 7) + 0.5*(nmes == 7)*(txdo==0.3) - 6*(nmes == 5), 0, 1)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    txevent <- sum(!apply(dropout, 2, prod))/n
  }
  
  if (missing == "threshMAR"){
    prob_Yobs_linpred <- -2.75 + 3*(dt_noisy > 30) - 1*(txdo == 0.3) + 1*(nmes == 7) - 0.25*(nmes==7)*(txdo == 0.3) + 1.5*(nmes==5) - 0.25*(nmes==5)*(txdo == 0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    prob_Yobs <- rbind(matrix(0, ncol = n, nrow = 1), prob_Yobs[-nmes,])
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    txevent <- sum(!apply(dropout, 2, prod))/n
  }
  
  if (missing == "incrMAR"){
    prob_Yobs_linpred <- 0.2*dt_noisy - 6 - 1.25*(txdo == 0.3) + 1*(nmes==7) + 1.5*(nmes == 5) + 0.25*(nmes==5)*(txdo==0.3)
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    prob_Yobs <- rbind(matrix(0, ncol = n, nrow = 1), prob_Yobs[-nmes,])
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    txevent <- sum(!apply(dropout, 2, prod))/n
  }

  if (missing == "threshMNAR"){
    prob_Yobs_linpred <- -3 + 3*(dt_noisy > 30) - 1.25*(txdo == 0.3) + 0.75*(nmes == 7) + 1*(nmes==5) +0.25*(nmes==5)*(txdo == 0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    txevent <- sum(!apply(dropout, 2, prod))/n
  }
  
  if (missing == "incrMNAR"){
    prob_Yobs_linpred <- 0.2*dt_noisy - 6.5 - 1*(txdo == 0.3) + 0.75*(nmes==7) - 0.25*(txdo == 0.3)*(nmes==7) + 1*(nmes == 5) - 0.25*(nmes==5)*(txdo==0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    txevent <- sum(!apply(dropout, 2, prod))/n
  }
  
  if (missing == "MCAR"){
    mat_Lt <- do.call(cbind, dt_Lt)
    prob_Yobs_linpred <- -7 + 0.5*mat_Lt + 1.25*(txdo == 0.3) - 0.05*(nmes==7) +1*(nmes==5) - 2.5*(nmes==5)*(txdo==0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    txevent <- sum(!apply(dropout, 2, prod))/n
  }
  
  # add dropout to noisy data
  dt_dropout <- ifelse(dropout == 0, NA, dt_noisy)
  #matplot(do.call(cbind, dt_Lt), dt_dropout, type="l", ylim = c(0,40))
  #matplot(do.call(cbind, dt_Lt), dt_noisy, type="l", ylim = c(0,40))
  
  # re-estimate FPCA on FPCA generated noisy data ==============================
  
  # preprocessing
  dt_Ly_noisy <- as.list(as.data.frame(dt_noisy))
  
  # estimation
  # fcpa_est_fve09_noisy <- FPCA(dt_Ly_noisy, dt_Lt, list(dataType = "Sparse", methodSelectK = 2, usergrid = FALSE))
  fcpa_est_fve09_noisy <- FPCA(dt_Ly_noisy, dt_Lt, list(dataType = "Sparse", methodSelectK = 2, usergrid = FALSE, lean = TRUE, imputeScores = TRUE))
  
  # extraction of FPC & scores
  fpc_noisy <- fcpa_est_fve09_noisy$phi
  scores_noisy <- fcpa_est_fve09_noisy$xiEst
  mu_noisy <- fcpa_est_fve09_noisy$mu
  times_noisy <- fcpa_est_fve09_noisy$workGrid
  
  # re-estimate FPCA on FPCA generated noisy data with dropout =================
  
  # preprocessing
  dt_Ly_dropout <- as.list(as.data.frame(dt_dropout))
  
  # estimation
  # fcpa_est_fve09_dropout <- FPCA(dt_Ly_dropout, dt_Lt, list(dataType = "Sparse", methodSelectK = 2, usergrid = FALSE))
  fcpa_est_fve09_dropout <- FPCA(dt_Ly_dropout, dt_Lt, list(dataType = "Sparse", methodSelectK = 2, usergrid = FALSE, lean = TRUE, imputeScores = TRUE))
  
  # extraction of FPC & scores
  fpc_dropout <- fcpa_est_fve09_dropout$phi
  scores_dropout <- fcpa_est_fve09_dropout$xiEst
  mu_dropout <- fcpa_est_fve09_dropout$mu
  times_dropout <- fcpa_est_fve09_dropout$workGrid
  
  return(list(#list(fpc_grid, mu_grid, times_grid), 
              # list(fpc_grid, mu_grid, xi, times_grid), 
              list(fpc_noisy, mu_noisy, times_noisy, apply(scores_noisy, 2, var)), 
              list(fpc_dropout, mu_dropout, times_dropout, apply(scores_dropout, 2, var))))
              # list(fpc_dropout, mu_dropout, scores_dropout, times_dropout)))
  
}

# simulation run and saving output =============================================

# loading data obtained after FPCA
n <- 200
load("dt_preds.Rdata")
if (nmes == 5){
  mu <- dt_preds[[1]][[1]]
  fpc <- dt_preds[[1]][[2]]
  scores <- dt_preds[[1]][[3]]
  times <- dt_preds[[1]][[4]]
  dt_Lt <- dt_preds[[1]][[5]]
  mu_grid <- dt_preds[[1]][[6]]
  fpc_grid <- dt_preds[[1]][[7]]
  times_grid <- dt_preds[[1]][[8]]
} else if (nmes == 7){
  mu <- dt_preds[[2]][[1]]
  fpc <- dt_preds[[2]][[2]]
  scores <- dt_preds[[2]][[3]]
  times <- dt_preds[[2]][[4]]
  dt_Lt <- dt_preds[[2]][[5]]
  mu_grid <- dt_preds[[2]][[6]]
  fpc_grid <- dt_preds[[2]][[7]]
  times_grid <- dt_preds[[2]][[8]]
} else if (nmes == 13){
  mu <- dt_preds[[3]][[1]]
  fpc <- dt_preds[[3]][[2]]
  scores <- dt_preds[[3]][[3]]
  times <- dt_preds[[3]][[4]]
  dt_Lt <- dt_preds[[3]][[5]]
  mu_grid <- dt_preds[[3]][[6]]
  fpc_grid <- dt_preds[[3]][[7]]
  times_grid <- dt_preds[[3]][[8]]
}
rm(dt_preds)

# sampling new scores 
xi <- rmvnorm(n, c(0,0), diag(apply(scores,2,var)))
rm(scores)

# re-create FPCA generated data with new generated scores
dt_pred <- matrix(NA, nrow = nmes, ncol = n)
for (i in seq(n)){
  ti <- dt_Lt[[i]]
  idxti <- which(round(times,7) %in% round(ti,7))
  dt_pred[,i] <- mu[idxti] + xi[i,] %*% t(fpc[idxti,])
}
rm(mu, xi, fpc)

# create noisy data from FPCA generated data ===================================
epsilon <- matrix(rnorm(nmes*n, 0, 3), nrow = nmes, ncol = n)
dt_noisy <- dt_pred + epsilon
rm(dt_pred, epsilon)

# computing FPCA on noisy data and on noisy data with dropout 
idx <- 1
out <- NULL

for (missing in c("fixed", "threshMAR", "incrMAR", "threshMNAR", "incrMNAR", "MCAR")){
  print(missing)
  out[[idx]] <- list(list(fpc_grid, mu_grid, times_grid), simu_fcpa(n, nmes, txdo, missing))
  idx <- idx+1
}

print("run done")

# saving output
if (dir.exists(paste0("output/n",n,"_nmes",nmes,"_txdo",txdo))){
  save(out, file = paste0("output/n",n,"_nmes",nmes,"_txdo",txdo,"/",job,"_",rep,".Rdata"))
} else {
  dir.create(paste0("output/n",n,"_nmes",nmes,"_txdo",txdo), recursive = TRUE)
  save(out, file = paste0("output/n",n,"_nmes",nmes,"_txdo",txdo,"/",job,"_",rep,".Rdata"))
}

print("output done")

# end of script ================================================================
q("no")