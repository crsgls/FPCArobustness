# load packages ================================================================
library(mvtnorm)
library(splines)
library(lme4)
library(fdapace)
library(MFPCA)
library(lcmm)
library(tidyverse)
library(plyr)
library(nlraa)
library(JM)

# functions for simulations ====================================================

# data generating function
gene_data <- function(n, nmes, tmax, scen){
  
  # initialization
  y <- matrix(NA, nrow = nmes, ncol = n)
  t <- matrix(0, nrow = nmes, ncol = n)
  tdef <- seq(0, tmax, length.out = nmes)

  for (i in seq(n)){
    
    # noise around planned time measurement
    t[,i] <- tdef + c(0, runif(nmes-1,-tmax/(2*(nmes-1)), tmax/(2*(nmes-1))))
    
  }
  
  for (i in seq(n)){
    # spline basis
    bst <- bs(t, df = 3)
    
    if (scen == 1){
      # random effects
      bi <- rmvnorm(1, mean = rep(0,4), sigma = diag(c(3,5,1,1)))
      # outcome generation
      y[, i] <- SSlogis5(t[,i], 10+bi[1], 30+bi[2], 6+bi[3], 5+bi[4], 1 + bi[4]/3) + rnorm(nmes, 0, 3)
    }
    
    if (scen == 2){
      # random effects
      bi <- rmvnorm(1, mean = rep(0,4), sigma = diag(c(3,5,2,2)))
      # outcome generation
      y[, i] <- SSlogis5(t[,i], 10+bi[1], 35+bi[2], 10+bi[3], 5+bi[4], 1 + bi[4]/3) + rnorm(nmes, 0, 3) 
    }
  }
  
  return(list(t, y))
  
}

# simulation function for comparing RMSE between lmm and fpca
simuRMSE <- function(dt, nmes, tmax, trainsize, testsize, missing, txdo, scen, est_compl){
  
  # generating dropout =========================================================
  
  # dropout fixed (MAR)
  if (missing == "fixed"){
    dtpast <- rbind(rep(10, n), dt[[2]][-nmes,])
    dropout <- matrix(NA, ncol = n, nrow = nmes)
    #if (scen == 1) dropout <- ifelse(dtpast > 29.5 + 3*(txdo==0.3) - 2*(nmes == 7) + 0.5*(nmes == 7)*(txdo==0.3) - 3.5*(nmes == 5), 0, 1) # 1
    if (scen == 2) dropout <- ifelse(dtpast > 24 + 4*(txdo==0.3) - 3.5*(nmes == 7) + 0.5*(nmes == 7)*(txdo==0.3) - 6*(nmes == 5), 0, 1) # 2
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    dt[[3]] <- ifelse(dropout == 0, NA, dt[[2]]) 
    
    txevent <- sum(!apply(dropout, 2, prod))/n
    meanNA <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x)))))
    meanNAdo <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))[which(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))>0)])
  }
  
  # dropout linked to threshold (MAR)
  if (missing == "threshMAR"){
    # if (scen == 1) prob_Yobs_linpred <- -3.25 + 3*(dt[[2]] > 30) - 1*(txdo == 0.3) + 1.25*(nmes == 7) - 0.25*(nmes==7)*(txdo == 0.3) + 1.5*(nmes==5) - 0.25*(nmes==5)*(txdo == 0.3) # 1
    if (scen == 2) prob_Yobs_linpred <- -2.75 + 3*(dt[[2]] > 30) - 1*(txdo == 0.3) + 1*(nmes == 7) - 0.25*(nmes==7)*(txdo == 0.3) + 1.5*(nmes==5) - 0.25*(nmes==5)*(txdo == 0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    prob_Yobs <- rbind(matrix(0, ncol = n, nrow = 1), prob_Yobs[-nmes,])
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dt[[3]] <- ifelse(dropout == 0, NA, dt[[2]]) 

    txevent <- sum(!apply(dropout, 2, prod))/n
    meanNA <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x)))))
    meanNAdo <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))[which(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))>0)])
  }
  
  # dropout whose probability increase with y (MAR)
  if (missing == "incrMAR"){
    # if (scen == 1) prob_Yobs_linpred <- 0.2*dt[[2]] - 7 - 1.25*(txdo == 0.3) + 0.75*(nmes==7) + 1.5*(nmes == 5) - 0.25*(nmes==5)*(txdo==0.3) # 1
    if (scen == 2) prob_Yobs_linpred <- 0.2*dt[[2]] - 6 - 1.25*(txdo == 0.3) + 1*(nmes==7) + 1.5*(nmes == 5) + 0.25*(nmes==5)*(txdo==0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    prob_Yobs <- rbind(matrix(0, ncol = n, nrow = 1), prob_Yobs[-nmes,])
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dt[[3]] <- ifelse(dropout == 0, NA, dt[[2]]) 
    
    txevent <- sum(!apply(dropout, 2, prod))/n
    meanNA <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x)))))
    meanNAdo <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))[which(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))>0)])
  }
  
  # dropout linked to threshold (MNAR)
  if (missing == "threshMNAR"){
    # if (scen == 1) prob_Yobs_linpred <-  -3.25 + 3*(dt[[2]] > 30) - 1.25*(txdo == 0.3) + 0.75*(nmes == 7) - 0.25*(nmes==7)*(txdo == 0.3) + 1*(nmes==5) +0.25*(nmes==5)*(txdo == 0.3) # 1
    if (scen == 2) prob_Yobs_linpred <- -3 + 3*(dt[[2]] > 30) - 1.25*(txdo == 0.3) + 0.75*(nmes == 7) + 1*(nmes==5) +0.25*(nmes==5)*(txdo == 0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    dt[[3]] <- ifelse(dropout == 0, NA, dt[[2]]) 
    
    txevent <- sum(!apply(dropout, 2, prod))/n
    meanNA <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x)))))
    meanNAdo <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))[which(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))>0)])
  }
  
  # dropout whose probability increase with y (MNAR)
  if (missing == "incrMNAR"){
    # if (scen == 1) prob_Yobs_linpred <- 0.2*dt[[2]] - 7.5 - 1*(txdo == 0.3) + 0.75*(nmes==7) - 0.25*(txdo == 0.3)*(nmes==7) + 1*(nmes == 5) # 1
    if (scen == 2) prob_Yobs_linpred <- 0.2*dt[[2]] - 6.5 - 1*(txdo == 0.3) + 0.75*(nmes==7) - 0.25*(txdo == 0.3)*(nmes==7) + 1*(nmes == 5) - 0.25*(nmes==5)*(txdo==0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    dt[[3]] <- ifelse(dropout == 0, NA, dt[[2]]) 
    
    txevent <- sum(!apply(dropout, 2, prod))/n
    meanNA <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x)))))
    meanNAdo <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))[which(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))>0)])
  }
  
  # dropout (MCAR)
  if (missing == "MCAR"){
    # if (scen == 1) prob_Yobs_linpred <- -7 + 0.5*dt[[1]]  - 1*(txdo == 0.3) + 0.5*(nmes==7) + 1*(nmes==5) # 1
    if (scen == 2) prob_Yobs_linpred <- -7 + 0.5*dt[[1]] + 1.25*(txdo == 0.3) - 0.05*(nmes==7) +1*(nmes==5) - 2.5*(nmes==5)*(txdo==0.3) # 2
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    dt[[3]] <- ifelse(dropout == 0, NA, dt[[2]]) 

    txevent <- sum(!apply(dropout, 2, prod))/n
    meanNA <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x)))))
    meanNAdo <- mean(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))[which(apply(dt[[3]], 2, function(x) return(sum(is.na(x))))>0)])
  }
  
  # discretizing time ==========================================================
  dt[[4]] <- round_any(dt[[1]], 0.5, f = ceiling)
  
  # train/test dataset =========================================================
  ids_train <- seq(as.numeric(trainsize))
  ids_test <- seq(as.numeric(testsize)) + as.numeric(trainsize)
  
  # lmm preprocessing estimation and prediction ================================
  
  # preprocessing
  dt_long <- suppressWarnings(data.frame("id" = sort(rep(seq(trainsize+testsize), nmes)),
                              "t" = as.numeric(paste(dt[[1]])), 
                              "t10" = as.numeric(paste(dt[[1]]))/10,
                              "yfull" = as.numeric(paste(dt[[2]])),
                              "y" = as.numeric(paste(dt[[3]]))))
  dt_long_train <- dt_long[dt_long$id %in% ids_train,]
  dt_long_test <- dt_long[!(dt_long$id %in% ids_train),]
  
  # estimation
  lmm_est_quad <- lcmm(fixed = y ~ 1 + t10 + I(t10^2), random = y ~ t10 + I(t10^2), subject = "id", ng = 1, data = dt_long_train)
  lmm_est_cub <- lcmm(fixed = y ~ 1 + t10 + I(t10^2) + I(t10^3), random = y ~ t10 + I(t10^2) + I(t10^3), subject = "id", ng = 1, data = dt_long_train)
  lmm_est_spl <- lcmm(fixed = y ~ 1 + ns(t, knots = quantile(t, probs = c(0.33,0.66))), random = y ~ 1 + ns(t, knots = quantile(t, probs = c(0.33,0.66))), subject = "id", ng = 1, data = dt_long_train)
  lmm_est_spl2 <- lcmm(fixed = y ~ 1 + ns(t, knots = c(4,8)), random = y ~ 1 + ns(t, knots = c(4,8)), subject = "id", ng = 1, data = dt_long_train)

  est_JM <- NA
  if ((missing == "threshMNAR") | (missing == "incrMNAR")){
    
    # JM
    # lmmObject
    lmmObject <- lme(y ~ 1 + t10 + I(t10^2), random = ~ 1 + t10 + I(t10^2) | id, data = dt_long_train, na.action = na.omit, control = list(returnObject = TRUE))
    # survObject
    dropout_train <- dropout[,ids_train]
    eventdropout_train <- ifelse(apply(dropout_train,2,sum)<nmes,1,0)
    idxtimetoevent_train <- pmin(nmes, apply(dropout_train,2,sum)+1)
    timetoevent_train <- dt[[1]][cbind(idxtimetoevent_train, ids_train)]
    dtsurv <- cbind(ids_train, timetoevent_train/10, eventdropout_train)
    dtsurv <- as.data.frame(dtsurv)
    colnames(dtsurv) <- c("id", "timetoevent", "event")
    survObject <- coxph(Surv(time = dtsurv$timetoevent+0.001, event = dtsurv$event)~1, data = dtsurv, x = TRUE, 
                        control = coxph.control(timefix = FALSE))
    est_JM <- jointModel(lmmObject, survObject, timeVar = "t10", method = "spline-PH-GH",
                         control = list(verbose = TRUE, nproc = 1, lng.in.kn = 3, optimizer = "optim"))
    
  }

  # prediction of complete data model on test data
  predRE_compl <- predictRE(est_compl, newdata = dt_long_test)
  pred_compl <- NULL
  bst <- bs(dt_long_test$t, df = 3)
  for (i in unique(dt_long_test$id)){
    bsti <- predict(bst, dt_long_test[dt_long_test$id == i,"t"] )
    pred_compl <- c(pred_compl,
                   as.numeric(est_compl$best)[14] + as.numeric(est_compl$best)[15]*
                     (predRE_compl[predRE_compl$id == i,2] + (as.numeric(est_compl$best)[1] + predRE_compl[predRE_compl$id == i,3]) * bsti[,1] + 
                        (as.numeric(est_compl$best)[2] + predRE_compl[predRE_compl$id == i,4])* bsti[,2]  + 
                        (as.numeric(est_compl$best)[3] + predRE_compl[predRE_compl$id == i,5])* bsti[,3])) 
  }
  
  # prediction on test data
  lmm_predRE_quad <- predictRE(lmm_est_quad, newdata = dt_long_test)
  lmm_pred_quad <- numeric(testsize*nmes); j <- 1
  for (i in unique(dt_long_test$id)){
    lmm_pred_quad[(nmes*(j-1)+1):(nmes*j)] <- as.numeric(lmm_est_quad$best)[9] + as.numeric(lmm_est_quad$best)[10]*
                         (lmm_predRE_quad[lmm_predRE_quad$id == i,2] + (as.numeric(lmm_est_quad$best)[1] + lmm_predRE_quad[lmm_predRE_quad$id == i,3]) * dt_long_test[dt_long_test$id == i,"t10"] + (as.numeric(lmm_est_quad$best)[2] + lmm_predRE_quad[lmm_predRE_quad$id == i,4])* dt_long_test[dt_long_test$id == i,"t10"]^2)
    j <- j + 1
  }
  
  lmm_predRE_cub <- predictRE(lmm_est_cub, newdata = dt_long_test)
  lmm_pred_cub <- numeric(testsize*nmes); j <- 1
  for (i in unique(dt_long_test$id)){
    lmm_pred_cub[(nmes*(j-1)+1):(nmes*j)] <- as.numeric(lmm_est_cub$best)[14] + as.numeric(lmm_est_cub$best)[15]*
                         (lmm_predRE_cub[lmm_predRE_cub$id == i,2] + (as.numeric(lmm_est_cub$best)[1] + lmm_predRE_cub[lmm_predRE_cub$id == i,3]) * dt_long_test[dt_long_test$id == i,"t10"] + 
                            (as.numeric(lmm_est_cub$best)[2] + lmm_predRE_cub[lmm_predRE_cub$id == i,4])* dt_long_test[dt_long_test$id == i,"t10"]^2 +
                            (as.numeric(lmm_est_cub$best)[3] + lmm_predRE_cub[lmm_predRE_cub$id == i,5])* dt_long_test[dt_long_test$id == i,"t10"]^3)
    j <- j + 1
  }
  
  lmm_predRE_spl <- predictRE(lmm_est_spl, newdata = dt_long_test)
  lmm_pred_spl <- numeric(testsize*nmes); j <- 1
  for (i in unique(dt_long_test$id)){
    lmm_pred_spl[(nmes*(j-1)+1):(nmes*j)] <- as.numeric(lmm_est_spl$best)[14] + as.numeric(lmm_est_spl$best)[15]*(lmm_predRE_spl[lmm_predRE_spl$id == i,2] + (as.numeric(lmm_est_spl$best[1:3] + lmm_predRE_spl[lmm_predRE_spl$id == i,3:5]) %*% t(as.matrix(predict(ns(dt_long_test$t, knots = quantile(dt_long_test$t, probs = c(0.33,0.66))), seq(0,tmax,length.out=nmes))))))
    j <- j + 1
  }
  
  lmm_predRE_spl2 <- predictRE(lmm_est_spl2, newdata = dt_long_test)
  lmm_pred_spl2 <- numeric(testsize*nmes); j <- 1
  for (i in unique(dt_long_test$id)){
    lmm_pred_spl2[(nmes*(j-1)+1):(nmes*j)] <- as.numeric(lmm_est_spl2$best)[14] + as.numeric(lmm_est_spl2$best)[15]*(lmm_predRE_spl2[lmm_predRE_spl2$id == i,2] + (as.numeric(lmm_est_spl2$best[1:3] + lmm_predRE_spl2[lmm_predRE_spl2$id == i,3:5]) %*% t(as.matrix(predict(ns(dt_long_test$t, knots = c(4,8)), seq(0,tmax,length.out=nmes))))))
    j <- j + 1
  }
  
  JM_pred <- NA
  if ((missing == "threshMNAR") | (missing == "incrMNAR")){
    
    JM_pred_mat <- matrix(NA, nrow = testsize, ncol = nmes); j <- 1;
    for (i in unique(dt_long_test$id)){
      JM_pred_mat[j,] <- c(predict.jointModel(est_JM, newdata = na.omit(dt_long_test[dt_long_test$id == i, ]), FtTimes = dt_long_test[dt_long_test$id == i, "t10"], type = "Subject"))
      j = j+1
    }
    JM_pred <- c(t(JM_pred_mat))
    
  }
  
  # fpca preprocessing estimation and prediction ===============================
  
  # preprocessing
  dt_Lt_train <- as.list(as.data.frame(dt[[4]][,ids_train]))
  dt_Ly_train <- as.list(as.data.frame(dt[[3]][,ids_train]))
  dt_Lt_test <- as.list(as.data.frame(dt[[4]][,ids_test]))
  dt_Ly_test <- as.list(as.data.frame(dt[[3]][,ids_test]))
  
  # estimation
  fpca_est_fve09 <- FPCA(dt_Ly_train, dt_Lt_train, list(dataType = "Sparse", FVEthreshold = 0.9, usergrid = TRUE))
  fpca_est_fve099 <- FPCA(dt_Ly_train, dt_Lt_train, list(dataType = "Sparse", FVEthreshold = 0.99, usergrid = TRUE))
  fpca_est_K4 <- FPCA(dt_Ly_train, dt_Lt_train, list(dataType = "Sparse", methodSelectK = as.integer(4), usergrid = TRUE))
  
  # prediction on test data
  fpca_pred_obj_fve09 <- predict(fpca_est_fve09, dt_Ly_test, dt_Lt_test)
  fpca_pred_fve09 <- as.numeric(t(fpca_pred_obj_fve09$predCurves))
  fpca_pred_fve09 <- fpca_pred_fve09[unlist(lapply(dt_Lt_test, function(x) return(fpca_pred_obj_fve09$predGrid %in% x)))]
  fpca_pred_obj_fve099 <- predict(fpca_est_fve099, dt_Ly_test, dt_Lt_test)
  fpca_pred_fve099 <- as.numeric(t(fpca_pred_obj_fve099$predCurves))
  fpca_pred_fve099 <- fpca_pred_fve099[unlist(lapply(dt_Lt_test, function(x) return(fpca_pred_obj_fve099$predGrid %in% x)))]
  fpca_pred_obj_K4 <- predict(fpca_est_K4, dt_Ly_test, dt_Lt_test)
  fpca_pred_K4 <- as.numeric(t(fpca_pred_obj_K4$predCurves))
  fpca_pred_K4 <- fpca_pred_K4[unlist(lapply(dt_Lt_test, function(x) return(fpca_pred_obj_K4$predGrid %in% x)))]
  
  # RMSE ======================================================================= HERE TO CHANGE !!!!!
  
  n_suj_NA <- testsize - sum(apply(dropout[,ids_test], 2, prod))
  
  # rmse for model learn on complete data
  lmm_predtab_complete <- data.frame("diff" = pred_compl - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
  lmm_predtab_complete_NA <- lmm_predtab_complete[lmm_predtab_complete$miss == 1,]
  lmm_predtab_complete_obs <- lmm_predtab_complete[lmm_predtab_complete$miss == 0,]
  lmm_RMSE_complete_NA <- sqrt(mean(as.numeric(by(lmm_predtab_complete_NA, lmm_predtab_complete_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
  lmm_RMSE_complete_obs <- sqrt(mean(as.numeric(by(lmm_predtab_complete_obs, lmm_predtab_complete_obs$id, function(x) return(t(x$diff) %*% x$diff)))))
  
  # rmse for lmm quad
  lmm_predtab_quad <- data.frame("diff" = lmm_pred_quad - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
  lmm_predtab_quad_NA <- lmm_predtab_quad[lmm_predtab_quad$miss == 1,]
  lmm_predtab_quad_obs <- lmm_predtab_quad[lmm_predtab_quad$miss == 0,]
  lmm_RMSE_quad_NA <- sqrt(mean(as.numeric(by(lmm_predtab_quad_NA, lmm_predtab_quad_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
  lmm_RMSE_quad_obs <- sqrt(mean(as.numeric(by(lmm_predtab_quad_obs, lmm_predtab_quad_obs$id, function(x) return(t(x$diff) %*% x$diff)))))
  
  # rmse for lmm cub
  lmm_predtab_cub <- data.frame("diff" = lmm_pred_cub - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
  lmm_predtab_cub_NA <- lmm_predtab_cub[lmm_predtab_cub$miss == 1,]
  lmm_predtab_cub_obs <- lmm_predtab_cub[lmm_predtab_cub$miss == 0,]
  lmm_RMSE_cub_NA <- sqrt(mean(as.numeric(by(lmm_predtab_cub_NA, lmm_predtab_cub_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
  lmm_RMSE_cub_obs <- sqrt(mean(as.numeric(by(lmm_predtab_cub_obs, lmm_predtab_cub_obs$id, function(x) return(t(x$diff) %*% x$diff)))))
  
  # rmse for lmm spl quant
  lmm_predtab_spl <- data.frame("diff" = lmm_pred_spl - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
  lmm_predtab_spl_NA <- lmm_predtab_spl[lmm_predtab_spl$miss == 1,]
  lmm_predtab_spl_obs <- lmm_predtab_spl[lmm_predtab_spl$miss == 0,]
  lmm_RMSE_spl_NA <- sqrt(mean(as.numeric(by(lmm_predtab_spl_NA, lmm_predtab_spl_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
  lmm_RMSE_spl_obs <- sqrt(mean(as.numeric(by(lmm_predtab_spl_obs, lmm_predtab_spl_obs$id, function(x) return(t(x$diff) %*% x$diff)))))
  
  # rmse for lmm spl equi
  lmm_predtab_spl2 <- data.frame("diff" = lmm_pred_spl2 - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
  lmm_predtab_spl2_NA <- lmm_predtab_spl2[lmm_predtab_spl2$miss == 1,]
  lmm_predtab_spl2_obs <- lmm_predtab_spl2[lmm_predtab_spl2$miss == 0,]
  lmm_RMSE_spl2_NA <- sqrt(mean(as.numeric(by(lmm_predtab_spl2_NA, lmm_predtab_spl2_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
  lmm_RMSE_spl2_obs <- sqrt(mean(as.numeric(by(lmm_predtab_spl2_obs, lmm_predtab_spl2_obs$id, function(x) return(t(x$diff) %*% x$diff)))))
  
  # rmse for JM
  JM_RMSE_NA <- NA; JM_RMSE_obs <- NA;
  if ((missing == "threshMNAR") | (missing == "incrMNAR")){
    JM_predtab <- data.frame("diff" = JM_pred - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
    JM_predtab_NA <- JM_predtab[JM_predtab$miss == 1,]
    JM_predtab_obs <- JM_predtab[JM_predtab$miss == 0,]
    JM_RMSE_NA <- sqrt(mean(as.numeric(by(JM_predtab_NA, JM_predtab_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
    JM_RMSE_obs <- sqrt(mean(as.numeric(by(JM_predtab_obs, JM_predtab_obs$id, function(x) return(t(x$diff) %*% x$diff)))))
  }
  
  # rmse for FPCA with fve = 0.9
  fpca_predtab_fve09 <- data.frame("diff" = fpca_pred_fve09 - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
  fpca_predtab_fve09_NA <- fpca_predtab_fve09[fpca_predtab_fve09$miss == 1,]
  fpca_predtab_fve09_obs <- fpca_predtab_fve09[fpca_predtab_fve09$miss == 0,]
  fpca_RMSE_fve09_NA <- sqrt(mean(as.numeric(by(fpca_predtab_fve09_NA, fpca_predtab_fve09_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
  fpca_RMSE_fve09_obs <- sqrt(mean(as.numeric(by(fpca_predtab_fve09_obs, fpca_predtab_fve09_obs$id, function(x) return(t(x$diff) %*% x$diff)))))
  
  # rmse for FPCA with fve = 0.99
  fpca_predtab_fve099 <- data.frame("diff" = fpca_pred_fve099 - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
  fpca_predtab_fve099_NA <- fpca_predtab_fve099[fpca_predtab_fve099$miss == 1,]
  fpca_predtab_fve099_obs <- fpca_predtab_fve099[fpca_predtab_fve099$miss == 0,]
  fpca_RMSE_fve099_NA <- sqrt(mean(as.numeric(by(fpca_predtab_fve099_NA, fpca_predtab_fve099_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
  fpca_RMSE_fve099_obs <- sqrt(mean(as.numeric(by(fpca_predtab_fve099_obs, fpca_predtab_fve099_obs$id, function(x) return(t(x$diff) %*% x$diff)))))
  
  # rmse for FPCA with K = 4
  fpca_predtab_K4 <- data.frame("diff" = fpca_pred_K4 - dt_long_test$yfull, "miss" = is.na(dt_long_test$y)*1, "id" = dt_long_test$id)
  fpca_predtab_K4_NA <- fpca_predtab_K4[fpca_predtab_K4$miss == 1,]
  fpca_predtab_K4_obs <- fpca_predtab_K4[fpca_predtab_K4$miss == 0,]
  fpca_RMSE_K4_NA <- sqrt(mean(as.numeric(by(fpca_predtab_K4_NA, fpca_predtab_K4_NA$id, function(x) return(t(x$diff) %*% x$diff)))))
  fpca_RMSE_K4_obs <- sqrt(mean(as.numeric(by(fpca_predtab_K4_obs, fpca_predtab_K4_obs$id, function(x) return(t(x$diff) %*% x$diff)))))

  res = c(txevent, meanNA, meanNAdo, n_suj_NA, fpca_est_fve09$selectK, fpca_est_fve099$selectK, lmm_RMSE_complete_NA, lmm_RMSE_complete_obs, 
          lmm_RMSE_quad_NA, lmm_RMSE_quad_obs, lmm_RMSE_cub_NA, lmm_RMSE_cub_obs, lmm_RMSE_spl_NA, lmm_RMSE_spl_obs, lmm_RMSE_spl2_NA, lmm_RMSE_spl2_obs, 
          JM_RMSE_NA, JM_RMSE_obs, fpca_RMSE_fve09_NA, fpca_RMSE_fve09_obs, fpca_RMSE_fve099_NA, fpca_RMSE_fve099_obs, fpca_RMSE_K4_NA, fpca_RMSE_K4_obs)
  
  return(res)

}

# simulation run and saving output =============================================

trainsize <- 200
nmes <- 5 # 5, 7 or 13
txdo <- 0.3 # 0.3 or 0.6
scen <- 2
testsize = 500
tmax = 12

error = TRUE
attempts <- 0

while (error == TRUE){
  
  # generating full data
  n <- trainsize + testsize
  dt <- gene_data(n, nmes, 12, scen)
  
  # estimation of flexible model on complete train data
  ids_train <- seq(as.numeric(trainsize))
  dt_long_train_full <- suppressWarnings(data.frame("id" = sort(rep(ids_train, nmes)),
                                                    "t" = as.numeric(paste(dt[[1]][, ids_train])),
                                                    "yfull" = as.numeric(paste(dt[[2]][, ids_train]))))
  est_compl <- lcmm(fixed = yfull ~ bs(t, df = 3), random = yfull ~ bs(t, df = 3), subject = "id", ng = 1, data = dt_long_train_full)
  
  # prediction of all models on all missing data schemes
  
  out <- matrix(NA, ncol = 24, nrow = 6); idx = 1
  for (missing in c("fixed", "threshMAR", "incrMAR", "threshMNAR", "incrMNAR", "MCAR")){
  # out <- matrix(NA, ncol = 24, nrow = 4); idx = 1
  # for (missing in c("fixed", "threshMAR", "incrMAR", "MCAR")){
    print(missing)
    outmiss <- try(simuRMSE(dt, nmes, tmax, trainsize, testsize, missing, txdo, scen, est_compl), silent = TRUE)
    error = inherits(outmiss, "try-error")
    
    if (error == TRUE){
      break;
    }
    
    out[idx,] <- outmiss
    idx <- idx+1
  }
  attempts <- attempts + 1
}

# save number of attempts
out <- cbind(out, rep(attempts, dim(out)[1]))

# saving output

if (dir.exists(paste0("output/trainsize",trainsize,"_scen",scen,"_nmes",nmes,"_txdo",txdo))){
  save(out, file = paste0("output/trainsize",trainsize,"_scen",scen,"_nmes",nmes,"_txdo",txdo,"/",job,"_",rep,".Rdata"))
} else {
  dir.create(paste0("output/trainsize",trainsize,"_scen",scen,"_nmes",nmes,"_txdo",txdo), recursive = TRUE)
  save(out, file = paste0("output/trainsize",trainsize,"_scen",scen,"_nmes",nmes,"_txdo",txdo,"/",job,"_",rep,".Rdata"))
}

# end of script ================================================================
q("no")