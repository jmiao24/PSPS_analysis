library(knockoff)
library(doParallel)
source("./Fun.R")
source("./knockoff-filter/stats_debiased_Lasso.R")

psps_knockoff <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, fdr = 0.05, knockoffs=create.second_order, family = "gaussian", parallel = T, boot = 200, offset=1){

  knockoffs=create.second_order
  # Input parameter
  n_lab = nrow(X_lab)
  n_unlab = nrow(X_unlab)
  p_lab = ncol(X_lab)
  p_unlab = ncol(X_unlab)

  # Create knockoff variables
  knock_lab = knockoffs(X_lab)
  knock_unlab = knockoffs(X_unlab)

  Xk_lab <- knock_lab
  Xk_unlab <- knock_unlab
  rm(knock_lab)
  rm(knock_unlab)

  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X_lab),1,0.5)
  swap.M.lab = matrix(swap,nrow=nrow(X_lab),ncol=length(swap),byrow=TRUE)
  X_lab.swap  = X_lab * (1-swap.M.lab) + Xk_lab * swap.M.lab
  Xk_lab.swap = X_lab * swap.M.lab + Xk_lab * (1-swap.M.lab)
  swap.M.unlab = matrix(swap,nrow=nrow(X_unlab),ncol=length(swap),byrow=TRUE)
  X_unlab.swap  = X_unlab * (1-swap.M.unlab) + Xk_unlab * swap.M.unlab
  Xk_unlab.swap = X_unlab * swap.M.unlab + Xk_unlab * (1-swap.M.unlab)

  X_lab_combined <- cbind(X_lab.swap, Xk_lab.swap)
  X_unlab_combined <- cbind(X_unlab.swap, Xk_unlab.swap)


  # Compute Statistics
  parallel <- T
  boot <- 200
  W <- psps_dlasso_stats(X_lab_combined, X_unlab_combined, Y_lab, Yhat_lab, Yhat_unlab, parallel = parallel, swap = swap, boot = boot)

  # Run the knockoff filter
  offset=1
  fdr <- 0.1
  t = knockoff.threshold(W, fdr=fdr, offset=offset)
  selected = sort(which(W >= t))
  X.names = colnames(X)
  names(selected) = X.names[selected]
  fdp(selected)
  power(selected)

  W1 <- dlasso_stats(X_lab_combined, X_unlab_combined, Y_lab, Yhat_lab, Yhat_unlab)
  # Run the knockoff filter
  offset=1
  fdr <- 0.1
  t = knockoff.threshold(W1, fdr=fdr, offset=offset)
  selected = sort(which(W1 >= t))
  X.names = colnames(X)
  names(selected) = X.names[selected]
  fdp(selected)
  power(selected)
  
  # Package up the results.
  structure(list(call = match.call(),
                 X = X_lab,
                 Xk = Xk_lab,
                 y = Y_lab,
                 statistic = W,
                 threshold = t,
                 selected = selected),
                 class = 'knockoff.result')
}

# False 
fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
# Power
power <- function(selected) sum(beta[selected] != 0) / sum(beta != 0)


psps_dlasso_stats <- function(X_lab_combined, X_unlab_combined, Y_lab, Yhat_lab, Yhat_unlab, parallel = T, boot = 200, swap){
  # boot <- 5
  dlasso_fit <- psps_dlasso(X_lab_combined, X_unlab_combined, Y_lab, Yhat_lab, Yhat_unlab, parallel = parallel, boot = boot)
  coef <- dlasso_fit$coef
  p <- length(coef)/2
  orig <- 1:p
  W = abs(coef[orig]) - abs(coef[orig+p])

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
  return(W)
}

dlasso_stats <- function(X_lab_combined, Y_lab, swap){
  coef <- dlasso(as.matrix(X_lab_combined), as.matrix(Y_lab))$unb.coef
  p <- length(coef)/2
  orig <- 1:p
  W = abs(coef[orig]) - abs(coef[orig+p])

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
  return(W)
}

psps_dlasso <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, parallel = T, boot = 200){
  source("./Fun.R")
  lab <- data.frame(Y = Y_lab, Y_hat = Yhat_lab, X_lab)
  unlab <- data.frame(Y_hat = Yhat_unlab, X_unlab)
  p <- ncol(X_lab)
  n <- nrow(X_lab)
  N <- nrow(X_unlab)

  # Change to perellel version
  if (parallel == T){
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)

    boot_result <- foreach(
    i = 1:boot, .combine = rbind, .packages = c("POPInf"),
    .errorhandling = "pass"
    ) %dopar% {
        source("./Fun.R")
        boot_lab <- lab[sample(1:n, n, replace = T), ]
        boot_unlab <- unlab[sample(1:N, N, replace = T), ]
        boot_X_lab <- boot_lab[, 3:ncol(boot_lab)]
        boot_X_unlab <- boot_unlab[, 2:ncol(boot_unlab)]
        boot_Y_lab <- boot_lab[, 1]
        boot_Yhat_lab <- boot_lab[, 2]
        boot_Yhat_unlab <- boot_unlab[, 1]
        est_lab_y_boot <- dlasso(as.matrix(boot_X_lab), as.matrix(boot_Y_lab))$unb.coef
        est_lab_yhat_boot <- dlasso(as.matrix(boot_X_lab), as.matrix(boot_Yhat_lab))$unb.coef
        est_unlab_yhat_boot <- dlasso(as.matrix(boot_X_unlab), as.matrix(boot_Yhat_unlab))$unb.coef
        df <- data.frame(est_lab_y_boot,
                        est_lab_yhat_boot,
                        est_unlab_yhat_boot
                      )
        t(df)
    }
    # print(head(boot_result))
    lab_y <- boot_result[seq(1, 3*boot, 3), ]
    lab_yhat <- boot_result[seq(2, 3*boot, 3), ]
    unlab_yhat <- boot_result[seq(3, 3*boot, 3), ]
    on.exit(stopCluster(cl))
  } else{
    est_lab_y_boot <- matrix(rep(0, boot), nrow = boot, ncol = p)
    est_lab_yhat_boot <- matrix(rep(0, boot), nrow = boot, ncol = p)
    est_unlab_yhat_boot <- matrix(rep(0, boot), nrow = boot, ncol = p)
    for (i in 1:boot){
        boot_lab <- lab[sample(1:n, n, replace = T), ]
        boot_unlab <- unlab[sample(1:N, N, replace = T), ]
        boot_X_lab <- boot_lab[, 3:ncol(boot_lab)]
        boot_X_unlab <- boot_unlab[, 2:ncol(boot_unlab)]
        boot_Y_lab <- boot_lab[, 1]
        boot_Yhat_lab <- boot_lab[, 2]
        boot_Yhat_unlab <- boot_unlab[, 1]
        est_lab_y_boot[i, ] <- dlasso(as.matrix(boot_X_lab), as.matrix(boot_Y_lab))$unb.coef
        est_lab_yhat_boot[i, ] <- dlasso(as.matrix(boot_X_lab), as.matrix(boot_Yhat_lab))$unb.coef
        est_unlab_yhat_boot[i, ] <- dlasso(as.matrix(boot_X_unlab), as.matrix(boot_Yhat_unlab))$unb.coef
    }
    lab_y <- est_lab_y_boot
    lab_yhat <- est_lab_yhat_boot
    unlab_yhat <- est_unlab_yhat_boot
  }

  est_lab_y <- colMeans(lab_y)
  est_lab_yhat <- colMeans(lab_yhat)
  est_unlab_yhat <- colMeans(unlab_yhat)

  # Apply PSPS
  est_psps <- c()
  se_psps <- c()
  for (j in 1:p){
      lab_y_tmp <- lab_y[, j]
      lab_yhat_tmp <- lab_yhat[, j]
      unlab_yhat_tmp <- unlab_yhat[, j]
      df_est <- data.frame(lab_y_tmp, lab_yhat_tmp)
      Sigma <- matrix(rep(0, 9), nrow = 3, ncol = 3)
      Sigma[1:2, 1:2] <- cov(df_est)
      Sigma[3, 3] <- var(unlab_yhat_tmp)
      sqrt(diag(Sigma))
      fit_psps <- PSPS(est_lab_y[j], est_lab_yhat[j], est_unlab_yhat[j], Sigma, alpha = 0.05)
      est_psps <- c(est_psps, fit_psps$Estimate)
      se_psps <- c(se_psps, fit_psps$Std.Error)
  }
  # P-value
  p_psps <-  2*pnorm(abs(est_psps/se_psps), lower.tail = F)
  return(list(coef = est_psps, se = se_psps, pval = p_psps))
}

knockoff_own <- function(X_lab_combined, Y_lab, fdr = 0.05, knockoffs=create.second_order, family = "gaussian", parallel = T, boot = 200, offset=1, swap){

  knockoffs=create.second_order
  # Input parameter
  n = nrow(X_lab_combined)
  p <- ncol(X_lab_combined)/2

  # Compute Statistics
  parallel <- T
  boot <- 200
  W <- dlasso_stats(X_lab_combined, Y_lab, swap)

  # Run the knockoff filter
  # offset=1
  # fdr <- 0.1
  t = knockoff.threshold(W, fdr=fdr, offset=offset)
  selected = sort(which(W >= t))
  X.names = colnames(X)
  names(selected) = X.names[selected]
  fdp(selected)
  power(selected)
  
  # Package up the results.
  return(list(selected = selected, W = W, t = t, fdr = fdr))
}

coef_to_W <- function(coef, swap){
  p <- length(coef)/2
  orig <- 1:p
  W = abs(coef[orig]) - abs(coef[orig+p])

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
  return(W)
}

W_to_variable <- function(W, fdr = 0.1, offset = 1){
  t = knockoff.threshold(W, fdr=fdr, offset=offset)
  selected = sort(which(W >= t))
  return(selected)
}

W_to_t <- function(W, fdr = 0.1, offset = 1){
  t = knockoff.threshold(W, fdr=fdr, offset=offset)
  return(t)
}
