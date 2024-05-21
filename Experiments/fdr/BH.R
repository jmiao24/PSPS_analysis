rm(list = ls())
library(ranger)
library(data.table)
source("./Fun.R")
source("./sslasso_code/lasso_inference.r")
source("./Fun_PSPS-knockoff.R")

# False 
fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
# Power
power <- function(selected) sum(beta[selected] != 0) / sum(beta != 0)

set.seed(42)
sim.times <- 200
n = 10^4        # number of observations
p = 150        # number of variables
k = 30            # number of variables with nonzero coefficients
amplitude = 1.5   # signal amplitude (for noise level = 1)

ratio <- 10
r <- 0.9

#-------------- Known Model ------
cl <- makeCluster(detectCores())
registerDoParallel(cl)
result <- foreach(
i = 1:sim.times, .combine = rbind, .packages = c("POPInf", "IPD", "ranger", "data.table"),
.errorhandling = "pass"
) %dopar% {

    # False 
    fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
    # Power
    power <- function(selected) sum(beta[selected] != 0) / sum(beta != 0)

    n = 10^4 
    mu = rep(0,p)
    rho = 0.25
    Sigma = toeplitz(rho^(0:(p-1)))
    X = matrix(rnorm(n*p),n) %*% chol(Sigma)
    nonzero = sample(1:p, k)
    beta = amplitude * (1:p %in% nonzero) / sqrt(100)
    y.sample = function(X) X %*% beta + rnorm(n)
    y = y.sample(X)
    theta_0 <- rep(beta[1], p)
    Z <- r*y + X %*% theta_0  + rnorm(n)
    n_train <- 500
    n_lab <- 500
    n_unlab <- ratio * n_lab
    full <- data.frame(y = y, Z, X)
    dat_train <- full[1:n_train, ]
    dat_lab <- full[(n_train + 1):(n_train + n_lab), ]
    dat_unlab <- full[(n_train + n_lab + 1):(n_train + n_lab + n_unlab), ]

    # Train random forest
    rf <- ranger::ranger(y ~ Z, data = dat_train, num.trees = 100)
    Yhat_lab <- predict(rf, data = data.frame(dat_lab))$predictions
    Yhat_unlab <- predict(rf, data = data.frame(dat_unlab))$predictions

    # Generate
    X_lab <- dat_lab[, -c(1, 2)]
    X_unlab <- dat_unlab[, -c(1, 2)]

    dat_lab$Yhat <- Yhat_lab
    dat_unlab$Yhat <- Yhat_unlab
    Yhat_lab <- dat_lab$Yhat
    Yhat_unlab <- dat_unlab$Yhat
    Y_lab <- dat_lab$y

    lab_Y <- data.frame(Y_lab, X_lab)
    colnames(lab_Y)[1] <- c("Y")
    lab_Yhat <- data.frame(Yhat_lab, X_lab)
    colnames(lab_Yhat)[1] <- c("Yhat")
    unlab <- data.frame(Yhat_unlab, X_unlab)
    colnames(unlab)[1] <- c("Yhat")

    N <- nrow(dat_unlab)
    n <- nrow(dat_lab)

    B <- 200
    est_lab_y_boot <- matrix(rep(0, B), nrow = B, ncol = p)
    est_lab_yhat_boot <- matrix(rep(0, B), nrow = B, ncol = p)
    est_unlab_yhat_boot <- matrix(rep(0, B), nrow = B, ncol = p)
    for (i in 1:B){
        index_lab <- sample(1:n, n, replace = T)
        index_unlab <- sample(1:N, N, replace = T)
        boot_lab_Y <- lab_Y[index_lab, ]
        boot_lab_Yhat <- lab_Yhat[index_lab, ]
        boot_unlab <- unlab[index_unlab, ]
        est_lab_y_boot[i, ] <- lm(Y ~ ., data = boot_lab_Y)$coefficients[-1]
        est_lab_yhat_boot[i, ] <- lm(Yhat ~ ., data = boot_lab_Yhat)$coefficients[-1]
        est_unlab_yhat_boot[i, ] <- lm(Yhat ~ ., data = boot_unlab)$coefficients[-1]
    }

    est_lab_y <- colMeans(est_lab_y_boot)
    se_lab_y <- apply(est_lab_y_boot, 2, sd)
    est_lab_yhat <- colMeans(est_lab_yhat_boot)
    est_unlab_yhat <- colMeans(est_unlab_yhat_boot)

    est_lab_y_boot <- as.data.frame(est_lab_y_boot)
    est_lab_yhat_boot <- as.data.frame(est_lab_yhat_boot)
    est_unlab_yhat_boot <- as.data.frame(est_unlab_yhat_boot)

    est_psps <- c()
    se_psps <- c()
    for (q in 1:p){
        lab_y_tmp <- est_lab_y_boot[, q]
        lab_yhat_tmp <- est_lab_yhat_boot[, q]
        unlab_yhat_tmp <- est_unlab_yhat_boot[, q]
        df_est <- data.frame(lab_y_tmp, lab_yhat_tmp)
        Sigma <- matrix(rep(0, 9), nrow = 3, ncol = 3)
        Sigma[1:2, 1:2] <- cov(df_est)
        Sigma[3, 3] <- var(unlab_yhat_tmp)
        sqrt(diag(Sigma))
        fit_psps <- PSPS(est_lab_y[q], est_lab_yhat[q], est_unlab_yhat[q], Sigma, alpha = 0.05)
        est_psps <- c(est_psps, fit_psps$Estimate)
        se_psps <- c(se_psps, fit_psps$Std.Error)
    }

    p_psps <- 2 * pnorm(-abs(est_psps / se_psps))

    # Fit classic linear regression
    fit_classic <- as.data.frame(summary(lm(Y_lab ~ as.matrix(X_lab)))$coefficient)
    fit_classic <- fit_classic[-1, ]
    p_classic <- fit_classic[, 4]
    fit_classic$fdr <- p.adjust(p_classic, method = "fdr")

    # Get the FDR results for 0.05, 0.1, 0.15
    fdr_classic <- p.adjust(p_classic, method = "fdr")
    selected_classic_05 <- which(fdr_classic < 0.05)
    selected_classic_10 <- which(fdr_classic < 0.1)
    selected_classic_15 <- which(fdr_classic < 0.15)
    # Power and FDP
    power_classic_05 <- power(selected_classic_05)
    power_classic_10 <- power(selected_classic_10)
    power_classic_15 <- power(selected_classic_15)
    fdp_classic_05 <- fdp(selected_classic_05)
    fdp_classic_10 <- fdp(selected_classic_10)
    fdp_classic_15 <- fdp(selected_classic_15)

    # Fit imputed linear regression
    fit_imputed <- as.data.frame(summary(lm(Yhat ~ ., data = unlab))$coefficient)
    fit_imputed <- fit_imputed[-1, ]
    p_imputed <- fit_imputed[, 4]
    fit_imputed$fdr <- p.adjust(p_imputed, method = "fdr")
    selected_imputed_05 <- which(fit_imputed$fdr < 0.05)
    selected_imputed_10 <- which(fit_imputed$fdr < 0.1)
    selected_imputed_15 <- which(fit_imputed$fdr < 0.15)
    power_imputed_05 <- power(selected_imputed_05)
    power_imputed_10 <- power(selected_imputed_10)
    power_imputed_15 <- power(selected_imputed_15)
    fdp_imputed_05 <- fdp(selected_imputed_05)
    fdp_imputed_10 <- fdp(selected_imputed_10)
    fdp_imputed_15 <- fdp(selected_imputed_15)

    # PSPS
    fdr_psps <- p.adjust(p_psps, method = "fdr")
    selected_psps_05 <- which(fdr_psps < 0.05)
    selected_psps_10 <- which(fdr_psps < 0.1)
    selected_psps_15 <- which(fdr_psps < 0.15)
    power_psps_05 <- power(selected_psps_05)
    power_psps_10 <- power(selected_psps_10)
    power_psps_15 <- power(selected_psps_15)
    fdp_psps_05 <- fdp(selected_psps_05)
    fdp_psps_10 <- fdp(selected_psps_10)
    fdp_psps_15 <- fdp(selected_psps_15)

    # Return the results
    df <- data.frame(
        power_psps_05, fdp_psps_05, power_classic_05, fdp_classic_05, power_imputed_05, fdp_imputed_05,
        power_psps_10, fdp_psps_10, power_classic_10, fdp_classic_10, power_imputed_10, fdp_imputed_10,
        power_psps_15, fdp_psps_15, power_classic_15, fdp_classic_15, power_imputed_15, fdp_imputed_15
    )
    df

}
stopCluster(cl)

fwrite(as.data.frame(result), paste0("./results/BH.txt"),
sep = "\t", quote = F, row.names = F, col.names = T)