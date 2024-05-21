require(POPInf)
require(IPD)
require(sandwich)
require(lmtest)
library(doParallel)
require(data.table)
source("./Fun.R")
source("./sslasso_code/lasso_inference.r")

set.seed(42)
sim.times <- 1000

ratio_vec <- c(15, 20, 25, 30)

for (ratio in ratio_vec){
    #-------------- Known Model ------
cl <- makeCluster(detectCores())
registerDoParallel(cl)
result <- foreach(
i = 1:sim.times, .combine = rbind, .packages = c("POPInf", "IPD", "sandwich", "lmtest"),
.errorhandling = "pass"
) %dopar% {
    r <- 0.9
    data <- sim_data_dlasso(r = r)
    X_lab = data$X_lab
    X_unlab = data$X_unlab
    Y_lab = data$Y_lab
    Yhat_lab = data$Yhat_lab
    Yhat_unlab = data$Yhat_unlab
    Y_unlab = data$Y_unlab

    lab <- cbind(cbind(Y_lab, Yhat_lab), X_lab)
    unlab <- cbind(cbind(Y_unlab, Yhat_unlab), X_unlab)
    colnames(unlab)[2] <- "Y_hat"
    full <- rbind(lab, unlab)

    n_lab <- nrow(Y_lab)
    n_unlab <- n_lab * ratio
    
    Yhat_unlab <- as.data.frame(Yhat_unlab)[1:n_unlab, ]
    X_unlab <- X_unlab[1:n_unlab, ]

    # Produce summary statistics
    lab <- data.frame(Y_lab, Yhat_lab, X_lab)
    colnames(lab) <- c("Y", "Yhat", "X")
    unlab <- data.frame(Yhat_unlab, X_unlab)
    colnames(unlab) <- c("Yhat", "X")
    fit_lab_y <- SSLasso(as.matrix(X_lab), as.matrix(Y_lab), verbose = TRUE)
    fit_lab_yhat <- SSLasso(as.matrix(X_lab), as.matrix(Yhat_lab), verbose = TRUE)
    fit_unlab_yhat <- SSLasso(as.matrix(X_unlab), as.matrix(Yhat_unlab), verbose = TRUE)

    est_lab_y <- fit_lab_y$unb.coef
    est_lab_yhat <- fit_lab_yhat$unb.coef
    est_unlab_yhat <- fit_unlab_yhat$unb.coef
    n <- nrow(lab)
    N <- nrow(unlab)


    # Calcualte the covariance between est_lab_y and est_lab_yhat using bootstrap
    B <- 200
    est_lab_y_boot <- matrix(rep(0, B), nrow = B, ncol = 200)
    est_lab_yhat_boot <- matrix(rep(0, B), nrow = B, ncol = 200)
    est_unlab_yhat_boot <- matrix(rep(0, B), nrow = B, ncol = 200)
    for (i in 1:B){
        print(i)
        boot_lab <- lab[sample(1:n, n, replace = T), ]
        boot_unlab <- unlab[sample(1:N, N, replace = T), ]
        boot_X_lab <- boot_lab[, 3:ncol(boot_lab)]
        boot_X_unlab <- boot_unlab[, 2:ncol(boot_unlab)]
        boot_Y_lab <- boot_lab[, 1]
        boot_Yhat_lab <- boot_lab[, 2]
        boot_Yhat_unlab <- boot_unlab[, 1]
        est_lab_y_boot[i, ] <- SSLasso(as.matrix(boot_X_lab), as.matrix(boot_Y_lab), verbose = FALSE)$unb.coef
        est_lab_yhat_boot[i, ] <- SSLasso(as.matrix(boot_X_lab), as.matrix(boot_Yhat_lab), verbose = FALSE)$unb.coef
        est_unlab_yhat_boot[i, ] <- SSLasso(as.matrix(boot_X_unlab), as.matrix(boot_Yhat_unlab), verbose = FALSE)$unb.coef
    }

    # For each parameter, calculate the covariance
    est_psps <- c()
    se_psps <- c()
    for (j in 1:200){
        est_lab_y_boot_tmp <- est_lab_y_boot[, j]
        est_lab_yhat_boot_tmp <- est_lab_yhat_boot[, j]
        est_unlab_yhat_boot_tmp <- est_unlab_yhat_boot[, j]
        df_est <- data.frame(est_lab_y_boot_tmp, est_lab_yhat_boot_tmp)
        Sigma <- matrix(rep(0, 9), nrow = 3, ncol = 3)
        Sigma[1:2, 1:2] <- cov(df_est)
        Sigma[3, 3] <- var(est_unlab_yhat_boot_tmp)
        sqrt(diag(Sigma))

        fit_psps <- PSPS(mean(est_lab_y_boot_tmp), mean(est_lab_yhat_boot_tmp), mean(est_unlab_yhat_boot_tmp), Sigma, alpha = 0.05)
        est_psps <- c(est_psps, fit_psps$Estimate)
        se_psps <- c(se_psps, fit_psps$Std.Error)
    }
    # Ground truth for beta
    true_beta <- lm(Y ~ X1, data = full)$coefficients[2]
    true_beta <- c(rep(true_beta, 5), rep(0, 200 - 5))
    
    width_psps <- mean(se_psps)
    ind_psps <- est_psps - 1.96 * se_psps < true_beta & true_beta < est_psps + 1.96 * se_psps
    cov_psps <-  sum(ind_psps)/length(ind_psps)


    # Classic method
    est_classic <- est_lab_y
    width_classic <- mean((fit_lab_y$up.lim - fit_lab_y$low.lim)/4)
    ind <- fit_lab_yhat$low.lim < true_beta & true_beta < fit_lab_y$up.lim
    cov_classic <- sum(ind)/length(ind)

    # Imputed method
    est_imputed <- est_unlab_yhat
    width_imputed <- mean((fit_unlab_yhat$up.lim - fit_unlab_yhat$low.lim)/4)
    ind <- fit_unlab_yhat$low.lim < true_beta & true_beta < fit_unlab_yhat$up.lim
    cov_imputed <- sum(ind)/length(ind)

    df <- data.frame(width_classic, cov_classic,
                    width_imputed, cov_imputed,
                    width_psps, cov_psps
                    )
    df
}
stopCluster(cl)

fwrite(as.data.frame(result), paste0("./results/dlasso_", ratio, ".txt"),
sep = "\t", quote = F, row.names = F, col.names = T)

}
