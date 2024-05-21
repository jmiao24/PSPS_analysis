require(POPInf)
require(quantreg)
library(doParallel)
require(data.table)
source("./Fun.R")

sim.times <- 1000
ratio_vec <- c(2, 5, 10, 20)

for (ratio in ratio_vec){
#-------------- Known Model ------
cl <- makeCluster(detectCores())
registerDoParallel(cl)
result <- foreach(
i = 1:sim.times, .combine = rbind, .packages = c("POPInf","quantreg", "sandwich", "lmtest"),
.errorhandling = "pass"
) %dopar% {
    r <- 0.9
    data <- sim_data_quantile(r = r)
    X_lab = data$X_lab
    X_unlab = data$X_unlab
    Y_lab = data$Y_lab
    Yhat_lab = data$Yhat_lab
    Yhat_unlab = data$Yhat_unlab
    Y_unlab = data$Y_unlab

    full <- data.frame(X = c(unlist(X_lab), unlist(X_unlab)), Y = c(unlist(Y_lab), unlist(Y_unlab)))

    n_lab <- nrow(X_lab)
    n_unlab <- ratio * n_lab
    Yhat_unlab <- Yhat_unlab[1:n_unlab, ]
    X_unlab <- X_unlab[1:n_unlab, ]

    lab <- data.frame(Y_lab, Yhat_lab, X_lab)
    colnames(lab) <- c("Y", "Yhat", "X")
    unlab <- data.frame(Yhat_unlab, X_unlab)
    colnames(unlab) <- c("Yhat", "X")
    
    est_lab_y <- summary(rq(Y ~ X, data = lab, tau = 0.75))$coefficients[2, 1]      
    est_lab_yhat <- summary(rq(Yhat ~ X, data = lab, tau = 0.75))$coefficients[2, 1] 
    est_unlab_yhat <- summary(rq(Yhat ~ X, data = unlab, tau = 0.75))$coefficients[2, 1] 
    n <- nrow(lab)
    N <- nrow(unlab)

    # Calcualte the covariance between est_lab_y and est_lab_yhat using bootstrap
    B <- 200
    est_lab_y_boot <- matrix(rep(0, B), nrow = B, ncol = 1)
    est_lab_yhat_boot <- matrix(rep(0, B), nrow = B, ncol = 1)
    est_unlab_yhat_boot <- matrix(rep(0, B), nrow = B, ncol = 1)
    for (i in 1:B){
        boot_lab <- lab[sample(1:n, n, replace = T), ]
        boot_unlab <- unlab[sample(1:N, N, replace = T), ]
        est_lab_y_boot[i, ] <- rq(Y ~ X, data = boot_lab, tau = 0.75)$coefficients[2]
        est_lab_yhat_boot[i, ] <- rq(Yhat ~ X, data = boot_lab, tau = 0.75)$coefficients[2]
        est_unlab_yhat_boot[i, ] <- rq(Yhat ~ X, data = boot_unlab, tau = 0.75)$coefficients[2]
    }
    df_est <- data.frame(est_lab_y_boot, est_lab_yhat_boot)

    Sigma <- matrix(rep(0, 9), nrow = 3, ncol = 3)
    Sigma[1:2, 1:2] <- cov(df_est)
    Sigma[3, 3] <- var(est_unlab_yhat_boot)
    sqrt(diag(Sigma))

    # Gold-standard results
    true_beta <- summary(rq(Y ~ X, data = full, tau = 0.75))$coefficients[2, 1]

    # Classic method
    est_classic <- est_lab_y
    se_classic <- sqrt(Sigma[1, 1])
    ind <- est_classic - 1.96 * se_classic < true_beta & true_beta < est_classic + 1.96 * se_classic
    cov_classic <- ifelse(ind, 1, 0)

    # Imputed method
    est_imputed <- summary(rq(Yhat ~ X, data = unlab, tau = 0.75))$coefficients[2, 1]
    se_imputed <- sqrt(Sigma[3, 3])
    ind <- est_imputed - 1.96 * se_imputed < true_beta & true_beta < est_imputed + 1.96 * se_imputed
    cov_imputed <- ifelse(ind, 1, 0)

    # PSPS method
    output <- PSPS(est_lab_y, est_lab_yhat, est_unlab_yhat, Sigma, alpha = 0.05)
    est_psps <- output$Estimate
    se_psps <- output$Std.Error
    ind <- est_psps - 1.96 * se_psps < true_beta & true_beta < est_psps + 1.96 * se_psps
    cov_psps <- ifelse(ind, 1, 0)

    df <- data.frame(est_classic, se_classic, cov_classic,
                     est_imputed, se_imputed, cov_imputed,
                     est_psps, se_psps, cov_psps
                    )
    df
}
stopCluster(cl)

fwrite(as.data.frame(result), paste0("./results/qr_", ratio, ".txt"),
sep = "\t", quote = F, row.names = F, col.names = T)
}


