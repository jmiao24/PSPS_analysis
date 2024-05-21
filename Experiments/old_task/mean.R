require(POPInf)
require(IPD)
require(sandwich)
require(lmtest)
library(doParallel)
require(data.table)
source("./Fun.R")

set.seed(42)
sim.times <- 1000

ratio_vec <- c(2, 5, 10, 20)

for (ratio in ratio_vec){
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    result <- foreach(
i = 1:sim.times, .combine = rbind, .packages = c("POPInf", "IPD", "sandwich", "lmtest"),
.errorhandling = "pass"
) %dopar% {
    r <- 0.9
    data <- sim_data(r = r)
    X_lab = data$X_lab
    X_unlab = data$X_unlab
    Y_lab = data$Y_lab
    Yhat_lab = data$Yhat_lab
    Yhat_unlab = data$Yhat_unlab
    Y_unlab = data$Y_unlab

    full <- data.frame(X = c(unlist(X_lab), unlist(X_unlab)), Y = c(unlist(Y_lab), unlist(Y_unlab)))

    n_lab <- nrow(Y_lab)
    n_unlab <- n_lab * ratio

    Yhat_unlab <- Yhat_unlab[1:n_unlab, ]
    X_unlab <- X_unlab[1:n_unlab, ]

    # Produce summary statistics
    lab <- data.frame(Y_lab, Yhat_lab, X_lab)
    colnames(lab) <- c("Y", "Yhat", "X")
    unlab <- data.frame(Yhat_unlab, X_unlab)
    colnames(unlab) <- c("Yhat", "X")
    est_lab_y <- mean(lab$Y)
    est_lab_yhat <- mean(lab$Yhat)
    est_unlab_yhat <- mean(unlab$Yhat)
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
        est_lab_y_boot[i, ] <- mean(boot_lab$Y)
        est_lab_yhat_boot[i, ] <- mean(boot_lab$Yhat)
        est_unlab_yhat_boot[i, ] <- mean(boot_unlab$Yhat)
    }
    df_est <- data.frame(est_lab_y_boot, est_lab_yhat_boot)

    Sigma <- matrix(rep(0, 9), nrow = 3, ncol = 3)
    Sigma[1:2, 1:2] <- cov(df_est)
    Sigma[3, 3] <- var(est_unlab_yhat_boot)
    sqrt(diag(Sigma))

    # Ground truth for beta
    true_beta <- mean(full$Y)

    # Classic method
    est_classic <- est_lab_y
    se_classic <- sqrt(Sigma[1, 1])
    ind <- est_classic - 1.96 * se_classic < true_beta & true_beta < est_classic + 1.96 * se_classic
    cov_classic <- ifelse(ind, 1, 0)

    # Imputed method
    est_imputed <- est_unlab_yhat
    se_imputed <- sqrt(Sigma[3, 3])
    ind <- est_imputed - 1.96 * se_imputed < true_beta & true_beta < est_imputed + 1.96 * se_imputed
    cov_imputed <- ifelse(ind, 1, 0)

    # PPI
    Y_lab <- as.matrix(lab$Y)
    Yhat_lab <- as.matrix(lab$Yhat)
    Yhat_unlab <- as.matrix(unlab$Yhat)

    fit_ppi <- IPD::ppi_plusplus_mean(Y_lab, Yhat_lab, Yhat_unlab, lhat = 1)
    est_ppi <- mean(fit_ppi)
    se_ppi <- (fit_ppi[2] - fit_ppi[1])/1.96/2
    ind <- est_ppi - 1.96 * se_ppi < true_beta & true_beta < est_ppi + 1.96 * se_ppi
    cov_ppi <- ifelse(ind, 1, 0)

    # PPI++
    fit_ppi_plus_plus <- IPD::ppi_plusplus_mean(Y_lab, Yhat_lab, Yhat_unlab)
    est_ppi_plus_plus <- mean(fit_ppi)
    se_ppi_plus_plus <- (fit_ppi[2] - fit_ppi[1])/1.96/2
    ind <- est_ppi_plus_plus - 1.96 * se_ppi_plus_plus < true_beta & true_beta < est_ppi_plus_plus + 1.96 * se_ppi_plus_plus
    cov_ppi_plus_plus <- ifelse(ind, 1, 0)

    # POP-Inf
    fit_popinf <- pop_M(Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
                  alpha = 0.05, method = "mean")
    est_popinf <- fit_popinf$Estimate
    se_popinf <- fit_popinf$Std.Error
    ind <- est_popinf - 1.96 * se_popinf < true_beta & true_beta < est_popinf + 1.96 * se_popinf
    cov_popinf <-  ifelse(ind, 1, 0)

    # PSPS
    fit_psps <- PSPS(est_lab_y, est_lab_yhat, est_unlab_yhat, Sigma, alpha = 0.05)
    est_psps <- fit_psps$Estimate
    se_psps <-fit_psps$Std.Error
    ind <- est_psps - 1.96 * se_psps < true_beta & true_beta < est_psps + 1.96 * se_psps
    cov_psps <- ifelse(ind, 1, 0)
    
    df <- data.frame(est_classic, se_classic, cov_classic,
                    est_imputed, se_imputed, cov_imputed,
                    est_popinf, se_popinf, cov_popinf,
                    est_psps, se_psps, cov_psps,
                    est_ppi, se_ppi, cov_ppi,
                    est_ppi_plus_plus = est_popinf, se_ppi_plus_plus = se_popinf, cov_ppi_plus_plus = cov_popinf)
    colnames(df) <- c("est_classic", "se_classic", "cov_classic",
                    "est_imputed", "se_imputed", "cov_imputed",
                    "est_popinf", "se_popinf", "cov_popinf",
                    "est_psps", "se_psps", "cov_psps",
                    "est_ppi", "se_ppi", "cov_ppi",
                    "est_ppi_plus_plus", "se_ppi_plus_plus", "cov_ppi_plus_plus")
    df
}
stopCluster(cl)

fwrite(as.data.frame(result), paste0("./results/mean_", ratio, ".txt"),
sep = "\t", quote = F, row.names = F, col.names = T)

}
