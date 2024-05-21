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
    data <- sim_data_wilcox_test(r = r, effect = F)
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

    # Measure the imputation quality
    corr <- cor(Yhat_lab, Y_lab)

    # Create the data frame
    lab_Y_0 <- lab[lab$X == 0, "Y"]
    lab_Y_1 <- lab[lab$X == 1, "Y"]
    lab_Yhat_0 <- lab[lab$X == 0, "Yhat"]
    lab_Yhat_1 <- lab[lab$X == 1, "Yhat"]
    unlab_Yhat_0 <- unlab[unlab$X == 0, "Yhat"]
    unlab_Yhat_1 <- unlab[unlab$X == 1, "Yhat"]

    est_lab_y <- wilcox_test(lab_Y_0, lab_Y_1)$Estimate
    est_lab_yhat <- wilcox_test(lab_Yhat_0, lab_Yhat_1)$Estimate
    est_unlab_yhat <- wilcox_test(unlab_Yhat_0, unlab_Yhat_1)$Estimate

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

        boot_lab_Y_0 <- boot_lab[boot_lab$X == 0, "Y"]
        boot_lab_Y_1 <- boot_lab[boot_lab$X == 1, "Y"]
        boot_lab_Yhat_0 <- boot_lab[boot_lab$X == 0, "Yhat"]
        boot_lab_Yhat_1 <- boot_lab[boot_lab$X == 1, "Yhat"]
        boot_unlab_Yhat_0 <- boot_unlab[boot_unlab$X == 0, "Yhat"]
        boot_unlab_Yhat_1 <- boot_unlab[boot_unlab$X == 1, "Yhat"]

        est_lab_y_boot[i, ] <- wilcox_test(boot_lab_Y_0, boot_lab_Y_1)$Estimate
        est_lab_yhat_boot[i, ] <- wilcox_test(boot_lab_Yhat_0, boot_lab_Yhat_1)$Estimate
        est_unlab_yhat_boot[i, ] <- wilcox_test(boot_unlab_Yhat_0, boot_unlab_Yhat_1)$Estimate
    }
    df_est <- data.frame(est_lab_y_boot, est_lab_yhat_boot)

    Sigma <- matrix(rep(0, 9), nrow = 3, ncol = 3)
    Sigma[1:2, 1:2] <- cov(df_est)
    Sigma[3, 3] <- var(est_unlab_yhat_boot)
    sqrt(diag(Sigma))

    # Classic method
    est_classic <- mean(est_lab_y_boot)
    se_classic <- sqrt(Sigma[1, 1])
    p_classic <- 2 * pnorm(-abs(est_classic / se_classic))

    # Imputed method
    est_imputed <- mean(est_unlab_yhat_boot)
    se_imputed <- sqrt(Sigma[3, 3])
    p_imputed <- 2 * pnorm(-abs(est_imputed / se_imputed))

    # PSPS
    est_lab_y <- mean(est_lab_y_boot)
    est_lab_yhat <- mean(est_lab_yhat_boot)
    est_unlab_yhat <- mean(est_unlab_yhat_boot)
    fit_psps <- PSPS(est_lab_y, est_lab_yhat, est_unlab_yhat, Sigma)
    est_psps <- fit_psps$Estimate
    se_psps <- fit_psps$Std.Error
    p_psps <- fit_psps$P.value
    
    df <- data.frame(est_classic, se_classic, p_classic, 
                     est_imputed, se_imputed, p_imputed, 
                     est_psps, se_psps, p_psps,
                     corr)
    colnames(df) <- c("est_classic", "se_classic", "p_classic", 
                      "est_imputed", "se_imputed", "p_imputed", 
                      "est_psps", "se_psps", "p_psps",
                      "corr")
    df
}
stopCluster(cl)

fwrite(as.data.frame(result), paste0("./results/rank_sum_t1e_", ratio, ".txt"),
sep = "\t", quote = F, row.names = F, col.names = T)

}
