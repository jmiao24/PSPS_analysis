rm(list = ls())
require(data.table)
source("./Fun.R")
source("./sslasso_code/lasso_inference.r")
source("./Fun_PSPS-knockoff.R")

sim.times <- 1000
df_out <- c()
for (j in 1:sim.times){
    print(j)
    y_lab_tmp <- fread(paste0("tmp_results/est_lab_y_sim_", j, "_boot.txt"))
    y_hat_lab <- fread(paste0("tmp_results/est_lab_yhat_sim_", j, "_boot.txt"))
    y_hat_unlab <- fread(paste0("tmp_results/est_unlab_yhat_sim_", j, "_boot.txt"))

    # Try PSPS
    # For each parameter, calculate the covariance
    est_lab_y <- colMeans(lab_y)
    est_lab_yhat <- colMeans(lab_yhat)
    est_unlab_yhat <- colMeans(unlab_yhat)

    lab_y <- as.data.frame(lab_y)
    lab_yhat <- as.data.frame(lab_yhat)
    unlab_yhat <- as.data.frame(unlab_yhat)

    est_psps <- c()
    se_psps <- c()
    for (q in 1:200){
        lab_y_tmp <- lab_y[, q]
        lab_yhat_tmp <- lab_yhat[, q]
        unlab_yhat_tmp <- unlab_yhat[, q]
        df_est <- data.frame(lab_y_tmp, lab_yhat_tmp)
        Sigma <- matrix(rep(0, 9), nrow = 3, ncol = 3)
        Sigma[1:2, 1:2] <- cov(df_est)
        Sigma[3, 3] <- var(unlab_yhat_tmp)
        sqrt(diag(Sigma))
        fit_psps <- PSPS(est_lab_y[q], est_lab_yhat[q], est_unlab_yhat[q], Sigma, alpha = 0.05)
        est_psps <- c(est_psps, fit_psps$Estimate)
        se_psps <- c(se_psps, fit_psps$Std.Error)
    }
    k <- 30
    p <- 150
    amplitude <- 5
    nonzero <- 1:k
    beta <- amplitude * (1:p %in% nonzero) / sqrt(100)

    swap <- as.numeric(unlist(fread(paste0(".data/swap_", j, ".txt"))[, 1])) 

    # PSPS results
    coef <- est_psps
    W <- coef_to_W(coef, swap)
    selected_05 <- W_to_variable(W, fdr = 0.05)
    selected_10 <- W_to_variable(W, fdr = 0.1)
    selected_15 <- W_to_variable(W, fdr = 0.15)
    power_psps_05 <- power(selected_05)
    fdp_psps_05 <- fdp(selected_05)
    power_psps_10 <- power(selected_10)
    fdp_psps_10 <- fdp(selected_10)
    power_psps_15 <- power(selected_15)
    fdp_psps_15 <- fdp(selected_15)

    # Classic results
    coef <- est_lab_y
    W <- coef_to_W(coef, swap)
    selected_05 <- W_to_variable(W, fdr = 0.05)
    selected_10 <- W_to_variable(W, fdr = 0.1)
    selected_15 <- W_to_variable(W, fdr = 0.15)
    power_classic_05 <- power(selected_05)
    fdp_classic_05 <- fdp(selected_05)
    power_classic_10 <- power(selected_10)
    fdp_classic_10 <- fdp(selected_10)
    power_classic_15 <- power(selected_15)
    fdp_classic_15 <- fdp(selected_15)

    # Imputed results
    coef <- est_unlab_yhat
    W <- coef_to_W(coef, swap)
    selected_05 <- W_to_variable(W, fdr = 0.05)
    selected_10 <- W_to_variable(W, fdr = 0.1)
    selected_15 <- W_to_variable(W, fdr = 0.15)
    power_imputed_05 <- power(selected_05)
    fdp_imputed_05 <- fdp(selected_05)
    power_imputed_10 <- power(selected_10)
    fdp_imputed_10 <- fdp(selected_10)
    power_imputed_15 <- power(selected_15)
    fdp_imputed_15 <- fdp(selected_15)
    

    # Compare with dlasso
    X_lab_combined <- fread(paste0("./data/X_lab_combined_", j, ".txt"))
    X_unlab_combined <- fread(paste0("./data/X_unlab_combined_", j, ".txt"))
    swap <- as.numeric(unlist(fread(paste0("./data/swap_", j, ".txt"))[, 1]))
    lab <- fread(paste0("./data/lab_", j, ".txt"))
    unlab <- fread(paste0("./data/unlab_", j, ".txt"))
    Y_lab <- lab[, 1]
    Yhat_lab <- lab[, ncol(lab)]
    Yhat_unlab <- unlab[, ncol(unlab)]
    coef <- dlasso(as.matrix(X_lab_combined), as.matrix(Y_lab))$unb.coef
    W <- coef_to_W(coef, swap)
    selected_05 <- W_to_variable(W, fdr = 0.05)
    selected_10 <- W_to_variable(W, fdr = 0.1)
    selected_15 <- W_to_variable(W, fdr = 0.15)
    power_classic_old_05 <- power(selected_05)
    fdp_classic_old_05 <- fdp(selected_05)
    power_classic_old_10 <- power(selected_10)
    fdp_classic_old_10 <- fdp(selected_10)
    power_classic_old_15 <- power(selected_15)
    fdp_classic_old_15 <- fdp(selected_15)


    df_tmp <- data.frame(power_psps_05, fdp_psps_05, power_psps_10, fdp_psps_10, power_psps_15, fdp_psps_15, power_classic_05, fdp_classic_05, power_classic_10, fdp_classic_10, power_classic_15, fdp_classic_15, power_imputed_05, fdp_imputed_05, power_imputed_10, fdp_imputed_10, power_imputed_15, fdp_imputed_15, power_classic_old_05, fdp_classic_old_05, power_classic_old_10, fdp_classic_old_10, power_classic_old_15, fdp_classic_old_15)
    df_out <- rbind(df_out, df_tmp)
}
colMeans(df_out)
fwrite(df_out, "./results/knockoff_results.txt",
sep = "\t", quote = F, row.names = F, col.names = T)
