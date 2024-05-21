require(data.table)

source("./Fun.R")
source("./Fun_PSPS-knockoff.R")

for (j in 1:1000){
    X_lab_combined <- fread(paste0("./data/X_lab_combined_", j, ".txt"))
    X_unlab_combined <- fread(paste0("./data/X_unlab_combined_", j, ".txt"))
    swap <- as.numeric(unlist(fread(paste0("./data/swap_", j, ".txt"))[, 1]))
    lab <- fread(paste0("./data/lab_", j, ".txt"))
    unlab <- fread(paste0("./data/unlab_", j, ".txt"))
    Y_lab <- lab[, 1, with = F]
    Yhat_lab <- lab[, ncol(lab), with = F]
    Yhat_unlab <- unlab[, ncol(unlab), with = F]

    # Run bootstrapping
    B <- 200
    n <- nrow(X_lab_combined)
    N <- nrow(X_unlab_combined)
    p <- ncol(X_lab_combined)
    est_lab_y_boot <- matrix(rep(0, B), nrow = B, ncol = p)
    est_lab_yhat_boot <- matrix(rep(0, B), nrow = B, ncol = p)
    est_unlab_yhat_boot <- matrix(rep(0, B), nrow = B, ncol = p)
    for (i in 1:B){
        print(i)
        index_lab <- sample(1:n, n, replace = T)
        index_unlab <- sample(1:N, N, replace = T)
        boot_X_lab_combined <- X_lab_combined[index_lab, ]
        boot_X_unlab_combined <- X_unlab_combined[index_unlab, ]
        boot_Y_lab <- Y_lab[index_lab]
        boot_Yhat_lab <- Yhat_lab[index_lab]
        boot_Yhat_unlab <- Yhat_unlab[index_unlab]
        est_lab_y_boot[i, ] <- dlasso(as.matrix(boot_X_lab_combined), as.matrix(boot_Y_lab))$unb.coef
        est_lab_yhat_boot[i, ] <- dlasso(as.matrix(boot_X_lab_combined), as.matrix(boot_Yhat_lab))$unb.coef
        est_unlab_yhat_boot[i, ] <- dlasso(as.matrix(boot_X_unlab_combined), as.matrix(boot_Yhat_unlab))$unb.coef
    }

    fwrite(est_lab_y_boot, 
        paste0("./tmp_results/est_lab_y_sim_", j, "_boot.txt"), 
        sep = "\t", col.names = T)
    fwrite(est_lab_yhat_boot, 
        paste0("./tmp_results/est_lab_yhat_sim_", j, "_boot.txt"), 
        sep = "\t", col.names = T)
    fwrite(est_unlab_yhat_boot, 
        paste0("./tmp_results/est_unlab_yhat_sim_", j, "_boot.txt"), 
        sep = "\t", col.names = T)
}
