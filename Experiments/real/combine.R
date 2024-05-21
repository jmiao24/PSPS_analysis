require(data.table)
for (chr in 1:22){
  est_lab_y <- fread(paste0("./results/gwas/Total.lab.y_chr.", chr, ".int_rank_score.glm.linear"))
  est_unlab_y_hat <- fread(paste0("./results/gwas/Total.unlab.y_hat_chr.", chr, ".int_rank_score.glm.linear"))
  est_lab_y_hat <- fread(paste0("./results/gwas/Total.lab.y_hat_chr.", chr, ".int_rank_score.glm.linear"))

  # P-value to
  est_lab_y$SE <- 1/sqrt(OBS_CT)
  est_lab_y_hat$SE <- 1/sqrt(OBS_CT)
  est_unlab_y_hat$SE <- 1/sqrt(OBS_CT)

  est_lab_y$BETA <- est_lab_y$T_STAT * est_lab_y$SE
  est_lab_y_hat$BETA <- est_lab_y_hat$T_STAT * est_lab_y_hat$SE
  est_unlab_y_hat$BETA <- est_unlab_y_hat$T_STAT * est_unlab_y_hat$SE

  Cov <- est_lab_y$SE * est_lab_y$SE * 0.2488

  # Apply PSPS
  psps_est <- c()
  psps_se <- c()
  psps_p <- c()
  for (i in 1:nrow(est_lab_y)){
    Sigma <- matrix(c(est_lab_y$SE[i]^2, Cov, 0, Cov, est_lab_y$SE[i]^2, 0, 0, 0, est_unlab_y_hat$SE[i]^2), nrow=2)
    psps <- PSPS(est_lab_y$BETA[i], est_lab_y_hat$BETA[i], est_unlab_y_hat$BETA[i], Sigma)
    est_psps <- psps$Esimate
    se_psps <- psps$Std.Error
    p_psps <- 2*pnorm(-abs(est_psps/se_psps))
    psps_est <- c(psps_est, est_psps)
    psps_se <- c(psps_se, se_psps)
    psps_p <- c(psps_p, p_psps)
  }
  fwrite(data.table(est_lab_y$SNP, est_lab_y$CHR, est_lab_y$BP, est_lab_y$A1, est_lab_y$A2, est_lab_y$BETA, est_lab_y$SE, est_lab_y$P, psps_est, psps_se, psps_p), paste0("./results/gwas/Total.psps.txt"),
  col.names = T, row.names = F, sep = "\t")
}