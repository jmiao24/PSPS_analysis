# est_lab_y is a Kx1 vector of the true labels for the labeled data
# est_lab_yhat is a Kx1 vector of the predicted labels for the labeled data
# est_unlab_yhat is a Kx1 vector of the predicted labels for the unlabeled data
# Sigma is a 3K x 3K matrix of the covariance of the predicted labels for the labeled and unlabeled data
PSPS <- function(est_lab_y, est_lab_yhat, est_unlab_yhat, Sigma, alpha = 0.05) {
  # Prepare inputs
  q <- length(est_lab_y)
  Sigma <- as.matrix(Sigma)
  v_est <- Sigma[1:q, 1:q]
  r <- Sigma[(q + 1):(2 * q), 1:q]
  v_eta_lab <- Sigma[(q + 1):(2 * q), (q + 1):(2 * q)]
  v_eta_unlab <- Sigma[(2 * q + 1):(3 * q), (2 * q + 1):(3 * q)]
  V <- v_eta_lab + v_eta_unlab
  omega_0 <- solve(V) %*% r

  # Get the results
  est <- est_lab_y + omega_0 %*% (est_unlab_yhat - est_lab_yhat)
  standard_errors <- sqrt(diag(v_est - t(omega_0) %*% r))
  lower_ci <- est - qnorm(1 - alpha / 2) * standard_errors
  upper_ci <- est + qnorm(1 - alpha / 2) * standard_errors
  p_values <- 2 * pnorm(abs(est / standard_errors), 0, 1, lower.tail = F)
  output_table <- data.frame(
    Estimate = est, Std.Error = standard_errors,
    Lower.CI = lower_ci, Upper.CI = upper_ci, P.value = p_values
  )
  colnames(output_table) <- c("Estimate", "Std.Error", "Lower.CI", "Upper.CI", "P.value")
  return(output_table)
}