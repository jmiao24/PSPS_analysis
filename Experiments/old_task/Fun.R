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


sim_data <- function(r = 0.9, binary = FALSE) {
  # Input parameters
  n_train <- 500
  n_lab <- 500
  n_unlab <- 5 * 10^4
  sigma_Y <- sqrt(1)

  # Simulate the data
  mu <- c(0, 0) # Mean vector
  Sigma <- matrix(c(1, 0, 0, 1), 2, 2) # Covariance matrix
  n_data <- n_unlab + n_lab + n_train
  data <- as.data.frame(MASS::mvrnorm(n_data, mu, Sigma))
  colnames(data) <- c("X1", "X2")
  beta_1 <- sqrt(0.08) * sigma_Y
  focal <- data$X1 * beta_1
  beta_2 <- 1 * sigma_Y
  Y_tmp <- focal + r * scale(data$X2 * beta_2 + data$X1^2 * beta_2 + data$X2^2 * beta_2) * sigma_Y
  data$epsilon <- rnorm(n_data, 0, sqrt(1 - var(Y_tmp))) * sigma_Y
  data$Y <- Y_tmp + data$epsilon

  if (binary) {
    data$Y <- ifelse(data$Y > median(unlist(data$Y)), 1, 0)
    # Split the data
    train_data <- data[1:n_train, ]
    lab_data <- data[(n_train + 1):(n_lab + n_train), ]
    unlab_data <- data[(n_lab + n_train + 1):n_data, ]
    # Fit the machine learning model
    train_data$Y <- as.factor(train_data$Y)
    train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)
    lab_data$Y_hat <- predict(train_fit, newdata = lab_data)
    unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

    X_lab <- as.data.frame(lab_data$X1)
    X_unlab <- as.data.frame(unlab_data$X1)
    Y_lab <- as.data.frame(lab_data$Y)
    Yhat_lab <- as.data.frame(as.numeric(lab_data$Y_hat) - 1)
    Yhat_unlab <- as.data.frame(as.numeric(unlab_data$Y_hat) - 1)
    Y_unlab <- as.data.frame(unlab_data$Y)
    colnames(X_lab) <- "X1"
    colnames(X_unlab) <- "X1"
    colnames(Y_lab) <- "Y"
    colnames(Yhat_lab) <- "Y_hat"
    colnames(Yhat_unlab) <- "Y_hat"
    colnames(Y_unlab) <- "Y"
    out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab, Y_unlab = Y_unlab)
  } else {
    # Split the data
    train_data <- data[1:n_train, ]
    lab_data <- data[(n_train + 1):(n_lab + n_train), ]
    unlab_data <- data[(n_lab + n_train + 1):n_data, ]

    # Fit the machine learning model
    train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)
    lab_data$Y_hat <- predict(train_fit, newdata = lab_data)
    unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

    X_lab <- as.data.frame(lab_data$X1)
    X_unlab <- as.data.frame(unlab_data$X1)
    Y_lab <- as.data.frame(lab_data$Y)
    Yhat_lab <- as.data.frame(lab_data$Y_hat)
    Yhat_unlab <- as.data.frame(unlab_data$Y_hat)
    Y_unlab <- as.data.frame(unlab_data$Y)
    colnames(X_lab) <- "X1"
    colnames(X_unlab) <- "X1"
    colnames(Y_lab) <- "Y"
    colnames(Yhat_lab) <- "Y_hat"
    colnames(Yhat_unlab) <- "Y_hat"
    colnames(Y_unlab) <- "Y"
    out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab, Y_unlab = Y_unlab)
  }
  return(out)
}

wilcox_test <- function(y1, y2) {
  # Perform the Mann-Whitney-Wilcoxon test
  test_result <- wilcox.test(y1, y2)

  # Calculation of U statistic (Mann-Whitney statistic)
  rank_sum <- sum(rank(c(y1, y2))[1:length(y1)]) # Rank sum for sample1
  n1 <- length(y1)
  n2 <- length(y2)
  U <- rank_sum - n1 * (n1 + 1) / 2 # Mann-Whitney U statistic

  # Expected value and variance under H0
  mu_U <- n1 * n2 / 2
  sigma_U <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)

  # Standardizing U
  Z <- (U - mu_U) / sigma_U
  P.value <- 2 * pnorm(-abs(Z)) # Using symmetry of the normal distribution

  df <- data.frame(Estimate = (U - mu_U) / mu_U, Std.Error = sigma_U / mu_U, Test.Stat = Z, P.value = 2 * pnorm(-abs(Z)))
  return(df)
}

sim_data_wilcox_test <- function(r = 0.9, effect = T) {
  # Input parameters
  n_train <- 500
  n_lab <- 500
  n_unlab <- 5 * 10^4
  sigma_Y <- sqrt(1)

  # Simulate the data
  mu <- c(0, 0) # Mean vector
  Sigma <- matrix(c(1, 0, 0, 1), 2, 2) # Covariance matrix
  n_data <- n_unlab + n_lab + n_train
  X1 <- rbinom(n_data, 1, 0.5)
  X2 <- rnorm(n_data, 0, 1)
  data <- data.frame(X1, X2)
  colnames(data) <- c("X1", "X2")
  if (effect == TRUE) {
    beta_1 <- sqrt(0.01) * sigma_Y
  } else {
    beta_1 <- 0
  }
  focal <- data$X1 * beta_1
  beta_2 <- 1 * sigma_Y
  Y_tmp <- focal + r * scale(data$X2 * beta_2 + data$X2^2 * beta_2) * sigma_Y
  data$epsilon <- rnorm(n_data, 0, sqrt(1 - var(Y_tmp))) * sigma_Y
  data$Y <- Y_tmp + data$epsilon

  # Split the data
  train_data <- data[1:n_train, ]
  lab_data <- data[(n_train + 1):(n_lab + n_train), ]
  unlab_data <- data[(n_lab + n_train + 1):n_data, ]

  # Fit the machine learning model
  train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)
  lab_data$Y_hat <- predict(train_fit, newdata = lab_data)
  unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

  X_lab <- as.data.frame(lab_data$X1)
  X_unlab <- as.data.frame(unlab_data$X1)
  Y_lab <- as.data.frame(lab_data$Y)
  Yhat_lab <- as.data.frame(lab_data$Y_hat)
  Yhat_unlab <- as.data.frame(unlab_data$Y_hat)
  Y_unlab <- as.data.frame(unlab_data$Y)
  colnames(X_lab) <- "X1"
  colnames(X_unlab) <- "X1"
  colnames(Y_lab) <- "Y"
  colnames(Yhat_lab) <- "Y_hat"
  colnames(Yhat_unlab) <- "Y_hat"
  colnames(Y_unlab) <- "Y"
  out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab, Y_unlab = Y_unlab)
  return(out)
}

sim_data_quantile <- function(r = 0.9) {
  # Input parameters
  n_train <- 500
  n_lab <- 500
  n_unlab <- 5 * 10^4
  sigma_Y <- sqrt(1)

  # Simulate the data
  mu <- c(0, 0) # Mean vector
  Sigma <- matrix(c(1, 0, 0, 1), 2, 2) # Covariance matrix
  n_data <- n_unlab + n_lab + n_train
  data <- as.data.frame(MASS::mvrnorm(n_data, mu, Sigma))
  colnames(data) <- c("X1", "X2")
  beta_1 <- sqrt(0.08) * sigma_Y
  focal <- data$X1 * beta_1
  beta_2 <- 1 * sigma_Y
  Y_tmp <- focal + r * scale(data$X1^2 * beta_2 + data$X2 * beta_2 + data$X2^2 * beta_2) * sigma_Y
  data$epsilon <- rnorm(n_data, 0, sqrt(1 - var(Y_tmp))) * sigma_Y
  data$Y <- Y_tmp + data$epsilon

  # Split the data
  train_data <- data[1:n_train, ]
  lab_data <- data[(n_train + 1):(n_lab + n_train), ]
  unlab_data <- data[(n_lab + n_train + 1):n_data, ]

  # Fit the machine learning model
  train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)
  lab_data$Y_hat <- predict(train_fit, newdata = lab_data)
  unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

  X_lab <- as.data.frame(lab_data$X1)
  X_unlab <- as.data.frame(unlab_data$X1)
  Y_lab <- as.data.frame(lab_data$Y)
  Yhat_lab <- as.data.frame(lab_data$Y_hat)
  Yhat_unlab <- as.data.frame(unlab_data$Y_hat)
  Y_unlab <- as.data.frame(unlab_data$Y)
  colnames(X_lab) <- "X1"
  colnames(X_unlab) <- "X1"
  colnames(Y_lab) <- "Y"
  colnames(Yhat_lab) <- "Y_hat"
  colnames(Yhat_unlab) <- "Y_hat"
  colnames(Y_unlab) <- "Y"
  out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab, Y_unlab = Y_unlab)
  return(out)
}

sim_data_nbinom <- function(r = 0.9) {
  # Input parameters
  n_train <- 500
  n_lab <- 500
  n_unlab <- 5 * 10^4
  sigma_Y <- sqrt(1)

  # Simulate the data
  mu <- c(0, 0) # Mean vector
  Sigma <- matrix(c(1, 0, 0, 1), 2, 2) # Covariance matrix
  n_data <- n_unlab + n_lab + n_train
  data <- as.data.frame(MASS::mvrnorm(n_data, mu, Sigma))
  colnames(data) <- c("X1", "X2")
  beta_1 <- sqrt(0.3) * sigma_Y
  focal <- data$X1 * beta_1
  beta_2 <- 1 * sigma_Y
  eta <- focal + r * scale(data$X2 * beta_2) * sigma_Y

  # Transform linear predictor to rate (mu) using exponential to ensure positivity
  mu <- exp(eta)
  k <- 1 # example value, adjust based on the desired variance

  data$Y <- rnbinom(n_data, s = k, mu = mu)

  # Split the data
  train_data <- data[1:n_train, ]
  lab_data <- data[(n_train + 1):(n_lab + n_train), ]
  unlab_data <- data[(n_lab + n_train + 1):n_data, ]

  # Fit the machine learning model
  train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)
  lab_data$Y_hat <- round(predict(train_fit, newdata = lab_data))
  unlab_data$Y_hat <- round(predict(train_fit, newdata = unlab_data))

  X_lab <- as.data.frame(lab_data$X1)
  X_unlab <- as.data.frame(unlab_data$X1)
  Y_lab <- as.data.frame(lab_data$Y)
  Yhat_lab <- as.data.frame(lab_data$Y_hat)
  Yhat_unlab <- as.data.frame(unlab_data$Y_hat)
  Y_unlab <- as.data.frame(unlab_data$Y)
  colnames(X_lab) <- "X1"
  colnames(X_unlab) <- "X1"
  colnames(Y_lab) <- "Y"
  colnames(Yhat_lab) <- "Y_hat"
  colnames(Yhat_unlab) <- "Y_hat"
  colnames(Y_unlab) <- "Y"
  out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab, Y_unlab = Y_unlab)

  return(out)
}

sim_data_iv <- function(r = 0.9) {
  # Input parameters
  n_train <- 500
  n_lab <- 500
  n_unlab <- 5 * 10^4
  sigma_Y <- sqrt(1)

  # Simulate the data
  mu <- c(0, 0) # Mean vector
  Sigma <- matrix(c(1, 0, 0, 1), 2, 2) # Covariance matrix
  n_data <- n_unlab + n_lab + n_train
  Z <- rnorm(n_data)
  eff <- 0.4
  X1 <- eff * Z + rnorm(n_data, 0, sqrt(1 - eff^2))
  data <- data.frame(X1, Z)
  colnames(data) <- c("X1", "Z")
  beta_1 <- sqrt(0.08) * sigma_Y
  focal <- data$X1 * beta_1
  beta_2 <- 1 * sigma_Y
  Y_tmp <- focal
  data$epsilon <- rnorm(n_data, 0, sqrt(1 - var(Y_tmp))) * sigma_Y
  data$Y <- Y_tmp + data$epsilon
  X2_tmp <- 0.3 * scale(Z) + r * scale(data$Y)
  data$X2 <- X2_tmp + rnorm(n_data, 0, sqrt(1 - var(X2_tmp)))

  # Split the data
  train_data <- data[1:n_train, ]
  lab_data <- data[(n_train + 1):(n_lab + n_train), ]
  unlab_data <- data[(n_lab + n_train + 1):n_data, ]

  # Fit the machine learning model
  train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)
  lab_data$Y_hat <- predict(train_fit, newdata = lab_data)
  unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

  X_lab <- as.data.frame(lab_data$X1)
  X_unlab <- as.data.frame(unlab_data$X1)
  Z_lab <- as.data.frame(lab_data$Z)
  Z_unlab <- as.data.frame(unlab_data$Z)
  Y_lab <- as.data.frame(lab_data$Y)
  Yhat_lab <- as.data.frame(lab_data$Y_hat)
  Yhat_unlab <- as.data.frame(unlab_data$Y_hat)
  Y_unlab <- as.data.frame(unlab_data$Y)
  colnames(X_lab) <- "X1"
  colnames(X_unlab) <- "X1"
  colnames(Z_lab) <- "Z"
  colnames(Z_unlab) <- "Z"
  colnames(Y_lab) <- "Y"
  colnames(Yhat_lab) <- "Y_hat"
  colnames(Yhat_unlab) <- "Y_hat"
  colnames(Y_unlab) <- "Y"
  out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab, Y_unlab = Y_unlab, Z_lab = Z_lab, Z_unlab = Z_unlab)
  return(out)
}

sim_data_dlasso <- function(r = 0.9, binary = FALSE) {
  # Input parameters
  n_train <- 500
  n_lab <- 100
  n_unlab <- 7500
  n_data <- n_train + n_lab + n_unlab
  sigma_Y <- sqrt(1)

  p <- 200
  s0 <- 15 # Number of non-zero coefficients
  b <- 1
  b0 <- 0
  sigma <- 0.5

  X <- rnorm(p * n_data)
  dim(X) <- c(n_data, p)
  theta0 <- c(rep(b, s0), rep(0, p - s0))
  Y_tmp <- r * scale(X %*% theta0)
  Y <- Y_tmp + rnorm(n_data, 0, sqrt(1 - var(Y_tmp)))

  data <- data.frame(Y, X)

  train_data <- data[1:n_train, ]
  lab_data <- data[(n_train + 1):(n_lab + n_train), ]
  unlab_data <- data[(n_lab + n_train + 1):n_data, ]

  # Fit the machine learning model
  train_fit <- ranger::ranger(Y ~ ., data = train_data[, 1:p], num.trees = 100)
  lab_data$Y_hat <- predict(train_fit, data = lab_data)$predictions
  unlab_data$Y_hat <- predict(train_fit, data = unlab_data)$predictions

  cor(lab_data$Y, lab_data$Y_hat)

  X_lab <- as.data.frame(lab_data[, paste0("X", 1:p)])
  X_unlab <- as.data.frame(unlab_data[, paste0("X", 1:p)])
  Y_lab <- as.data.frame(lab_data$Y)
  Yhat_lab <- as.data.frame(lab_data$Y_hat)
  Yhat_unlab <- as.data.frame(unlab_data$Y_hat)
  Y_unlab <- as.data.frame(unlab_data$Y)
  colnames(Y_lab) <- "Y"
  colnames(Yhat_lab) <- "Y_hat"
  colnames(Yhat_unlab) <- "Y_hat"
  colnames(Y_unlab) <- "Y"
  out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab, Y_unlab = Y_unlab)
  return(out)
}
