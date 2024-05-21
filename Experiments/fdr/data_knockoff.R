rm(list = ls())
library(knockoff)
library(ranger)
library(data.table)
n = 10^4        # number of observations
p = 150        # number of variables
k = 30            # number of variables with nonzero coefficients
amplitude = 5   # signal amplitude (for noise level = 1)

set.seed(42)
sim.times <- 1000

ratio <- 10
r <- 0.9
for (i in 1:sim.times){
    print(i)
    # Generate the variables from a multivariate normal distribution
    mu = rep(0,p)
    rho = 0.25
    Sigma = toeplitz(rho^(0:(p-1)))
    X = matrix(rnorm(n*p),n) %*% chol(Sigma)

    # Generate the response from a linear model
    nonzero = 1:k
    beta = amplitude * (1:p %in% nonzero) / sqrt(100)
    y.sample = function(X) X %*% beta + rnorm(n)
    y = y.sample(X)

    # Get Z as a function of 
    theta_0 <- rep(beta[1], p)
    Z <- r*y + X %*% theta_0  + rnorm(n)

    n_train <- 500
    n_lab <- 100
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


    # Create knockoff variables
    knockoffs=create.second_order
    knock_lab = knockoffs(as.matrix(X_lab))
    knock_unlab = knockoffs(as.matrix(X_unlab))

    Xk_lab <- knock_lab
    Xk_unlab <- knock_unlab

    # Randomly swap columns of X and Xk
    swap = rbinom(ncol(X_lab),1,0.5)
    swap.M.lab = matrix(swap,nrow=nrow(X_lab),ncol=length(swap),byrow=TRUE)
    X_lab.swap  = X_lab * (1-swap.M.lab) + Xk_lab * swap.M.lab
    Xk_lab.swap = X_lab * swap.M.lab + Xk_lab * (1-swap.M.lab)
    swap.M.unlab = matrix(swap,nrow=nrow(X_unlab),ncol=length(swap),byrow=TRUE)
    X_unlab.swap  = X_unlab * (1-swap.M.unlab) + Xk_unlab * swap.M.unlab
    Xk_unlab.swap = X_unlab * swap.M.unlab + Xk_unlab * (1-swap.M.unlab)

    X_lab_combined <- cbind(X_lab.swap, Xk_lab.swap)
    X_unlab_combined <- cbind(X_unlab.swap, Xk_unlab.swap)


    fwrite(dat_lab, paste0("./data/lab_", i, ".txt"), sep = "\t", col.names = T)
    fwrite(dat_unlab, paste0("./data/unlab_", i, ".txt"), sep = "\t", col.names = T)

    fwrite(X_lab_combined, paste0("./data/X_lab_combined_", i, ".txt"), sep = "\t", col.names = T)
    fwrite(X_unlab_combined, paste0("./data/X_unlab_combined_", i, ".txt"), sep = "\t", col.names = T)
    fwrite(data.frame(swap), paste0("./data/swap_", i, ".txt"), sep = "\t", col.names = T)
}
