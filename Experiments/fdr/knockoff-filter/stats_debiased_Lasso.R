source("./stats_glmnet_cv.R") # Under the same repo

stat.sslasso_coefdiff <- function(X, X_k, y, family='gaussian', cores=2, ...) {
  source("./sslasso_code/lasso_inference.r")
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  parallel=T
  if (!requireNamespace('doParallel', quietly=T)) {
    warning('doParallel is not installed. Without parallelization, the statistics will be slower to compute', call.=F,immediate.=T)
    parallel=F
  }
  if (!requireNamespace('parallel', quietly=T)) {
    warning('parallel is not installed. Without parallelization, the statistics will be slower to compute.', call.=F,immediate.=T)
    parallel=F
  }
    
  # Register cores for parallel computation
  if (parallel) {
    ncores = parallel::detectCores(all.tests = TRUE, logical = TRUE)
    if( cores==2 ) {
      cores = min(2,ncores)
    }
    else {
      if (cores > ncores ) {
        warning(paste("The requested number of cores is not available. Using instead",ncores,"cores"),immediate.=T)
        cores = ncores
      }
    }
    if (cores>1) {
      doParallel::registerDoParallel(cores=cores)
      parallel = TRUE
    }
    else {
      parallel = FALSE
    }
  }
  
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)

  p = ncol(X)

  # Compute statistics
  glmnet.lambda = cv_lambda_glmnet(cbind(X.swap, Xk.swap), y, family=family, parallel=parallel)
  glmnet.coefs = SSLasso(cbind(X.swap, Xk.swap), y, lambda = glmnet.lambda)$unb.coef
  # cv_coeffs_glmnet(cbind(X.swap, Xk.swap), y, family=family, parallel=parallel)

  # Lasso(X, y, lambda = glmnet.lambda)$coef
  Z <- glmnet.coefs[1:(2*p)]
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)

  # Stop the parallel cluster (if applicable)  
  if (parallel) {
    if (cores>1) {
      doParallel::stopImplicitCluster()
    }
  }
  return(W)
}
