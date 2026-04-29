library(mvtnorm)

gradBetaLik = function(B, Y, X, varphi){
  phi = exp(varphi)
  u = as.vector(exp(X %*% B))
  as.vector(t(X) %*% (phi * (Y - u) / (u + phi)))
}

gradBetaPrior = function(B, tau2){
  -B / tau2
}

gradBeta = function(B, Y, X, varphi, tau2, ...){
  -1 * (gradBetaLik(B, Y, X, varphi) + gradBetaPrior(B, tau2))
}

hess = function(B, Y, X, varphi, tau2){
  phi = exp(varphi)
  u = as.vector(exp(X %*% B))
  p = length(B)
  
  w = phi * u * (Y + phi) / (u + phi)^2
  H = -t(X) %*% diag(w) %*% X - diag(1 / tau2, p)
  
  return(-H)
}

logLik = function(B, Y, X, varphi){
  sum(dnbinom(
    Y,
    mu = as.vector(exp(X %*% B)),
    size = exp(varphi),
    log = TRUE
  ))
}

logPrior = function(B, varphi, tau2, a, b){
  sum(dnorm(B, mean = 0, sd = sqrt(tau2), log = TRUE)) +
    dgamma(exp(varphi), shape = a, rate = b, log = TRUE) +
    varphi
}

objectiveFn = function(B, Y, X, varphi, tau2, a, b){
  -1 * (logLik(B, Y, X, varphi) + logPrior(B, varphi, tau2, a, b))
}

argmaxBeta = function(B, Y, X, varphi, tau2, a, b){
  fit = optim(
    par = B,
    fn = objectiveFn,
    gr = gradBeta,
    Y = Y,
    X = X,
    varphi = varphi,
    tau2 = tau2,
    a = a,
    b = b,
    hessian = FALSE,
    method = "BFGS"
  )
  
  return(fit)
}

dVarphi = function(B, Y, X, varphi, tau2, a, b, Bcov){
  exp(-1 * objectiveFn(B, Y, X, varphi, tau2, a, b)) /
    dmvnorm(B, mean = B, sigma = Bcov)
}
# fit INLA


INLA = function(Y, X, varphi_grid, B_start, tau2, a, b, n_post = 1000){
  
  stopifnot(length(Y) == nrow(X))
  stopifnot(length(B_start) == ncol(X))
  
  laplace_results = vector("list", length(varphi_grid))
  
  for(k in seq_along(varphi_grid)){
    varphi = varphi_grid[k]
    
    fit = argmaxBeta(
      B = B_start,
      Y = Y,
      X = X,
      varphi = varphi,
      tau2 = tau2,
      a = a,
      b = b
    )
    
    Bhat = fit$par
    Q = hess(Bhat, Y, X, varphi, tau2)
    Bcov = solve(Q)
    
    log_dvarphi =
      -objectiveFn(Bhat, Y, X, varphi, tau2, a, b) -
      dmvnorm(Bhat, mean = Bhat, sigma = Bcov, log = TRUE)
    
    laplace_results[[k]] = list(
      varphi = varphi,
      Bhat = Bhat,
      Q = Q,
      Bcov = Bcov,
      log_dvarphi = log_dvarphi,
      convergence = fit$convergence
    )
    
    B_start = Bhat
  }
  
  logw = sapply(laplace_results, function(z) z$log_dvarphi)
  
  logw_shift = logw - max(logw)
  
  delta = c(diff(varphi_grid), tail(diff(varphi_grid), 1))
  
  w_unnorm = exp(logw_shift) * delta
  w = w_unnorm / sum(w_unnorm)
  
  Bmat = do.call(rbind, lapply(laplace_results, function(z) z$Bhat))
  
  E_beta = colSums(Bmat * w)
  E_phi = sum(exp(varphi_grid) * w)
  
  Blist = lapply(laplace_results, function(z) z$Bhat)
  Bcovlist = lapply(laplace_results, function(z) z$Bcov)
  
  k_samp = sample(
    seq_along(varphi_grid),
    size = n_post,
    replace = TRUE,
    prob = w
  )
  
  Bpost = matrix(NA, nrow = n_post, ncol = length(Blist[[1]]))
  
  for(s in 1:n_post){
    k = k_samp[s]
    Bpost[s, ] = rmvnorm(
      n = 1,
      mean = Blist[[k]],
      sigma = Bcovlist[[k]]
    )
  }
  
  varphi_post = varphi_grid[k_samp]
  phi_post = exp(varphi_post)
  
  posterior = cbind(Bpost, phi_post)
  colnames(posterior) = c(paste0("B", seq_len(ncol(X)) - 1), "phi")
  
  return(list(
    posterior = posterior,
    E_beta = E_beta,
    E_phi = E_phi,
    weights = w,
    laplace_results = laplace_results
  ))
}


# usage -------------------------------------------------------------------

# Y = readRDS('C:/Users/dcginger/Desktop/st502/project/data/data_simulator___sim1.rds')
# Y$params
# Y = Y$counts[, 5]
# case_idx = as.numeric(grepl("case", names(Y)))
# X = cbind(1, case_idx)
# varphi_grid = seq(-4, 4, length.out = 100)
# B_start = rep(0, ncol(X))
# tau2 = 1
# a = 1
# b = 2
# 
# 
# 
# fit = INLA(
#   Y = Y,
#   X = X,
#   varphi_grid = seq(-4, 4, length.out = 100),
#   B_start = rep(0, ncol(X)),
#   tau2 = tau2,
#   a = a,
#   b = b,
#   n_post = 1000
# )
# 
# posterior = fit$posterior
# plot(density(posterior[, 'B1'])
# )