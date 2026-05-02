library(mvtnorm)
library(dplyr)

bayesGaussian = function(X, Y, a, b, tau2){
  n = nrow(X)
  POSTERIOR = data.frame()
  for (i in 1:ncol(Y)){
    message(i, '/', ncol(Y))
    y = Y[, i]
    
    
    Ip = diag(1, nrow = ncol(X))
    XTX = t(X) %*% X
    XTY = t(X) %*% y
    YTY = t(y) %*% y
    
    # params of beta coefficient's posterior
    L0 = 1/tau2 * Ip
    V = solve(XTX + L0)
    Bhat = V %*% XTY
    
    # params of sigma2's posterior
    ahat = a + n/2
    ahat = as.numeric(ahat)
    
    bhat = b + 1/2 * (YTY - t(Bhat) %*% solve(V) %*% Bhat)
    bhat = as.numeric(bhat)
    # posterior of beta coefficients
    nu = 2 * ahat
    crit = qt(0.975, df = nu)
    beta_scale = sqrt(bhat / ahat * diag(V))
    ci_lower = Bhat - crit * beta_scale
    ci_upper = Bhat + crit * beta_scale
    
    posterior = data.frame(u = Bhat[2], df = nu, 
                           scale = beta_scale[2], ci_lower = ci_lower[2], ci_upper = ci_upper[2])
    POSTERIOR = rbind(POSTERIOR, posterior)
    
  }
  return(POSTERIOR)
}


# usage -------------------------------------------------------------------

# setwd('C:/Users/dcginger/Desktop/st502/project/data/')
# sim = readRDS('data_simulator___sim1.rds')
# Y = sim$counts
# params = sim$params
# dim(Y)
# case_idx = grepl('case', rownames(Y)) %>% as.numeric()
# X = cbind(1, case_idx)
# 
# 
# # log normalize
# lib_size = rowSums(Y)
# scale_factor = median(lib_size)
# Y = sweep(Y, MARGIN = 1, STATS = lib_size, FUN = '/')
# Y = Y * scale_factor
# Y = log1p(Y)
# 
# a = 2
# b = 0.5
# tau2 = 5
# results = bayesGaussian(X, Y, a, b, tau2)

# sanity check ------------------------------------------------------------
# RESULTS = data.frame()
# for ( i in 1:ncol(Y)){
#   message(i, '/', ncol(Y))
#   y = Y[, i]
#   h1 = glm(y ~ X) 
#   u = summary(h1)$coefficients[2, 1]
#   tval = summary(h1)$coefficients[2, 3]
#   results = data.frame(u = u, tval = tval)
#   RESULTS = rbind(RESULTS, results)
# }
# plot(RESULTS$u, POSTERIOR$u)
# plot(RESULTS$tval, POSTERIOR$u / POSTERIOR$scale)
