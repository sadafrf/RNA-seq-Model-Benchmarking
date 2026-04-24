library(tidyr)
library(dplyr)

# Frequentist Negative Binomial Model for ST502 Project
# ========================================================================================
# Input (simulated) data is assumed to be in the following form:
  # Response variable, Y, is a matrix with 100 columns (samples) and 1000 rows (genes)
    # there are n subjects and p genes
  # Predictor variables, X, is a matrix with two columns: intercept (all 1) and case/control indicator
    # Intercept gives baseline of 1 to all genes
    # Rows are samples

all_nb_genes <- function(X,Y) {

  p <- ncol(X)

  # Step 1: compute the negative log likelihood
  # -------------------------------------------------------
  neg_log_likelihood <- function(par, X, Y){
    beta <- par[1:p] #regression coefficient B
    phi <- par[p+1] #dispersion parameter phi
    mu <- exp(X %*% beta)   # = exp(x_i^T beta) --> linear predictor (XB). This is the mean of a negative binomial GLM with log link
    theta <- exp(phi)       # = e^phi --> ensures dispersion parameter is positive
    #now compute negative binomial log likelihood
    ll <- sum(lgamma(Y+theta) -lgamma(theta) -lgamma(Y+1) + theta*log(theta) + Y*log(mu) - (Y+theta)*log(mu+theta))
    return(-ll) #converts max log-likelihood to min log-likelihood
  }

  init_par <- c(rep(0, p), 0) #starting point for the optimizer (beta,phi)=(0,0)
  results <- vector("list", nrow(Y))

  # Optimize Beta_hat because it can't be solved in closed form
  # -------------------------------------------------------
  for (gene in 1:nrow(Y)) {
    y <- as.numeric(Y[gene, ])
    fit <- optim(par = init_par, fn = neg_log_likelihood, X = X, y = y, method = "BFGS")
    #optimization algorithm. It approximates the Hessian (second derivatives) -->
    # uses gradient information (or approximates it) --> updates parameters iteratively

    # Extract estimates
    # -------------------------------------------------------
    beta_hat <- fit$par[1:ncol(X)]
    phi_hat  <- fit$par[ncol(X) + 1]
    # compute fitted values
    mu_hat <- exp(X %*% beta_hat)

    # Hessian (replaces XtX in gaussian)
    # -------------------------------------------------------
    h <- optimHess(fit$par, neg_loglik, X = X, y = y)
    # get the variance-covariance matrix
    vcov_mat <- solve(h)

    # Get standard errors and test Beta1
    # -------------------------------------------------------
    se_beta <- sqrt(diag(vcov_mat))[1:ncol(X)]

    beta1 <- beta_hat[2]
    se_beta1 <- se_beta[2]
    t_stat <- beta1 / se_beta1
    p_val <- 2 * pnorm(-abs(t_stat))

    # calculate CIs
    # -------------------------------------------------------
    ci_lower <- beta1 - 1.96 * se_beta1
    ci_upper <- beta1 + 1.96 * se_beta1

    # save results
    # -------------------------------------------------------
    results[[gene]] <- data.frame(
      Gene_id = rownames(Y)[gene],
      beta0 = beta_hat[1],
      beta1 = beta1,
      SE_beta1 = se1,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      T_Statistic = t_stat,
      P_Value = p_val
    )
  }
  do.call(rbind, results)
}
