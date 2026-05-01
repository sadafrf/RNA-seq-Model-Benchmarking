library(tidyr)
library(dplyr)

# Frequentist Negative Binomial Model for ST502 Project
# ========================================================================================
# Input (simulated) data is assumed to be in the following form:
  # Response variable, Y, is a matrix with 1000 columns (genes) and 100 rows (samples)
    # there are n subjects and p genes
  # Predictor variables, X, is a matrix with two columns: intercept (all 1) and case/control indicator
    # Intercept gives baseline of 1 to all genes
    # Rows are samples

all_nb_genes <- function(X,Y) {
  stopifnot(all(rownames(Y) == rownames(X)))
  Y <- t(Y)
  stopifnot(all(colnames(Y) == rownames(X)))
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
    fit <- optim(par = init_par, fn = neg_log_likelihood, X = X, Y = y, method = "BFGS")
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
    h <- optimHess(fit$par, neg_log_likelihood, X = X, Y = y)

    if (fit$convergence != 0) {
      warning(paste("Gene", gene, "did not converge"))
      next
    }

    h <- optimHess(fit$par, neg_log_likelihood, X = X, Y = y)
    if (any(is.nan(h)) || det(h) <= 0) {
      warning(paste("Gene", gene, "has invalid Hessian"))
      next
    }

    vcov_mat <- solve(h)
    # get the variance-covariance matrix

    # Get standard errors and test Beta1
    # -------------------------------------------------------
    se_beta <- sqrt(diag(vcov_mat))[1:ncol(X)]

    beta1 <- beta_hat[2]
    se_beta1 <- se_beta[2]
    z_stat <- beta1 / se_beta1
    p_val <- 2 * pnorm(-abs(z_stat))

    # calculate CIs
    # -------------------------------------------------------
    z <- qnorm(0.975)
    ci_lower <- beta1 - z * se_beta1
    ci_upper <- beta1 + z * se_beta1

    # save results
    # -------------------------------------------------------
    results[[gene]] <- data.frame(
      Gene_id = rownames(Y)[gene],
      beta0 = beta_hat[1],
      beta1 = beta1,
      SE_beta1 = se_beta1,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      Z_Statistic = z_stat,
      P_Value = p_val
    )
  }
  results <- do.call(rbind, results)
  # False Discovery Rate Adjustment for P-values
  results$FDR_Adj_P<- p.adjust(results$P_Value, method = "BH")
  return(results)
}

# ============================================================================
# FREQUENTIST NEGATIVE BINOMIAL - Evaluation Criteria (same as Frequentist Gaussian)
# called de, power, FDP, AUC
# ============================================================================
evaluate_freq_nb <- function(results, true_de) {
  # DE call
  #-------------------------------------
  # CI does not include 0 and p-value is significant
  called_de <- ((results$CI_Upper < 0 | results$CI_Lower > 0) & results$FDR_Adj_P < 0.05)
  # Returns TRUE or FALSE

  # Power
  # power closed form solution?
  #-------------------------------------
  # in order to get the true DE, we need to compare the test results with the original simulated data that has a call for DE or not
  # de_p is a vector from the simulated data with the calls for DE genes
  true_de_vector <- true_de[results$Gene_id]
  #power <- sum(called_de & true_de_vector == 1) / sum(true_de_vector == 1)
  # True positives: called_de=TRUE  & true_de_vector=1
  # Divide true positives called by the true number of DEs


  # at my fixed threshold (FDR < 0.05 + CI):
  true_pos  <- sum(called_de & true_de_vector == 1)   # same numerator as power
  false_pos <- sum(called_de & true_de_vector == 0)   # already computed above
  false_neg <- sum(!called_de & true_de_vector == 1)  # missed DE genes
  true_neg  <- sum(!called_de & true_de_vector == 0)

  precision <- ifelse((true_pos + false_pos) == 0, 0, true_pos / (true_pos + false_pos))
  n_de <- sum(true_de_vector == 1)
  recall <- ifelse(n_de == 0, NA, true_pos / n_de)
  f1        <- ifelse((precision + recall) == 0, 0,
                      2 * precision * recall / (precision + recall))

  # FDP and Type I Error
  #False discovery proportion: false pos as a proportion of the total called genes
  total_called <- sum(called_de)
  fdp <- ifelse(total_called == 0, 0, false_pos / total_called)
  non_de_genes <- true_de_vector == 0
  n_null <- sum(non_de_genes)
  type_i_error <- ifelse(n_null == 0, NA, false_pos / n_null)

  list(
    Type_I_error = type_i_error,
    called_de    = called_de,
    #power        = power,          # = recall so I'm removing
    fdp          = fdp,
    #auc          = auc,
    true_pos     = true_pos,
    false_pos     = false_pos,
    false_neg = false_neg,
    true_neg = true_neg,
    precision    = precision,
    recall       = recall,         # same as power
    f1           = f1
  )
}
