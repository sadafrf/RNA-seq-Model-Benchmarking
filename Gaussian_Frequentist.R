library(tidyr)
library(dplyr)


# Frequentist Models for ST502 Project
# ========================================================================================
# Input (simulated) data is assumed to be in the following form:
  # Response variable, Y, is a matrix with 1000 columns (genes) and 100 rows (samples)
      # there are n subjects and p genes
  # Predictor variables, X, is a matrix with two columns: intercept (all 1) and case/control indicator
      # Intercept gives baseline of 1 to all genes
      # Rows are samples


#==========================================================================================
# FREQUENTIST GAUSSIAN
#==========================================================================================
all_gaussian_genes <- function(X, Y) {
  stopifnot(all(rownames(Y) == rownames(X)))
  Y <- t(Y)
  stopifnot(all(colnames(Y) == rownames(X)))
  # Updated to divide by library size then multiple by median counts, then log transform
  #--> this more accurately captures RNA seq data variability
  lib_sizes <- colSums(Y) #library sizes (per sample)
  median_lib <- median(lib_sizes) #median library size
  Y_norm <- sweep(Y, 2, lib_sizes, FUN = "/")   # divide each column
  Y_scaled <- Y_norm * median_lib              # rescale
  Y_log <- log2(Y_scaled + 1)      #log transform
  XtX_inv <- solve(t(X) %*% X)
  XtX_inv_times_X <- XtX_inv %*% t(X)
  p <- ncol(X)

  results <- vector("list", nrow(Y))
  for (gene in 1:nrow(Y)) {
    y_log <- as.numeric(Y_log[gene, ])

    # β̂ is a vector: it contains β₀ and β₁
    # β0 is intercept
    # β1 is difference between case and controls -> if β1=0 then no DE. If β1 not 0 then DE (after checking CI)
    # %*% = matrix multiplication
    beta_hat <- XtX_inv_times_X %*% y_log # MLE Closed Form Solution is: β_hat = (X^T X)^(−1) X^T Y

    # Beta hat coefficients
    beta0 = beta_hat[1]
    beta1 = beta_hat[2] #B1 is log fold change on log2 scale (due to normalization step)

    # residuals and sigma ^2
    df <- length(y_log) - p
    residuals <- y_log - X %*% beta_hat
    sigma2_hat <- sum(residuals^2) / df

    # Standard error of beta1
    # SE = sqrt( sigma^2 * (X^T X)^{-1}[2,2] )
    se_beta1 <- sqrt(sigma2_hat * XtX_inv[2,2])

    # Confidence intervals
    #ci_lower <- beta1 - qnorm(0.975) * se_beta1
    #ci_upper <- beta1 + qnorm(0.975) * se_beta1
    # need to use t distribution instead of normal because we estimated variance (sigma^2)
    t_crit <- qt(0.975, df)
    ci_lower <- beta1 - t_crit * se_beta1
    ci_upper <- beta1 + t_crit * se_beta1

    # FREQUENTIST GAUSSIAN- HYPOTHESIS TESTING: Is δ = 0?
    t_stat <- beta1/se_beta1
    p_val <- 2 * pt(-abs(t_stat), df = df)

    # results
    results[[gene]] <- data.frame(
      Gene_id = rownames(Y)[gene],
      beta0 = beta0,
      beta1 = beta1,
      sigma2 = sigma2_hat,
      SE_beta1 = se_beta1,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      T_Statistic = t_stat,
      P_Value = p_val
    )
  }
  results <- do.call(rbind, results)

  # False Discovery Rate Adjustment for P-values
  results$FDR_Adj_P<- p.adjust(results$P_Value, method = "BH")
  return(results)
}


# ============================================================================
# FREQUENTIST GAUSSIAN- Evaluation Criteria
# called de, power, FDP, AUC
# ============================================================================
evaluate_freq_gaussian <- function(results, true_de) {
  # DE call
  #-------------------------------------
  called_de <- results$FDR_Adj_P < 0.05
  called_de[is.na(called_de)] <- FALSE

  # True DEs
  true_de_vector <- true_de[as.character(results$Gene_id)]
  true_de_vector[is.na(true_de_vector)] <- 0


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

  # --- β1 metrics ---
  df <- results[!is.na(results$beta1) & !is.na(results$true_beta1), ]

  bias   <- mean(df$beta1 - df$true_beta1)
  mse    <- mean((df$beta1 - df$true_beta1)^2)
  #corval is the correlation between the estimated B1 and the true B1
  corval <- cor(df$beta1, df$true_beta1)
  #coverage is the % of time the true B1 falls in the 95% CI
  cover  <- mean(df$true_beta1 >= df$CI_Lower &
                   df$true_beta1 <= df$CI_Upper)

  list(
    Type_I_error = type_i_error,
    called_de    = called_de,
    fdp          = fdp,
    true_pos     = true_pos,
    false_pos     = false_pos,
    false_neg = false_neg,
    true_neg = true_neg,
    precision    = precision,
    recall       = recall,         # same as power
    f1           = f1,
    beta_bias    = bias,
    beta_mse     = mse,
    beta_cor     = corval,
    beta_coverage = cover
  )
}
