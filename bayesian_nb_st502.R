library(ggplot2)
library(dplyr)
library(tidyr)
library(rjags)
library(coda)


counts <- matrix(nrow=1000,ncol=100)
counts <- as.data.frame(counts)
for (i in (1:nrow(counts))){
  rownames(counts)[i] <- paste0("Gene", i) 
}

#Initialize alpha dispersion parameter with gamma prior. This prior is uninformative.
#The names for the dispersion list are set to be the rownames for the count matrix
alpha <- rgamma(n=1000, shape=2,rate=0.5)
names(alpha) <- rownames(counts)

#Initalize bernoulli (1 or 0) if a gene is DE expressed with probability 0.1. 
#The names of the probability list are set to be the rownames for the count matrix
de_p <- rbinom(n=1000, size=1, prob=0.1)
names(de_p) <- rownames(counts)

#Initialize delta as a normal distribution with mean 0 and variance 2
delta <- rnorm(n=length(rownames(counts[de_p == 1, ])), mean = 1, sd = sqrt(2))
names(delta) <- rownames(counts[de_p == 1, ])

#Initialize gene average (mu) as a poisson prior with sampling rate = 10 
mu <- rgamma(n=1000, shape=10, rate=1)
names(mu) <- rownames(counts)


colnames(counts)[1:50] <- "Control"
colnames(counts)[51:100] <- "Case" 

# Generate count data with error handling
for (g in 1:nrow(counts)) {
  gene_name <- rownames(counts)[g]
  # Generate control samples
  counts[g, 1:50] <- rnbinom(n = 50, mu = mu[g], size = alpha[g])
  # Generate case samples
  if (de_p[g] == 1) {
    # Differentially expressed gene
    mu_case <- mu[g] * exp(delta[gene_name])
    
    # Check if mu_case is valid (not too large, not NA, not Inf)
    if (is.finite(mu_case) && mu_case < 1e6) {
      counts[g, 51:100] <- rnbinom(n = 50, mu = mu_case, size = alpha[g])
    } else {
      # If overflow, cap the fold change
      warning(sprintf("Gene %d: mu_case too large, capping delta", g))
      delta[gene_name] <- log(1000 / mu[g])  # Cap at 1000x fold change
      mu_case <- mu[g] * exp(delta[gene_name])
      counts[g, 51:100] <- rnbinom(n = 50, mu = mu_case, size = alpha[g])
    }
  } else {
    # Not differentially expressed
    counts[g, 51:100] <- rnbinom(n = 50, mu = mu[g], size = alpha[g])
  }
}

print("Data is successfully simulated")
# Summary of simulated data
cat("\n=== SIMULATION SUMMARY ===\n")
cat(sprintf("Total genes: %d\n", nrow(counts)))
cat(sprintf("DE genes: %d (%.1f%%)\n", sum(de_p), 100 * sum(de_p) / length(de_p)))
cat(sprintf("Non-DE genes: %d (%.1f%%)\n", sum(de_p == 0), 100 * sum(de_p == 0) / length(de_p)))
cat(sprintf("Mean count (control): %.2f\n", mean(as.matrix(counts[, 1:50]))))
cat(sprintf("Mean count (case): %.2f\n", mean(as.matrix(counts[, 51:100]))))
cat(sprintf("Delta range: [%.2f, %.2f]\n", min(delta), max(delta)))

##Bayesian Inference##
# Log-likelihood function for Negative Binomial
log_likelihood_nb <- function(y, mu, alpha) {
  if (mu <= 0 || alpha <= 0) return(-Inf)
  sum(dnbinom(y, mu = mu, size = alpha, log = TRUE))
}

# Log-prior functions
log_prior_delta <- function(delta, mean = 0, sd = sqrt(2)) {
  dnorm(delta, mean = mean, sd = sd, log = TRUE)
}

log_prior_alpha <- function(alpha, shape = 2, rate = 0.5) {
  if (alpha <= 0) return(-Inf)
  dgamma(alpha, shape = shape, rate = rate, log = TRUE)
}

#Use Gamma prior for mu (continuous parameter)
log_prior_mu <- function(mu, shape = 10, rate = 1) {
  if (mu <= 0) return(-Inf)
  dgamma(mu, shape = shape, rate = rate, log = TRUE)
}

# Log-posterior function for a single gene
log_posterior <- function(params, y_control, y_case) {
  delta <- params[1]
  alpha <- params[2]
  mu <- params[3]
  
  # Constraints
  if (alpha <= 0 || mu <= 0) return(-Inf)
  
  # Likelihood for control group
  ll_control <- log_likelihood_nb(y_control, mu, alpha)
  if (!is.finite(ll_control)) return(-Inf)
  
  # Likelihood for case group
  mu_case <- mu * exp(delta)
  if (mu_case <= 0 || !is.finite(mu_case)) return(-Inf)
  ll_case <- log_likelihood_nb(y_case, mu_case, alpha)
  if (!is.finite(ll_case)) return(-Inf)
  
  # Priors
  lp_delta <- log_prior_delta(delta)
  lp_alpha <- log_prior_alpha(alpha)
  lp_mu <- log_prior_mu(mu)  # NOW USING GAMMA PRIOR
  
  # Log-posterior
  log_post <- ll_control + ll_case + lp_delta + lp_alpha + lp_mu
  
  if (!is.finite(log_post)) return(-Inf)
  
  return(log_post)
}



#Metropolis-Hastings sampler
metropolis_hastings <- function(y_control, y_case, 
                                n_iter = 10000, 
                                burn_in = 2000,
                                proposal_sd = c(0.1, 0.5, 1)) {
  delta <- rnorm(1, 0, 0.1)
  alpha <- rgamma(1, shape = 2, rate = 0.5)
  mu <- max(mean(y_control), 1)
  current_params <- c(delta, alpha, mu)
  current_log_post <- log_posterior(current_params, y_control, y_case)
  attempts <- 0
  while (!is.finite(current_log_post) && attempts < 10) {
    delta <- rnorm(1, 0, 0.1)
    alpha <- rgamma(1, shape = 5, rate = 1)
    mu <- max(mean(y_control), mean(y_case), 1)
    current_params <- c(delta, alpha, mu)
    current_log_post <- log_posterior(current_params, y_control, y_case)
    attempts <- attempts + 1
  }
  if (!is.finite(current_log_post)) {
    warning("Could not find valid starting values")
    return(NULL)
  }
  samples <- matrix(NA, nrow = n_iter, ncol = 3)
  colnames(samples) <- c("delta", "alpha", "mu")
  accepted <- 0
  for (i in 1:n_iter) {
    # Propose new parameters
    proposal <- current_params + rnorm(3, 0, proposal_sd)
    if (proposal[2] <= 0) proposal[2] <- abs(proposal[2])
    if (proposal[3] <= 0) proposal[3] <- abs(proposal[3])
    proposal_log_post <- log_posterior(proposal, y_control, y_case)
    # Handle -Inf case
    if (!is.finite(proposal_log_post)) {
      # Reject proposal
      samples[i, ] <- current_params
      next
    }
    # Acceptance ratio (log scale)
    log_alpha <- proposal_log_post - current_log_post
    # Accept or reject
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      current_params <- proposal
      current_log_post <- proposal_log_post
      accepted <- accepted + 1
    }
    # Store samples
    samples[i, ] <- current_params
  }
  # Remove burn-in
  samples_post_burnin <- samples[(burn_in + 1):n_iter, ]
  # Acceptance rate
  accept_rate <- accepted / n_iter
  return(list(
    samples = samples_post_burnin,
    all_samples = samples,
    acceptance_rate = accept_rate
  ))
}


analyze_gene <- function(gene_idx, counts_data, n_iter = 10000, burn_in = 2000) {
  y_control <- as.numeric(counts_data[gene_idx, 1:50])
  y_case <- as.numeric(counts_data[gene_idx, 51:100])
  
  # Run MH
  mh_results <- metropolis_hastings(y_control, y_case, 
                                    n_iter = n_iter, 
                                    burn_in = burn_in,
                                    proposal_sd = c(0.15, 0.3, 0.5))
  
  # Check if MH failed or returned NULL
  if (is.null(mh_results)) {
    warning(sprintf("MH failed for gene %d", gene_idx))
    return(NULL)
  }
  
  # Check if samples exist and have sufficient rows
  if (is.null(mh_results$samples) || nrow(mh_results$samples) == 0) {
    warning(sprintf("No samples generated for gene %d", gene_idx))
    return(NULL)
  }
  
  # Posterior summaries
  posterior_summary <- data.frame(
    gene = rownames(counts_data)[gene_idx],
    delta_mean = mean(mh_results$samples[, "delta"]),
    delta_median = median(mh_results$samples[, "delta"]),
    delta_sd = sd(mh_results$samples[, "delta"]),
    delta_ci_lower = quantile(mh_results$samples[, "delta"], 0.025),
    delta_ci_upper = quantile(mh_results$samples[, "delta"], 0.975),
    alpha_mean = mean(mh_results$samples[, "alpha"]),
    mu_mean = mean(mh_results$samples[, "mu"]),
    acceptance_rate = mh_results$acceptance_rate,
    prob_delta_positive = mean(mh_results$samples[, "delta"] > 0),
    prob_delta_zero = mean(abs(mh_results$samples[, "delta"]) < 0.1)  # ROPE
  )
  
  return(list(
    summary = posterior_summary,
    samples = mh_results$samples
  ))
}

# ============================================================================
# RUN ANALYSIS WITH IMPROVED ERROR HANDLING
# ============================================================================

cat("Running Metropolis-Hastings...\n")

n_genes_to_analyze <- 1000 # Change to 1000 for all genes
results_list <- list()
posterior_summaries <- data.frame()
failed_genes <- c()

for (g in 1:n_genes_to_analyze) {
  if (g %% 10 == 0) cat(sprintf("Processing gene %d/%d\n", g, n_genes_to_analyze))
  
  result <- tryCatch({
    analyze_gene(g, counts, n_iter = 5000, burn_in = 1000)
  }, error = function(e) {
    cat(sprintf("Error processing gene %d: %s\n", g, e$message))
    return(NULL)
  })
  
  if (!is.null(result)) {
    results_list[[g]] <- result
    posterior_summaries <- rbind(posterior_summaries, result$summary)
  } else {
    failed_genes <- c(failed_genes, g)
  }
}

cat(sprintf("\nMetropolis-Hastings complete!\n"))
cat(sprintf("Successfully analyzed: %d/%d genes\n", 
            nrow(posterior_summaries), n_genes_to_analyze))

if (length(failed_genes) > 0) {
  cat(sprintf("Failed genes: %s\n", paste(failed_genes, collapse = ", ")))
}

# ============================================================================
# HYPOTHESIS TESTING: Is δ = 0?
# ============================================================================

if (nrow(posterior_summaries) > 0) {
  posterior_summaries <- posterior_summaries %>%
    mutate(
      # Method 1: Credible interval excludes 0
      reject_H0_CI = (delta_ci_lower > 0 | delta_ci_upper < 0),
      
      # Method 2: Posterior probability that delta is practically zero
      reject_H0_ROPE = prob_delta_zero < 0.05,
      
      # Method 3: Posterior probability delta > 0 is very high or very low
      reject_H0_prob = (prob_delta_positive > 0.95 | prob_delta_positive < 0.05)
    )
  
  print(posterior_summaries)
  
  # ============================================================================
  # VISUALIZATION FOR FIRST SUCCESSFULLY ANALYZED GENE
  # ============================================================================
  
  # Find first non-NULL result
  first_valid_idx <- which(!sapply(results_list, is.null))[9]
  
  if (!is.na(first_valid_idx) && !is.null(results_list[[first_valid_idx]])) {
    gene_samples <- results_list[[first_valid_idx]]$samples
    gene_df <- as.data.frame(gene_samples)
    gene_name <- posterior_summaries$gene[9]
    
    # Trace plot for delta
    p1 <- ggplot(gene_df, aes(x = 1:nrow(gene_df), y = delta)) +
      geom_line(alpha = 0.7) +
      labs(title = sprintf("Trace Plot for δ (%s)", gene_name), 
           x = "Iteration", y = "δ") +
      theme_minimal()
    
    # Posterior distribution of delta
    p2 <- ggplot(gene_df, aes(x = delta)) +
      geom_histogram(aes(y = ..density..), bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_density(color = "red", size = 1) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
      geom_vline(xintercept = mean(gene_df$delta), linetype = "dashed", 
                 color = "red", size = 0.8) +
      labs(title = sprintf("Posterior Distribution of δ (%s)", gene_name), 
           x = "δ", y = "Density") +
      theme_minimal()
    
    # Trace plots for all parameters
    gene_df_long <- gene_df %>%
      mutate(iteration = 1:n()) %>%
      pivot_longer(cols = c(delta, alpha, mu), 
                   names_to = "parameter", 
                   values_to = "value")
    
    p3 <- ggplot(gene_df_long, aes(x = iteration, y = value)) +
      geom_line(alpha = 0.7) +
      facet_wrap(~parameter, scales = "free_y", ncol = 1) +
      labs(title = sprintf("Trace Plots for All Parameters (%s)", gene_name),
           x = "Iteration", y = "Value") +
      theme_minimal()
    
    print(p1)
    print(p2)
    print(p3)
  }
  
  # Summary of hypothesis tests
  cat(sprintf("Genes with δ ≠ 0 (CI method): %d/%d\n", 
              sum(posterior_summaries$reject_H0_CI), 
              nrow(posterior_summaries)))
  cat(sprintf("Genes with δ ≠ 0 (ROPE method): %d/%d\n", 
              sum(posterior_summaries$reject_H0_ROPE), 
              nrow(posterior_summaries)))
  
  # Compare with true DE status
  if (n_genes_to_analyze <= length(de_p)) {
    gene_indices <- as.numeric(gsub("Gene", "", posterior_summaries$gene))
    true_de <- de_p[gene_indices]
    posterior_summaries$true_DE <- true_de
    posterior_summaries$true_delta <- ifelse(true_de == 1, 
                                             delta[names(delta) %in% posterior_summaries$gene],
                                             0)
    
    cat("\n=== COMPARISON WITH TRUE DE STATUS ===\n")
    confusion_matrix <- table(
      Predicted = posterior_summaries$reject_H0_CI,
      True = posterior_summaries$true_DE
    )
    print(confusion_matrix)
    
    # Calculate accuracy metrics
    if (sum(confusion_matrix) > 0) {
      accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
      cat(sprintf("\nAccuracy: %.2f%%\n", accuracy * 100))
      
      # Sensitivity and Specificity
      if (sum(confusion_matrix[, "1"]) > 0) {
        sensitivity <- confusion_matrix["TRUE", "1"] / sum(confusion_matrix[, "1"])
        cat(sprintf("Sensitivity (True Positive Rate): %.2f%%\n", sensitivity * 100))
      }
      if (sum(confusion_matrix[, "0"]) > 0) {
        specificity <- confusion_matrix["FALSE", "0"] / sum(confusion_matrix[, "0"])
        cat(sprintf("Specificity (True Negative Rate): %.2f%%\n", specificity * 100))
      }
    }
  }
  
  # Display acceptance rates
  cat("\n=== ACCEPTANCE RATE SUMMARY ===\n")
  cat(sprintf("Mean acceptance rate: %.2f%%\n", 
              mean(posterior_summaries$acceptance_rate) * 100))
  cat(sprintf("Min acceptance rate: %.2f%%\n", 
              min(posterior_summaries$acceptance_rate) * 100))
  cat(sprintf("Max acceptance rate: %.2f%%\n", 
              max(posterior_summaries$acceptance_rate) * 100))
  
} else {
  cat("\nNo genes were successfully analyzed. Check warnings for details.\n")
}

# Display warnings
cat("\n=== WARNINGS ===\n")
warnings()
