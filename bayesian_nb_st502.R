library(ggplot2)
library(dplyr)
library(tidyr)


counts <- matrix(nrow=1000,ncol=100)
counts <- as.data.frame(counts)
for (i in (1:nrow(counts))){
  rownames(counts)[i] <- paste0("Gene", i) 
}
all_posterior_summaries <- data.frame()   # collects posterior_summaries from every iteration
all_conv_df <- data.frame()   # collects conv_df from every iteration
all_results_list <- list()         # optional: keeps all per-gene MH samples
iteration_metadata <- data.frame()  # tracks simulation-level stats per iteration
all_results <- data.frame()


for(i in 1:100) {
  #Initialize alpha dispersion parameter with gamma prior. This prior is uninformative.
  #The names for the dispersion lit are set to be the rownames for the count matrix
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
    mh_results <- metropolis_hastings(y_control, y_case, 
                                      n_iter = n_iter, 
                                      burn_in = burn_in,
                                      proposal_sd = c(0.15, 0.3, 0.5))
    
    
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
    
    posterior_summary <- posterior_summary %>% mutate(z = abs(delta_mean / delta_sd))
    posterior_summary <- posterior_summary[order(posterior_summary$z, decreasing = TRUE), ]
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
  thresholds_list <- list()
  results <- data.frame()
  
  for (g in 1:n_genes_to_analyze) {
    if (g %% 500 == 0) cat(sprintf("Processing gene %d/%d\n", g, n_genes_to_analyze))
    result <- analyze_gene(g, counts, n_iter = 5000, burn_in = 1000)
    results_list[[g]] <- result
    posterior_summaries <- rbind(posterior_summaries, result$summary)
  }
  
  gene_indices <- as.numeric(gsub("Gene", "", posterior_summaries$gene))
  posterior_summaries$true_DE <- de_p[gene_indices]
  
  for (j in 1:nrow(posterior_summaries)) {
      posterior_summaries$predicted_DE <- FALSE
      posterior_summaries$predicted_DE[1:j] <- TRUE
      true <- posterior_summaries$true_DE
      pred <- posterior_summaries$predicted_DE
      TP <- sum(pred & true)
      FP <- sum(pred & !true)
      TN <- sum(!pred & !true)
      FN <- sum(!pred & true)
      precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
      recall    <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)  # TPR
      TPR       <- recall
      FPR       <- ifelse((FP + TN) > 0, FP / (FP + TN), NA)
      results <- rbind(results, data.frame(
        threshold_index = j,
        z_threshold = posterior_summaries$z[j],
        precision = precision,
        recall = recall,
        TPR = TPR,
        FPR = FPR
      ))
  }
  cat(sprintf("\nMetropolis-Hastings complete!\n"))
  
  # ============================================================================
  # HYPOTHESIS TESTING: Is δ = 0?
  # ============================================================================
  
  posterior_summaries <- posterior_summaries %>%
      mutate(
        reject_H0_CI = (delta_ci_lower > 0 | delta_ci_upper < 0),
        reject_H0_ROPE = prob_delta_zero < 0.05,
        reject_H0_prob = (prob_delta_positive > 0.95 | prob_delta_positive < 0.05)
      )
  posterior_summaries$iteration <- i
  results$iteration <- i
  all_posterior_summaries <- rbind(all_posterior_summaries, posterior_summaries)
  all_results <- rbind(all_results, results)
    
    # ============================================================================
    # VISUALIZATION FOR FIRST SUCCESSFULLY ANALYZED GENE
    # ============================================================================
    
    # Find first non-NULL result
  first_valid_idx <- which(!sapply(results_list, is.null))[9]
  if (!is.na(first_valid_idx) && !is.null(results_list[[first_valid_idx]])) {
      gene_samples <- results_list[[first_valid_idx]]$samples
      gene_df <- as.data.frame(gene_samples)
      gene_name <- posterior_summaries$gene[9]
    
    # Summary of hypothesis tests
  cat(sprintf("Genes with δ ≠ 0 (CI method): %d/%d\n", 
                sum(posterior_summaries$reject_H0_CI), 
                nrow(posterior_summaries)))
  cat(sprintf("Genes with δ ≠ 0 (ROPE method): %d/%d\n", 
              sum(posterior_summaries$reject_H0_ROPE), 
                nrow(posterior_summaries)))
    
    # Compare with true DE status
  if (n_genes_to_analyze <= length(de_p)) {
  # Map true delta values onto every gene (0 for non-DE genes)
    true_delta_vec <- rep(0, nrow(posterior_summaries))
    names(true_delta_vec) <- posterior_summaries$gene
    matched <- intersect(names(delta), names(true_delta_vec))
    true_delta_vec[matched] <- delta[matched]
    posterior_summaries$true_delta <- ifelse(posterior_summaries$true_DE == 1,true_delta_vec,0)
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
  
  iteration_metadata <- rbind(iteration_metadata, data.frame(
    iteration   = i,
    n_de_genes  = sum(de_p),
    n_analyzed  = nrow(posterior_summaries),
    n_failed    = length(failed_genes),
    mean_accept = mean(posterior_summaries$acceptance_rate),
    sensitivity = confusion_matrix["TRUE", "1"] / sum(confusion_matrix[, "1"]),
    specificity = confusion_matrix["FALSE", "0"] / sum(confusion_matrix[, "0"]),
    accuracy = sum(diag(confusion_matrix)) / sum(confusion_matrix),
    acceptance_rate = posterior_summaries$acceptance_rate,
    true_negative = confusion_matrix[1,1],
    true_positive = confusion_matrix[2,2],
    false_negative = confusion_matrix[1,2],
    false_positive = confusion_matrix[2,1]
  ))
  cat("\n=== WARNINGS ===\n")
  warnings()
}


write.table(iteration_metadata, "iteration_metadata.txt", sep = "\t")
write.table(all_results, "z_thresholding_for_auc_curve.txt", sep = "\t")

metrics_long <- iteration_metadata %>%
  select(iteration, sensitivity, specificity, accuracy) %>%
  distinct() %>%                          # drop duplicated rows from acceptance_rate expansion
  pivot_longer(cols = c(sensitivity, specificity, accuracy),
               names_to  = "metric",
               values_to = "value")

pdf("MCMC_metrics.pdf")
ggplot(metrics_long, aes(x = iteration, y = value, color = metric)) +
  geom_line(alpha = 0.7) +
  geom_smooth(se = TRUE, method = "loess", linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_color_manual(values = c(accuracy    = "#2196F3",
                                sensitivity = "#4CAF50",
                                specificity = "#F44336")) +
  labs(title    = "Classification Performance Across 100 Simulations",
       x        = "Simulation Iteration",
       y        = "Value",
       color    = "Metric") +
  theme_minimal()

# ── 2. Confusion matrix counts across iterations ──────────────────────────────
cm_long <- iteration_metadata %>%
  select(iteration, true_positive, true_negative, false_positive, false_negative) %>%
  distinct() %>%
  pivot_longer(cols = c(true_positive, true_negative, false_positive, false_negative),
               names_to  = "outcome",
               values_to = "count")

ggplot(cm_long, aes(x = iteration, y = count, color = outcome)) +
  geom_line(alpha = 0.6) +
  geom_smooth(se = FALSE, method = "loess", linewidth = 0.8) +
  scale_color_manual(values = c(true_positive  = "#4CAF50",
                                true_negative  = "#2196F3",
                                false_positive = "#FF9800",
                                false_negative = "#F44336")) +
  labs(title  = "Confusion Matrix Counts Across 100 Simulations",
       x      = "Simulation Iteration",
       y      = "Gene Count",
       color  = "Outcome") +
  theme_minimal()

# ── 3. Number of DE genes simulated vs detected ───────────────────────────────
de_detect <- iteration_metadata %>%
  select(iteration, n_de_genes, true_positive) %>%
  distinct() %>%
  pivot_longer(cols = c(n_de_genes, true_positive),
               names_to  = "type",
               values_to = "count")

ggplot(de_detect, aes(x = iteration, y = count, color = type)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values  = c(n_de_genes   = "steelblue",
                                 true_positive = "darkgreen"),
                     labels  = c(n_de_genes    = "True DE Genes (simulated)",
                                 true_positive = "True Positives (detected)")) +
  labs(title  = "DE Genes Simulated vs Correctly Detected",
       x      = "Simulation Iteration",
       y      = "Gene Count",
       color  = NULL) +
  theme_minimal()

# ── 4. Mean acceptance rate across iterations ─────────────────────────────────
accept_df <- iteration_metadata %>%
  select(iteration, mean_accept) %>%
  distinct()

ggplot(accept_df, aes(x = iteration, y = mean_accept)) +
  geom_line(color = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = c(0.2, 0.5), linetype = "dashed",
             color = c("red", "green"), linewidth = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  annotate("text", x = 5, y = 0.21, label = "20% (min ideal)", size = 3, color = "red") +
  annotate("text", x = 5, y = 0.51, label = "50% (max ideal)", size = 3, color = "darkgreen") +
  labs(title  = "Mean MH Acceptance Rate Across 100 Simulations",
       x      = "Simulation Iteration",
       y      = "Mean Acceptance Rate") +
  theme_minimal()

# ── 5. Distributions of key metrics (boxplots / violin) ───────────────────────
dist_df <- iteration_metadata %>%
  select(iteration, sensitivity, specificity, accuracy) %>%
  distinct() %>%
  pivot_longer(cols = -iteration, names_to = "metric", values_to = "value")

ggplot(dist_df, aes(x = metric, y = value, fill = metric)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = 21) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_fill_manual(values = c(accuracy    = "#2196F3",
                               sensitivity = "#4CAF50",
                               specificity = "#F44336")) +
  labs(title = "Distribution of Performance Metrics Across 100 Simulations",
       x     = NULL,
       y     = "Value") +
  theme_minimal() +
  theme(legend.position = "none")
  

roc_avg <- all_results %>%
  group_by(threshold_index) %>%
  summarise(
    mean_TPR = mean(TPR, na.rm = TRUE),
    mean_FPR = mean(FPR, na.rm = TRUE)
  )


ggplot(roc_avg, aes(x = mean_FPR, y = mean_TPR)) +
  geom_line(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Average ROC Curve (100 Simulations)",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal()

dev.off()
