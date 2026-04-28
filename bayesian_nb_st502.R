set.seed(42)

library(ggplot2)
library(dplyr)
library(tidyr)
library(pROC)

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
# 100-SIMULATION LOOP
# ============================================================================

n_sims <- 100
n_genes <- 1000
all_sim_results <- list()

for (sim in 1:n_sims) {
  if (sim %% 10 == 0) cat(sprintf("Simulation %d/100\n", sim))

  # --- Simulate data ---
  counts <- matrix(nrow = n_genes, ncol = 100)
  counts <- as.data.frame(counts)
  for (i in seq_len(nrow(counts))) {
    rownames(counts)[i] <- paste0("Gene", i)
  }
  colnames(counts)[1:50]  <- "Control"
  colnames(counts)[51:100] <- "Case"

  alpha <- rgamma(n = n_genes, shape = 2, rate = 0.5)
  names(alpha) <- rownames(counts)

  de_p <- rbinom(n = n_genes, size = 1, prob = 0.1)
  names(de_p) <- rownames(counts)

  delta <- rnorm(n = sum(de_p == 1), mean = 1, sd = sqrt(2))
  names(delta) <- rownames(counts)[de_p == 1]

  mu <- rgamma(n = n_genes, shape = 10, rate = 1)
  names(mu) <- rownames(counts)

  for (g in 1:n_genes) {
    gene_name <- rownames(counts)[g]
    counts[g, 1:50] <- rnbinom(n = 50, mu = mu[g], size = alpha[g])
    if (de_p[g] == 1) {
      mu_case <- mu[g] * exp(delta[gene_name])
      if (is.finite(mu_case) && mu_case < 1e6) {
        counts[g, 51:100] <- rnbinom(n = 50, mu = mu_case, size = alpha[g])
      } else {
        delta[gene_name] <- log(1000 / mu[g])
        mu_case <- mu[g] * exp(delta[gene_name])
        counts[g, 51:100] <- rnbinom(n = 50, mu = mu_case, size = alpha[g])
      }
    } else {
      counts[g, 51:100] <- rnbinom(n = 50, mu = mu[g], size = alpha[g])
    }
  }

  # --- Run MH analysis on all genes ---
  results_list <- list()
  posterior_summaries <- data.frame()

  for (g in 1:n_genes) {
    result <- tryCatch({
      analyze_gene(g, counts, n_iter = 5000, burn_in = 1000)
    }, error = function(e) NULL)

    if (!is.null(result)) {
      # Only store MCMC samples for the first 5 genes per simulation
      if (g <= 5) {
        results_list[[g]] <- result
      } else {
        results_list[[g]] <- list(summary = result$summary, samples = NULL)
      }
      posterior_summaries <- rbind(posterior_summaries, result$summary)
    }
  }

  # Build true_de and true_delta aligned to posterior_summaries
  gene_indices <- as.numeric(gsub("Gene", "", posterior_summaries$gene))
  true_de_vec  <- de_p[gene_indices]
  true_delta_vec <- ifelse(
    true_de_vec == 1,
    delta[posterior_summaries$gene],
    0
  )
  true_delta_vec[is.na(true_delta_vec)] <- 0

  all_sim_results[[sim]] <- list(
    posterior_summaries = posterior_summaries,
    true_de             = true_de_vec,
    true_delta          = true_delta_vec,
    results_list        = results_list,
    sim_id              = sim
  )
}

cat("All simulations complete.\n")

# ============================================================================
# POOL RESULTS ACROSS SIMULATIONS
# ============================================================================

pooled_list <- lapply(seq_along(all_sim_results), function(sim) {
  res <- all_sim_results[[sim]]
  df  <- res$posterior_summaries
  df$true_DE    <- res$true_de
  df$true_delta <- res$true_delta
  df$sim_id     <- sim
  df
})

pooled_results <- bind_rows(pooled_list)

# Add MLE estimate of delta
# Raw counts are not retained across simulations; approximate using posterior mu_mean:
#   mean_control ≈ mu_mean,  mean_case ≈ mu_mean * exp(delta_mean)
pooled_results <- pooled_results %>%
  mutate(
    delta_MLE = log((mu_mean * exp(delta_mean) + 0.5) / (mu_mean + 0.5))
  )

# Hypothesis testing columns on pooled data
pooled_results <- pooled_results %>%
  mutate(
    reject_H0_CI   = (delta_ci_lower > 0 | delta_ci_upper < 0),
    reject_H0_ROPE = prob_delta_zero < 0.05,
    reject_H0_prob = (prob_delta_positive > 0.95 | prob_delta_positive < 0.05)
  )

# ============================================================================
# POOLED SUMMARY STATISTICS
# ============================================================================

roc_pooled <- roc(
  response  = pooled_results$true_DE,
  predictor = pooled_results$prob_delta_positive,
  levels    = c(0, 1),
  direction = "<",
  quiet     = TRUE
)
mean_auc <- as.numeric(auc(roc_pooled))

sim_aucs <- sapply(seq_len(n_sims), function(s) {
  sub_df <- pooled_results[pooled_results$sim_id == s, ]
  if (length(unique(sub_df$true_DE)) < 2) return(NA_real_)
  as.numeric(auc(roc(sub_df$true_DE, sub_df$prob_delta_positive,
                     levels = c(0, 1), direction = "<", quiet = TRUE)))
})

# Per-simulation sensitivity and specificity (CI method)
sim_sens <- sapply(seq_len(n_sims), function(s) {
  sub_df <- pooled_results[pooled_results$sim_id == s, ]
  tp <- sum(sub_df$reject_H0_CI & sub_df$true_DE == 1)
  fn <- sum(!sub_df$reject_H0_CI & sub_df$true_DE == 1)
  if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
})
sim_spec <- sapply(seq_len(n_sims), function(s) {
  sub_df <- pooled_results[pooled_results$sim_id == s, ]
  tn <- sum(!sub_df$reject_H0_CI & sub_df$true_DE == 0)
  fp <- sum(sub_df$reject_H0_CI  & sub_df$true_DE == 0)
  if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
})

cat(sprintf("\n=== POOLED RESULTS SUMMARY (%d simulations) ===\n", n_sims))
cat(sprintf("Total gene-simulation pairs analyzed: %d\n", nrow(pooled_results)))
cat(sprintf("Mean AUC (prob_delta_positive): %.4f\n", mean(sim_aucs, na.rm = TRUE)))
cat(sprintf("Mean acceptance rate: %.2f%%\n",
            mean(pooled_results$acceptance_rate) * 100))
cat(sprintf("Mean sensitivity: %.2f%%\n", mean(sim_sens, na.rm = TRUE) * 100))
cat(sprintf("Mean specificity: %.2f%%\n", mean(sim_spec, na.rm = TRUE) * 100))

# ============================================================================
# VISUALIZATIONS — save to simulation_results.pdf
# ============================================================================

pdf("simulation_results.pdf", width = 10, height = 8)

# --- A. ROC Curve (pooled) ---
roc_prob <- roc(pooled_results$true_DE, pooled_results$prob_delta_positive,
                levels = c(0, 1), direction = "<", quiet = TRUE)
roc_rope <- roc(pooled_results$true_DE, 1 - pooled_results$prob_delta_zero,
                levels = c(0, 1), direction = "<", quiet = TRUE)
roc_mle  <- roc(pooled_results$true_DE, abs(pooled_results$delta_MLE),
                levels = c(0, 1), direction = "<", quiet = TRUE)

roc_df <- bind_rows(
  data.frame(FPR = 1 - roc_prob$specificities, TPR = roc_prob$sensitivities,
             Method = sprintf("P(\u03b4>0), AUC=%.3f", as.numeric(auc(roc_prob)))),
  data.frame(FPR = 1 - roc_rope$specificities, TPR = roc_rope$sensitivities,
             Method = sprintf("ROPE, AUC=%.3f", as.numeric(auc(roc_rope)))),
  data.frame(FPR = 1 - roc_mle$specificities,  TPR = roc_mle$sensitivities,
             Method = sprintf("MLE, AUC=%.3f",  as.numeric(auc(roc_mle))))
)

pA <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Method)) +
  geom_line(linewidth = 1.1) +
  geom_abline(linetype = "dashed", color = "grey50") +
  labs(title = "ROC Curve (Pooled Across All Simulations)",
       x = "False Positive Rate", y = "True Positive Rate",
       color = "Method") +
  theme_minimal()
print(pA)

# --- B. Posterior distribution of delta — DE vs non-DE ---
pooled_results$de_label <- factor(pooled_results$true_DE,
                                  levels = c(0, 1),
                                  labels = c("Non-DE", "DE"))

pB <- ggplot(pooled_results, aes(x = delta_mean, fill = de_label, color = de_label)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Pooled Posterior Mean of \u03b4: DE vs Non-DE Genes",
       x = "Posterior Mean \u03b4", y = "Density",
       fill = "Gene Type", color = "Gene Type") +
  theme_minimal()
print(pB)

# --- C. True delta vs estimated delta (DE genes only) ---
de_only <- pooled_results %>% filter(true_DE == 1)

pC <- ggplot(de_only, aes(x = true_delta, y = delta_mean, color = factor(sim_id))) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "True \u03b4 vs Posterior Mean \u03b4 (DE Genes, Pooled)",
       x = "True \u03b4", y = "Posterior Mean \u03b4",
       color = "Simulation") +
  theme_minimal() +
  theme(legend.position = "none")
print(pC)

# --- D. MLE vs Bayesian delta comparison (DE genes only) ---
cor_val <- cor(de_only$delta_MLE, de_only$delta_mean, use = "complete.obs")

pD <- ggplot(de_only, aes(x = delta_MLE, y = delta_mean)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("r = %.3f", cor_val), size = 4) +
  labs(title = "MLE vs Bayesian Posterior Mean \u03b4",
       x = "MLE \u03b4", y = "Posterior Mean \u03b4") +
  theme_minimal()
print(pD)

# --- E. Acceptance rate distribution ---
mean_ar <- mean(pooled_results$acceptance_rate)

pE <- ggplot(pooled_results, aes(x = acceptance_rate)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean_ar, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = mean_ar, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Mean = %.3f", mean_ar), color = "red") +
  labs(title = "Distribution of MH Acceptance Rates (Pooled)",
       x = "Acceptance Rate", y = "Count") +
  theme_minimal()
print(pE)

# --- F. Trace plots — representative DE gene across simulations ---
# Only Gene1-Gene5 have MCMC samples stored; find one that is DE in >= 10 sims
de_gene_counts <- table(
  pooled_results$gene[pooled_results$true_DE == 1 &
                        as.numeric(gsub("Gene", "", pooled_results$gene)) <= 5]
)
rep_gene_candidates <- names(de_gene_counts)[de_gene_counts >= 10]

trace_plotted <- FALSE
if (length(rep_gene_candidates) > 0) {
  rep_gene <- rep_gene_candidates[1]
  gene_num  <- as.numeric(gsub("Gene", "", rep_gene))

  # Collect trace data from simulations where this gene was DE and has samples
  trace_list <- list()
  for (s in seq_len(n_sims)) {
    if (length(trace_list) >= 9) break
    res_s <- all_sim_results[[s]]
    if (is.null(res_s$results_list[[gene_num]])) next
    if (is.null(res_s$results_list[[gene_num]]$samples)) next
    samps <- res_s$results_list[[gene_num]]$samples
    trace_list[[length(trace_list) + 1]] <- data.frame(
      iteration = seq_len(nrow(samps)),
      delta     = samps[, "delta"],
      sim_id    = s
    )
  }

  if (length(trace_list) > 0) {
    trace_df <- bind_rows(trace_list)
    pF <- ggplot(trace_df, aes(x = iteration, y = delta)) +
      geom_line(alpha = 0.7, color = "steelblue") +
      facet_wrap(~sim_id, ncol = 3) +
      labs(title = sprintf("Trace Plots for \u03b4 \u2014 %s Across Simulations", rep_gene),
           x = "Iteration", y = "\u03b4") +
      theme_minimal()
    print(pF)
    trace_plotted <- TRUE
  }
}

if (!trace_plotted) {
  # Fallback: trace for first gene with samples across first 9 simulations
  trace_list <- list()
  for (s in seq_len(min(9, n_sims))) {
    res_s <- all_sim_results[[s]]
    for (gnum in seq_len(5)) {
      if (!is.null(res_s$results_list[[gnum]]) &&
          !is.null(res_s$results_list[[gnum]]$samples)) {
        samps   <- res_s$results_list[[gnum]]$samples
        gene_nm <- paste0("Gene", gnum)
        trace_list[[length(trace_list) + 1]] <- data.frame(
          iteration = seq_len(nrow(samps)),
          delta     = samps[, "delta"],
          sim_id    = s,
          gene      = gene_nm
        )
        break
      }
    }
  }
  if (length(trace_list) > 0) {
    trace_df <- bind_rows(trace_list)
    pF <- ggplot(trace_df, aes(x = iteration, y = delta)) +
      geom_line(alpha = 0.7, color = "steelblue") +
      facet_wrap(~sim_id, ncol = 3) +
      labs(title = "Trace Plots for \u03b4 \u2014 Representative Gene Across Simulations",
           x = "Iteration", y = "\u03b4") +
      theme_minimal()
    print(pF)
  }
}

# --- G. Simulation-level AUC distribution ---
sim_auc_df <- data.frame(sim_id = seq_len(n_sims), auc = sim_aucs)
mean_sim_auc <- mean(sim_auc_df$auc, na.rm = TRUE)

pG <- ggplot(sim_auc_df, aes(x = auc)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean_sim_auc, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = mean_sim_auc, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Mean = %.3f", mean_sim_auc), color = "red") +
  labs(title = "Distribution of AUC Across 100 Simulations",
       x = "AUC", y = "Count") +
  theme_minimal()
print(pG)

# --- H. Confusion matrix heatmap (pooled, CI method) ---
cm_df <- pooled_results %>%
  mutate(
    Predicted = ifelse(reject_H0_CI, "Reject H0", "Fail to Reject"),
    True      = ifelse(true_DE == 1, "DE", "Non-DE")
  ) %>%
  count(Predicted, True)

pH <- ggplot(cm_df, aes(x = True, y = Predicted, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Pooled Confusion Matrix (CI Method)",
       x = "True DE Status", y = "Prediction",
       fill = "Count") +
  theme_minimal()
print(pH)

dev.off()
cat("Plots saved to simulation_results.pdf\n")
