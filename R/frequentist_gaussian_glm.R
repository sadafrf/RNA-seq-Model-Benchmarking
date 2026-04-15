# Frequentist Gaussian GLM for differential expression
#
# Model for one gene:
#   Y = X beta + epsilon
#   epsilon ~ Normal(0, sigma^2 I)
#
# With two groups:
#   beta0 = control-group mean
#   beta1 = case mean - control mean
#
# This implementation uses the closed-form least-squares estimate
#   beta_hat = (X'X)^(-1) X'Y
# and estimates standard errors using the observed residual variance.

library(ggplot2)
library(dplyr)
library(tidyr)

build_design_matrix <- function(group) {
  cbind(
    intercept = 1,
    case = as.integer(group == 1 | group == "Case" | group == "case")
  )
}

fit_gaussian_gene <- function(y, X) {
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(X)

  if (n <= p) {
    stop("Need more observations than parameters.")
  }

  xtx <- crossprod(X)
  xty <- crossprod(X, y)

  # Cholesky-based solve for beta_hat = (X'X)^(-1) X'Y.
  chol_xtx <- chol(xtx)
  beta_hat <- solve(chol_xtx, solve(t(chol_xtx), xty))
  beta_hat <- as.vector(beta_hat)
  names(beta_hat) <- colnames(X)

  fitted_values <- as.vector(X %*% beta_hat)
  residuals <- y - fitted_values
  sigma2_hat <- sum(residuals^2) / (n - p)

  xtx_inv <- chol2inv(chol_xtx)
  beta_var <- sigma2_hat * xtx_inv
  beta_se <- sqrt(diag(beta_var))
  names(beta_se) <- colnames(X)

  t_stat <- beta_hat / beta_se
  p_value <- 2 * pt(abs(t_stat), df = n - p, lower.tail = FALSE)

  list(
    beta_hat = beta_hat,
    beta_se = beta_se,
    beta_var = beta_var,
    sigma2_hat = sigma2_hat,
    fitted_values = fitted_values,
    residuals = residuals,
    t_stat = t_stat,
    p_value = p_value,
    df = n - p
  )
}

analyze_gaussian_genes <- function(counts_data, group, alpha = 0.05) {
  X <- build_design_matrix(group)

  summaries <- lapply(seq_len(nrow(counts_data)), function(gene_idx) {
    fit <- fit_gaussian_gene(counts_data[gene_idx, ], X)

    beta1 <- fit$beta_hat["case"]
    beta1_se <- fit$beta_se["case"]
    critical_value <- qt(1 - alpha / 2, df = fit$df)

    data.frame(
      gene = rownames(counts_data)[gene_idx],
      beta0_hat = fit$beta_hat["intercept"],
      beta1_hat = beta1,
      beta1_se = beta1_se,
      beta1_t = fit$t_stat["case"],
      beta1_p_value = fit$p_value["case"],
      beta1_ci_lower = beta1 - critical_value * beta1_se,
      beta1_ci_upper = beta1 + critical_value * beta1_se,
      sigma2_hat = fit$sigma2_hat,
      reject_H0 = fit$p_value["case"] < alpha
    )
  })

  bind_rows(summaries) %>%
    mutate(
      beta1_p_adj_BH = p.adjust(beta1_p_value, method = "BH"),
      reject_H0_BH = beta1_p_adj_BH < alpha
    )
}

simulate_gaussian_gene_data <- function(
    n_genes = 1000,
    n_control = 50,
    n_case = 50,
    de_prob = 0.10,
    seed = 1
) {
  set.seed(seed)

  n_samples <- n_control + n_case
  gene_names <- paste0("Gene", seq_len(n_genes))
  group <- c(rep(0, n_control), rep(1, n_case))

  beta0_true <- rgamma(n_genes, shape = 10, rate = 1)
  beta1_true <- rep(0, n_genes)
  true_DE <- rbinom(n_genes, size = 1, prob = de_prob)
  beta1_true[true_DE == 1] <- rnorm(sum(true_DE), mean = 2, sd = 1)
  sigma_true <- rgamma(n_genes, shape = 3, rate = 1)

  counts <- matrix(NA_real_, nrow = n_genes, ncol = n_samples)
  rownames(counts) <- gene_names
  colnames(counts) <- c(
    paste0("Control", seq_len(n_control)),
    paste0("Case", seq_len(n_case))
  )

  for (g in seq_len(n_genes)) {
    mu <- beta0_true[g] + beta1_true[g] * group
    counts[g, ] <- rnorm(n_samples, mean = mu, sd = sigma_true[g])
  }

  truth <- data.frame(
    gene = gene_names,
    true_DE = true_DE,
    beta0_true = beta0_true,
    beta1_true = beta1_true,
    sigma_true = sigma_true
  )

  list(
    counts = as.data.frame(counts),
    group = group,
    truth = truth
  )
}

plot_gaussian_results <- function(results) {
  p_values <- ggplot(results, aes(x = beta1_p_value)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8) +
    labs(
      title = "Gaussian GLM p-value distribution",
      x = "p-value for beta1",
      y = "Number of genes"
    ) +
    theme_minimal()

  effects <- ggplot(results, aes(x = beta1_hat, y = -log10(beta1_p_value))) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Gaussian GLM effect sizes",
      x = "beta1 estimate: case mean - control mean",
      y = "-log10(p-value)"
    ) +
    theme_minimal()

  list(
    p_values = p_values,
    effects = effects
  )
}

# Example run.
# If your teammate's simulation already creates `counts`, you can skip this
# simulation block and run:
#   group <- c(rep(0, 50), rep(1, 50))
#   gaussian_results <- analyze_gaussian_genes(counts, group)

sim <- simulate_gaussian_gene_data(n_genes = 1000, seed = 123)

gaussian_results <- analyze_gaussian_genes(
  counts_data = sim$counts,
  group = sim$group,
  alpha = 0.05
) %>%
  left_join(sim$truth, by = "gene")

print(head(gaussian_results, 20))

cat("\n=== FREQUENTIST GAUSSIAN GLM SUMMARY ===\n")
cat(sprintf("Total genes analyzed: %d\n", nrow(gaussian_results)))
cat(sprintf("Rejected H0 at alpha = 0.05: %d/%d\n",
            sum(gaussian_results$reject_H0),
            nrow(gaussian_results)))
cat(sprintf("Rejected H0 after BH correction: %d/%d\n",
            sum(gaussian_results$reject_H0_BH),
            nrow(gaussian_results)))

confusion_matrix <- table(
  Predicted = gaussian_results$reject_H0_BH,
  True = gaussian_results$true_DE
)

cat("\n=== COMPARISON WITH TRUE DE STATUS ===\n")
print(confusion_matrix)

plots <- plot_gaussian_results(gaussian_results)
print(plots$p_values)
print(plots$effects)

