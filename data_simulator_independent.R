# Independent negative-binomial count simulator

validate_independent_simulator_inputs <- function(n,
                                                  p,
                                                  alpha,
                                                  lambda,
                                                  de_count) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer giving the number of subjects per group.")
  }

  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p <= 0 || p != as.integer(p)) {
    stop("p must be a positive integer giving the number of genes.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a positive scalar.")
  }

  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("lambda must be a positive scalar.")
  }

  if (!is.numeric(de_count) || length(de_count) != 1L || !is.finite(de_count) ||
      de_count < 0 || de_count != as.integer(de_count)) {
    stop("de_count must be a non-negative integer.")
  }

  if (de_count > p) {
    stop("de_count cannot exceed p.")
  }
}

simulate_independent_nb_data <- function(n,
                                         p,
                                         alpha,
                                         lambda,
                                         de_count = 100L,
                                         seed = NULL,
                                         include_dimnames = TRUE) {
  validate_independent_simulator_inputs(
    n = n,
    p = p,
    alpha = alpha,
    lambda = lambda,
    de_count = de_count
  )

  n <- as.integer(n)
  p <- as.integer(p)
  de_count <- as.integer(de_count)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  I_DE <- integer(p)
  de_idx <- if (de_count > 0L) sort(sample.int(p, size = de_count, replace = FALSE)) else integer(0)
  I_DE[de_idx] <- 1L

  beta0 <- rnorm(p, mean = 0, sd = 1)
  beta1 <- rnorm(p, mean = 0, sd = 1)

  mu0 <- exp(beta0)
  mu1 <- exp(beta0 + beta1 * I_DE)
  phi <- rgamma(p, shape = alpha, rate = lambda)

  Y0 <- vapply(
    seq_len(p),
    function(j) rnbinom(n = n, size = phi[j], mu = mu0[j]),
    numeric(n)
  )

  Y1 <- vapply(
    seq_len(p),
    function(j) rnbinom(n = n, size = phi[j], mu = mu1[j]),
    numeric(n)
  )

  Y <- rbind(Y0, Y1)
  group <- factor(
    c(rep("control", n), rep("case", n)),
    levels = c("control", "case")
  )

  if (include_dimnames) {
    gene_names <- paste0("gene_", seq_len(p))
    rownames(Y0) <- paste0("control_", seq_len(n))
    rownames(Y1) <- paste0("case_", seq_len(n))
    rownames(Y) <- c(rownames(Y0), rownames(Y1))
    colnames(Y0) <- gene_names
    colnames(Y1) <- gene_names
    colnames(Y) <- gene_names
    names(I_DE) <- gene_names
    names(beta0) <- gene_names
    names(beta1) <- gene_names
    names(mu0) <- gene_names
    names(mu1) <- gene_names
    names(phi) <- gene_names
  }

  gene_parameters <- data.frame(
    gene = if (include_dimnames) colnames(Y) else seq_len(p),
    I_DE = I_DE,
    beta0 = beta0,
    beta1 = beta1,
    mu0 = mu0,
    mu1 = mu1,
    phi = phi,
    stringsAsFactors = FALSE
  )

  list(
    Y = Y,
    Y0 = Y0,
    Y1 = Y1,
    group = group,
    ground_truth_I_DE = I_DE,
    I_DE = I_DE,
    de_idx = de_idx,
    gene_parameters = gene_parameters,
    beta0 = beta0,
    beta1 = beta1,
    mu0 = mu0,
    mu1 = mu1,
    phi = phi,
    n_per_group = n,
    n_total = 2L * n,
    p = p,
    alpha = alpha,
    lambda = lambda,
    de_count = de_count,
    seed = seed
  )
}

print_independent_nb_simulation <- function(sim) {
  if (!is.list(sim) || is.null(sim$Y) || is.null(sim$group) || is.null(sim$phi)) {
    stop("sim must be the result of simulate_independent_nb_data().")
  }

  cat("Independent negative-binomial simulation\n")
  cat("  subjects per group =", sim$n_per_group, "\n")
  cat("  total subjects     =", sim$n_total, "\n")
  cat("  genes              =", sim$p, "\n")
  cat("  DE genes           =", sim$de_count, "\n")
  cat("  mean(phi)          =", signif(mean(sim$phi), 6), "\n")
  invisible(sim)
}


# simulate counts for 1000 genes ---------------------------------------------------------


sim <- simulate_independent_nb_data(
  n = 1000,
  p = 1000,
  alpha = 2,
  lambda = 1,
  de_count = 100,
  seed = 1
)
dim(sim$Y)
head(sim$gene_parameters)
sim$Y[1:10, 1:10]

setwd('C:/Users/dcginger/Desktop/st502/project/data')
saveRDS(sim, 'data_simulator_independent___simulated_counts.rds')
