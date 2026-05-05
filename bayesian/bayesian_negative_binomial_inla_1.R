# negative-binomial INLA-style posterior approximation

logspace_add <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

validate_nb_inla_inputs <- function(y, X, tau2, a, b, varphi_grid, beta_init) {
  if (!is.numeric(y) || any(!is.finite(y)) || any(y < 0)) {
    stop("y must be a finite numeric vector of non-negative counts.")
  }

  X <- as.matrix(X)
  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix.")
  }

  n <- length(y)
  if (nrow(X) != n) {
    stop("nrow(X) must equal length(y).")
  }

  if (!is.numeric(tau2) || length(tau2) != 1L || !is.finite(tau2) || tau2 <= 0) {
    stop("tau2 must be a positive scalar.")
  }

  if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a <= 0) {
    stop("a must be a positive scalar.")
  }

  if (!is.numeric(b) || length(b) != 1L || !is.finite(b) || b <= 0) {
    stop("b must be a positive scalar.")
  }

  if (!is.numeric(varphi_grid) || any(!is.finite(varphi_grid))) {
    stop("varphi_grid must be a finite numeric vector.")
  }

  if (length(unique(varphi_grid)) != length(varphi_grid)) {
    stop("varphi_grid must not contain duplicate values.")
  }

  p <- ncol(X)
  if (is.null(beta_init)) {
    beta_init <- rep(0, p)
  } else {
    beta_init <- as.numeric(beta_init)
    if (length(beta_init) != p || any(!is.finite(beta_init))) {
      stop("beta_init must be NULL or a finite numeric vector of length ncol(X).")
    }
  }

  list(
    y = as.numeric(y),
    X = X,
    tau2 = tau2,
    a = a,
    b = b,
    varphi_grid = as.numeric(varphi_grid),
    beta_init = beta_init
  )
}

nb_eta_mu_phi <- function(beta, varphi, X) {
  eta <- drop(X %*% beta)
  mu <- exp(eta)
  phi <- exp(varphi)
  list(eta = eta, mu = mu, phi = phi)
}

nb_log_likelihood <- function(beta, varphi, y, X) {
  pieces <- nb_eta_mu_phi(beta, varphi, X)
  eta <- pieces$eta
  mu <- pieces$mu
  phi <- pieces$phi
  log_denom <- logspace_add(eta, varphi)

  sum(
    lgamma(y + phi) -
      lgamma(phi) -
      lgamma(y + 1) +
      phi * varphi +
      y * eta -
      (y + phi) * log_denom
  )
}

nb_log_prior <- function(beta, varphi, tau2, a, b) {
  p <- length(beta)

  -(p / 2) * log(2 * pi * tau2) -
    sum(beta * beta) / (2 * tau2) +
    a * log(b) -
    lgamma(a) +
    a * varphi -
    b * exp(varphi)
}

nb_log_posterior <- function(beta, varphi, y, X, tau2, a, b) {
  nb_log_likelihood(beta, varphi, y, X) +
    nb_log_prior(beta, varphi, tau2, a, b)
}

nb_gradient_log_likelihood <- function(beta, varphi, y, X) {
  pieces <- nb_eta_mu_phi(beta, varphi, X)
  mu <- pieces$mu
  phi <- pieces$phi

  beta_grad <- drop(crossprod(X, phi * (y - mu) / (mu + phi)))

  g_i <- digamma(y + phi) -
    digamma(phi) +
    varphi + 1 -
    log(mu + phi) -
    (y + phi) / (mu + phi)

  varphi_grad <- sum(phi * g_i)

  list(beta = beta_grad, varphi = varphi_grad)
}

nb_gradient_log_prior <- function(beta, varphi, tau2, a, b) {
  list(
    beta = -beta / tau2,
    varphi = a - b * exp(varphi)
  )
}

nb_gradient_log_posterior <- function(beta, varphi, y, X, tau2, a, b) {
  grad_lik <- nb_gradient_log_likelihood(beta, varphi, y, X)
  grad_prior <- nb_gradient_log_prior(beta, varphi, tau2, a, b)

  c(
    grad_lik$beta + grad_prior$beta,
    grad_lik$varphi + grad_prior$varphi
  )
}

nb_hessian_log_likelihood <- function(beta, varphi, y, X) {
  pieces <- nb_eta_mu_phi(beta, varphi, X)
  mu <- pieces$mu
  phi <- pieces$phi

  w <- phi * mu * (y + phi) / (mu + phi)^2
  h_bb <- -crossprod(X, X * w)

  c_vec <- phi * mu * (y - mu) / (mu + phi)^2
  h_bv <- matrix(drop(crossprod(X, c_vec)), ncol = 1L)

  g_i <- digamma(y + phi) -
    digamma(phi) +
    varphi + 1 -
    log(mu + phi) -
    (y + phi) / (mu + phi)

  h_i <- trigamma(y + phi) -
    trigamma(phi) +
    1 / phi -
    1 / (mu + phi) -
    (mu - y) / (mu + phi)^2

  h_vv <- sum(phi * g_i + phi^2 * h_i)

  list(
    beta_beta = h_bb,
    beta_varphi = h_bv,
    varphi_varphi = h_vv
  )
}

nb_hessian_log_prior <- function(p, varphi, tau2, b) {
  list(
    beta_beta = -diag(1 / tau2, p),
    beta_varphi = matrix(0, nrow = p, ncol = 1L),
    varphi_varphi = -b * exp(varphi)
  )
}

nb_hessian_log_posterior <- function(beta, varphi, y, X, tau2, a, b) {
  p <- length(beta)
  h_lik <- nb_hessian_log_likelihood(beta, varphi, y, X)
  h_prior <- nb_hessian_log_prior(p, varphi, tau2, b)

  h_bb <- h_lik$beta_beta + h_prior$beta_beta
  h_bv <- h_lik$beta_varphi + h_prior$beta_varphi
  h_vv <- h_lik$varphi_varphi + h_prior$varphi_varphi

  rbind(
    cbind(h_bb, h_bv),
    cbind(t(h_bv), h_vv)
  )
}

nb_beta_gradient <- function(beta, varphi, y, X, tau2) {
  grad_lik <- nb_gradient_log_likelihood(beta, varphi, y, X)$beta
  grad_lik - beta / tau2
}

nb_beta_hessian <- function(beta, varphi, y, X, tau2) {
  h_lik <- nb_hessian_log_likelihood(beta, varphi, y, X)$beta_beta
  h_lik - diag(1 / tau2, length(beta))
}

nb_beta_log_posterior <- function(beta, varphi, y, X, tau2) {
  nb_log_likelihood(beta, varphi, y, X) - sum(beta * beta) / (2 * tau2)
}

safe_chol <- function(M, base_jitter = 1e-10, max_tries = 8L) {
  stopifnot(nrow(M) == ncol(M))

  jitter <- 0
  eye <- diag(nrow(M))

  for (attempt in seq_len(max_tries)) {
    candidate <- if (jitter == 0) M else M + jitter * eye
    chol_fit <- tryCatch(chol(candidate), error = function(e) NULL)
    if (!is.null(chol_fit)) {
      return(list(chol = chol_fit, jitter = jitter))
    }
    jitter <- if (jitter == 0) base_jitter else jitter * 10
  }

  stop("Cholesky factorization failed even after diagonal jitter.")
}

solve_chol_system <- function(chol_upper, rhs) {
  z <- forwardsolve(t(chol_upper), rhs)
  backsolve(chol_upper, z)
}

grid_cell_widths <- function(grid) {
  k <- length(grid)
  if (k == 1L) {
    return(1)
  }

  widths <- numeric(k)
  widths[1] <- (grid[2] - grid[1]) / 2
  widths[k] <- (grid[k] - grid[k - 1]) / 2

  if (k > 2L) {
    widths[2:(k - 1)] <- (grid[3:k] - grid[1:(k - 2)]) / 2
  }

  widths
}

normalize_log_grid_density <- function(log_values, grid) {
  delta <- grid_cell_widths(grid)
  max_log <- max(log_values)
  scaled <- exp(log_values - max_log)
  normalizer <- sum(scaled * delta)

  if (!is.finite(normalizer) || normalizer <= 0) {
    stop("Failed to normalize the varphi grid posterior.")
  }

  log_norm_const <- max_log + log(normalizer)
  density <- exp(log_values - log_norm_const)
  weight <- density * delta

  list(
    delta = delta,
    log_norm_const = log_norm_const,
    density = density,
    weight = weight / sum(weight)
  )
}

nb_beta_mode_fixed_varphi <- function(varphi,
                                      y,
                                      X,
                                      tau2,
                                      beta_init = NULL,
                                      tol = 1e-8,
                                      maxit = 100L,
                                      line_search_halvings = 25L,
                                      jitter = 1e-10) {
  p <- ncol(X)
  beta <- if (is.null(beta_init)) rep(0, p) else as.numeric(beta_init)

  current_lp <- nb_beta_log_posterior(beta, varphi, y, X, tau2)

  converged <- FALSE
  iterations <- 0L

  for (iter in seq_len(maxit)) {
    iterations <- iter
    grad <- nb_beta_gradient(beta, varphi, y, X, tau2)
    grad_max <- max(abs(grad))

    if (!is.finite(grad_max)) {
      stop("Encountered a non-finite beta gradient during optimization.")
    }

    if (grad_max < tol) {
      converged <- TRUE
      break
    }

    q_mat <- -nb_beta_hessian(beta, varphi, y, X, tau2)
    chol_fit <- safe_chol(q_mat, base_jitter = jitter)
    step <- solve_chol_system(chol_fit$chol, grad)

    accepted <- FALSE
    step_scale <- 1
    beta_candidate <- beta
    candidate_lp <- current_lp

    for (ls_iter in seq_len(line_search_halvings)) {
      beta_trial <- beta + step_scale * step
      lp_trial <- nb_beta_log_posterior(beta_trial, varphi, y, X, tau2)

      if (is.finite(lp_trial) && lp_trial >= current_lp) {
        beta_candidate <- beta_trial
        candidate_lp <- lp_trial
        accepted <- TRUE
        break
      }

      step_scale <- step_scale / 2
    }

    if (!accepted) {
      break
    }

    beta <- beta_candidate
    current_lp <- candidate_lp

    if (max(abs(step_scale * step)) < tol) {
      converged <- TRUE
      break
    }
  }

  h_bb <- nb_beta_hessian(beta, varphi, y, X, tau2)
  q_star <- -h_bb
  chol_q <- safe_chol(q_star, base_jitter = jitter)

  list(
    beta = beta,
    log_posterior = current_lp,
    gradient = nb_beta_gradient(beta, varphi, y, X, tau2),
    hessian = h_bb,
    q_star = q_star,
    chol_q = chol_q$chol,
    q_jitter = chol_q$jitter,
    converged = converged,
    iterations = iterations
  )
}

negative_binomial_inla <- function(y,
                                   X,
                                   tau2,
                                   a,
                                   b,
                                   varphi_grid,
                                   beta_init = NULL,
                                   newton_tol = 1e-8,
                                   newton_maxit = 100L,
                                   line_search_halvings = 25L,
                                   jitter = 1e-10,
                                   return_covariance = TRUE) {
  checked <- validate_nb_inla_inputs(y, X, tau2, a, b, varphi_grid, beta_init)

  y <- checked$y
  X <- checked$X
  tau2 <- checked$tau2
  a <- checked$a
  b <- checked$b
  beta_init <- checked$beta_init

  sort_idx <- order(checked$varphi_grid)
  varphi_grid_sorted <- checked$varphi_grid[sort_idx]

  n <- length(y)
  p <- ncol(X)
  k <- length(varphi_grid_sorted)

  beta_modes <- matrix(NA_real_, nrow = k, ncol = p)
  log_tilde_varphi <- numeric(k)
  q_list <- vector("list", k)
  covariance_list <- vector("list", k)
  convergence <- logical(k)
  iterations <- integer(k)
  q_jitter <- numeric(k)

  beta_start <- beta_init

  for (idx in seq_len(k)) {
    varphi_k <- varphi_grid_sorted[idx]

    beta_fit <- nb_beta_mode_fixed_varphi(
      varphi = varphi_k,
      y = y,
      X = X,
      tau2 = tau2,
      beta_init = beta_start,
      tol = newton_tol,
      maxit = newton_maxit,
      line_search_halvings = line_search_halvings,
      jitter = jitter
    )

    beta_hat <- beta_fit$beta
    beta_modes[idx, ] <- beta_hat
    q_list[[idx]] <- beta_fit$q_star
    q_jitter[idx] <- beta_fit$q_jitter
    convergence[idx] <- beta_fit$converged
    iterations[idx] <- beta_fit$iterations

    log_det_q <- 2 * sum(log(diag(beta_fit$chol_q)))
    log_tilde_varphi[idx] <-
      nb_log_likelihood(beta_hat, varphi_k, y, X) -
      sum(beta_hat * beta_hat) / (2 * tau2) +
      a * varphi_k -
      b * exp(varphi_k) -
      0.5 * log_det_q

    if (return_covariance) {
      covariance_list[[idx]] <- chol2inv(beta_fit$chol_q)
    }

    beta_start <- beta_hat
  }

  grid_norm <- normalize_log_grid_density(log_tilde_varphi, varphi_grid_sorted)
  weight <- grid_norm$weight

  posterior_mean_beta <- drop(crossprod(weight, beta_modes))
  posterior_mean_varphi <- sum(varphi_grid_sorted * weight)
  posterior_mean_phi <- sum(exp(varphi_grid_sorted) * weight)

  posterior_cov_beta <- NULL
  posterior_sd_beta <- NULL

  if (return_covariance) {
    second_moment <- matrix(0, nrow = p, ncol = p)

    for (idx in seq_len(k)) {
      beta_k <- beta_modes[idx, ]
      second_moment <- second_moment +
        weight[idx] * (covariance_list[[idx]] + tcrossprod(beta_k))
    }

    posterior_cov_beta <- second_moment - tcrossprod(posterior_mean_beta)
    posterior_sd_beta <- sqrt(pmax(diag(posterior_cov_beta), 0))
  }

  beta_marginals <- lapply(seq_len(p), function(j) {
    data.frame(
      varphi = varphi_grid_sorted,
      phi = exp(varphi_grid_sorted),
      weight = weight,
      mean = beta_modes[, j],
      sd = if (return_covariance) {
        sqrt(pmax(vapply(covariance_list, function(s) s[j, j], numeric(1)), .Machine$double.eps))
      } else {
        NA_real_
      }
    )
  })

  varphi_posterior <- data.frame(
    varphi = varphi_grid_sorted,
    phi = exp(varphi_grid_sorted),
    delta = grid_norm$delta,
    log_unnormalized = log_tilde_varphi,
    density = grid_norm$density,
    weight = weight,
    beta_converged = convergence,
    beta_iterations = iterations,
    q_jitter = q_jitter
  )

  structure(
    list(
      call = match.call(),
      n = n,
      p = p,
      y = y,
      X = X,
      hyperparameters = list(tau2 = tau2, a = a, b = b),
      varphi_grid = varphi_grid_sorted,
      phi_grid = exp(varphi_grid_sorted),
      beta_modes = beta_modes,
      q_matrices = q_list,
      covariance_matrices = covariance_list,
      varphi_posterior = varphi_posterior,
      beta_marginals = beta_marginals,
      posterior_mean_beta = posterior_mean_beta,
      posterior_mean_varphi = posterior_mean_varphi,
      posterior_mean_phi = posterior_mean_phi,
      posterior_cov_beta = posterior_cov_beta,
      posterior_sd_beta = posterior_sd_beta,
      beta_mode_converged = convergence,
      beta_mode_iterations = iterations
    ),
    class = "negative_binomial_inla_fit"
  )
}

evaluate_beta_marginal <- function(fit, j, x) {
  if (!inherits(fit, "negative_binomial_inla_fit")) {
    stop("fit must be the result of negative_binomial_inla().")
  }

  if (!is.numeric(j) || length(j) != 1L || j < 1 || j > fit$p) {
    stop("j must be a valid coefficient index.")
  }

  components <- fit$beta_marginals[[j]]
  if (anyNA(components$sd)) {
    stop("beta marginal standard deviations are unavailable because covariance matrices were not stored.")
  }

  component_sd <- pmax(components$sd, .Machine$double.eps)

  vapply(
    x,
    function(x0) sum(components$weight * dnorm(x0, mean = components$mean, sd = component_sd)),
    numeric(1)
  )
}

print.negative_binomial_inla_fit <- function(x, ...) {
  cat("Negative-binomial INLA-style approximation\n")
  cat("  n =", x$n, "\n")
  cat("  p =", x$p, "\n")
  cat("  varphi grid points =", length(x$varphi_grid), "\n")
  cat("  Posterior mean phi =", signif(x$posterior_mean_phi, 6), "\n")
  cat("  Beta MAP convergence =", sum(x$beta_mode_converged), "/", length(x$beta_mode_converged), "\n")
  invisible(x)
}

build_case_control_design <- function(group) {
  group <- as.character(group)

  if (!all(group %in% c("control", "case"))) {
    stop("group must contain only 'control' and 'case'.")
  }

  X <- cbind(
    beta_control = as.numeric(group == "control"),
    beta_case = as.numeric(group == "case")
  )

  storage.mode(X) <- "double"
  X
}

empirical_case_control_beta_init <- function(y, group, pseudocount = 0.1) {
  group <- as.character(group)

  c(
    log(mean(y[group == "control"]) + pseudocount),
    log(mean(y[group == "case"]) + pseudocount)
  )
}

phi_posterior_summary <- function(fit) {
  phi_grid <- fit$varphi_posterior$phi
  weight <- fit$varphi_posterior$weight
  phi_mean <- sum(phi_grid * weight)
  phi_var <- sum((phi_grid^2) * weight) - phi_mean^2

  c(mean = phi_mean, sd = sqrt(max(phi_var, 0)))
}

summarize_negative_binomial_inla_fit <- function(fit,
                                                 gene_index = NA_integer_,
                                                 gene_name = NA_character_,
                                                 is_de = NA_integer_,
                                                 true_beta_control = NA_real_,
                                                 true_beta_case = NA_real_,
                                                 true_phi = NA_real_) {
  beta_sd <- fit$posterior_sd_beta
  if (is.null(beta_sd)) {
    beta_sd <- rep(NA_real_, fit$p)
  }

  phi_summary <- phi_posterior_summary(fit)

  data.frame(
    gene_index = gene_index,
    gene = gene_name,
    is_de = is_de,
    posterior_mean_beta_control = fit$posterior_mean_beta[1],
    posterior_mean_beta_case = fit$posterior_mean_beta[2],
    posterior_sd_beta_control = beta_sd[1],
    posterior_sd_beta_case = beta_sd[2],
    posterior_mean_phi = phi_summary["mean"],
    posterior_sd_phi = phi_summary["sd"],
    all_beta_modes_converged = all(fit$beta_mode_converged),
    max_beta_mode_iterations = max(fit$beta_mode_iterations),
    true_beta_control = true_beta_control,
    true_beta_case = true_beta_case,
    true_phi = true_phi,
    stringsAsFactors = FALSE
  )
}

plot_case_control_posteriors <- function(fit,
                                         gene_label = NULL,
                                         true_beta = NULL,
                                         true_phi = NULL) {
  if (fit$p != 2L) {
    stop("plot_case_control_posteriors() expects a two-column case/control design.")
  }

  beta_names <- colnames(fit$X)
  if (is.null(beta_names)) {
    beta_names <- c("beta_control", "beta_case")
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 3))

  for (j in seq_len(2)) {
    components <- fit$beta_marginals[[j]]
    x_lower <- min(components$mean - 4 * components$sd)
    x_upper <- max(components$mean + 4 * components$sd)
    x_grid <- seq(x_lower, x_upper, length.out = 400)
    y_grid <- evaluate_beta_marginal(fit, j, x_grid)

    plot(
      x_grid, y_grid,
      type = "l",
      lwd = 2,
      col = "steelblue4",
      xlab = beta_names[j],
      ylab = "Posterior density",
      main = if (is.null(gene_label)) beta_names[j] else paste(gene_label, beta_names[j])
    )

    if (!is.null(true_beta) && length(true_beta) >= j && is.finite(true_beta[j])) {
      abline(v = true_beta[j], col = "firebrick3", lty = 2, lwd = 2)
      legend("topright", legend = "True value", col = "firebrick3", lty = 2, lwd = 2, bty = "n")
    }
  }

  phi_df <- fit$varphi_posterior
  phi_density <- phi_df$density / phi_df$phi

  plot(
    phi_df$phi, phi_density,
    type = "l",
    lwd = 2,
    col = "darkgreen",
    xlab = expression(phi == e^varphi),
    ylab = "Posterior density",
    main = if (is.null(gene_label)) expression(phi) else paste(gene_label, "phi")
  )

  if (!is.null(true_phi) && is.finite(true_phi)) {
    abline(v = true_phi, col = "firebrick3", lty = 2, lwd = 2)
    legend("topright", legend = "True value", col = "firebrick3", lty = 2, lwd = 2, bty = "n")
  }
}


# implementation ----------------------------------------------------------

simulated_data_path <- file.path(
  "../data",
  "data_simulator_independent___simulated_counts.rds"
)

if (file.exists(simulated_data_path)) {
  sim <- readRDS(simulated_data_path)

  ground_truth_I_DE <- if (!is.null(sim$ground_truth_I_DE)) sim$ground_truth_I_DE else sim$I_DE
  if (is.null(ground_truth_I_DE)) {
    stop("The simulated data object must contain ground_truth_I_DE or I_DE.")
  }

  gene_names <- colnames(sim$Y)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(ncol(sim$Y)))
  }

  X_case_control <- build_case_control_design(sim$group)

  tau2_prior <- 100
  a_prior <- if (!is.null(sim$alpha)) sim$alpha else 2
  b_prior <- if (!is.null(sim$lambda)) sim$lambda else 1
  varphi_grid <- seq(log(0.05), log(20), length.out = 121)

  de_gene_idx <- which(ground_truth_I_DE == 1L)[1]
  if (is.na(de_gene_idx)) {
    stop("No ground-truth DE gene was found in the simulated data.")
  }

  de_gene_name <- gene_names[de_gene_idx]
  de_gene_truth <- if (!is.null(sim$gene_parameters)) sim$gene_parameters[de_gene_idx, ] else NULL
  de_gene_y <- sim$Y[, de_gene_idx]

  de_gene_fit <- negative_binomial_inla(
    y = de_gene_y,
    X = X_case_control,
    tau2 = tau2_prior,
    a = a_prior,
    b = b_prior,
    varphi_grid = varphi_grid,
    beta_init = empirical_case_control_beta_init(de_gene_y, sim$group)
  )

  true_beta_de <- if (!is.null(de_gene_truth)) {
    c(
      de_gene_truth$beta0,
      de_gene_truth$beta0 + de_gene_truth$beta1 * de_gene_truth$I_DE
    )
  } else {
    NULL
  }

  true_phi_de <- if (!is.null(de_gene_truth)) de_gene_truth$phi else NULL

  plot_case_control_posteriors(
    fit = de_gene_fit,
    gene_label = de_gene_name,
    true_beta = true_beta_de,
    true_phi = true_phi_de
  )

  de_gene_summary <- summarize_negative_binomial_inla_fit(
    fit = de_gene_fit,
    gene_index = de_gene_idx,
    gene_name = de_gene_name,
    is_de = ground_truth_I_DE[de_gene_idx],
    true_beta_control = if (!is.null(de_gene_truth)) de_gene_truth$beta0 else NA_real_,
    true_beta_case = if (!is.null(de_gene_truth)) de_gene_truth$beta0 + de_gene_truth$beta1 * de_gene_truth$I_DE else NA_real_,
    true_phi = if (!is.null(de_gene_truth)) de_gene_truth$phi else NA_real_
  )

  print(de_gene_summary)

  n_genes <- ncol(sim$Y)
  all_gene_laplace_summary <- vector("list", n_genes)

  for (j in seq_len(n_genes)) {
    if (j %% 50L == 0L) {
      message("Completed ", j, " / ", n_genes, " genes")
    }

    y_j <- sim$Y[, j]
    truth_j <- if (!is.null(sim$gene_parameters)) sim$gene_parameters[j, ] else NULL

    fit_j <- negative_binomial_inla(
      y = y_j,
      X = X_case_control,
      tau2 = tau2_prior,
      a = a_prior,
      b = b_prior,
      varphi_grid = varphi_grid,
      beta_init = empirical_case_control_beta_init(y_j, sim$group)
    )

    all_gene_laplace_summary[[j]] <- summarize_negative_binomial_inla_fit(
      fit = fit_j,
      gene_index = j,
      gene_name = gene_names[j],
      is_de = ground_truth_I_DE[j],
      true_beta_control = if (!is.null(truth_j)) truth_j$beta0 else NA_real_,
      true_beta_case = if (!is.null(truth_j)) truth_j$beta0 + truth_j$beta1 * truth_j$I_DE else NA_real_,
      true_phi = if (!is.null(truth_j)) truth_j$phi else NA_real_
    )
  }

  all_gene_laplace_summary <- do.call(rbind, all_gene_laplace_summary)
  rownames(all_gene_laplace_summary) <- NULL

  head(all_gene_laplace_summary)
} else {
  message("Simulated data file not found: ", simulated_data_path)
}

