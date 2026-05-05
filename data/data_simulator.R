simulateData <- function(b0, b1, phi, n_control, n_case, n_genes, de_idx) {
  
  Y <- list()
  
  for (i in seq_len(n_genes)) {
    control_obs <- rnbinom(n_control, mu = exp(b0[i]), size = phi[i])
    case_obs <- rnbinom(n_case, mu = exp(b0[i] + b1[i]), size = phi[i])
    Y[[i]] <- c(control_obs, case_obs)
  }
  
  Y <- do.call(cbind, Y)
  colnames(Y) <- paste0("gene", seq_len(ncol(Y)))
  rownames(Y) <- c(
    paste0("control", seq_len(n_control)),
    paste0("case", seq_len(n_case))
  )
  
  params <- data.frame(
    de_idx = de_idx,
    b0 = b0,
    b1 = b1,
    phi = phi
  )
  
  sim <- list(counts = Y, params = params)
  return(sim)
}


N_CONTROL <- 50
N_CASE <- 50
N_GENES <- 1000

DE_IDX <- rep(0, N_GENES)
DE_IDX[sample(seq_len(N_GENES), 100)] <- 1

B0 <- rnorm(N_GENES)
B1 <- rnorm(N_GENES) * DE_IDX
PHI <- rgamma(N_GENES, shape = 1, rate = 2)

out_dir <- "."

for (k in 1:100) {
  message(k, "/", 100)
  
  sim <- simulateData(
    b0 = B0,
    b1 = B1,
    phi = PHI,
    n_control = N_CONTROL,
    n_case = N_CASE,
    n_genes = N_GENES,
    de_idx = DE_IDX
  )
  
  saveRDS(sim, file.path(out_dir, paste0("data_simulator___sim", k, ".rds")))
}

