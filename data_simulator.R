simulateData = function(b0, b1, phi, n_control, n_case, n_genes, de_idx){
  
  Y = list()
  for (i in 1:n_genes){
    control_obs = rnbinom(n_control, mu = exp(b0[i]), size = phi[i])
    case_obs = rnbinom(n_case, mu = exp(b0[i] + b1[i]), size = phi[i])
    Y[[i]] = c(control_obs, case_obs)
  }
  Y = do.call(what = cbind, args = Y)
  colnames(Y) = paste0('gene', 1:ncol(Y))
  rownames(Y) = c(paste0('control', 1:n_control), paste0('case', 1:n_case))
  
  params = data.frame(
    de_idx = de_idx, 
    b0 = b0,
    b1 = b1
  )
  
  sim = list(counts = Y, params = params)
}



B0 = rnorm(n_genes)
B1 = rnorm(n_genes) * de_idx
PHI = rgamma(n_genes, shape = 1, rate = 2)
N_CONTROL = 50
N_CASE = 50
N_GENES = 1000
DE_IDX = rep(0, 1000)
DE_IDX[sample(1:1000, 100)] = 1

for (k in 1:100){
  message(k, '/', 100)
  sim = simulateData(b0 = B0, b1 = B1, phi = PHI, 
                     n_control = N_CONTROL, n_case = N_CASE, 
                     n_genes = N_GENES, de_idx = DE_IDX)
  setwd('C:/Users/dcginger/Desktop/st502/project/data/')
  saveRDS(sim, paste0('data_simulator___sim', k, '.rds'))
}



