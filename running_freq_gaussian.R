# Open up the 100 simulated datasets
#each has 100 gene columns
#each has 50 control sample rows, and 50 case sample rows

library(pROC)
library(PRROC)

# Load in the Y
#-------------------------------------------------
in_dir <- "/Users/kateeverly/Documents/GradGeneral/ST502/ST502_Project"
sim_list <- vector("list", 100)
for (k in 1:100) {
  message(k, "/", 100)
  sim_list[[k]] <- readRDS(
    file.path(in_dir, paste0("data_simulator___sim", k, ".rds"))
  )
}

# X is the same for ALL DATASETS
#-------------------------------------------------
n_control <- 50
n_case <- 50
X <- cbind(
  intercept = rep(1, n_control + n_case),
  indicator = c(rep(0, n_control), rep(1, n_case))
)
x_y_rownames <- row.names(sim_list[[1]]$counts)
row.names(X) <- x_y_rownames

# I'm gonna try running!
#------------------------------------------------------------
all_freq_gaussian_results <- vector("list", length(sim_list))
all_freq_gaussian_evals   <- vector("list", length(sim_list))

for (i in seq_along(sim_list)) {
  message(i, "/", length(sim_list))
  sim <- sim_list[[i]]
  X <- X
  Y <- as.matrix(sim$counts)
  mode(Y) <- "numeric"

  # Name true_de so evaluate_gaussian can index by gene name
  true_de <- sim$params$de_idx
  names(true_de) <- paste0("gene", seq_len(length(true_de)))

  # MODEL FUNCTION
  res <- all_gaussian_genes(X, Y)
  # Add simulation ID
  res$sim_id <- i
  all_freq_gaussian_results[[i]] <- res

  # EVALUATION FUNCTION
  eval_res <- evaluate_gaussian(res, true_de)
  eval_res$sim_id <- i
  all_freq_gaussian_evals[[i]] <- eval_res
}

# Combine into one big dataframe
final_freq_gaussian_df <- dplyr::bind_rows(all_freq_gaussian_results)

# Evaluation metrics (one row per simulation)
final_freq_gaussian_evals <- dplyr::bind_rows(
  lapply(all_freq_gaussian_evals, function(e) {
    data.frame(
      sim_id      = e$sim_id,
      type_i_error = e$Type_I_error,
      fdp         = e$fdp,
      #auc         = e$auc,
      true_pos    = e$true_pos,
      false_pos   = e$false_pos,
      false_neg   = e$false_neg,
      true_neg    = e$true_neg,
      precision   = e$precision,
      recall      = e$recall,
      f1          = e$f1
    )
  })
)


# MAKE ROC AND PROC OBJECTS
# ----------------------------------------------------
mean_fpr <- seq(0, 1, length.out = 100) #false positive rate
tprs_gaussian <- vector("list", 100) #true positive rate
prcs_gaussian <- vector("list", 100)

for (i in seq_along(sim_list)) {
  sim_res  <- final_freq_gaussian_df[final_freq_gaussian_df$sim_id == i, ]

  true_de <- sim_list[[i]]$params$de_idx
  names(true_de) <- paste0("gene", seq_len(length(true_de)))
  true_de_vector <- true_de[sim_res$Gene_id]

  # ROC
  roc_obj <- pROC::roc(true_de_vector, sim_res$beta1, quiet = TRUE)
  tprs_gaussian[[i]] <- approx(1 - roc_obj$specificities, roc_obj$sensitivities,
                               xout = mean_fpr, rule = 2)$y
  # PRC
  prec_obj <- PRROC::pr.curve(
    scores.class0 = sim_res$beta1[true_de_vector == 1],
    scores.class1 = sim_res$beta1[true_de_vector == 0],
    curve = TRUE
  )
  prcs_gaussian[[i]] <- approx(prec_obj$curve[, 1], prec_obj$curve[, 2],
                               xout = mean_fpr, rule = 2)$y
}

mean_tpr <- colMeans(do.call(rbind, tprs_gaussian))
mean_prc <- colMeans(do.call(rbind, prcs_gaussian))

