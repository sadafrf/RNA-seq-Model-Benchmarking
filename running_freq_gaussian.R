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
  Y <- as.matrix(sim$counts)
  mode(Y) <- "numeric"

  # Name true_de so evaluate_gaussian can index by gene name
  true_de <- sim$params$de_idx
  names(true_de) <- paste0("gene", seq_len(length(true_de)))

  # --- Model ---
  res <- all_gaussian_genes(X, Y)
  res$sim_id <- i

  # --- Attach true beta1 ---
  true_beta1 <- sim$params$b1
  names(true_beta1) <- paste0("gene", seq_len(length(true_beta1)))
  res$true_beta1 <- true_beta1[res$Gene_id]

  # --- DE evaluation ---
  eval_row <- evaluate_freq_gaussian(res, true_de)
  eval_row$sim_id <- i

  all_freq_gaussian_results[[i]] <- res
  all_freq_gaussian_evals[[i]]   <- eval_row
}

# Combine into one big dataframe
final_freq_gaussian_df <- dplyr::bind_rows(all_freq_gaussian_results)
final_freq_gaussian_eval_df <- dplyr::bind_rows(all_freq_gaussian_evals)


final_freq_gaussian_evals <- final_freq_gaussian_eval_df %>%
  group_by(sim_id) %>%
  summarise(
    type_i_error = first(Type_I_error),
    fdp          = first(fdp),
    true_pos     = first(true_pos),
    false_pos    = first(false_pos),
    false_neg    = first(false_neg),
    true_neg     = first(true_neg),
    precision    = first(precision),
    recall       = first(recall),
    f1           = first(f1),
    beta_bias    = first(beta_bias),
    beta_mse     = first(beta_mse),
    beta_cor     = first(beta_cor),
    beta_coverage = first(beta_coverage),
    .groups = "drop"
  )


# MAKE ROC AND PROC OBJECTS
# ----------------------------------------------------
mean_fpr_gaussian <- seq(0, 1, length.out = 100) #false positive rate
mean_recall_gaussian <- seq(0, 1, length.out = 100)
tprs_gaussian <- vector("list", 100) #true positive rate
prcs_gaussian <- vector("list", 100)

for (i in seq_along(sim_list)) {
  sim_res  <- final_freq_gaussian_df[final_freq_gaussian_df$sim_id == i, ]

  true_de <- sim_list[[i]]$params$de_idx
  names(true_de) <- paste0("gene", seq_len(length(true_de)))
  true_de_vector <- true_de[as.character(sim_res$Gene_id)]
  true_de_vector[is.na(true_de_vector)] <- 0
  true_de_vector <- as.numeric(true_de_vector)

  # clean p-values
  pvals <- sim_res$P_Value
  pvals[is.na(pvals)] <- 1
  scores <- -log10(pmax(pvals, 1e-300))
  idx_pos <- which(true_de_vector == 1)
  idx_neg <- which(true_de_vector == 0)
  if (length(idx_pos) == 0 || length(idx_neg) == 0) next

  # ROC
  roc_obj <- pROC::roc(true_de_vector, scores, quiet = TRUE)
  tprs_gaussian[[i]] <- approx(1 - roc_obj$specificities, roc_obj$sensitivities,
                               xout = mean_fpr_gaussian, rule = 2)$y
  # PRC
  prec_obj <- PRROC::pr.curve(
    scores.class0 = scores[idx_pos],
    scores.class1 = scores[idx_neg],
    curve = TRUE
  )

  prcs_gaussian[[i]] <- approx(prec_obj$curve[, 1], prec_obj$curve[, 2],
                               xout = mean_recall_gaussian, rule = 2)$y
}

mean_tpr_gaussian <- colMeans(do.call(rbind, tprs_gaussian))
mean_prc_gaussian <- colMeans(do.call(rbind, prcs_gaussian))

