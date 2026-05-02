# Open up the 100 simulated datasets
#each has 100 gene columns
#each has 50 control sample rows, and 50 case sample rows

library(pROC)
library(PRROC)

# X and Y already loaded from the Frequentist Gaussian run

#------------------------------------------------------------
all_freq_nb_results <- vector("list", length(sim_list))
all_freq_nb_evals   <- vector("list", length(sim_list))

for (i in seq_along(sim_list)) {
  message(i, "/", length(sim_list))

  sim <- sim_list[[i]]
  Y <- as.matrix(sim$counts)
  mode(Y) <- "numeric"

  true_de <- sim$params$de_idx
  names(true_de) <- paste0("gene", seq_len(length(true_de)))

  # --- Model ---
  res <- all_nb_genes(X, Y)

  # --- Attach true beta1 ---
  true_beta1 <- sim$params$b1
  names(true_beta1) <- paste0("gene", seq_len(length(true_beta1)))
  res$true_beta1 <- true_beta1[res$Gene_id]

  # --- DE evaluation ---
  eval_row <- evaluate_freq_nb(res, true_de)
  eval_row$sim_id <- i

  all_freq_nb_results[[i]] <- res
  all_freq_nb_evals[[i]]   <- eval_row

}

# Combine into one big dataframe
final_freq_nb_df <- dplyr::bind_rows(all_freq_nb_results)
final_freq_nb_eval_df <- dplyr::bind_rows(all_freq_nb_evals)

final_freq_nb_evals <- final_freq_nb_eval_df %>%
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
mean_fpr_nb <- seq(0, 1, length.out = 100) #false positive rate
mean_recall_nb <- seq(0, 1, length.out = 100)
tprs_nb <- vector("list", 100) #true positive rate
prcs_nb <- vector("list", 100)
aucs_nb <- numeric(100)

for (i in seq_along(sim_list)) {
  sim_res  <- final_freq_nb_df[final_freq_nb_df$sim_id == i, ]

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
  aucs_nb[i] <- as.numeric(roc_obj$auc)            # add this line
  tprs_nb[[i]] <- approx(1 - roc_obj$specificities, roc_obj$sensitivities,
                               xout = mean_fpr_nb, rule = 2)$y
  # PRC

  prec_obj <- PRROC::pr.curve(
    scores.class0 = scores[idx_pos],
    scores.class1 = scores[idx_neg],
    curve = TRUE
  )

  prcs_nb[[i]] <- approx(prec_obj$curve[, 1],  # recall
                         prec_obj$curve[, 2],  # precision
                         xout = mean_recall_nb,
                         rule = 2
                         )$y
}

mean_tpr_nb <- colMeans(do.call(rbind, tprs_nb))
mean_prc_nb <- colMeans(do.call(rbind, prcs_nb), na.rm = TRUE)
mean_auc_nb <- mean(aucs_nb)
sd_auc_nb   <- sd(aucs_nb)

