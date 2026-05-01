# Open up the 100 simulated datasets
#each has 100 gene columns
#each has 50 control sample rows, and 50 case sample rows

library(pROC)
library(PRROC)

# X and Y already loaded from the Frequentist Gaussian run

# I'm gonna try running!
#------------------------------------------------------------
all_freq_nb_results <- vector("list", length(sim_list))
all_freq_nb_evals   <- vector("list", length(sim_list))

for (i in seq_along(sim_list)) {
  message(i, "/", length(sim_list))
  sim <- sim_list[[i]]
  X <- X
  Y <- as.matrix(sim$counts)
  mode(Y) <- "numeric"

  # Name true_de so evaluate_freq_nb can index by gene name
  true_de <- sim$params$de_idx
  names(true_de) <- paste0("gene", seq_len(length(true_de)))

  # MODEL FUNCTION
  res <- all_nb_genes(X, Y)
  # Add simulation ID
  res$sim_id <- i
  all_freq_nb_results[[i]] <- res

  # EVALUATION FUNCTION
  eval_res <- evaluate_freq_nb(res, true_de)
  eval_res$sim_id <- i
  all_freq_nb_evals[[i]] <- eval_res
}

# Combine into one big dataframe
final_freq_nb_df <- dplyr::bind_rows(all_freq_nb_results)

# Evaluation metrics (one row per simulation)
final_freq_nb_evals <- dplyr::bind_rows(
  lapply(all_freq_nb_evals, function(e) {
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
mean_fpr_nb <- seq(0, 1, length.out = 100) #false positive rate
mean_recall_nb <- seq(0, 1, length.out = 100)
tprs_nb <- vector("list", 100) #true positive rate
prcs_nb <- vector("list", 100)

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

