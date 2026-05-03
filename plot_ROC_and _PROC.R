# Plot ROC curves (true positive on y-axis vs. false positive on x-axis)
# Plot PROC curves (precision on y-axis vs. recall on x-axis)
#"Precision measures accuracy—out of all items labeled "positive," how many were actually correct.
# Recall measures completeness—out of all actual positive items, how many did the model find"

gaus_de_p <- final_freq_gaussian_df %>%
  filter(P_Value <0.05)

gaus_de_adj_p <- final_freq_gaussian_df %>%
  filter(FDR_Adj_P <0.05)

nb_de_p <- final_freq_nb_df %>%
  filter(P_Value <0.05)

nb_de_adj_p <- final_freq_nb_df %>%
  filter(FDR_Adj_P <0.05)

# From Jonathan: Bayesian NB
# ----------------------------------------------------
# ----------------------------------------------------
#bayes_nb_means <- read.table("/Users/kateeverly/Documents/GradGeneral/ST502/bayes_nb_input_for_roc.txt")
#bayes_nb_all_data <- read.table("/Users/kateeverly/Documents/GradGeneral/ST502/z_thresholding_for_auc_curve.txt")
bayes_nb_metadata <- read.table("/Users/kateeverly/Documents/GradGeneral/ST502/iteration_metadata.txt", header=TRUE)

bayes_nb_roc_data <- read.table("/Users/kateeverly/Documents/GradGeneral/ST502/roc_input_negative_binomial_bayesian.txt")
bayes_nb_roc_obj <- roc(response = bayes_nb_roc_data$True, predictor = bayes_nb_roc_data$Predicted)
bayes_nb_roc_df <- data.frame(
  fpr = 1 - bayes_nb_roc_obj$specificities,
  tpr = bayes_nb_roc_obj$sensitivities
)

bayes_nb_pr <- pr.curve(
  scores.class0 = bayes_nb_roc_data$Predicted[bayes_nb_roc_data$True == 1],  # positives
  scores.class1 = bayes_nb_roc_data$Predicted[bayes_nb_roc_data$True == 0],  # negatives
  curve = TRUE
)

bayes_nb_pr_df <- data.frame(
  recall = bayes_nb_pr$curve[,1],
  precision = bayes_nb_pr$curve[,2]
)

# From Dan: Bayesian NB (INLA)
# ---------------------------------------------
bayes_nb_inla_roc_data <- readRDS("/Users/kateeverly/Documents/GradGeneral/ST502/auc___auroc_nb_inla.rds")
bayes_nb_inla_pr_data <- readRDS("/Users/kateeverly/Documents/GradGeneral/ST502/auc___auprc_nb_inla.rds")

b_nb_inla_tpr <- bayes_nb_inla_roc_data [["sensitivities"]]
b_nb_inla_fpr <- 1 - bayes_nb_inla_roc_data [["specificities"]]
B_NB_INLA_df <- data.frame(fpr = b_nb_inla_fpr, tpr = b_nb_inla_tpr )

B_NB_INLA_PR_df <- data.frame(
  recall = bayes_nb_inla_pr_data$curve[, 1],
  precision = bayes_nb_inla_pr_data$curve[, 2],
  threshold = bayes_nb_inla_pr_data$curve[, 3]
)

mean_recall <- seq(0, 1, length.out = 100)
# sort by recall first (CRITICAL)
B_NB_INLA_PR_df <- B_NB_INLA_PR_df[order(B_NB_INLA_PR_df$recall), ]
interp_precision <- approx(
  x = B_NB_INLA_PR_df$recall,
  y = B_NB_INLA_PR_df$precision,
  xout = mean_recall,
  ties = "ordered"
)$y
B_NB_INLA_PR_df_smooth <- data.frame(
  recall = mean_recall,
  precision = interp_precision
)

# From Dan: Bayesian Gaussian
# ----------------------------------------------------
# ----------------------------------------------------
bayesian_gaussian_auroc <- readRDS("/Users/kateeverly/Documents/GradGeneral/ST502/benchmarking___auroc_bayesian_gaussian.rds")
bayesian_gaussian_pr <- readRDS("/Users/kateeverly/Documents/GradGeneral/ST502/benchmarking___auprc_bayesian_gaussian.rds")
bayesian_gaussian_results <- readRDS("/Users/kateeverly/Documents/GradGeneral/ST502/benchmarking___results_bayesian_gaussian.rds")

bayes_gaus_tpr <- bayesian_gaussian_auroc[["sensitivities"]]
bayes_gaus_fpr <- 1 - bayesian_gaussian_auroc[["specificities"]]
bayes_gaus_roc_df <- data.frame(fpr = bayes_gaus_fpr, tpr = bayes_gaus_tpr)

bayes_gaus_pr_df <- data.frame(
  recall = bayesian_gaussian_pr$curve[, 1],
  precision = bayesian_gaussian_pr$curve[, 2],
  threshold = bayesian_gaussian_pr$curve[, 3]
)

# PLOT CURVES
# ----------------------------------------------------
library(ggplot2)

# Build dataframes for plotting
roc_df <- rbind (
  data.frame(
    fpr  = mean_fpr_gaussian,
    tpr  = mean_tpr_gaussian,
    model = "Frequentist\nGaussian"
  ),
  data.frame(
    fpr  = mean_fpr_nb,
    tpr  = mean_tpr_nb,
    model = "Frequentist Negative\nBinomial"
  ),
  data.frame(
    fpr  = bayes_nb_roc_df$fpr,
    tpr  = bayes_nb_roc_df$tpr,
    model = "Bayesian Negative\nBinomial (Metro)"
  ),
  data.frame(
    fpr  = bayes_gaus_roc_df$fpr,
    tpr  = bayes_gaus_roc_df$tpr,
    model = "Bayesian\nGaussian"
  ),
  data.frame(
    fpr  = B_NB_INLA_df$fpr,
    tpr  = B_NB_INLA_df$tpr,
    model = "Bayesian Negative\nBinomial (Laplace)"
  )
)

# ROC curve
ggplot(roc_df, aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray") +
  labs(
    title = "Mean ROC Curve",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Model"
  ) +
  #scale_color_manual(values = c(
  #  "Frequentist Gaussian" = "#4ecf70",
  #  "Frequentist Negative\nBinomial" = "#7d44be",
  #  "Bayesian Gaussian" = "#f55f74",
  #  "Bayesian Negative\nBinomial" = "#009cae"
  #)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, hjust=0.5),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  )


prc_df <- rbind(
  data.frame(
    recall = mean_recall_gaussian,
    precision = mean_prc_gaussian,
    model = "Frequentist\nGaussian"
  ),
  data.frame(
    recall = mean_recall_nb,
    precision = mean_prc_nb,
    model = "Frequentist Negative\nBinomial"
  ),
  data.frame(
    recall = bayes_nb_pr_df$recall,
    precision = bayes_nb_pr_df$precision,
    model = "Bayesian Negative\nBinomial (Metro)"
  ),
  data.frame(
    recall = bayes_gaus_pr_df$recall,
    precision = bayes_gaus_pr_df$precision,
    model = "Bayesian\nGaussian"
  ),
  data.frame(
    recall = B_NB_INLA_PR_df_smooth$recall,
    precision = B_NB_INLA_PR_df_smooth$precision,
    model = "Bayesian Negative\nBinomial (Laplace)"
  )
)

# PRC curve
ggplot(prc_df, aes(x = recall, y = precision, color = model)) +
  #geom_line(linewidth = 1) +
  geom_smooth(se=FALSE, method="loess", span=0.6, linewidth=1) +
  labs(
    title = "Mean Precision-Recall Curve",
    x = "Recall",
    y = "Precision",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, hjust=0.5),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  )
