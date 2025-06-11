#   ---------------------------
#   
# Script name: predpickr 
# 
# Purpose of script: Picks the predictors that are relevant for modelling from a dataset.
# 
# Last changed by: Eisnecker, Philipp
# 
# Date Created: 2025-06-11
# 
# Email: eisnecker.p@gmail.com
# 
# ---------------------------



predpickr <- function(X, 
                      y, 
                      cor_thresh = 0.8,
                      vif_thresh = 5,
                      method = c("rf", "boruta", "aic"),
                      keep_top_n = 3,
                      verbose = TRUE) {
  
  # Required packages
  suppressPackageStartupMessages({
    require(randomForest)
    require(Boruta)
    require(leaps)
    require(Hmisc)
    require(car)
    require(dplyr)
  })
  
  # Match method
  method <- match.arg(method)
  
  if (verbose) message("Step 1: Starting feature ranking...")
  
  # --- Step 1: Feature ranking ---
  if (method == "rf") {
    rf_model <- randomForest(X, y)
    importance_vals <- importance(rf_model)[, 1]
    rank_df <- data.frame(Variable = names(importance_vals), Importance = importance_vals)
    rank_df <- rank_df[order(-rank_df$Importance), ]
    
  } else if (method == "boruta") {
    boruta_model <- Boruta(X, y, doTrace = 0)
    importance_vals <- attStats(boruta_model)
    importance_vals <- importance_vals[importance_vals$decision == "Confirmed", ]
    rank_df <- data.frame(Variable = rownames(importance_vals), 
                          Importance = importance_vals$meanImp)
    rank_df <- rank_df[order(-rank_df$Importance), ]
    
  } else if (method == "aic") {
    df <- data.frame(y = y, X)
    null_model <- lm(y ~ 1, data = df)
    full_model <- lm(y ~ ., data = df)
    step_model <- step(null_model, scope = list(lower = null_model, upper = full_model),
                       direction = "forward", trace = 0)
    vars <- names(coef(step_model))[-1]
    rank_df <- data.frame(Variable = vars, Importance = 1:length(vars))
  }
  
  # Ensure only ranked variables are kept
  X_ranked <- X[, rank_df$Variable, drop = FALSE]
  
  if (verbose) message("Step 2: Correlation filtering...")
  
  # --- Step 2: Remove highly correlated variables ---
  cor_mat <- cor(X_ranked, use = "pairwise.complete.obs")
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  high_cor <- which(abs(cor_mat) > cor_thresh, arr.ind = TRUE)
  
  remove_cor <- character()
  if (nrow(high_cor) > 0) {
    for (i in seq_len(nrow(high_cor))) {
      var1 <- colnames(cor_mat)[high_cor[i, 1]]
      var2 <- colnames(cor_mat)[high_cor[i, 2]]
      imp1 <- rank_df$Importance[rank_df$Variable == var1]
      imp2 <- rank_df$Importance[rank_df$Variable == var2]
      drop_var <- ifelse(imp1 < imp2, var1, var2)
      remove_cor <- c(remove_cor, drop_var)
    }
    remove_cor <- unique(remove_cor)
  }
  X_filtered <- X_ranked[, !(colnames(X_ranked) %in% remove_cor), drop = FALSE]
  
  if (verbose) message("Step 3: Clustering via varclus...")
  
  # --- Step 3: Clustering redundant variables (Varclus) ---
  varclus_fit <- varclus(as.matrix(X_filtered))
  cluster_groups <- cutree(as.hclust(varclus_fit$hclust), h = 1 - cor_thresh)
  var_groups <- split(names(cluster_groups), cluster_groups)
  
  selected_vars <- character()
  for (grp in var_groups) {
    if (length(grp) == 1) {
      selected_vars <- c(selected_vars, grp)
    } else {
      # Select highest ranked variable per cluster
      imp <- rank_df$Importance[match(grp, rank_df$Variable)]
      selected <- grp[which.max(imp)]
      selected_vars <- c(selected_vars, selected)
    }
  }
  X_clustered <- X_filtered[, selected_vars, drop = FALSE]
  
  if (verbose) message("Step 4: VIF filtering...")
  
  # --- Step 4: Remove high VIF variables ---
  vif_filtered <- X_clustered
  repeat {
    df <- data.frame(y = y, vif_filtered)
    mod <- lm(y ~ ., data = df)
    vifs <- vif(mod)
    if (all(vifs < vif_thresh)) break
    drop_var <- names(which.max(vifs))
    vif_filtered <- vif_filtered[, !(names(vif_filtered) %in% drop_var), drop = FALSE]
  }
  X_final <- vif_filtered
  
  if (verbose) message("Step 5: Reinserting top variables if lost...")
  
  # --- Step 5: Reinsert top-N variables if lost ---
  top_vars <- head(rank_df$Variable, keep_top_n)
  lost_vars <- setdiff(top_vars, colnames(X_final))
  
  if (length(lost_vars) > 0) {
    reinserts <- X[, lost_vars, drop = FALSE]
    X_final <- cbind(X_final, reinserts)
    if (verbose) message("Reinserted: ", paste(lost_vars, collapse = ", "))
  }
  
  # Final dataset
  if (verbose) message("Selection complete: ", ncol(X_final), " features selected.")
  
  return(list(
    X_selected = X_final,
    ranking = rank_df
  ))
}

