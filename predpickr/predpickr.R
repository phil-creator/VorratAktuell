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




# Kapitel 0: Laden der benötigten Pakete
requiredPackages <- c("randomForest", "Boruta", "car", "Hmisc", "utils")
for(pkg in requiredPackages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# Setze einen zufälligen Seed für Reproduzierbarkeit (z.B. Random Forest, Boruta)
set.seed(123)

# Kapitel 1: Einlesen und Bereinigung des Datensatzes
# Pfad zum Datensatz anpassen!
data_raw <- read.csv("path/to/dataset.csv", stringsAsFactors = FALSE)
# Beispiel: Zielvariable festlegen (bitte anpassen)
response_var <- "Zielvariable"
# In Faktor umwandeln, falls notwendig (bei Klassifikation)
# data_raw[[response_var]] <- as.factor(data_raw[[response_var]]) 

# Fehlende Werte entfernen
data_clean <- na.omit(data_raw)

# Spaltennamen auf gültige R-Namen prüfen
names(data_clean) <- make.names(names(data_clean))
# Prädiktorvariablen (ohne Zielvariable)
predictor_vars <- setdiff(names(data_clean), response_var)


# Kapitel 2: Berechnung des kombinierten Variablenrankings

# (a) AIC-basiertes Ranking durch drop1
lm_full <- lm(as.formula(paste(response_var, "~ .")), data = data_clean)
drop1_res <- drop1(lm_full, test = "Chisq")  # Standardtest AIC
# AIC-Werte für Modelle ohne jeweilige Variable
full_AIC <- AIC(lm_full)
# Drop1-Ergebnis enthält Zeile "<none>" für das volle Modell; diese entfernen
drop1_table <- as.data.frame(drop1_res)[-1, , drop = FALSE]
# Wichtigkeit definieren als AIC(ohne Var) - AIC(full); höhere Werte = wichtiger
aic_diff <- drop1_table$AIC - full_AIC
names(aic_diff) <- rownames(drop1_table)
aic_rank <- rank(-aic_diff, ties.method = "min")  # Negativ, damit größte Diff Rang 1

# (b) Random Forest-Importance
rf_model <- randomForest(as.formula(paste(response_var, "~ .")), data = data_clean, importance = TRUE)
rf_imp <- importance(rf_model)
# Extrahiere Gini (IncNodePurity) und Permutation (%IncMSE)
gini_imp <- rf_imp[, "IncNodePurity"]
perm_imp <- rf_imp[, "%IncMSE"]
names(gini_imp) <- rownames(rf_imp)
names(perm_imp) <- rownames(rf_imp)
gini_rank <- rank(-gini_imp, ties.method = "min")
perm_rank <- rank(-perm_imp, ties.method = "min")

# (c) Boruta-Feature-Selection
boruta_model <- Boruta(as.formula(paste(response_var, "~ .")), data = data_clean, maxRuns = 100, doTrace = 0)
boruta_stats <- attStats(boruta_model)
# Median-Importance extrahieren
boruta_imp <- boruta_stats$medianImp
names(boruta_imp) <- rownames(boruta_stats)
boruta_rank <- rank(-boruta_imp, ties.method = "min")

# (d) Gesamtranking zusammenführen
all_vars <- sort(unique(c(names(aic_diff), names(gini_imp), names(boruta_imp))))
ranking_df <- data.frame(Variable = all_vars,
                         AIC_diff = as.numeric(aic_diff[all_vars]),
                         Rank_AIC = as.numeric(aic_rank[all_vars]),
                         RF_Gini = as.numeric(gini_imp[all_vars]),
                         Rank_RF_Gini = as.numeric(gini_rank[all_vars]),
                         RF_Perm = as.numeric(perm_imp[all_vars]),
                         Rank_RF_Perm = as.numeric(perm_rank[all_vars]),
                         Boruta_imp = as.numeric(boruta_imp[all_vars]),
                         Rank_Boruta = as.numeric(boruta_rank[all_vars]),
                         stringsAsFactors = FALSE)
# Fehlende Werte aus Boruta, RF (falls variable z.B. konstant war) auf NA setzen
ranking_df[is.na(ranking_df)] <- 0

# Gesamt-Rang (z.B. durchschnittlicher Rang)
ranking_df$MeanRank <- rowMeans(ranking_df[, c("Rank_AIC","Rank_RF_Gini","Rank_RF_Perm","Rank_Boruta")])
ranking_df <- ranking_df[order(ranking_df$MeanRank), ]
# CSV-Ausgabe des Rankings
write.csv(ranking_df, "output/variable_ranking.csv", row.names = FALSE)


# Kapitel 3: Funktion zur Multikollinearitätsbehandlung

handle_multicollinearity <- function(data, response, ranking_df, method = c("tolerant","strict","very_strict")) {
  method <- match.arg(method)
  # Schwellenwerte je nach Methode festlegen
  thresholds <- list(
    tolerant    = list(cor = 0.8, vif = 10, varclus = 0.75),
    strict      = list(cor = 0.7, vif = 5,  varclus = 0.5),
    very_strict = list(cor = 0.5, vif = 3,  varclus = 0.3)
  )
  th <- thresholds[[method]]
  
  vars <- setdiff(names(data), response)
  ranks <- setNames(ranking_df$MeanRank, ranking_df$Variable)
  # Stelle sicher, dass alle zu prüfenden Variablen einen Rang haben
  ranks <- ranks[intersect(names(ranks), vars)]
  
  # Schritt 1: Paarweise Korrelation
  cor_matrix <- abs(cor(data[, vars], use = "pairwise.complete.obs", method = "spearman"))
  diag(cor_matrix) <- 0
  repeat {
    # Finde alle Paare über Schwellenwert
    corr_pairs <- which(cor_matrix > th$cor, arr.ind = TRUE)
    if(nrow(corr_pairs) == 0) break
    # Betrachte nur i<j (Doppelzählung vermeiden)
    corr_pairs <- corr_pairs[corr_pairs[,1] < corr_pairs[,2], ]
    if(nrow(corr_pairs) == 0) break
    # Nehme Paar mit maximaler Korrelation
    max_pair <- corr_pairs[which.max(cor_matrix[corr_pairs]), ]
    var1 <- colnames(cor_matrix)[max_pair[2]]
    var2 <- rownames(cor_matrix)[max_pair[1]]
    # Wähle die weniger wichtige Variable (höherer Rang-Wert)
    drop_var <- if(ranks[var1] > ranks[var2]) var1 else var2
    # Entferne Variable
    vars <- setdiff(vars, drop_var)
    cor_matrix <- cor_matrix[vars, vars, drop = FALSE]
  }
  
  # Schritt 2: VIF-Schritt
  repeat {
    if(length(vars) == 0) break
    lm_temp <- lm(as.formula(paste(response, "~", paste(vars, collapse = "+"))), data = data)
    vif_vals <- try(car::vif(lm_temp), silent = TRUE)
    if(inherits(vif_vals, "try-error")) break
    vif_vals <- vif_vals[!is.na(vif_vals)]
    if(max(vif_vals, na.rm = TRUE) <= th$vif) break
    # Variable mit größtem VIF entfernen (wenn gleich, niedrigere Rang entfernen)
    max_vif_var <- names(which.max(vif_vals))
    # Falls mehrere, wähle nach Rang
    if(length(max_vif_var) > 1) {
      max_vif_var <- max_vif_var[which.max(sapply(max_vif_var, function(x) -ranks[x]))]
    }
    vars <- setdiff(vars, max_vif_var)
  }
  
  # Schritt 3: Varclus (Cluster-Analyse)
  if(length(vars) > 1) {
    # Hmisc::varclus erfordert keine NA
    vc <- Hmisc::varclus(as.matrix(data[, vars]), similarity = "spearman")
    # Cut the tree at height = 1 - threshold (Spearman Ähnlichkeit)
    hc <- vc$hclust
    # Höhe schneiden (0 = perfekt korreliert, 1 = unkorreliert)
    cluster_assign <- cutree(hc, h = 1 - th$varclus)
    # In jedem Cluster nur die wichtigste Variable behalten
    to_keep <- tapply(vars, cluster_assign, function(cluster_vars) {
      if(length(cluster_vars) == 1) return(cluster_vars)
      # Wähle Variable mit dem niedrigsten Rang (wichtigste)
      cluster_vars[which.min(ranks[cluster_vars])]
    })
    vars <- unname(to_keep)
  }
  
  return(vars)
}


# Kapitel 4: Anwenden der Multikollinearitäts-Funktion
selected_tolerant    <- handle_multicollinearity(data_clean, response_var, ranking_df, method = "tolerant")
selected_strict      <- handle_multicollinearity(data_clean, response_var, ranking_df, method = "strict")
selected_very_strict <- handle_multicollinearity(data_clean, response_var, ranking_df, method = "very_strict")

# Ausgabe der ausgewählten Variablen (optional)
print(selected_tolerant)
print(selected_strict)
print(selected_very_strict)



# Kapitel 5: Export der finalen Datensätze
final_tolerant <- data_clean[, c(response_var, selected_tolerant)]
final_strict   <- data_clean[, c(response_var, selected_strict)]
final_very     <- data_clean[, c(response_var, selected_very_strict)]

# Pfade anpassen!
write.csv(final_tolerant,    "output/final_dataset_tolerant.csv",    row.names = FALSE)
write.csv(final_strict,      "output/final_dataset_strict.csv",      row.names = FALSE)
write.csv(final_very,        "output/final_dataset_very_strict.csv", row.names = FALSE)

