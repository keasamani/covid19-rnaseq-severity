
# 04_ml_model.R
# Purpose: Random Forest classifier - ICU vs non-ICU COVID-19 prediction
# Author: Elvis Kwabena Asamani
# Date: June 2026

library(randomForest)
library(caret)
library(pROC)

# ── 1. Load data ──────────────────────────────────────────────────────────────
counts <- readRDS("data/counts_clean.rds")
metadata <- readRDS("data/metadata_clean.rds")
res_sig_df <- read.csv("results/deseq2_significant_DEGs.csv")

# ── 2. Prepare feature matrix ─────────────────────────────────────────────────
top_genes <- head(res_sig_df$gene, 100)
counts_top <- counts[rownames(counts) %in% top_genes, ]
X <- as.data.frame(t(counts_top))
y <- factor(metadata$`icu:ch1`, levels = c("no", "yes"))

# ── 3. Train/test split ───────────────────────────────────────────────────────
set.seed(42)
train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]

# ── 4. Train Random Forest with 5-fold CV ─────────────────────────────────────
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

rf_model <- train(
  x = X_train,
  y = y_train,
  method = "rf",
  metric = "ROC",
  trControl = ctrl,
  ntree = 500
)

# ── 5. Evaluate on test set ───────────────────────────────────────────────────
y_pred <- predict(rf_model, X_test)
y_prob <- predict(rf_model, X_test, type = "prob")[, "yes"]

cm <- confusionMatrix(y_pred, y_test, positive = "yes")
print(cm)

roc_obj <- roc(y_test, y_prob)
cat("Test AUC:", round(auc(roc_obj), 3), "
")

# ── 6. Save plots ─────────────────────────────────────────────────────────────
png("results/roc_curve.png", width = 700, height = 600)
plot(roc_obj,
     col = "steelblue",
     lwd = 2,
     main = paste0("Random Forest ROC Curve (AUC = ", round(auc(roc_obj), 3), ")"),
     print.auc = TRUE,
     print.auc.y = 0.4)
dev.off()

imp <- varImp(rf_model)
png("results/feature_importance.png", width = 800, height = 600)
plot(imp, top = 20, main = "Top 20 Important Genes - Random Forest")
dev.off()

cat("Done.
")

