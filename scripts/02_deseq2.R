
# 02_deseq2.R
# Purpose: Differential expression analysis - ICU vs non-ICU COVID-19 patients
# Author: Elvis Kwabena Asamani
# Date: June 2026

library(DESeq2)
library(ggplot2)

# ── 1. Load cleaned data ──────────────────────────────────────────────────────
counts <- readRDS("data/counts_clean.rds")
metadata <- readRDS("data/metadata_clean.rds")

# ── 2. Build DESeq2 object ────────────────────────────────────────────────────
counts_int <- round(counts)
condition <- factor(metadata$`icu:ch1`, levels = c("no", "yes"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData = data.frame(condition = condition, row.names = colnames(counts_int)),
  design = ~ condition
)

# ── 3. Filter low-count genes ─────────────────────────────────────────────────
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "
")

# ── 4. Run DESeq2 ─────────────────────────────────────────────────────────────
dds <- DESeq(dds)

# ── 5. Extract results ────────────────────────────────────────────────────────
res <- results(dds, contrast = c("condition", "yes", "no"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Significant DEGs: padj < 0.05, |LFC| > 1
res_sig_df <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
res_sig_df <- res_sig_df[order(res_sig_df$padj), ]
cat("Significant DEGs:", nrow(res_sig_df), "
")

# ── 6. Save results ───────────────────────────────────────────────────────────
write.csv(res_df, "results/deseq2_all_results.csv", row.names = FALSE)
write.csv(res_sig_df, "results/deseq2_significant_DEGs.csv", row.names = FALSE)

# ── 7. Volcano plot ───────────────────────────────────────────────────────────
res_df$significance <- "Not significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up in ICU"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down in ICU"

png("results/volcano_plot.png", width = 800, height = 600)
ggplot(res_df[!is.na(res_df$padj), ],
       aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Up in ICU" = "red",
                                "Down in ICU" = "blue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "COVID-19 ICU vs Non-ICU: Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Significance") +
  theme_classic()
dev.off()

cat("Done.
")

