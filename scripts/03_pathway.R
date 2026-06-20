
# 03_pathway.R
# Purpose: KEGG pathway enrichment analysis of DESeq2 results
# Author: Elvis Kwabena Asamani
# Date: June 2026

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ── 1. Load DESeq2 significant DEGs ──────────────────────────────────────────
res_sig_df <- read.csv("results/deseq2_significant_DEGs.csv")
cat("DEGs loaded:", nrow(res_sig_df), "
")

# ── 2. Convert gene symbols to Entrez IDs ─────────────────────────────────────
gene_entrez <- bitr(
  res_sig_df$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
cat("Genes mapped:", nrow(gene_entrez), "
")

# ── 3. KEGG enrichment ────────────────────────────────────────────────────────
kegg_res <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
cat("Enriched pathways:", nrow(kegg_res), "
")

# ── 4. Save results ───────────────────────────────────────────────────────────
write.csv(kegg_res@result, "results/kegg_enrichment_results.csv", row.names = FALSE)

# ── 5. Dotplot ────────────────────────────────────────────────────────────────
png("results/kegg_dotplot.png", width = 900, height = 600)
dotplot(kegg_res, showCategory = 13, 
        title = "KEGG Pathway Enrichment - ICU vs Non-ICU COVID-19") +
  theme(axis.text.y = element_text(size = 10))
dev.off()

cat("Done.
")

