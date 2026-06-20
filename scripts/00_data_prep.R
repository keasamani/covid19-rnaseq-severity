
# 00_data_prep.R
# Purpose: Download and prepare GSE157103 count matrix and metadata
# Author: Elvis Kwabena Asamani
# Date: June 2026

library(GEOquery)

# ── 1. Download metadata ──────────────────────────────────────────────────────
gse <- getGEO("GSE157103", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])

# ── 2. Download count matrix ──────────────────────────────────────────────────
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157103/suppl/GSE157103_genes.ec.tsv.gz",
  destfile = "data/GSE157103_genes.ec.tsv.gz",
  method = "curl"
)

counts <- read.table(
  gzfile("data/GSE157103_genes.ec.tsv.gz"),
  header = TRUE,
  sep = "	",
  row.names = 1
)

# ── 3. Align metadata and count matrix positionally ───────────────────────────
colnames(counts) <- rownames(metadata)

# ── 4. Extract ICU labels ─────────────────────────────────────────────────────
metadata$icu_status <- metadata$`icu:ch1`

# ── 5. Save cleaned objects ───────────────────────────────────────────────────
saveRDS(counts, "data/counts_clean.rds")
saveRDS(metadata, "data/metadata_clean.rds")

cat("Done. Dimensions:", dim(counts), "
")
cat("ICU distribution:
")
print(table(metadata$icu_status))

