library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# create ranked gene lists per gene
create_ranks_df <- function(gene_name = NULL,
                            counts_df = NULL) {
  estimator <- counts_df[[gene_name]]

  ranks_df <-
    data.frame(
      rho = apply(counts_df, 2, function(x) {
        cor.test(x, estimator, method = "spearman")$estimate
      }),
      p.val = apply(counts_df, 2, function(x) {
        cor.test(x, estimator, method = "spearman")$p.value
      })
    ) %>%
    mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    mutate(rank = -log10(p.adj) * rho) %>%
    arrange(desc(rank))
  ranks_df <-
    ranks_df %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    mutate(
      ENTREZID = AnnotationDbi::mapIds(org.Hs.eg.db,
                                       keys = .$gene_symbol,
                                       keytype = "SYMBOL", column = "ENTREZID"),
      ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys = .$gene_symbol,
                                      keytype = "SYMBOL", column = "ENSEMBL")
    )
  ranks_df <- ranks_df[ranks_df$gene_symbol != gene_name, ]

  return(ranks_df)
}

# save enrichement analysis results
write_enrichments_cytoscape <- function(enrichment, gene_name, path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
  enrichment_df <- as.data.frame(enrichment)
  enrichment_file <- data.frame(
    GO.ID = enrichment_df$ID,
    Description = enrichment_df$Description,
    p.Val = enrichment_df$pvalue,
    FDR = enrichment_df$qvalue
  ) %>%
    mutate(Phenotype = case_when(
      enrichment_df$NES > 0 ~ "+1",
      enrichment_df$NES < 0 ~ "-1"
    ))

  write_delim(enrichment_file, paste0(path, gene_name, ".txt"), delim = "\t")
}