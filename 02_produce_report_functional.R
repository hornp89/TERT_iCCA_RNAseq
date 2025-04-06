library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)

# generate gene set for enrichment maps ####
x <- org.Hs.egGO
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
if (length(xx) > 0) {
  # Try the first one
  got <- xx[[1]]
  got[[1]][["GOID"]]
  got[[1]][["Ontology"]]
  got[[1]][["Evidence"]]
}

xx <- as.list(org.Hs.egGO2ALLEGS)
if (length(xx) > 0) {
  # Gets the Entrez Gene identifiers for the top 2nd and 3nd GO identifiers
  goids <- xx[2:3]
  # Gets all the Entrez Gene identifiers for the first element of goids
  goids[[1]]
  # Evidence code for the mappings
  names(goids[[1]])
}


library(pathwayPCA)
library(GO.db)
gmt <- list(
  pathways = xx,
  TERMS = names(xx),
  description = names(xx)
)

GO <- as.list(GOTERM)

gmt$description <- mapIds(GO.db,
  keys = gmt$TERMS,
  column = "TERM",
  keytype = "GOID"
)

dir.create("reports/coexpression")
dir.create("reports/coexpression/enrichment")

write_gmt(gmt, file = "reports/coexpression/enrichment/GO_BP.gmt")

# generate ranked gene lists per gene ####
source("src/coexpression_analysis.R")

counts_df <- read_csv("data/interim/gene_counts.csv")

genes <- c("TERT")

ranks_df <- list()
for (gene in genes) {
  ranks_df[[gene]] <- create_ranks_df(gene_name = gene, counts_df = counts_df)
  write_csv(
    ranks_df[[gene]] %>% as.data.frame(),
    paste0("reports/coexpression/", gene, "_ranks.csv")
  )
}

# perform enrichment analysis ####
gse_results <- list()
for (gene in genes) {
  ranks <- ranks_df[[gene]] %>%
    pull(rank) %>%
    setNames(ranks_df[[gene]]$ENTREZID)
  ranks <- ranks[!is.na(ranks)]

  gse_results[[gene]] <- gseGO(ranks,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  ) %>%
    setReadable(org.Hs.eg.db, keyType = "ENTREZID") %>%
    simplify()

  # save results
  write_csv(
    gse_results[[gene]] %>% as.data.frame(),
    paste0("reports/coexpression/enrichment/", gene, "_enrichment.csv")
  )
}

# save enrichement analysis results
for (gene in genes) {
  write_enrichments_cytoscape(gse_results[[gene]], gene,
                              path = "reports/coexpression/enrichment/")
}
