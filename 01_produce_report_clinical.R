source("src/survival_analysis.R")

# Imports ####
# import rna expression and clinical data
data_clinical <- read_csv("data/interim/clinical_data.csv")
counts_df <- read_csv("data/interim/gene_counts.csv")


# filter clinical data
data_clinical %>% filter(RNA_seq == "Yes")

data_clinical <- data_clinical %>% filter(
  data_clinical$patient_id %in% counts_df$patient_id
)

# Survival analysis ####
genes <- c("TERT")

# generate survival analysis results for each gene
for (gene in genes) {
  write_survival_results(gene, "OS_event_2yrs", "OS_time_2yrs",
    counts_df, data_clinical,
    path = paste0("reports/survival_analysis/univariate/")
  )
}

# generate survival plots for each gene
plots <- list()
for (gene in genes) {
  plots[[gene]] <- plot_survival(
    gene, "OS_event_2yrs", "OS_time_2yrs",
    counts_df, data_clinical
  )
}

plots$TERT <- plots$TERT +
  annotate("text", x = 21, y = 0.85, color = "#9A9A02",
           label = expression(italic("TERT")~"low")) +
  annotate("text", x = 19, y = 0.62, color = "#465FAF",
           label = expression(italic("TERT")~"high")) +
  annotate("text", x = 2, y = 0.05, color = "black",
           label = "log-rank\np=0.118") +
  theme(plot.title = element_blank(),
        plot.legend.position = "none")

plots$TERT

if(!dir.exists("reports/survival_analysis/univariate/figures")) {
  dir.create("reports/survival_analysis/univariate/figures", recursive = TRUE)
}

library(svglite)
for (gene in genes) {
  svglite(filename  = paste0("reports/survival_analysis/univariate/figures/", gene, "_survival_plot.svg"),
          width = 5, height = 4, fix_text_size = FALSE)
  print(plots[gene])
  dev.off()
}

for (gene in genes) {
  pdf(file = paste0("reports/survival_analysis/univariate/figures/", gene, "_survival_plot.pdf"),
      width = 5, height = 4)
  print(plots[gene])
  dev.off()
}


# multivariate analysis ####
data_clinical$log_bilirubin <- log2(data_clinical$TB_total_bilirubin_)

mult_surv_fit <-
  lapply(genes, function(gene) {
    multivariate_survival_analysis(
    gene = gene,
    covariates_num = c("Age"),
    covariates_cat = c("Sex",
                       "TNM_stage"),
    survival_est = "OS_event_2yrs",
    survival_time = "OS_time_2yrs",
    counts_df = counts_df,
    data_clinical = data_clinical
  )
  }) %>% setNames(genes)

lapply(genes, function(x) summary(mult_surv_fit[[x]]))

# save multivariate analysis results
for (gene in genes) {
  write_multivariate_survival_results(
    mult_surv_fit[[gene]],
    gene,
    path = "reports/survival_analysis/multivariate/"
  )
}
