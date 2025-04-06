library(tidyverse)
library(dplyr)
library(openxlsx)
library(cutpointr)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(survival)
library(survminer)

# define function for dichotomizing continuous variables
dichotomise <- function(x, cutpoint) {
  case_when(
    x >= cutpoint ~ 1,
    x < cutpoint ~ 0
  )
}

calculate_cutpoint <- function(gene = NULL, survival_est = NULL,
                               counts_df = NULL,
                               data_clinical = NULL) {
  # first we calculate the youden indices for all possible cutpoints
  youden_gene <- lapply(counts_df[[gene]], function(x) {
    prediction <- dichotomise(counts_df[[gene]], x)
    freq <- table(prediction, data_clinical[[survival_est]]) %>%
      as.data.frame() %>%
      pull("Freq")
    youden(freq[4], freq[2], freq[1], freq[3])
  }) %>%
    unlist()

  # then we select the cutpoint that maximizes the youden index
  cutpoint_gene <- counts_df[[gene]][which.max(youden_gene)]
  return(cutpoint_gene)
}

# a function to add the dichotomized variable, based on the youden_cutpoint,
# to the clinical data
add_dichotomized <- function(gene = NULL, survival_est = NULL,
                             counts_df = NULL,
                             data_clinical = NULL) {
  cutpoint <- calculate_cutpoint(
    gene, survival_est,
    counts_df, data_clinical
  )

  data_clinical[[gene]] <- dichotomise(counts_df[[gene]], cutpoint)
  return(data_clinical)
}

# functions for survival analysis ####

# helper function for plotting survival curves
plot_survival <- function(gene = NULL,
                          survival_est = NULL,
                          survival_time = NULL,
                          counts_df = NULL,
                          data_clinical = NULL) {
  data_df <- add_dichotomized(gene, survival_est, counts_df, data_clinical)

  km_fit <- survfit(as.formula(
    paste0("Surv(", survival_time, " , ", survival_est, ") ~ ", gene)
  ), data = data_df)


  autoplot(km_fit, conf.int = FALSE, censor.shape = '|', censor.size = 2) +
    ggtitle(gene) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    labs(x = "OS (months)", y = "Cumulative Survival") +
    scale_color_discrete(labels = c("low", "high"), type = c("#9A9A02", "#465FAF")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 8.8, color = "black"),
          axis.ticks = element_line(size = 0.5, color = "black")
          )
}

# helper function to generate statistics from survival analysis
survival_stats <- function(gene = NULL,
                           survival_est = NULL,
                           survival_time = NULL,
                           counts_df = NULL,
                           data_clinical = NULL) {
  data_df <- add_dichotomized(gene, survival_est, counts_df, data_clinical)

  km_fit <- survfit(as.formula(
    paste0("Surv(", survival_time, " , ", survival_est, ") ~ ", gene)
  ), data = data_df)

  surv_table <- broom::tidy(km_fit)

  sum_survival <- summary(km_fit)$table %>%
    as.data.frame()

  surv_stats <- survdiff(
    as.formula(
      paste0("Surv(", survival_time, " , ", survival_est, ") ~ ", gene)
    ),
    data = data_df
  )

  return(list("surv_table" = surv_table,
              "sum_survival" = sum_survival,
              "logrank_stats" = broom::tidy(surv_stats),
              "logrank_pval" = broom::glance(surv_stats)))
}

# function to save survival analysis results
write_survival_results <- function(gene = NULL,
                                   survival_est = NULL,
                                   survival_time = NULL,
                                   counts_df = NULL,
                                   data_clinical = NULL,
                                   path = NULL) {
  test_res <- survival_stats(
    gene, survival_est, survival_time,
    counts_df, data_clinical
  )

  wb <- createWorkbook()

  if(!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  addWorksheet(wb, "surv_table")
  writeData(wb, "surv_table", test_res$surv_table)

  addWorksheet(wb, "sum_survival")
  writeData(wb, "sum_survival", test_res$sum_survival)

  addWorksheet(wb, "logrank_stats")
  writeData(wb, "logrank_stats", test_res$logrank_stats)

  addWorksheet(wb, "logrank_pval")
  writeData(wb, "logrank_pval", test_res$logrank_pval)

  saveWorkbook(wb, paste0(path, gene, "_univariate_survival.xlsx"), overwrite = TRUE)
}

# function to perform univariate survival analysis to test single variables
univariate_survival_analysis <- function(variable = NULL,
                                         type = c("continuous", "categorical"),
                                         survival_est = NULL,
                                         survival_time = NULL,
                                         counts_df = NULL,
                                         data_clinical = NULL) {
  data_df <- add_dichotomized(gene, survival_est, counts_df, data_clinical)

  if (type == "categorical"){
    data_df[[variable]] <- as.factor(data_df[[variable]])
  }

  fit <- coxph(
    as.formula(
      paste0("Surv(", survival_time, " , ", survival_est, ") ~ ", variable)
    ),
    data = data_df
  )

  return(summary(fit))
}

# function to perform multivariate survival analysis
multivariate_survival_analysis <- function(gene = NULL,
                                           covariates_num = NULL,
                                           covariates_cat = NULL,
                                           survival_est = NULL,
                                           survival_time = NULL,
                                           counts_df = NULL,
                                           data_clinical = NULL) {
  data_df <- data_clinical

  data_df <- add_dichotomized(gene, survival_est, counts_df, data_df)

  for (covariate in covariates_cat) {
    data_df[[covariate]] <- as.factor(data_df[[covariate]])
  }

  fit <- coxph(
    as.formula(
      paste0("Surv(", survival_time, " , ", survival_est, ") ~ ",
             paste(gene, collapse = " + "), " + ", 
             paste(covariates_num, collapse = " + "),
             " + ",
             paste(covariates_cat, collapse = " + "))
    ),
    data = data_df
  )

  return(fit)
}

# save multivariate survival analysis results
write_multivariate_survival_results <- function(coxph_model = NULL,
                                                gene = NULL,
                                                path = NULL) {
  wb <- createWorkbook()

  addWorksheet(wb, "coxph_model")
  writeData(wb, "coxph_model", summary(coxph_model))

  addWorksheet(wb, "residuals")
  writeData(wb, "residuals", broom::augment(coxph_model))

  addWorksheet(wb, "statistics")
  writeData(wb, "statistics", broom::glance(coxph_model))

  if(!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  saveWorkbook(wb, paste0(path, gene, "multivariate_survival.xlsx"), overwrite = TRUE)

}

