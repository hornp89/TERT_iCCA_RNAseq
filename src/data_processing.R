library(tidyverse)
library(downloader)
library(readxl)
library(matrixStats)

# gene expression data ####
expression_mrna <-
  read_xlsx("data/external/mmc2.xlsx", sheet = 4, skip = 1) %>%
  column_to_rownames(colnames(.)[1])

# transform gene expression data
counts_df <- expression_mrna %>%
  t() %>%
  as.data.frame()

no_expressed_samples <- rowSums(expression_mrna > 0)

keep <- no_expressed_samples > 0
counts_df <- counts_df[, keep] %>% as.data.frame()

counts_df[1:10] %>% head()

counts_df <- counts_df %>% rownames_to_column("patient_id")

counts_df[1:10] %>% head()

write_csv(counts_df, "data/interim/gene_counts.csv", col_names = TRUE)

# protein expression data
expression_protein <- read_xlsx("mmc2.xlsx", sheet = 5, skip = 1) %>%
  column_to_rownames(colnames(.)[1])

# phosphosite expression data
expression_phosphosites <- read_xlsx("mmc2.xlsx", sheet = 6, skip = 1) %>%
  column_to_rownames(colnames(.)[1])

# phosphoprotein expression data
expression_phosphoprotein <- read_xlsx("mmc2.xlsx", sheet = 7, skip = 1) %>%
  column_to_rownames(colnames(.)[1])


### transform clinical data ####
data_clinical <- read_xlsx("data/external/mmc2.xlsx", sheet = 2, skip = 1)

data_clinical <- data_clinical %>%
  rename(
    OS_day = `OS, overall survival (day)`
  )

colnames(data_clinical) <- str_replace_all(colnames(data_clinical), " ", "_")
colnames(data_clinical) <- str_replace_all(colnames(data_clinical), "-", "_")
for (i in c(",", ";", ":", "(", ")", "U/L", "U/mL", "μg/L", "g/L", "μmol/L", "ng/mL")) {
  colnames(data_clinical) <- str_remove(colnames(data_clinical), fixed(i))
}
for (i in c(",", ";", ":", "(", ")", "U/L", "U/mL", "μg/L", "g/L", "μmol/L")) {
  colnames(data_clinical) <- str_remove(colnames(data_clinical), fixed(i))
}
colnames(data_clinical)

data_clinical$OS_event_2yrs <- case_when(
  data_clinical$OS_day >= 365 * 2 ~ 0,
  data_clinical$OS_day < 365 * 2 ~ data_clinical$OS_event
)

data_clinical$OS_time_2yrs <- case_when(
  data_clinical$OS_day >= 365 * 2 ~ (365 * 2) / 30,
  data_clinical$OS_day < 365 * 2 ~ data_clinical$OS_day / 30
)

data_clinical <- data_clinical %>% rename(patient_id = Patient_ID)

write_csv(data_clinical, "data/interim/clinical_data.csv", col_names = TRUE)
