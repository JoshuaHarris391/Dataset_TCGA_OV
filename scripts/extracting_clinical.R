# Importing clinical info
library(tidyverse)
clinical <- read.csv(file = "data/clinical/clinical.cases_selection.2021-12-02/clinical.tsv", sep = "\t", header = T)
# Converting '-- to na
clinical <- lapply(clinical, function(x){gsub(pattern = "\'--", replacement = NA, x)} ) %>% as.data.frame()

# Converting variables to numeric
integer_vars <- c("age_at_index", "days_to_death", "year_of_birth", "age_at_diagnosis", "year_of_death", 
                  "days_to_diagnosis", "days_to_last_follow_up", "year_of_diagnosis")
for (i in integer_vars) {
  clinical[, i] <- clinical[, i] %>% as.integer()
}

# Converting variables to factor
factor_vars <- c("ethnicity", "gender", "race", "figo_stage", "icd_10_code", "last_known_disease_status", "morphology",
                 "primary_diagnosis", "prior_treatment", "progression_or_recurrence", "site_of_resection_or_biopsy", "tissue_or_organ_of_origin",
                 "tumor_grade", "treatment_or_therapy", "treatment_type")
for (i in factor_vars) {
  clinical[, i] <- clinical[, i] %>% as.factor()
}

# Converting age of diagnosis to years
clinical$age_at_diagnosis <- (clinical$age_at_diagnosis/365) %>% round(digits = 2)

# Adding overall survival
clinical$os_event <- ifelse(clinical$vital_status == "Alive", 0, 
                            ifelse(clinical$vital_status == "Dead", 1, NA)) %>% as.integer()
clinical$os_time <- clinical$days_to_last_follow_up
