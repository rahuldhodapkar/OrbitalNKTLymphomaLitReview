#!/usr/bin/env Rscript
## LoadData.R
# Simple pre-processing script for loading raw data from CSV file into
# cleaned R object for downstream analysis
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

library(stringr)
library(dplyr)
library(puddlr)
library(ComplexHeatmap)
library(viridis)
library(cowplot)
library(ggrepel)

######## Create output scaffolding

s <- ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)
s <- ifelse(!dir.exists("./calc/general"), dir.create("./calc/general"), FALSE)

################################################################################
## LOAD DATA
################################################################################

data <- read.csv('./data/SystematicReviewMachineReadableDatasheet.csv')

stage.numeric <- sapply(str_match(data$Dx_AnnArborStage, '^([IV]*)')[,1], function(x) {
  if (x == 'I') {
    return(1)
  } else if (x == 'II') {
    return(2)
  } else if (x == 'III') {
    return(3)
  } else if (x == 'IV') {
    return(4)
  } else{ 
    stop('Invalid Ann Arbor Stage provided in column "AnnArborStageInferred"')
  }
}, USE.NAMES=FALSE)

data$Dx_AnnArborStageNumeric <- stage.numeric
data$Dx_IsLateStage <- stage.numeric >= 3

data$status <- (!str_starts(data$Surv_TimeMonths, ">"))
data$time <- data$Surv_TimeMonths %>%
  str_replace(">", "") %>%
  as.integer()

data$Demo_IsAdvancedAge <- as.numeric(data$Demo_Age) > 65
data$Demo_IsMale <- ifelse(data$Demo_Sex == 'M', 1, 0)

chemo.tx.cols <- c('Tx_Chemotherapy_CHOP',
                   'Tx_Chemotherapy_Methotrexate',
                   'Tx_Chemotherapy_SMILE',
                   'Tx_Chemotherapy_DHAP',
                   'Tx_Chemotherapy_CEOP',
                   'Tx_Chemotherapy_DeVIC',
                   'Tx_Chemotherapy_GELOX')

for (c in chemo.tx.cols) {
  data[,c] <- ifelse(data[,c] == 'YES', 1, 0)
}

data$Tx_Chemotherapy_Non_L_Asp <- (data$Tx_Chemotherapy_CHOP 
                                   | data$Tx_Chemotherapy_CEOP 
                                   | data$Tx_Chemotherapy_DeVIC
                                   | data$Tx_Chemotherapy_DHAP)

data$Tx_Chemotherapy_L_Asp <- (data$Tx_Chemotherapy_SMILE
                               | data$Tx_Chemotherapy_GELOX
                               | data$Tx_Chemotherapy_Other == 'L_ASPARAGINASE')

data$Loc_IsRecurrence <- data$Loc_TumorPrimaryVRecurrence == 'RECURRENCE'
data$Dx_GTE_3Mo_Misdiagnosis <- suppressWarnings(as.numeric(data$Dx_Misdiagnosis_TimeToDiagnosis) > 3)

data$Dx_IHC_Ki67_GTE_80 <- as.numeric(as.numeric(data$Dx_IHC_Ki67) >= 0.8)

meta.data <- data[,c(
  'status', 'time',
  'Demo_Age',
  'Demo_Sex'
)]

saveRDS(data, file='./calc/general/data.rds')
saveRDS(meta.data, file='./calc/general/meta_data.rds')

print('All done!')

 