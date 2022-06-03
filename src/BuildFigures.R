#!/usr/bin/env Rscript
#
## BuildFigures.R
#
# Generate simple figures and summary statistics for systematic review data.
#

library(ggplot2)
library(survival) # must import survival before survminer
library(survminer)
library(SurvMI)
library(stringr)
library(cowplot)
library(MASS)
library(reshape2)

######## Create output scaffolding
s <- ifelse(!dir.exists("./fig"), dir.create("./fig"), FALSE)
s <- ifelse(!dir.exists("./fig/general"), dir.create("./fig/general"), FALSE)

s <- ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)
s <- ifelse(!dir.exists("./calc/general"), dir.create("./calc/general"), FALSE)

################################################################################
## LOAD DATA
################################################################################

data <- read.csv('./data/SystematicReviewMachineReadableDatasheet.csv')

################################################################################
## CLEAN DATA
################################################################################

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

data$status <- (!str_starts(data$Surv_TimeMonths, ">"))
data$time <- data$Surv_TimeMonths %>%
  str_replace(">", "") %>%
  as.integer()


numeric.cols <- c(
  'Dx_Misdiagnosis_TimeToDiagnosis'
  , 'Tx_MaxRadiotherapyDoseGy'
  , 'Demo_Age'
  , 'Meta_ReportYear'
)

for (c in numeric.cols) {
  data[,c] <- as.numeric(data[,c])
}

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

for (c in chemo.tx.cols) {
  print(c)
  print(var(data[,c]))
}

data$Tx_Chemotherapy_Non_L_Asp <- (data$Tx_Chemotherapy_CHOP 
                                            | data$Tx_Chemotherapy_CEOP 
                                            | data$Tx_Chemotherapy_DeVIC
                                            | data$Tx_Chemotherapy_DHAP)

data$Tx_Chemotherapy_L_Asp <- (data$Tx_Chemotherapy_SMILE
                                        | data$Tx_Chemotherapy_GELOX
                                        | data$Tx_Chemotherapy_Other == 'L_ASPARAGINASE')




################################################################################
## Run Cox Proportional Hazards Regression
################################################################################

tx.surv.data <- data
tx.surv.data$L_Asp_Only <- (tx.surv.data$Tx_Chemotherapy_L_Asp
                           & !tx.surv.data$Tx_Chemotherapy_Non_L_Asp)

tx.surv.data <- subset(tx.surv.data, tx.surv.data$Tx_Chemotherapy_L_Asp
                                     | tx.surv.data$Tx_Chemotherapy_Non_L_Asp)

km.plots <- ggsurvplot(
  fit = surv_fit(Surv(time, status) ~ Tx_Surgical, data=data),
  xlab = "Months",
  xlim = c(1,23),
  break.x.by = 3,
  ylab = "Overall Survival Probability",
  risk.table=TRUE,
  conf.int=FALSE,
  pval=TRUE,
  pval.method=TRUE)

km.plots[[1]] + background_grid()



res.cox <- coxph(
  formula=Surv(time, status) ~
    (Dx_AnnArborStageNumeric >= 3) +
    (Demo_Age >= 76) +
    Demo_Sex,
  data = data,
  method = 'exact')
summary(res.cox)



res.cox <- coxph(
  formula=Surv(time, status) ~
    strata(Dx_AnnArborStageNumeric >= 3) +
    Tx_Chemotherapy_Non_L_Asp +
    Tx_Chemotherapy_L_Asp +
    Tx_Chemotherapy_Methotrexate,
  data = data,
  method = 'exact')
summary(res.cox)


res.cox <- coxph(
  formula=Surv(time, status) ~
    strata(Dx_AnnArborStageNumeric >= 3) +
    Tx_Chemotherapy +
    Tx_Radiotherapy +
    Tx_Surgical,
  data = data,
  method = 'exact')
summary(res.cox)



res.cox <- coxph(
  formula=Surv(time, status) ~
    Sx_VisionLoss	+
    Sx_VisionWorseThan20200 +
    Sx_RestrictedEOM +
    Sx_PeriorbitalSwelling +
    Sx_Chemosis +
    Sx_Proptosis +
    Sx_EyelidSwelling +
    Sx_Ptosis,
  data = data,
  method = 'exact')
summary(res.cox)


# no sx significant after stratification for ann arbor staging
res.cox <- coxph(
  formula=Surv(time, status) ~ 
    strata(Dx_AnnArborStageNumeric >= 3) +
    Sx_VisionLoss	+
    Sx_VisionWorseThan20200 +
    Sx_RestrictedEOM +
    Sx_PeriorbitalSwelling +
    Sx_Chemosis +
    Sx_Proptosis +
    Sx_EyelidSwelling +
    Sx_Ptosis,
  data = data,
  method = 'exact')
summary(res.cox)


# no locations significant after stratification for ann arbor staging
res.cox <- coxph(
  formula=Surv(time, status) ~ 
    strata(Dx_AnnArborStageNumeric >= 3) +
    Loc_TumorPrimaryVRecurrence +
    Loc_OrbitalLocationOrbital +
    Loc_OrbitalLocationLacrimalGland +
    Loc_OrbitalLocationLacrimalDrainage +
    Loc_OrbitalLocationConjunctival +
    Loc_OrbitalLocationUveal +
    Loc_OrbitalIntraocularExtension +
    Loc_OrbitalLocationEyelid +
    Loc_LocationNasosinus +
    Loc_LocationCNS,
  data = data,
  method = 'exact')
summary(res.cox)


res.cox <- coxph(
  formula=Surv(time, status) ~
    strata(Dx_AnnArborStageNumeric >= 3) +
    Demo_Age +
    Demo_Sex +
    Tx_Chemotherapy_Non_L_Asp +
    Tx_Chemotherapy_L_Asp +
    Tx_Chemotherapy_Methotrexate +
    Tx_Radiotherapy +
    Tx_Surgical +
    Sx_VisionLoss	+
    Sx_VisionWorseThan20200 +
    Sx_RestrictedEOM +
    Sx_PeriorbitalSwelling +
    Sx_Chemosis +
    Sx_Proptosis +
    Sx_EyelidSwelling +
    Sx_Ptosis +
    Loc_TumorPrimaryVRecurrence +
    Loc_OrbitalLocationOrbital +
    Loc_OrbitalLocationLacrimalGland +
    Loc_OrbitalLocationLacrimalDrainage +
    Loc_OrbitalLocationConjunctival +
    Loc_OrbitalLocationUveal +
    Loc_OrbitalIntraocularExtension +
    Loc_OrbitalLocationEyelid +
    Loc_LocationNasosinus +
    Loc_LocationCNS,
  data = data,
  method = 'exact')
summary(res.cox)

################################################################################
## Run Cox Proportional Hazards Regression with Multiple Imputation
################################################################################

#
# Cox Regression with Multiple Imputation for IHC data.
#

ihc.cols <- colnames(data)[
    str_detect(colnames(data), "^Dx_IHC")
    & !(colnames(data) %in% c('Dx_IHC_LCA', 'Dx_IHC_CD79a', 'Dx_IHC_CD30', 'Dx_IHC_CD2') )]
ihc.data <- data[,ihc.cols]
for (c in colnames(ihc.data)) {
  if (c == 'Dx_IHC_Ki67') {
    next
  }
  
  ihc.data[,c] <- factor(ihc.data[,c], levels=c('YES', 'NO'), ordered=T)
  ihc.data[,c] <- ifelse(ihc.data[,c] == 'YES', 1, 0)
}

ihc.data$time <- data$time
ihc.data$status <- data$status

dat <- mice(ihc.data)

models <- with(dat, coxph(formula=Surv(time, status) ~
                            Dx_IHC_CD3 +
                            Dx_IHC_CD4 +
                            Dx_IHC_CD5 +
                            Dx_IHC_CD7 + 
                            Dx_IHC_CD8 + 
                            Dx_IHC_CD20 +
                            Dx_IHC_CD45 +
                            Dx_IHC_CD56 + 
                            Dx_IHC_TIA1 +
                            Dx_IHC_LMP1 +
                            Dx_IHC_GrB +
                            Dx_IHC_Perforin +
                            Dx_IHC_EBER +
                            Dx_IHC_Ki67))
summary(pool(models))





res.cox <- coxph(
  formula=Surv(time, status) ~ 
    Dx_AnnArborStageNumeric +
    Demo_Age +
    Demo_Sex +
    Demo_Ethnicity +
    Sx_VisionLoss	+
    Sx_VisionWorseThan20200 +
    Sx_RestrictedEOM +
    Sx_PeriorbitalSwelling +
    Sx_Chemosis +
    Sx_Proptosis +
    Sx_EyelidSwelling +
    Sx_Ptosis +
    Loc_OrbitalLocationOrbital +
    Loc_OrbitalLocationLacrimalGland +
    Loc_OrbitalLocationLacrimalDrainage +
    Loc_OrbitalLocationConjunctival +
    Loc_OrbitalLocationUveal +
    Loc_OrbitalIntraocularExtension +
    Loc_OrbitalLocationEyelid +
    Loc_LocationNasosinus +
    Loc_LocationCNS +
    Tx_Surgical +
    Tx_Surgical_FESS +
    Tx_Surgical_Orbitotomy +
    Tx_Chemotherapy +
    Tx_Chemotherapy_CHOP +
    Tx_Chemotherapy_Methotrexate +
    Tx_Chemotherapy_SMILE +
    Tx_Chemotherapy_DHAP +
    Tx_Chemotherapy_CEOP +
    Tx_Chemotherapy_DeVIC +
    Tx_Chemotherapy_GELOX +
    Tx_Radiotherapy +
    Tx_Radiotherapy_Localized +
    Tx_MaxRadiotherapyDoseGy +
    Tx_Radiotherapy_WholeBody,
  data = data,
  method = 'exact')


print('All done!')