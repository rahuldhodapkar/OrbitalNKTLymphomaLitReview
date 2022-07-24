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
library(mice)

######## Create output scaffolding
s <- ifelse(!dir.exists("./fig"), dir.create("./fig"), FALSE)
s <- ifelse(!dir.exists("./fig/general"), dir.create("./fig/general"), FALSE)

s <- ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)
s <- ifelse(!dir.exists("./calc/general"), dir.create("./calc/general"), FALSE)

################################################################################
## Generate Screening Statistics
################################################################################

search.metadata <- read.csv2('./data/search_metadata.csv', quote='', sep=',')

search.data <- NULL
for (fn in search.metadata$filename) {
  x <- read.csv(paste0('./data/searches/', fn))
  
  if (is.null(search.data)) {
    search.data <- x
  } else {
    search.data <- rbind(search.data,x)
  }
}

dedup.data <- read.csv('./data/deduplicated_searches.csv')

print(paste0("Initial Search Results: ", dim(search.data)[[1]]))
print(paste0("Duplicate Records Removed: ", 
             nrow(search.data) - nrow(dedup.data)))

table(dedup.data$ExclusionReason)

################################################################################
## Run Cox Proportional Hazards Regression
################################################################################

data <- readRDS('./calc/general/data.rds')
meta.data <- readRDS('./calc/general/meta_data.rds')

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

data$Dx_IsLateStage <- factor(data$Dx_IsLateStage, levels=c(TRUE, FALSE), ordered=T)

# plot survival by Ann Arbor Stage
km.plots <- ggsurvplot(
  fit = surv_fit(Surv(time, status) ~ Dx_IsLateStage, data=data),
  xlab = "Months",
  xlim = c(1,23),
  break.x.by = 3,
  ylab = "Overall Survival Probability",
  risk.table=TRUE,
  conf.int=FALSE,
  pval=TRUE,
  pval.method=TRUE,
  palette='npg')

km.plots[[1]] + background_grid()
ggsave('./fig/general/ann_arbor_kaplan_meier_labeled.png', width=8, height=8)

# plot for manuscript
km.plots <- ggsurvplot(
  fit = surv_fit(Surv(time, status) ~ Dx_IsLateStage, data=data),
  xlab = "Months",
  xlim = c(1,23),
  break.x.by = 3,
  ylab = "Overall Survival Probability",
  risk.table=FALSE,
  conf.int=FALSE,
  pval=FALSE,
  pval.method=FALSE,
  palette='npg')

km.plots[[1]] + background_grid()
ggsave('./fig/general/ann_arbor_kaplan_meier_unlabeled.png', width=8, height=8)

survfit(Surv(time, status) ~ Dx_AnnArborStageNumeric >= 3, data = data) 



print('===== BEGIN Model 1 =====')

res.cox <- coxph(
  formula=Surv(time, status) ~
    (Dx_AnnArborStageNumeric >= 3) +
    (Demo_Age >= 76) +
    Demo_Sex,
  data = data,
  method = 'exact')
summary(res.cox)

print('===== End Model 1 =====')

print('===== BEGIN Model 2 =====')

res.cox <- coxph(
  formula=Surv(time, status) ~
    strata(Dx_AnnArborStageNumeric >= 3) +
    Tx_Chemotherapy_Non_L_Asp +
    Tx_Chemotherapy_L_Asp +
    Tx_Chemotherapy_Methotrexate,
  data = data,
  method = 'exact')
summary(res.cox)

print('===== End Model 2 =====')

print('===== BEGIN Model 3 =====')

res.cox <- coxph(
  formula=Surv(time, status) ~
    strata(Dx_AnnArborStageNumeric >= 3) +
    Tx_Chemotherapy +
    Tx_Radiotherapy +
    Tx_Surgical,
  data = data,
  method = 'exact')
summary(res.cox)

print('===== END Model 3 =====')

print('===== BEGIN Model 4 =====')

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

print('===== END Model 4 =====')


print('===== BEGIN Model 5 =====')

# no ptosis significantly negative association after
# stratificiation by ann arbor stage
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

print('===== END Model 5 =====')

print('===== BEGIN Model 6 =====')

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

print('===== END Model 6 =====')

print('===== BEGIN Model 7 =====')

# no ptosis + radiotherapy are significant
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
print('===== END Model 6 =====')


################################################################################
## Run Cox Proportional Hazards Regression with Multiple Imputation
################################################################################

data <- readRDS('./calc/general/data.rds')
meta.data <- readRDS('./calc/general/meta_data.rds')

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

print('===== BEGIN Pooled Models with Multiple Imputation =====')
summary(pool(models))
print('===== END Pooled Models with Multiple Imputation =====')

################################################################################
## Run Basic Statistics
################################################################################

data <- readRDS('./calc/general/data.rds')
meta.data <- readRDS('./calc/general/meta_data.rds')


survfit(Surv(time, status) ~ 1, data = data) 
mean(data$time, na.rm=T)
sd(data$time, na.rm=T)

# misdiagnosis portions
sum(data$Dx_Misdiagnosis_InitialDiagnosisCellulitis == 'YES')
sum(data$Dx_Misdiagnosis_InitialDiagnosisPseudotumor == 'YES')
sum(data$Dx_Misdiagnosis_InitialDiagnosisSinusitis == 'YES')
sum(data$Dx_Misdiagnosis_InitialDiagnosisUveitis == 'YES')
sum(data$Dx_Misdiagnosis_TimeToDiagnosis == 'Postmortem') +
  sum(as.numeric(data$Dx_Misdiagnosis_TimeToDiagnosis) > 3, na.rm=T)

# Age
sum(data$Demo_Age < 40)
sum(data$Demo_Age >= 40 & data$Demo_Age < 65 )
sum(data$Demo_Age >= 65)

mean(data$Demo_Age)
sd(data$Demo_Age)

# test if male age of onset different from female age of onset
t.test(data$Demo_Age[data$Demo_Sex == 'M'],
       data$Demo_Age[data$Demo_Sex == 'F'])

# test if number of males > number of females.
binom.test(x=sum(data$Demo_Sex == 'M'), n=nrow(data), p=.5)

# Sex
table(data$Demo_Sex)
table(data$Demo_Sex) / sum(table(data$Demo_Sex))

# Geographic Location
table(data$Meta_ReportGeographicLocation)

# AAS
table(data$Dx_AnnArborStageNumeric)
table(data$Dx_AnnArborStageNumeric) / table(data$Dx_AnnArborStageNumeric) %>% sum()

# Ki-67
sum(data$Dx_IHC_Ki67 < 0.5, na.rm=T)
sum(data$Dx_IHC_Ki67 >= 0.5 & data$Dx_IHC_Ki67 < 0.8, na.rm=T)
sum(data$Dx_IHC_Ki67 >= 0.8, na.rm=T)


# Collected Features
for (c in colnames(data)[str_detect(colnames(data), '^Sx')]) {
  print(paste0('===== BEGIN ', c, ' ====='))
  print(table(data[,c]))
  print(table(data[,c]) / sum(!is.na(data[,c])))
  print(paste0('===== END ', c, ' ====='))
}

for (c in colnames(data)[str_detect(colnames(data), '^Loc')]) {
  print(paste0('===== BEGIN ', c, ' ====='))
  print(table(data[,c]))
  print(table(data[,c]) / sum(!is.na(data[,c])))
  print(paste0('===== END ', c, ' ====='))
}


for (c in colnames(data)[str_detect(colnames(data), '^Dx')]) {
  print(paste0('===== BEGIN ', c, ' ====='))
  print(table(data[,c]))
  print(table(data[,c]) / sum(!is.na(data[,c])))
  print(paste0('===== END ', c, ' ====='))
}

ihc.cols <- colnames(data)[
  str_detect(colnames(data), "^Dx_IHC")]
ihc.data <- data[,ihc.cols]
num.nonempty.rows <- 0
for (i in 1:nrow(ihc.data)) {
  for (j in 1:ncol(ihc.data)) {
    if (! (is.na(ihc.data[i,j]) || ihc.data[i,j] == '') ) {
      num.nonempty.rows <- num.nonempty.rows + 1
      break
    }
  }
}

# number of pts with at least partial immunophenotyping
sum(rowSums(ihc.data[,1:17] != '') > 0)

for (c in ihc.cols) {
  print(paste0('===== BEGIN ', c, ' ====='))
  print(table(ihc.data[,c]))
  print(table(ihc.data[,c]) / sum(!(ihc.data[,c] == '')))
  print(96 -  sum(ihc.data[,c] == ''))
  print( (96 -  sum(ihc.data[,c] == '')) / 96 )
  print(paste0('===== END ', c, ' ====='))
}


for (c in colnames(data)[str_detect(colnames(data), '^Tx_Surgical')]) {
  print(paste0('===== BEGIN ', c, ' ====='))
  print(table(data[,c]))
  print(table(data[,c]) / sum(!is.na(data[,c])))
  print(paste0('===== END ', c, ' ====='))
}

for (c in colnames(data)[str_detect(colnames(data), '^Tx_Radiotherapy')]) {
  print(paste0('===== BEGIN ', c, ' ====='))
  print(table(data[,c]))
  print(table(data[,c]) / sum(!is.na(data[,c])))
  print(paste0('===== END ', c, ' ====='))
}

for (c in colnames(data)[str_detect(colnames(data), '^Tx_Chemotherapy')]) {
  print(paste0('===== BEGIN ', c, ' ====='))
  print(table(data[,c]))
  print(table(data[,c]) / sum(!is.na(data[,c])))
  print(paste0('===== END ', c, ' ====='))
}

sum(data$Tx_Chemotherapy_Other != 'NO')

# survival numbers

sum(data$time < 3, na.rm=T)
sum(data$time >= 3 & data$time < 6, na.rm=T)
sum(data$time >= 6 & data$time < 12, na.rm=T)
sum(data$time >= 12, na.rm=T)

sum(data$time >= 3, na.rm=T)
sum(data$time >= 6, na.rm=T)
sum(data$time >= 12, na.rm=T)

# misdiagnosis
misdiag.cols <- colnames(data)[
  str_detect(colnames(data), "^Dx_Misdiagnosis")]
misdiag.data <- data[,misdiag.cols]
num.nonempty.rows <- 0
for (i in 1:nrow(misdiag.data)) {
  for (j in 1:ncol(misdiag.data)) {
    if (! (is.na(misdiag.data[i,j]) || misdiag.data[i,j] == '') ) {
      num.nonempty.rows <- num.nonempty.rows + 1
      break
    }
  }
}

sum(data$Dx_Misdiagnosis_InitialDiagnosisCellulitis == 'YES')
sum(data$Dx_Misdiagnosis_InitialDiagnosisSinusitis == 'YES')
sum(data$Dx_Misdiagnosis_InitialDiagnosisUveitis == 'YES')
sum(data$Dx_Misdiagnosis_InitialDiagnosisPseudotumor== 'YES')

################################################################################
## CHOP vs L-Asparaginase Chemotherapy
################################################################################

data <- readRDS('./calc/general/data.rds')
meta.data <- readRDS('./calc/general/meta_data.rds')

data <- data[xor(data$Tx_Chemotherapy_CHOP, data$Tx_Chemotherapy_L_Asp),]

km.plots <- ggsurvplot(
  fit = surv_fit(Surv(time, status) ~ strata(Dx_AnnArborStageNumeric >= 3) + Tx_Chemotherapy_CHOP, data=data),
  xlab = "Months",
  xlim = c(1,23),
  break.x.by = 3,
  ylab = "Overall Survival Probability",
  risk.table=TRUE,
  conf.int=FALSE,
  pval=TRUE,
  pval.method=TRUE)

km.plots[[1]] + background_grid()

survfit(Surv(time, status) ~ Tx_Chemotherapy_CHOP, data = data) 

res.cox <- coxph(
  formula=Surv(time, status) ~
    Tx_Chemotherapy_L_Asp,
  data = data,
  method = 'exact')
summary(res.cox)



print('All done!')