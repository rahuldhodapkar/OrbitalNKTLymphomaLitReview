#!/usr/bin/env Rscript
#
## DimensionalityReductionCoxRegression.R
#
# Generate simple figures and summary statistics for systematic review data.
#

library(stringr)
library(dplyr)
library(puddlr)
library(ComplexHeatmap)
library(viridis)
library(cowplot)
library(ggrepel)
library(RColorBrewer)
library(ggnewscale)

######## Create output scaffolding
s <- ifelse(!dir.exists("./fig"), dir.create("./fig"), FALSE)
s <- ifelse(!dir.exists("./fig/dimreduc"), dir.create("./fig/dimreduc"), FALSE)

s <- ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)
s <- ifelse(!dir.exists("./calc/dimreduc"), dir.create("./calc/dimreduc"), FALSE)

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
data$Dx_GTE_3Mo_Misdiagnosis <- as.numeric(data$Dx_Misdiagnosis_TimeToDiagnosis) > 3
data$Dx_GTE_3Mo_Misdiagnosis[is.na(data$Dx_GTE_3Mo_Misdiagnosis)] <- FALSE


!is.na(data$time)
meta.data <- data[,c(
  'status', 'time',
  'Demo_Age',
  'Demo_Sex'
)]


# extract binary data from target cols
data <- data[,
  c('Dx_IsLateStage',
    'Tx_Chemotherapy_L_Asp',
    'Tx_Chemotherapy',
    'Tx_Radiotherapy',
    'Tx_Surgical',
    'Sx_VisionLoss',
    'Sx_RestrictedEOM',
    'Sx_PeriorbitalSwelling',
    'Sx_Chemosis',
    'Sx_Proptosis',
    'Sx_EyelidSwelling',
    'Sx_Ptosis',
    'Loc_IsRecurrence',
    'Loc_OrbitalLocationOrbital',
    'Loc_OrbitalLocationLacrimalGland',
    'Loc_OrbitalLocationLacrimalDrainage',
    'Loc_OrbitalLocationConjunctival',
    'Loc_OrbitalLocationUveal',
    'Loc_OrbitalIntraocularExtension',
    'Loc_OrbitalLocationEyelid',
    'Loc_LocationNasosinus',
    'Loc_LocationCNS')
]

for (c in colnames(data)) {
  if (typeof(data[,c]) == 'logical') {
    next
  } else if (typeof(data[,c]) == 'double') {
    next
  }
  data[,c] <- data[,c] == 'YES'
}

for (c in colnames(data)) {
  data[,c] <- as.numeric(data[,c])
}

data.puddlr <- data

################################################################################
## BUILD Complex Heatmap
################################################################################

data <- readRDS('./calc/general/data.rds')
meta.data <- readRDS('./calc/general/meta_data.rds')

data <- data[,
             c('Sx_VisionLoss',
               'Sx_RestrictedEOM',
               'Sx_PeriorbitalSwelling',
               'Sx_Chemosis',
               'Sx_Proptosis',
               'Sx_EyelidSwelling',
               'Sx_Ptosis',
               'Loc_OrbitalLocationOrbital',
               'Loc_OrbitalLocationLacrimalGland',
               'Loc_OrbitalLocationLacrimalDrainage',
               'Loc_OrbitalLocationConjunctival',
               'Loc_OrbitalLocationUveal',
               'Loc_OrbitalIntraocularExtension',
               'Loc_OrbitalLocationEyelid',
               'Loc_LocationNasosinus',
               'Loc_LocationCNS',
               'Dx_IHC_CD2',
               'Dx_IHC_CD3',
               'Dx_IHC_CD4',
               'Dx_IHC_CD5',
               'Dx_IHC_CD7',
               'Dx_IHC_CD8',
               'Dx_IHC_CD20',
               'Dx_IHC_CD30',
               'Dx_IHC_CD45',
               'Dx_IHC_CD56',
               'Dx_IHC_CD79a',
               'Dx_IHC_TIA1',
               'Dx_IHC_LCA',
               'Dx_IHC_LMP1',
               'Dx_IHC_GrB',
               'Dx_IHC_Perforin',
               'Dx_IHC_Ki67_GTE_80',
               'Dx_IHC_EBER',
               'Tx_Surgical',
               'Tx_Surgical_FESS',
               'Tx_Surgical_Orbitotomy',
               'Tx_Chemotherapy',
               'Tx_Chemotherapy_CHOP',
               'Tx_Chemotherapy_Methotrexate',
               'Tx_Chemotherapy_DHAP',
               'Tx_Chemotherapy_CEOP',
               'Tx_Chemotherapy_DeVIC',
               'Tx_Chemotherapy_SMILE',
               'Tx_Chemotherapy_GELOX',
               'Tx_Radiotherapy',
               'Tx_Radiotherapy_Localized')
]

for (c in colnames(data)) {
  if (typeof(data[,c]) == 'logical') {
    next
  } else if (typeof(data[,c]) == 'double') {
    next
  }
  yes.rows <- data[,c] == 'YES'
  no.rows <- data[,c] == 'NO'
  na.rows <- data[,c] == ''
  
  data[,c] <- NA
  data[yes.rows,c] <- 1
  data[no.rows,c] <- 0
}

hamming.dist <- function(m) {
  dist.mat <- matrix(0, nrow=nrow(m), ncol=nrow(m))
  for (i in 1:nrow(m)) {
    for (j in 1:nrow(m)) {
      ixs.to.compare <- !is.na(m[i,]) & !is.na(m[j,])
      dist.mat[i,j] <- sum(m[i,ixs.to.compare] != m[j,ixs.to.compare])
    }
  }
  return(as.dist(dist.mat))
}

jaccard <- function(a, b) {
  return (sum( a & b )/sum( a | b ))
}

jaccard.dist <- function(m) {
  dist.mat <- matrix(0, nrow=nrow(m), ncol=nrow(m))
  for (i in 1:nrow(m)) {
    for (j in 1:nrow(m)) {
      dist.mat[i,j] <- 1 - jaccard(m[i,], m[j,])
    }
  }
  return(as.dist(dist.mat))
}

data.matrix <- as.matrix(data) %>% t()
rownames(data.matrix) <- colnames(data)

ha <- HeatmapAnnotation(
  age = meta.data$Demo_Age,
  sex = meta.data$Demo_Sex,
  time = meta.data$time,
  status = meta.data$status,
  col = list(
    'age' = circlize::colorRamp2(
      breaks=seq(min(meta.data$Demo_Age, na.rm=T), 
                 max(meta.data$Demo_Age, na.rm=T),
                 length.out=100),
      colors=magma(100)
    ),
    'sex' = c(
      'M' = '#347DC1', 'F' = '#E6A6C7'
    ),
    'time' = circlize::colorRamp2(
      breaks=seq(min(meta.data$time, na.rm=T), 
                 max(meta.data$time, na.rm=T),
                 length.out=100),
      colors=viridis(100)
    ),
    'status' = c(
      'TRUE' = '#c43025', 'FALSE' = '#c99591'
    )
  )
)


paired.pal <- brewer.pal(n=8,name='Paired')

col1 = circlize::colorRamp2(c(0, 0.5, 1), c(paired.pal[[1]], 'white', paired.pal[[2]]))
col2 = circlize::colorRamp2(c(0, 0.5, 1), c(paired.pal[[3]],'white', paired.pal[[4]]))
col3 = function(x) {
  ifelse(is.na(x), '#FFFFFF', circlize::colorRamp2(c(0, 0.5, 1), c(paired.pal[[5]], 'white', paired.pal[[6]]))(x))
}
col4 = circlize::colorRamp2(c(0, 0.5, 1), c(paired.pal[[7]],'white', paired.pal[[8]]))


png('./fig/dimreduc/symptoms_heatmap_colclust.png',
    width=14,height=13.5,units="in",res=300)
htmp <- Heatmap(data.matrix,
                col = list('0' = '#DDDDDD', '1' = '#222244', 'NA'='#FFFFFF'),
                clustering_distance_columns = function(m) hamming.dist(m),
                cluster_columns=TRUE,
                cluster_rows=FALSE,
                top_annotation=ha,
                show_column_names = FALSE,
                show_row_names = TRUE,
                heatmap_legend_param = list(
                  direction = "horizontal"),
                rect_gp = gpar(type='none'),
                show_heatmap_legend = FALSE,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(str_starts(rownames(data.matrix)[i], "Sx_")) {
                    grid.rect(x, y, w, h, gp = gpar(fill = col1(data.matrix[i, j]), col = '#DDDDDD'))
                  } else if (str_starts(rownames(data.matrix)[i], "Loc_")) {
                    grid.rect(x, y, w, h, gp = gpar(fill = col2(data.matrix[i, j]), col = '#DDDDDD'))
                  } else if (str_starts(rownames(data.matrix)[i], "Dx_")){
                    grid.rect(x, y, w, h, gp = gpar(fill = col3(data.matrix[i, j]), col = '#DDDDDD'))
                  } else if (str_starts(rownames(data.matrix)[i], "Tx_")){
                    grid.rect(x, y, w, h, gp = gpar(fill = col4(data.matrix[i, j]), col = '#DDDDDD'))
                  }
                })
draw(htmp,
     heatmap_legend_side='bottom', 
     annotation_legend_side='bottom',
     legend_grouping='original',
     padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()

################################################################################
## PCA and Cox Regression
################################################################################

data <- readRDS('./calc/general/data.rds')
meta.data <- readRDS('./calc/general/meta_data.rds')

data <- data[,
             c('Dx_IsLateStage',
               'Tx_Chemotherapy_L_Asp',
               'Tx_Chemotherapy',
               'Tx_Radiotherapy',
               'Tx_Surgical',
               'Sx_VisionLoss',
               'Sx_RestrictedEOM',
               'Sx_PeriorbitalSwelling',
               'Sx_Chemosis',
               'Sx_Proptosis',
               'Sx_EyelidSwelling',
               'Sx_Ptosis',
               'Loc_IsRecurrence',
               'Loc_OrbitalLocationOrbital',
               'Loc_OrbitalLocationLacrimalGland',
               'Loc_OrbitalLocationLacrimalDrainage',
               'Loc_OrbitalLocationConjunctival',
               'Loc_OrbitalLocationUveal',
               'Loc_OrbitalIntraocularExtension',
               'Loc_OrbitalLocationEyelid',
               'Loc_LocationNasosinus',
               'Loc_LocationCNS')
]

for (c in colnames(data)) {
  if (typeof(data[,c]) == 'logical') {
    next
  } else if (typeof(data[,c]) == 'double') {
    next
  }
  data[,c] <- data[,c] == 'YES'
}

for (c in colnames(data)) {
  data[,c] <- as.numeric(data[,c])
}

# need to set full set of predictors

data.matrix <- as.matrix(data)
rownames(data.matrix) <- 1:nrow(data.matrix)

puddlr.obj <- CreatePuddlrObject(meta.data$time, data.matrix)
puddlr.obj <- NormalizePredictors(puddlr.obj)
puddlr.obj <- RunPCA(puddlr.obj)

k.cross <- 10
rand.seed <- 42
reduction <- 'pca'

avg.train.rsq <- c()
avg.rmse <- c()

for (n.components in 2:20) {
  print(paste0("Testing with [", 
               as.character(n.components) ,"] components"))
  folds <- cut(sample(1:length(puddlr.obj$response)),breaks=k.cross,labels=FALSE)
  
  rmse.values <- c()
  rsq.values <- c()
  for(k.i in 1:k.cross){
    #Segement data by fold
    test.ixs <- which(folds==k.i,arr.ind=TRUE)
    train.ixs <- which(folds!=k.i,arr.ind=TRUE)
    
    # train data
    train.puddlr <- puddlr.obj
    train.puddlr$reductions[[reduction]]$rotation <- (
      train.puddlr$reductions[[reduction]]$rotation)
    train.puddlr$reductions[[reduction]]$embedding <- (
      train.puddlr$reductions[[reduction]]$embedding[train.ixs,])
    
    train.puddlr$predictors <- (train.puddlr$predictors[train.ixs,])
    train.puddlr$response <- (train.puddlr$response[train.ixs])
    
    # test data
    test.puddlr <- puddlr.obj
    test.puddlr$reductions[[reduction]]$rotation <- (
      test.puddlr$reductions[[reduction]]$rotation)
    test.puddlr$reductions[[reduction]]$embedding <- (
      test.puddlr$reductions[[reduction]]$embedding[test.ixs,])
    
    test.puddlr$predictors <- (test.puddlr$predictors[test.ixs,])
    test.puddlr$response <- (test.puddlr$response[test.ixs])
    
    cox.df <- as.data.frame(train.puddlr$reductions$pca$embedding[,1:n.components])
    cox.df$time <- meta.data[rownames(cox.df), 'time']
    cox.df$status <- meta.data[rownames(cox.df), 'status']
    
    res.cox <- coxph(
      formula=Surv(time, status) ~ .,
      data = cox.df,
      method = 'exact')
    
    test.df <- data.frame(test.puddlr$reductions[[reduction]]$embedding)
    test.df$time <- meta.data[rownames(test.df), 'time']
    test.df$status <- meta.data[rownames(test.df), 'status']
    
    # Root mean square for the martingale residuals
    mresid <- (test.df$status-1) - predict(res.cox, test.df, type='expected')
    
    rmse.values <- c(rmse.values, 
                     sqrt(mean((mresid)^2)))
    rsq.values <- c(rsq.values, summary(res.cox)$concordance[[1]])
  }
  avg.rmse <- c(avg.rmse, mean(rmse.values))
  avg.train.rsq <- c(avg.train.rsq, mean(rsq.values))
}

plot.df <- data.frame(
  concordance = avg.train.rsq,
  rms.martingale.residual = avg.rmse,
  epoch = 1:length(avg.train.rsq) + 1
)

ggplot(plot.df, aes(x=epoch, y=concordance)) +
  geom_line() + 
  ylab('Concordance') +
  xlab('Number of Principal Components') +
  theme_cowplot() + background_grid()
ggsave('./fig/dimreduc/concordance_plot.png', width=8, height=4)

ggplot(plot.df, aes(x=epoch, y=log10(rms.martingale.residual))) +
  geom_line() +
  ylab('Log10(RMS Martingale Residual)') +
  xlab('Number of Principal Components') +
  theme_cowplot() + background_grid()
ggsave('./fig/dimreduc/rms_martingale_plot.png', width=8, height=4)

################################################################################
## PCA and Cox Regression
################################################################################

n.components <- 10

cox.df <- as.data.frame(puddlr.obj$reductions$pca$embedding[,1:n.components])
cox.df$time <- meta.data[rownames(cox.df), 'time']
cox.df$status <- meta.data[rownames(cox.df), 'status']

res.cox <- coxph(
  formula=Surv(time, status) ~ .,
  data = cox.df,
  method = 'exact')
summary(res.cox)


sum.feature.estimate <- (
  rowSums(t(
    t(puddlr.obj$reductions[['pca']]$rotation[,1:n.components])
    * res.cox$coefficients[1:n.components]
  ))
)


sum.feature.var <- c()
for (i in 1:ncol(puddlr.obj$predictors)) {
  vconv.unscaled <- vcov(res.cox)[
    1:n.components,
    1:n.components]
  v <- puddlr.obj$reductions[['pca']]$rotation[i,1:n.components]
  var.i <- sum(t(vconv.unscaled * v) * v)
  sum.feature.var <- c(sum.feature.var, var.i)
}
names(sum.feature.var) <- colnames(puddlr.obj$predictors)

stats.df <- data.frame(
  VariableName = names(sum.feature.var),
  GLMEstimate = sum.feature.estimate,
  StdErr = sqrt( sum.feature.var )
)
stats.df$ZScore <- stats.df$GLMEstimate / stats.df$StdErr
stats.df$PValue <- 2*pnorm(-abs(stats.df$ZScore))

# generate hazard ratios
stats.df$HazardRatioEstimate <- exp(stats.df$GLMEstimate)
stats.df$Low95CI <- exp(stats.df$GLMEstimate - 1.96 * stats.df$StdErr)
stats.df$Hi95CI <- exp(stats.df$GLMEstimate + 1.96 * stats.df$StdErr)

write.csv(stats.df, './calc/dimreduc/constrained_cox_model.csv', row.names=F)

# create color map for segments
paired.pal <- brewer.pal(n=8,name='Paired')

# create volcano plot
plot.df <- stats.df
plot.df <- plot.df[order(plot.df$PValue, decreasing = FALSE),]
plot.df$VariableName <- factor(plot.df$VariableName,
                               levels=plot.df$VariableName,
                               ordered=T)
plot.df$NegLog10P <- -log10(plot.df$PValue)
plot.df$VariableSegment <- str_match(plot.df$VariableName, '([A-Za-z]+)_')[,2]
plot.df$VariableSegment <- factor(plot.df$VariableSegment,
                                  levels=c('Sx', 'Loc', 'Dx', 'Tx'),
                                  ordered = T)

callout.colors <- paired.pal[1:4 * 2]
names(callout.colors) <- c('Sx', 'Loc', 'Dx', 'Tx')
grey.colors <- paired.pal[1:4 * 2 - 1]
names(grey.colors) <- c('Sx', 'Loc', 'Dx', 'Tx')

callout.df <- plot.df[plot.df$PValue < 0.05,]
grey.df <- plot.df[! (plot.df$VariableName %in% callout.df$VariableName),]

ggplot(grey.df,
       aes(x=log2(exp(GLMEstimate)), y=NegLog10P, label=VariableName)) +
  geom_point(color='#CCCCCC') +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') +
  geom_point(data = callout.df, color='#F8766D', size=4) +
  geom_text_repel(max.overlaps = Inf,
                  nudge_x = ifelse(callout.df$GLMEstimate > 0,  4, -4),
                  nudge_y = 1,
                  force = 10,
                  data=callout.df) +
  xlab('Log2(Hazard Ratio)') +
  ylab('-Log10(P Value)') +
  ylim(0, 6) +
  xlim(-0.75, 0.75) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  theme_minimal_grid()

ggsave('./fig/dimreduc/volcano_labeled.png', width=8, height=8)


ggplot(grey.df,
       aes(x=log2(exp(GLMEstimate)), y=NegLog10P, label=VariableName)) +
  geom_point(color='#CCCCCC') +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') +
  geom_point(data = callout.df, fill='#F8766D', 
             size=4, shape=21, stroke=1.15) +
  xlab('Log2(Hazard Ratio)') +
  ylab('-Log10(P Value)') +
  ylim(0, 6) +
  xlim(-0.5, 0.5) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  theme_minimal_grid()

ggsave('./fig/dimreduc/volcano_unlabeled.png', width=8, height=8)





