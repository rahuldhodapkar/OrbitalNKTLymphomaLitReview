#!/usr/bin/env Rscript
#
## BuildFigures.R
#
# Generate simple figures and summary statistics for systematic review data.
#

library(ggplot2)
library(survival) # must import survival before survminer
library(survminer)
library(stringr)
library(cowplot)
library(MASS)
library(reshape2)

######## Create output scaffolding
ifelse(!dir.exists("./fig"), dir.create("./fig"), FALSE)
ifelse(!dir.exists("./fig/general"), dir.create("./fig/general"), FALSE)

ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)
ifelse(!dir.exists("./calc/general"), dir.create("./calc/general"), FALSE)

#################################################
## Helper Functions
#################################################

hashmap <- function(keys, vals) {
    x <- new.env(hash=TRUE)
    for(i in 1:length(keys)) {
        x[[as.character(keys[[i]])]] <- vals[[i]]
    }
    return(x)
}

happly <- function(hm, vec, ...) {
    return(sapply(
      X = vec, 
      FUN = function(x) { hm[[as.character(x)]] },
      ...))
}

#################################################
## LOAD DATA
#################################################

data <- read.csv('./data/SystematicReview.csv')

#####
# Clean Data
#####
stage.numeric <- sapply(str_match(data$AnnArborStageInferred, '^([IV]*)')[,1], function(x) {
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

data$stage.numeric <- stage.numeric

# extract geographic region information
data$nktl.region.rate <- factor(
        happly(hashmap(
            c('EastAsia', 'SouthAsia', 'SouthEastAsia',
                'SouthAmerica',
                'Europe', 'MiddleEast', 'NorthAmerica', 'Australia'),
            c(rep('Asia', 3),
              'SouthAmerica',
              rep('Other', 4))
        ), data$GeographicLocation),
    levels=c('Other', 'Asia', 'SouthAmerica'),
    ordered=FALSE)


data$ebv.endemic.loc <- data$GeographicLocation %in% c('EastAsia', 'SouthAsia', 'SouthEastAsia')

# extract tumor location
has.orbital <- c()
has.intraocular <- c()
has.nasosinus <- c()
has.lacrimal <- c()
has.lacrimal.gland <- c()
has.lacrimal.drainage <- c()
has.conjunctival <- c()
has.cns <- c()

for (i in 1:nrow(data)) {
    loc.str <- paste0(
        data$TumorLocationOrbital[[i]],
        data$TumorLocationSystemic[[i]], sep=';')

    loc.toks <- str_to_lower(str_split(loc.str, pattern=';', simplify=FALSE))
    has.orbital <- c(has.orbital, 
        any(str_detect(unlist(loc.toks), pattern=c('orbital'))))
    has.intraocular <- c(has.intraocular,
        any(str_detect(unlist(loc.toks), pattern=c('intraocular'))))
    has.nasosinus <- c(has.nasosinus,
        any(str_detect(unlist(loc.toks), pattern=c('nose', 'nasopharynx', 'nasal'))))
    has.lacrimal <- c(has.lacrimal,
        any(str_detect(unlist(loc.toks), pattern=c('lacrimal'))))
    has.lacrimal.gland <- c(has.lacrimal.gland,
        any(str_detect(unlist(loc.toks), pattern=c('lacrimal gland'))))
    has.lacrimal.drainage <- c(has.lacrimal.drainage,
        any(str_detect(unlist(loc.toks), pattern=c('nasolacrimal', 'lacrimal sac'))))
    has.conjunctival <- c(has.conjunctival,
        any(str_detect(unlist(loc.toks), pattern=c('conjunctiva'))))
    has.cns <- c(has.cns,
        any(str_detect(unlist(loc.toks), pattern=c('cns'))))
}

data$has.orbital <- has.orbital
data$has.intraocular <- has.intraocular
data$has.nasosinus <- has.nasosinus
data$has.lacrimal <- has.lacrimal
data$has.lacrimal.gland <- has.lacrimal.gland
data$has.lacrimal.drainage <- has.lacrimal.drainage
data$has.conjunctival <- has.conjunctival
data$has.cns <- has.cns


# extract therapy type
has.surgical <- c()
has.orbitotomy.surgical <- c()

has.ct <- c()

has.gelox.ct <- c()
has.dhap.ct <- c()
has.methotrexate.ct <- c()
has.chop.ct <- c()
has.ceop.ct <- c()
has.smile.ct <- c()
has.devic.ct <- c()

has.rt <- c()
has.localized.rt <- c()
has.total.body.rt <- c()

for (i in 1:nrow(data)) {
    loc.str <- data$Treatment[[i]]
    loc.toks <- str_to_lower(str_split(loc.str, pattern=';', simplify=FALSE))
    has.surgical <- c(has.surgical, 
        any(str_detect(unlist(loc.toks), pattern=c('surgical'))))
    has.ct <- c(has.ct,
        any(str_detect(unlist(loc.toks), pattern=c(' ct', '^ct'))))
    has.rt <- c(has.rt,
        any(str_detect(unlist(loc.toks), pattern=c(' rt'))))

    has.gelox.ct <- c(has.gelox.ct,
        any(str_detect(unlist(loc.toks), pattern=c('gelox'))))
    has.dhap.ct <- c(has.dhap.ct,
        any(str_detect(unlist(loc.toks), pattern=c('dhap'))))
    has.methotrexate.ct <- c(has.methotrexate.ct,
        any(str_detect(unlist(loc.toks), pattern=c('methotrexate'))))
    has.chop.ct <- c(has.chop.ct,
        any(str_detect(unlist(loc.toks), pattern=c('chop'))))
    has.ceop.ct <- c(has.ceop.ct,
        any(str_detect(unlist(loc.toks), pattern=c('ceop'))))
    has.smile.ct <- c(has.smile.ct,
        any(str_detect(unlist(loc.toks), pattern=c('smile'))))
    has.devic.ct <- c(has.devic.ct,
        any(str_detect(unlist(loc.toks), pattern=c('2/3devic'))))

    has.orbitotomy.surgical <- c(has.orbitotomy.surgical, 
        any(str_detect(unlist(loc.toks), pattern=c('orbitotomy'))))

    has.localized.rt <- c(has.localized.rt, 
        any(str_detect(unlist(loc.toks), pattern=c('localized rt'))))
    has.total.body.rt <- c(has.total.body.rt, 
        any(str_detect(unlist(loc.toks), pattern=c('total body rt'))))
}
data$has.surgical <- has.surgical
data$has.orbitotomy.surgical <- has.orbitotomy.surgical
data$has.unspecified.surgical <- (
  (data$has.orbitotomy.surgical 
  - data$has.surgical) == -1
)

data$has.ct <- has.ct
data$has.gelox.ct <- has.gelox.ct
data$has.dhap.ct <- has.dhap.ct
data$has.methotrexate.ct <- has.methotrexate.ct
data$has.chop.ct <- has.chop.ct
data$has.ceop.ct <- has.ceop.ct
data$has.smile.ct <- has.smile.ct
data$has.devic.ct <- has.devic.ct
data$has.unspecified.ct <- (
  ((data$has.gelox.ct
    | data$has.dhap.ct
    | data$has.methotrexate.ct
    | data$has.chop.ct
    | data$has.ceop.ct
    | data$has.smile.ct
    | data$has.devic.ct)
   - data$has.ct) == -1
)

data$has.rt <- has.rt
data$has.localized.rt <- has.localized.rt
data$has.total.body.rt <- has.total.body.rt
data$has.unspecified.rt <- (
  ((data$has.localized.rt
    | data$has.total.body.rt)
   - data$has.rt) == -1
)

# extract ocular symptoms
has.vision.sx <- c()
has.eom.sx <- c()
has.lid.sx <- c()
has.proptosis <- c()
has.ptosis <- c()
for (i in 1:nrow(data)) {
    loc.str <- data$Presenting.Symptoms[[i]]
    loc.toks <- str_to_lower(str_split(loc.str, pattern=';', simplify=FALSE))
    has.vision.sx <- c(has.vision.sx, 
        any(str_detect(unlist(loc.toks), pattern=c('vision'))))
    has.eom.sx <- c(has.eom.sx,
        any(str_detect(unlist(loc.toks), pattern=c('eom', 'fixation'))))
    has.lid.sx <- c(has.lid.sx,
        any(str_detect(unlist(loc.toks), pattern=c('lid'))))
    has.proptosis <- c(has.proptosis,
        any(str_detect(unlist(loc.toks), pattern=c('proptosis'))))
    has.ptosis <- c(has.ptosis,
        any(str_detect(unlist(loc.toks), pattern=c('[^o]ptosis'))))
}
data$has.vision.sx <- has.vision.sx
data$has.eom.sx <- has.eom.sx
data$has.lid.sx <- has.lid.sx
data$has.proptosis <- has.proptosis
data$has.ptosis <- has.ptosis

#' extract markers from a vector of delimited strings
#' 
#' @param x A vector of delimited strings
#' @param markers.list a list of markers to extract
#' @param split.pattern The delimiter
#' @param pos.indicator.string the suffix indicating a positive entry
#' @param neg.indicator.string the suffix indicating a negative entry
#' @return Data frame of length equal to length(x), with all markers
#'
ExtractMarkers <- function(x, marker.list, 
    split.pattern=';',
    pos.indicator.string='\\+',
    neg.indicator.string='\\-'
    ) {

    df.entries <- lapply(
        marker.list,
        function(name) {
            rep(as.logical(NA), length(x))
        }
    )
    names(df.entries) <- marker.list
    df <- do.call(data.frame, df.entries)

    marker2name <- hashmap(marker.list, colnames(df))
    for (i in 1:length(x)) {
        toks <- str_split(tolower(x[[i]]), pattern=split.pattern, simplify=FALSE)
        for(m in marker.list) {
            if(any(str_detect(unlist(toks), 
                    pattern=c(paste0(tolower(m),pos.indicator.string))))) {
                df[[i,happly(marker2name, m)]] <- TRUE
            }
            if(any(str_detect(unlist(toks), 
                    pattern=c(paste0(tolower(m),neg.indicator.string))))) {
                if(!is.na(df[[i,m]])){
                    stop(paste0('Corrupt data at SystematicReview.csv: line ',
                        as.character(i), ' marker ',m,' is both + and -.'))
                }

                df[[i,happly(marker2name, m)]] <- FALSE
            }
        }
    }
    return(df)
}

ihc.markers <- ExtractMarkers(data$Immunohistochemistry,
    c('CD3','CD4','CD8','CD20','CD56','perforin','GrB','TIA-1','EBER'))

data <- cbind(data, ihc.markers)

sum(rowSums(!is.na(as.matrix(ihc.markers))) == 0)

# add EBV serology data
ebv.sero.status <- rep(as.logical(NA), nrow(data))
ebv.sero.status[grepl('Positive', data$SerologicalEBVStatus)] <- TRUE
ebv.sero.status[grepl('Negative', data$SerologicalEBVStatus)] <- FALSE
data$ebv.sero.status <- ebv.sero.status

################################################################################
## Persist cleaned data to CSV
################################################################################

write.csv(data, file='./calc/general/cleaned_data.csv', row.names=FALSE)
message('cleaned data persisted at `./calc/general/cleaned_data.csv`')

################################################################################
# PRISMA diagram
################################################################################

prisma.data <- read.csv('./data/deduplicated_searches.csv')
print('===== Data Exclusion Reason for PRISMA Diagram =====')
table(prisma.data$ExclusionReason)

################################################################################
# AGE OF ONSET HISTOGRAM
################################################################################

age.to.bin.map <- hashmap(
    c(
        seq(0,9),
        seq(10,19),
        seq(20,29),
        seq(30,39),
        seq(40,49),
        seq(50,59),
        seq(60,69),
        seq(70,79),
        seq(80,89),
        seq(90,99)
    ),
    c(
        rep("0-9", 10),
        rep("10-19", 10),
        rep("20-29", 10),
        rep("30-39", 10),
        rep("40-49", 10),
        rep("50-59", 10),
        rep("60-69", 10),
        rep("70-79", 10),
        rep("80-89", 10),
        rep("90-99", 10)
    )
)

data$binned.age <- happly(age.to.bin.map, data$Age)

ggplot(data, aes(x=binned.age, fill=Sex)) +
    geom_bar() +
    xlab('Age') + ylab('Number of Cases') +
    theme_cowplot()
ggsave('./fig/general/age_histogram.png', width=8.5, height=8)

ggplot(data, aes(x=Age, fill=Sex)) +
    geom_histogram(binwidth=10, position='stack')

ggplot(data, aes(x=Age, fill=Sex)) + 
    geom_density(alpha=0.5)

################################################################################
# KAPLAN-MEIER SURVIVAL ANALYSIS
################################################################################

# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
# status FALSE := censored, TRUE := dead
data$status <- (!str_starts(data$Survival.From.Diagnosis..mo., ">"))
data$time <- data$Survival.From.Diagnosis..mo. %>%
            str_replace(">", "") %>%
            as.integer()

km.plots <- ggsurvplot(
    fit = surv_fit(Surv(time, status) ~ 1, data=data),
    xlab = "Months",
    xlim = c(1,23),
    break.x.by = 3,
    ylab = "Overall Survival Probability",
    risk.table=TRUE)

km.plots[[1]] + background_grid()
ggsave('./fig/general/kaplan_meier_all.png', height=8.5, width=8)

print('===== Overall Median Survival =====')
surv_median(surv_fit(Surv(time, status) ~ 1, data=data))

################################################################################
# Additional Analysis for manuscript
################################################################################

## 
print("Number of cases surviving at time of report writing or for >12mo:")
sum(data$time >= 12 | !data$status)
sum(data$time >= 6 | !data$status)
sum(data$time >= 3 | !data$status)
print("total number of cases:")
nrow(data)

#####
# Collate Data for Table 1
#####

PrintDataForTable <- function(df, col.name) {
  print(paste0('Counts for [', col.name, ']'))
  print(table(data[,col.name]))
  print(paste0('Percentage for [', col.name, ']'))
  print(table(data[,col.name]) / sum(table(data[,col.name])))
}

# Age
#<44,45-64,65+
age.to.bin.map.tab1 <- hashmap(
    c(
        seq(0,39),
        seq(40,64),
        seq(65,99)
    ),
    c(
        rep("<40", length(seq(0,39))),
        rep("40-64", length(seq(40,64))),
        rep("65+", length(seq(65,99)))
    )
)
data$age.tab1 <- factor(
    happly(age.to.bin.map.tab1, data$Age),
    levels=c('<40','40-64', '65+'),
    ordered=TRUE)

print('===== Age Table Data =====')
PrintDataForTable(data, 'age.tab1')

# Sex
print('===== Sex Table Data =====')
PrintDataForTable(data, 'Sex')

print('t-test for sex: ')
t.test(data$Age[data$Sex == 'M'],
       data$Age[data$Sex == 'F'],
       paired=FALSE)

# Region
print('===== Region Table Data =====')
PrintDataForTable(data, 'GeographicLocation')

# Symptoms
print('===== Symptoms Table Data =====')
sx.vec <- c(
  "has.vision.sx",
  "has.eom.sx",
  "has.lid.sx",
  "has.ptosis",
  "has.proptosis"
)
for (s in sx.vec) {
  print('===== > ')
  PrintDataForTable(data, s)
}

# Tumor Location
print('===== Tumor Location Table Data =====')
loc.vec <- c(
  "has.orbital",
  "has.intraocular",
  "has.lacrimal.gland",
  "has.lacrimal.drainage",
  "has.conjunctival",
  "has.nasosinus",
  "has.cns"
)
for (l in loc.vec) {
  print('===== > ')
  PrintDataForTable(data, l)
}

# Ki67
print('===== Ki67 Table Data =====')
data$Ki67.binned <- sapply(data$Ki67, function(x) {
  if (is.na(x)) {
    return(NA)
  } else {
    if      (x < 0.5)             { return ('Ki67 < 0.5') }
    else if (x >= 0.5 && x < 0.8) { return ('0.5 <= Ki67 < 0.8') }
    else if (x >= 0.8)            { return ('0.8 <= Ki67') }
    else                          { stop('Invalid Ki67 Value!') }
  }
})
PrintDataForTable(data, 'Ki67.binned')

# Therapy
print('===== Therapy Table Data =====')

PrintDataForTable(data, 'has.surgical')
for (x in c('has.orbitotomy.surgical', 'has.unspecified.surgical')) {
  print('===== > ')
  PrintDataForTable(data, x)
}

PrintDataForTable(data, 'has.rt')
for (x in c('has.localized.rt',
            'has.total.body.rt',
            'has.unspecified.rt')) {
  print('===== > ')
  PrintDataForTable(data, x)
}

PrintDataForTable(data, 'has.ct')
for (x in c('has.gelox.ct',
            'has.dhap.ct',
            'has.methotrexate.ct',
            'has.chop.ct',
            'has.ceop.ct',
            'has.smile.ct',
            'has.devic.ct')) {
  print('===== > ')
  PrintDataForTable(data, x)
}

#######
## Multivariate Cox Regression
#######
print('===== Run Linear Regression =====')
res.cox <- coxph(
    formula=Surv(time, status) ~
        strata(age.tab1) + # Demo
        strata(Sex) +
        stage.numeric +
        has.surgical + # Tx
        has.ct +
        has.rt +
        has.vision.sx + # SX
        has.eom.sx +
        has.lid.sx +
        has.ptosis +
        has.proptosis +
        has.orbital + # Loc
        has.intraocular +
        has.lacrimal.gland +
        has.lacrimal.drainage +
        has.conjunctival +
        has.nasosinus +
        has.cns,
    data = data,
    method = 'exact')
summary(res.cox)

#
# For the correlation analsis between symptoms and tumor location, a
# backward stepwise linear regression using the `stepAIC` function
# from the `MASS` R package (https://www.stats.ox.ac.uk/pub/MASS4/)
#

print('===== Run Stepwise Regression and Symptom Correlation Analysis =====')
step.model <- stepAIC(res.cox, direction = "backward", 
                      trace = FALSE)
summary(step.model)

# consider only those symptoms that are included in the final model
sx.vec <- c(
  "has.eom.sx",
  "has.ptosis",
  "has.proptosis"
)

# consider only those locations included in the final model
loc.vec <- c(
  "has.orbital",
  "has.lacrimal.drainage",
  "has.conjunctival",
  "has.nasosinus"
)

# show significant correlation between tumor location and reported symptoms.
jaccard.sim <- matrix(0, nrow=length(sx.vec), ncol=length(loc.vec),
                      dimnames=list(sx.vec, loc.vec))
cor.val <- matrix(0, nrow=length(sx.vec), ncol=length(loc.vec),
                  dimnames=list(sx.vec, loc.vec))
fisher.p <- matrix(0, nrow=length(sx.vec), ncol=length(loc.vec),
                      dimnames=list(sx.vec, loc.vec))

for (i in 1:length(sx.vec)) {
  sx <- sx.vec[[i]]
  for (j in 1:length(loc.vec)) {
    loc <- loc.vec[[j]]
    c.tab <- table(data[,sx], data[,loc])
    
    jaccard.sim[[i,j]] <- (
      sum(data[,sx] & data[,loc]) / sum(data[,sx] | data[,loc])
    )
    cor.val[[i,j]] <- cor(data[,sx], data[,loc])
    fisher.p[[i,j]] <- fisher.test(table(data[,sx], data[,loc]))$p.value
  }
}

fisher.p.adj <- fisher.p %>% 
  as.numeric() %>%
  p.adjust(method='bonferroni') %>%
  matrix(nrow=length(sx.vec), ncol=length(loc.vec),
    dimnames = list(sx.vec, loc.vec), byrow = FALSE)
fisher.p.adj

cor.val.df <- melt(cor.val)
fisher.p.adj.df <- melt(fisher.p.adj)

# plot correlation, along with significance data.
ggplot(cor.val.df, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  geom_point(
    data=subset(fisher.p.adj.df, value < 0.05),
    aes(x=Var1, y=Var2, fill=value),
    shape=8
  ) +
  theme_cowplot() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  )
ggsave('./fig/general/symptom_location_correlation.png', width=8, height = 8)
print('p-value for eom <> orbital correlation = 0.02233017')

print('All done!')