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
library(hashmap)

######## Create output scaffolding
ifelse(!dir.exists("./fig"), dir.create("./fig"), FALSE)
ifelse(!dir.exists("./fig/general"), dir.create("./fig/general"), FALSE)

ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)
ifelse(!dir.exists("./calc/general"), dir.create("./calc/general"), FALSE)

#################################################
## LOAD DATA
#################################################

data <- read.csv('./data/SystematicReview.csv')


#####
# Clean Data
#####
stage.numeric <- data$AnnArborStaging
stage.numeric <- str_replace(stage.numeric, '[AE]$', '')
char2num <- hashmap(
    c('','I', 'II', 'III', 'IV'),
    c(NA,1,2,3,4)
)
data$stage.numeric <- char2num[[stage.numeric]]

# extract geographic region information
data$nktl.region.rate <- hashmap(
    c('EastAsia', 'SouthAsia', 'SoutEastAsia', 
        'SouthAmerica', 
        'Europe', 'MiddleEast', 'NorthAmerica', 'Australia'),
    c(rep('Asia', 3),
      'SouthAmerica',
      rep('Other', 4))
)[[data$GeographicLocation]]


data$ebv.endemic.loc <- data$GeographicLocation %in% c('EastAsia', 'SouthAsia', 'SoutEastAsia')

# extract tumor location
has.orbital <- c()
has.intraocular <- c()
has.nasosinus <- c()
has.lacrimal <- c()
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
        any(str_detect(unlist(loc.toks), pattern=c('nose', 'naso', 'nasal'))))
    has.lacrimal <- c(has.lacrimal,
        any(str_detect(unlist(loc.toks), pattern=c('lacrimal'))))
    has.conjunctival <- c(has.conjunctival,
        any(str_detect(unlist(loc.toks), pattern=c('conjunctiva'))))
    has.cns <- c(has.cns,
        any(str_detect(unlist(loc.toks), pattern=c('cns'))))
}

data$has.orbital <- has.orbital
data$has.intraocular <- has.intraocular
data$has.nasosinus <- has.nasosinus
data$has.lacrimal <- has.lacrimal
data$has.conjunctival <- has.conjunctival
data$has.cns <- has.cns


# extract therapy type
has.surgical <- c()
has.ct <- c()
has.rt <- c()

for (i in 1:nrow(data)) {
    loc.str <- data$Treatment[[i]]
    loc.toks <- str_to_lower(str_split(loc.str, pattern=';', simplify=FALSE))
    has.surgical <- c(has.surgical, 
        any(str_detect(unlist(loc.toks), pattern=c('surgical'))))
    has.ct <- c(has.ct,
        any(str_detect(unlist(loc.toks), pattern=c('ct'))))
    has.rt <- c(has.rt,
        any(str_detect(unlist(loc.toks), pattern=c('rt'))))
}
data$has.surgical <- has.surgical
data$has.ct <- has.ct
data$has.rt <- has.rt

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

head(data[,c('has.ptosis', 'has.proptosis', 'Presenting.Symptoms')])

#####
# PRISMA diagram
#####

prisma.data <- read.csv('./data/deduplicated_searches.csv')
table(prisma.data$ExclusionReason)

#####
# AGE OF ONSET HISTOGRAM
#####

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

data$binned.age <- age.to.bin.map[[data$Age]]

ggplot(data, aes(x=binned.age, fill=Sex)) +
    geom_bar() +
    xlab('Age') + ylab('Number of Cases') +
    theme_cowplot()
ggsave('./fig/general/age_histogram.png', width=8.5, height=8)

ggplot(data, aes(x=Age, fill=Sex)) +
    geom_histogram(binwidth=10, position='stack')

ggplot(data, aes(x=Age, fill=Sex)) + 
    geom_density(alpha=0.5)

#####
# KAPLAN-MEIER SURVIVAL ANALYSIS
#####

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

surv_median(surv_fit(Surv(time, status) ~ 1, data=data))
#####
# Additional Analysis for manuscript
#####

## 
print("Number of cases surviving at time of report writing or for >12mo:")
sum(data$time >= 12 | !data$status)
print("total number of cases:")
nrow(data)

## Cox-proportional hazards to identify factors assoc. with favorable prognosis

res.cox <- coxph(
    formula=Surv(time, status) ~ nktl.region.rate == 'Asia',
    data = data,
    method = 'exact')
summary(res.cox)

ggsurvplot(
    fit = surv_fit(Surv(time, status) ~ (
        nktl.region.rate == 'Asia'
    ), data=data),
    xlab = "Months",
    ylab = "Overall Survival Probability")

