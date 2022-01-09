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
data$nktl.region.rate <- factor(
        hashmap(
            c('EastAsia', 'SouthAsia', 'SouthEastAsia',
                'SouthAmerica',
                'Europe', 'MiddleEast', 'NorthAmerica', 'Australia'),
            c(rep('Asia', 3),
              'SouthAmerica',
              rep('Other', 4))
        )[[data$GeographicLocation]],
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
}
data$has.surgical <- has.surgical
data$has.orbitotomy.surgical <- has.orbitotomy.surgical

data$has.ct <- has.ct

data$has.gelox.ct <- has.gelox.ct
data$has.dhap.ct <- has.dhap.ct
data$has.methotrexate.ct <- has.methotrexate.ct
data$has.chop.ct <- has.chop.ct
data$has.ceop.ct <- has.ceop.ct
data$has.smile.ct <- has.smile.ct
data$has.devic.ct <- has.devic.ct

data$has.rt <- has.rt
data$has.localized.rt <- has.localized.rt

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
                df[[i,marker2name[[m]]]] <- TRUE
            }
            if(any(str_detect(unlist(toks), 
                    pattern=c(paste0(tolower(m),neg.indicator.string))))) {
                if(!is.na(df[[i,m]])){
                    stop(paste0('Corrupt data at SystematicReview.csv: line ',
                        as.character(i), ' marker ',m,' is both + and -.'))
                }

                df[[i,marker2name[[m]]]] <- FALSE
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

########
## Persist cleaned data to CSV
########

write.csv(data, file='./calc/general/cleaned_data.csv', row.names=FALSE)
message('cleaned data persisted at `./calc/general/cleaned_data.csv`')

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
sum(data$time >= 6 | !data$status)
sum(data$time >= 3 | !data$status)
print("total number of cases:")
nrow(data)

## Cox-proportional hazards to identify factors assoc. with favorable prognosis

res.cox <- coxph(
    formula=Surv(time, status) ~ nktl.region.rate == 'Asia',
    data = data,
    method = 'exact')
summary(res.cox)


res.cox <- coxph(
    formula=Surv(time, status) ~ Age + Sex + nktl.region.rate,
    data = data,
    method = 'exact')
summary(res.cox)

ggsurvplot(
    fit = surv_fit(Surv(time, status) ~ (
        nktl.region.rate == 'Asia'
    ), data=data),
    xlab = "Months",
    ylab = "Overall Survival Probability")

ggsurvplot(
    fit = surv_fit(Surv(time, status) ~ (
        has.smile.ct
    ), data=subset(data, xor(has.chop.ct, has.smile.ct))),
    xlab = "Months",
    ylab = "Overall Survival Probability")

km.plots <- ggsurvplot(
    fit = surv_fit(Surv(time, status) ~ has.smile.ct,
        data=subset(data, xor(has.chop.ct, has.smile.ct))),
    #pval = TRUE, conf.int = TRUE,
    xlab = "Months",
    xlim = c(1,23),
    break.x.by = 3,
    ylab = "Overall Survival Probability",
    risk.table=TRUE)

km.plots[[1]] + background_grid()
ggsave('./fig/general/kaplan_meier_chop_smile.png', height=8.5, width=8)

survdiff(Surv(time, status) ~ has.smile.ct,
    data=subset(data, xor(has.chop.ct, has.smile.ct)))

surv_median(surv_fit(Surv(time, status) ~ has.smile.ct,
    data=subset(data, xor(has.chop.ct, has.smile.ct))))

surv_median(
    surv_fit(Surv(time, status) ~ has.smile.ct,
    data=data))
surv_median(
    surv_fit(Surv(time, status) ~ has.chop.ct,
    data=data))


#####
# Collate Data for Table 1
#####

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
    age.to.bin.map.tab1[[data$Age]],
    levels=c('<40','40-64', '65+'),
    ordered=TRUE)


table(data$age.tab1)

table(data$age.tab1) / sum(table(data$age.tab1))


# Sex
table(data$Sex)
table(data$Sex) / sum(table(data$Sex))

t.test(data$Age[data$Sex == 'M'],
       data$Age[data$Sex == 'F'],
       paired=FALSE)

# Region

# Symptoms
table(data$has.vision.sx)
table(data$has.eom.sx)
table(data$has.lid.sx)
table(data$has.ptosis)
table(data$has.proptosis)

# Tumor Location
table(data$has.orbital)
table(data$has.intraocular)
table(data$has.lacrimal)
table(data$has.conjunctival)

table(data$has.nasosinus)
table(data$has.cns)

# Therapy
table(data$Ki67 < 0.2)

table(data$Ki67 >= 0.2 & data$Ki67 < 0.8)

table(data$Ki67 >= 0.8)

#######
## Multivariate Cox Regression
#######
data$is.asia <- data$nktl.region.rate == 'Asia'

res.cox <- coxph(
    formula=Surv(time, status) ~
        age.tab1 + # Demo
        Sex +
        is.asia +
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

print("===== all specific CT regimens =====")
res.cox <- coxph(
    formula=Surv(time, status) ~
        age.tab1 + # Demo
        Sex +
        is.asia +
        has.surgical + # Tx
        has.chop.ct +
        has.methotrexate.ct +
        has.smile.ct +
        has.dhap.ct +
        has.ceop.ct +
        has.devic.ct +
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

res.cox <- coxph(
    formula=Surv(time, status) ~
        age.tab1 + # Demo
        Sex +
        is.asia +
        has.surgical + # Tx
        strata(has.ct) +
        has.rt +
        has.vision.sx + # SX
        has.eom.sx +
        has.lid.sx +
        has.ptosis +
        has.proptosis +
        has.orbital + # Loc
        has.intraocular +
        has.lacrimal +
        has.conjunctival +
        has.nasosinus +
        has.cns,
    data = data,
    method = 'exact')
summary(res.cox)

# repeat with age-stratification
res.cox <- coxph(
    formula=Surv(time, status) ~
        strata(age.tab1) +
        Sex +
        is.asia +
        has.surgical +
        has.ct +
        has.rt
        ,
    data = data,
    method = 'exact')
summary(res.cox)

# repeat with region-stratification
res.cox <- coxph(
    formula=Surv(time, status) ~
        Age +
        Sex +
        strata(is.asia) +
        has.surgical +
        has.ct +
        has.rt
        ,
    data = data,
    method = 'exact')
summary(res.cox)


# plot age distribution, segmented by region
ggplot(data, aes(x=Age, fill=nktl.region.rate)) +
    geom_density(alpha = 0.5)

## Immunophenotyping data

for (m in colnames(ihc.markers)) {
    print("###############")
    print(m)
    print(table(ihc.markers[,m]))
    print(table(ihc.markers[,m])/sum(table(ihc.markers[,m])))
    print(sum(table(ihc.markers[,m])))
    print(sum(table(ihc.markers[,m]))/64)
    print("###############")
}

table(data$ebv.sero.status) / sum(table(data$ebv.sero.status))

## gather data for management

ct.subset.cols <- c(
    'has.gelox.ct', 'has.dhap.ct', 'has.methotrexate.ct',
    'has.chop.ct', 'has.smile.ct', 'has.devic.ct'
)

for (m in ct.subset.cols) {
    print("###############")
    print(m)
    print(table(data[,m]))
    print(table(data[,m])/sum(table(data[,m])))
    print(sum(table(data[,m])))
    print(sum(table(data[,m]))/64)
    print("###############")
}
