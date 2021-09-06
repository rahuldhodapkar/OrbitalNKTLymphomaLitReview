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

ggsurvplot(
    fit = surv_fit(Surv(time, status) ~ Sex, data=data),
    xlab = "Months",
    ylab = "Overall Survival Probability")

survdiff(Surv(time, status) ~ Sex, data=data)

summary(survfit(Surv(time, status) ~ 1, data = data), 
            times = c(1,6,12,24))


