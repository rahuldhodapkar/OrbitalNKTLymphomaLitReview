#!/usr/bin/env Rscript
#
## MergeSearches.R
#
# Merge data from searches to identify the number of duplicates
#

library(stringr)

######## Create output scaffolding
ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)

#################################################
## LOAD DATA
#################################################

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

print(paste0("Initial Search Results: ", dim(search.data)[[1]]))

dedup.data <- data.frame(
    PMID = unique(search.data$PMID)
)

dedup.data$Title <-            search.data$Title[match(dedup.data$PMID, search.data$PMID)]
dedup.data$Publication.Year <- search.data$Publication.Year[match(dedup.data$PMID, search.data$PMID)]
dedup.data$Authors <-          search.data$Authors[match(dedup.data$PMID, search.data$PMID)]

write.csv(dedup.data, './calc/deduplicated_searches.csv', row.names=F)

print(paste0("Deduplicated Search Results: ", dim(dedup.data)[[1]]))

print("All done!")
