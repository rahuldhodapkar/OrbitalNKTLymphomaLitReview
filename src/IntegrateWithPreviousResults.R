#!/usr/bin/env Rscript
#
## IntegrateWithPreviousResults.R
#
# Integrate deduplicated_searches with previous work
#

library(stringr)

######## Create output scaffolding
ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)

#################################################
## LOAD DATA
#################################################

searches <- read.csv('./calc/deduplicated_searches.csv')
previous.searches <- read.csv('./data/202110_deduplicated_searches.csv')


searches$ExclusionReason <- previous.searches$ExclusionReason[
    match(searches$PMID, previous.searches$PMID)]

# Print Merge Summary
print('Number of Overlap Searches from Previous:')
sum(searches$PMID %in% previous.searches$PMID, na.rm=T)

print('Number of Old Searches not in New List:')
sum(! (previous.searches$PMID %in% searches$PMID), na.rm=T)

print('PMIDs for searches not in New List:')
previous.searches$PMID[
    !(previous.searches$PMID %in% searches$PMID)]

write.csv(searches, './calc/annotated_deduplicated_searches.csv', row.names=F)

print('All done!')
