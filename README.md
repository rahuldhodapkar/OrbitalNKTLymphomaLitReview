# NK/T cell lymphoma of the orbit and adnexa
Literature review of NK/T cell lymphoma involving the orbit and adnexa

## Inclusion Criteria

- Pubmed indexed
- English language
- Case description of NK/T cell lymphoma involving the orbit and adnexa

## Exclusion Criteria

- Did not have per-patient reporting of age, sex, and particular orbital invovement
- Did not have description of treatment
- Did not have description of outcome

## Methods

Using the PRISMA guidelines, searches were performed on pubmed on 9/3/21 using a predetermined
set of search criteria. Articles were evaluated manually on the basis of title and abstract as to
whether they fit the inclusion and exclusion criteria. Only primary articles describing cases were
included, if articles were literature reviews or summaries of previously described cases, citations
for the original article with case descriptions were extracted from the review references and
subsequently evaluated using the inclusion and exclusion criteria above.

## Running Code
All scripts are designed to be run from the root directory of this repository. To install
dependencies, invoke:

    make setup

To perform search term deduplication (before manual literature review stage), invoke:

    make merge-searches

To build figures from systematic review data, invoke:

    make run
