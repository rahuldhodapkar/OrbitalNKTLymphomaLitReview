# Makefile
#
# 	Easy-to-use scaffolding to run all analysis in a repeatable
#	fashion, with minimal documentation.
#
# @author Rahul Dhodapkar <rahuldhodapkar>
# @version 2022.03.10
#

.PHONY: setup merge-searches run

setup:
	./Dependencies.R

merge-searches:
	./src/MergeSearches.R

run:
	mkdir -p ./calc/log
	./src/BuildFigures.R > ./calc/log/output.txt 2> ./calc/log/warnings.txt
