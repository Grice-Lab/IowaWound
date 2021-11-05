#!/bin/bash
# Amy Campbell
# Updated 11/2021
# Filter to included samples and controls; merge feature tables and sequence files for the two runs 

# Filter run 32 
###############

# filter samples
qiime feature-table filter-samples \
  --i-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table32.qza \
  --m-metadata-file /Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/samples-to-keep-32.tsv \
  --o-filtered-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/id-filtered-table32.qza

# filter sequences based on new	table
qiime feature-table filter-seqs \
  --i-data /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs32.qza \
  --i-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/id-filtered-table32.qza \
  --o-filtered-data /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/filtered-rep-seqs32.qza

# Filter run 35
###############

# filter samples
qiime feature-table filter-samples \
  --i-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table35.qza \
  --m-metadata-file /Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/samples-to-keep-35.tsv \
  --o-filtered-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/id-filtered-table35.qza

# filter sequences based on new table
qiime feature-table filter-seqs \
  --i-data /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs35.qza \
  --i-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/id-filtered-table35.qza \
  --o-filtered-data /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/filtered-rep-seqs35.qza

# Merge
########

qiime feature-table merge \
  --i-tables /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/id-filtered-table32.qza \
  --i-tables /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/id-filtered-table35.qza \
  --o-merged-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/mergedtable.qza

qiime feature-table merge-seqs \
  --i-data /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/filtered-rep-seqs32.qza \
  --i-data /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/filtered-rep-seqs35.qza \
  --o-merged-data /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs-merged.qza
