#!/bin/bash
# Amy Campbell
# 06/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# Filter run 32 
###############

# filter samples
qiime feature-table filter-samples \

  --i-table /home/acampbe/IowaWoundData/MiSeqV1V3_32/table32.qza \
  --m-metadata-file /home/acampbe/GitHub/IowaWound/mappings/samples-to-keep-32.tsv \
  --o-filtered-table /home/acampbe/IowaWoundData/MiSeqV1V3_32/id-filtered-table32.qza

# filter sequences based on new	table
qiime feature-table filter-seqs \
  --i-data /home/acampbe/IowaWoundData/MiSeqV1V3_32/rep-seqs32.qza \
  --i-table /home/acampbe/IowaWoundData/MiSeqV1V3_32/id-filtered-table32.qza \
  --o-filtered-data /home/acampbe/IowaWoundData/MiSeqV1V3_32/filtered-rep-seqs32.qza

# Filter run 35
###############

# filter samples
qiime feature-table filter-samples \
  --i-table /home/acampbe/IowaWoundData/MiSeqV1V3_35/table35.qza \
  --m-metadata-file /home/acampbe/GitHub/IowaWound/mappings/samples-to-keep-35.tsv \
  --o-filtered-table /home/acampbe/IowaWoundData/MiSeqV1V3_35/id-filtered-table35.qza

# filter sequences based on new table
qiime feature-table filter-seqs \
  --i-data /home/acampbe/IowaWoundData/MiSeqV1V3_35/rep-seqs35.qza \
  --i-table /home/acampbe/IowaWoundData/MiSeqV1V3_35/id-filtered-table35.qza \
  --o-filtered-data /home/acampbe/IowaWoundData/MiSeqV1V3_35/filtered-rep-seqs35.qza

# Merge
########

qiime feature-table merge \
  --i-tables /home/acampbe/IowaWoundData/MiSeqV1V3_32/id-filtered-table32.qza \
  --i-tables /home/acampbe/IowaWoundData/MiSeqV1V3_35/id-filtered-table35.qza \
  --o-merged-table /home/acampbe/IowaWoundData/MergedRuns/mergedtable.qza


qiime feature-table merge-seqs \
  --i-data /home/acampbe/IowaWoundData/MiSeqV1V3_32/filtered-rep-seqs32.qza \
  --i-data /home/acampbe/IowaWoundData/MiSeqV1V3_35/filtered-rep-seqs35.qza \
  --o-merged-data /home/acampbe/IowaWoundData/MergedRuns/rep-seqs-merged.qza
