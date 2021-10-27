#!/bin/bash
# Amy Campbell
# 07/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

qiime vsearch cluster-features-closed-reference \
  --i-table /home/acampbe/IowaWoundData/MergedRuns/mergedtable.qza \
  --i-sequences /home/acampbe/IowaWoundData/MergedRuns/rep-seqs-merged.qza \
  --i-reference-sequences /home/acampbe/DownloadedDatabases/Silva97_V1V3/silva-97-seqs-27f-534r-lca.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table /home/acampbe/IowaWoundData/MergedRuns/table-cr-97.qza \
  --o-clustered-sequences /home/acampbe/IowaWoundData/MergedRuns/rep-seqs-cr-97.qza \
  --o-unmatched-sequences /home/acampbe/IowaWoundData/MergedRuns/unmatched-cr-97.qza
