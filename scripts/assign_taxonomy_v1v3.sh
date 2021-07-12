#!/bin/bash
# Amy Campbell
# 06/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env
qiime feature-classifier classify-sklearn \
  --i-classifier /home/acampbe/DownloadedDatabases/Silva97_V1V3/silva-97-27f-534r-classifier.qza \
  --i-reads /home/acampbe/IowaWoundData/MergedRuns/rep-seqs-cr-97.qza \
  --o-classification /home/acampbe/IowaWoundData/MergedRuns/taxonomy_v1v3.qza
