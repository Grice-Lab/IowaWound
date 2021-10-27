#!/bin/bash
# Amy Campbell
# 06/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env
qiime feature-classifier classify-sklearn \
  --i-classifier /home/acampbe/DownloadedDatabases/silva-132-99-nb-classifier.qza \
  --i-reads /home/acampbe/IowaWoundData/MergedRuns/rep-seqs-merged.qza \
  --o-classification /home/acampbe/IowaWoundData/MergedRuns/taxonomy.qza
