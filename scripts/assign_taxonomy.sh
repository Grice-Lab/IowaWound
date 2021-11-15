#!/bin/bash
# Amy Campbell
# Upated 11/2021

qiime feature-classifier classify-sklearn \
  --i-classifier /Users/amycampbell/Documents/IowaWoundData2021/databases/silva-138-99-nb-classifier_1.qza \
  --i-reads /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs-merged.qza \
  --o-classification /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/taxonomy.qza

qiime metadata tabulate \
  --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/taxonomy.qza \
  --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/taxonomy.qzv
