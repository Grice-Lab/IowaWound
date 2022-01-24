#!/bin/bash
# Amy Campbell
# Updated 11/2021


# Make tree for downstream diversity analyses 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs-merged.qza \
  --o-alignment /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/aligned-rep-seqs-merged.qza \
  --o-masked-alignment /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/masked-aligned-rep-seqs-merged.qza \
  --o-tree /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/unrooted-tree-merged.qza \
  --o-rooted-tree /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rooted-tree-merged.qza

