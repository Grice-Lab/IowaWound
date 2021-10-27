#!/bin/bash
# Amy Campbell
# 06/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env


# Make tree for downstream diversity analyses 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /home/acampbe/IowaWoundData/MergedRuns/rep-seqs-merged.qza \
  --o-alignment /home/acampbe/IowaWoundData/MergedRuns/aligned-rep-seqs-merged.qza \
  --o-masked-alignment /home/acampbe/IowaWoundData/MergedRuns/masked-aligned-rep-seqs-merged.qza \
  --o-tree /home/acampbe/IowaWoundData/MergedRuns/unrooted-tree-merged.qza \
  --o-rooted-tree /home/acampbe/IowaWoundData/MergedRuns/rooted-tree-merged.qza

