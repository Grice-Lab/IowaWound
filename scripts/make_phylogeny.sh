#!/bin/bash
# Amy Campbell
# 06/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# Run 32 
#########

# Make tree for downstream diversity analyses 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /home/acampbe/IowaWoundData/MiSeqV1V3_32/rep-seqs32.qza \
  --o-alignment /home/acampbe/IowaWoundData/MiSeqV1V3_32/aligned-rep-seqs32.qza \
  --o-masked-alignment /home/acampbe/IowaWoundData/MiSeqV1V3_32/masked-aligned-rep-seqs32.qza \
  --o-tree /home/acampbe/IowaWoundData/MiSeqV1V3_32/unrooted-tree32.qza \
  --o-rooted-tree /home/acampbe/IowaWoundData/MiSeqV1V3_32/rooted-tree32.qza


# Run 35
#########
# Make tree for downstream diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /home/acampbe/IowaWoundData/MiSeqV1V3_35/rep-seqs35.qza \
  --o-alignment /home/acampbe/IowaWoundData/MiSeqV1V3_35/aligned-rep-seqs35.qza \
  --o-masked-alignment /home/acampbe/IowaWoundData/MiSeqV1V3_35/masked-aligned-rep-seqs35.qza \
  --o-tree /home/acampbe/IowaWoundData/MiSeqV1V3_35/unrooted-tree35.qza \
  --o-rooted-tree /home/acampbe/IowaWoundData/MiSeqV1V3_35/rooted-tree35.qza

