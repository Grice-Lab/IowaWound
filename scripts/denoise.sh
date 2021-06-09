#!/bin/bash
# Amy Campbell
# 05/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# denoise miseq32
#################

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /home/acampbe/IowaWoundData/MiSeqV1V3_32/paired-end-demux32.qza \
  --p-trunc-len-f 272 \
  --p-trunc-len-r 257 \
  --o-table /home/acampbe/IowaWoundData/MiSeqV1V3_32/table32.qza \
  --o-representative-sequences /home/acampbe/IowaWoundData/MiSeqV1V3_32/rep-seqs32.qza \
  --o-denoising-stats /home/acampbe/IowaWoundData/MiSeqV1V3_32/denoising-stats32.qza

# denoise miseq35
###################

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /home/acampbe/IowaWoundData/MiSeqV1V3_35/paired-end-demux35.qza \
  --p-trunc-len-f 262 \
  --p-trunc-len-r 232 \
  --o-table /home/acampbe/IowaWoundData/MiSeqV1V3_35/table35.qza \
  --o-representative-sequences /home/acampbe/IowaWoundData/MiSeqV1V3_35/rep-seqs35.qza \
  --o-denoising-stats /home/acampbe/IowaWoundData/MiSeqV1V3_35/denoising-stats35.qza
