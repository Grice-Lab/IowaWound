#!/bin/bash
# Amy Campbell
# Updated 11/2021

# dont need this because I'm on my local machine
# source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# denoise miseq32
#################

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux32.qza \
  --p-trunc-len-f 267 \
  --p-trunc-len-r 257 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table32.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs32.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32.qza

# denoise miseq35
###################

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux35.qza \
  --p-trunc-len-f 282 \
  --p-trunc-len-r 247 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table35.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs35.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35.qza
