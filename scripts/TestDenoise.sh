#!/bin/bash
# Amy Campbell
# Updated 11/2021

# dont need this because I'm on my local machine
# source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# denoise miseq32
######################

#Phred20 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux32.qza \
  --p-trunc-len-f 300 \
  --p-trunc-len-r 282 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table32_phred20.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs32_phred20.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred20.qza

qiime metadata tabulate \
   --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred20.qza \
   --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred20.qzv

#Phred23 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux32.qza \
  --p-trunc-len-f 282 \
  --p-trunc-len-r 262 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table32_phred23.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs32_phred23.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred23.qza

qiime metadata tabulate \
   --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred23.qza \
   --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred23.qzv

#Phred25 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux32.qza \
  --p-trunc-len-f 267 \
  --p-trunc-len-r 257 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table32_phred25.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs32_phred25.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred25.qza

qiime metadata tabulate \
   --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred25.qza \
   --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred25.qzv


#Phred27 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux32.qza \
  --p-trunc-len-f 257 \
  --p-trunc-len-r 242 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table32_phred27.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs32_phred27.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred27.qza

qiime metadata tabulate \
   --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred27.qza \
   --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats32_phred27.qzv


# denoise miseq35
###################

#Phred20 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux35.qza \
  --p-trunc-len-f 300 \
  --p-trunc-len-r 267 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table35_phred20.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs35_phred20.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred20.qza

qiime metadata tabulate \
   --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred20.qza \
   --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred20.qzv

# Phred23
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux35.qza \
  --p-trunc-len-f 282 \
  --p-trunc-len-r 247 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table35_phred23.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs35_phred23.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred23.qza

qiime metadata tabulate \
   --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred23.qza \
   --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred23.qzv


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux35.qza \
  --p-trunc-len-f 262 \
  --p-trunc-len-r 232 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table35_phred25.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs35_phred25.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred25.qza

qiime metadata tabulate \
   --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred25.qza \
   --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred25.qzv

# Phred27 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux35.qza \
  --p-trunc-len-f 257 \
  --p-trunc-len-r 217 \
  --o-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/table35_phred27.qza \
  --o-representative-sequences /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rep-seqs35_phred27.qza \
  --o-denoising-stats /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred27.qza

qiime metadata tabulate \
   --m-input-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred27.qza \
   --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/denoising-stats35_phred27.qzv