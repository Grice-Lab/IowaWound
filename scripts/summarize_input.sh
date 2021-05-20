#!/bin/bash
# Amy Campbell
# 05/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# summarize miseq 32
####################
qiime demux summarize \
  --i-data /home/acampbe/IowaWoundData/MiSeqV1V3_32/paired-end-demux32.qza \
  --o-visualization /home/acampbe/IowaWoundData/MiSeqV1V3_32/summary-paired-end-demux32.qzv

# summarize miseq35
###################
qiime demux summarize \
  --i-data /home/acampbe/IowaWoundData/MiSeqV1V3_35/paired-end-demux35.qza \
  --o-visualization /home/acampbe/IowaWoundData/MiSeqV1V3_35/summary-paired-end-demux35.qzv

