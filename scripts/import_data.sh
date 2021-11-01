#!/bin/bash
# Amy Campbell
# Updated 11/2021
# import demultiplexed data into qiime using manifest files 
source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# import data for MiSeqV1V3_32 
###############################

mkdir -p  /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/
# set manifest path 
manifest=/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/Manifest32.tsv
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $manifest \
  --output-path /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux32.qza \
  --input-format PairedEndFastqManifestPhred33V2 


# import data for MiSeqV1V3_35
##############################
# set manifest path
manifest=/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/Manifest35.tsv

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $manifest \
  --output-path /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/paired-end-demux35.qza \
  --input-format PairedEndFastqManifestPhred33V2

