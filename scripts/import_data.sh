#!/bin/bash
# Amy Campbell
# 05/2021
# import demultiplexed data into qiime using manifest files 
source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# import data for MiSeqV1V3_32 
###############################

# set manifest path 
manifest=/home/acampbe/GitHub/IowaWound/mappings/Manifest32.tsv

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $manifest \
  --output-path /home/acampbe/IowaWoundData/MiSeqV1V3_32/paired-end-demux32.qza \
  --input-format PairedEndFastqManifestPhred33V2 


# import data for MiSeqV1V3_35
##############################
# set manifest path
manifest=/home/acampbe/GitHub/IowaWound/mappings/Manifest35.tsv

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $manifest \
  --output-path /home/acampbe/IowaWoundData/MiSeqV1V3_35/paired-end-demux35.qza \
  --input-format PairedEndFastqManifestPhred33V2



