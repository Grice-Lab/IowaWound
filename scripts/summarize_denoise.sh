#!/bin/bash
# Amy Campbell
# 06/2021

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

# denoise miseq32
#################
  # Summary of feature table
qiime feature-table summarize \
  --i-table /home/acampbe/IowaWoundData/MiSeqV1V3_32/table32.qza \
  --o-visualization /home/acampbe/IowaWoundData/MiSeqV1V3_32/table32.qzv \
  --m-sample-metadata-file /home/acampbe/Club_Grice/mapping_files/run_maps/MiSeqV1V3_32.tsv

  # Summary of feature sequences
qiime feature-table tabulate-seqs \
  --i-data /home/acampbe/IowaWoundData/MiSeqV1V3_32/rep-seqs32.qza \
  --o-visualization /home/acampbe/IowaWoundData/MiSeqV1V3_32/rep-seqs32.qzv

  # Summary of denoising stats
qiime metadata tabulate \
   --m-input-file /home/acampbe/IowaWoundData/MiSeqV1V3_32/denoising-stats32.qza \
   --o-visualization /home/acampbe/IowaWoundData/MiSeqV1V3_32/denoising-stats32.qzv


# summarize denoising miseq35
#############################
  # Summary of feature table
qiime feature-table summarize \
  --i-table /home/acampbe/IowaWoundData/MiSeqV1V3_35/table35.qza \
  --o-visualization /home/acampbe/IowaWoundData/MiSeqV1V3_35/table35.qzv \
  --m-sample-metadata-file /home/acampbe/Club_Grice/mapping_files/run_maps/MiSeqV1V3_35.tsv

  # Summary of feature sequences
qiime feature-table tabulate-seqs \
  --i-data /home/acampbe/IowaWoundData/MiSeqV1V3_35/rep-seqs35.qza \
  --o-visualization /home/acampbe/IowaWoundData/MiSeqV1V3_35/rep-seqs35.qzv

  # Summary of denoising stats
qiime metadata tabulate \
   --m-input-file /home/acampbe/IowaWoundData/MiSeqV1V3_35/denoising-stats35.qza \
   --o-visualization /home/acampbe/IowaWoundData/MiSeqV1V3_35/denoising-stats35.qzv
