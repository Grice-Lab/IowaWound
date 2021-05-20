#!/bin/bash
# Amy Campbell
# 5/2021
# Checking demultiplexed reads prior to denoising with Dada2

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

r1string="_1.fastq.gz"
r2string="_2.fastq.gz"

READSPATH35=/home/acampbe/IowaWoundData/MiSeqV1V3_35/demultiplexed
OUTPUT35=/home/acampbe/IowaWoundData/MiSeqV1V3_32/fastqc

READSPATH32=/home/acampbe/IowaWoundData/MiSeqV1V3_32/demultiplexed
OUTPUT32=/home/acampbe/IowaWoundData/MiSeqV1V3_32/fastqc

# make output directories for the fastqc
# in each data path  
mkdir -p $OUTPUT32
mkdir -p $OUTPUT35

MiSeqV1V3_35
#############

for filename in $READSPATH35/*_1.fastq.gz; do
	run1=$filename	
	run2=${run1/$r1string/$r2string}	
	fastqc -o $OUTPUT35 $run1 -f fastq
	fastqc -o $OUTPUT35 $run2 -f fastq

done		
	
multiqc -n MultiQC_PostDemuxed35 $OUTPUT35

# MiSeqV1V3_32
###############

for filename in $READSPATH32/*_1.fastq.gz; do
        run1=$filename	
        run2=${run1/$r1string/$r2string}
        fastqc -o $OUTPUT32 $run1 -f fastq
        fastqc -o $OUTPUT32 $run2 -f fastq

done

multiqc -n MultiQC_PostDemuxed32 $OUTPUT32 

