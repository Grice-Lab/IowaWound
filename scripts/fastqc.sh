#!/bin/bash
# Amy Campbell
# 5/2021
# Checking demultiplexed reads prior to denoising with Dada2

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env

r1string="_1.fastq.gz"
r2string="_2.fastq.gz"

READSPATH35=/home/acampbe/IowaWoundData/MiSeqV1V3_35/demultiplexed

READSPATH32=/home/acampbe/IowaWoundData/MiSeqV1V3_32/demultiplexed
OUTPUT=/home/acampbe/IowaWoundData/fastqc

# make output directories for the fastqc
# in each data path  
mkdir -p $OUTPUT

MiSeqV1V3_35
#############

for filename in $READSPATH35/*_1.fastq.gz; do
	run1=$filename	
	run2=${run1/$r1string/$r2string}	
	fastqc -o $OUTPUT $run1 -f fastq
	fastqc -o $OUTPUT $run2 -f fastq

done		
	

# MiSeqV1V3_32
###############

for filename in $READSPATH32/*_1.fastq.gz; do
        run1=$filename	
        run2=${run1/$r1string/$r2string}
        fastqc -o $OUTPUT $run1 -f fastq
        fastqc -o $OUTPUT $run2 -f fastq

done

multiqc -n MultiQC_PostDemuxed $OUTPUT
mv MultiQC_PostDemuxed /home/acampbe/IowaWoundData/MultiQC_PostDemuxed
