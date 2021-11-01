#!/bin/bash
# Amy Campbell
# 5/2021
# Checking demultiplexed reads prior to denoising with Dada2


r1string="_1.fastq.gz"
r2string="_2.fastq.gz"

READSPATH35=/Users/amycampbell/Documents/IowaWoundData2021/MiSeqV1V3_35/demultiplexed

READSPATH32=/Users/amycampbell/Documents/IowaWoundData2021/MiSeqV1V3_32/demultiplexed
OUTPUT=/Users/amycampbell/Documents/IowaWoundData2021/fastqc

# make output directories for the fastqc
# in each data path  
mkdir -p $OUTPUT

#MiSeqV1V3_35
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
mv MultiQC_PostDemuxed /Users/amycampbell/Documents/IowaWoundData2021/MultiQC_PostDemuxed
