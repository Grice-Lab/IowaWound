#!/bin/bash
# Amy Campbell
# Iowa Wound Pipeline 2020 

#################################
# Step 1:Demultiplexing raw reads
#################################
# (Qiime 1, due to dual indexing)


conda activate qiime1

# Set path variables
####################
Raw32="/Users/amycampbell/Iowa/Iowa2020/MiSeqV1V3_32_raw"


# Preprocess & Demultiplex MiSeqV1V3_32
#######################################
extract_barcodes.py --input_type barcode_paired_end -f "${Raw32}Undetermined_S0_L001_I1_001.fastq.gz" -r "${Raw32}Undetermined_S0_L001_I2_001.fastq.gz" --bc1_len 12 --bc2_len 12 -o "${Raw32}parsed_barcodes/"

mv "${Raw32}parsed_barcodes/barcodes.fastq" "${Raw32}parsed_barcodes/MiSeqV1V3_32_barcodes.fastq"


