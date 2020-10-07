#!/bin/bash
# Amy Campbell
# Iowa Wound Pipeline 2020 

#################################
# Step 1:Demultiplexing raw reads
#################################
# (Qiime 1, due to dual indexing)

conda activate qiime1

# Set path variables with positional arguments
##############################################
# $1 should be path to folder containing raw data
# (ex: "/Users/amycampbell/Iowa/Iowa2020/MiSeqV1V3_32_raw/")

# $2 should be a string you want to use to name the barcodes file
# (ex: "MiSeqV1V3_32")
  
Raw=$1
String=$2

# Preprocess & Demultiplex MiSeqV1V3_32
#######################################
extract_barcodes.py --input_type barcode_paired_end -f "${Raw}Undetermined_S0_L001_I1_001.fastq.gz" -r "${Raw}Undetermined_S0_L001_I2_001.fastq.gz" --bc1_len 12 --bc2_len 12 -o "${Raw}parsed_barcodes/"

mv "${Raw}parsed_barcodes/barcodes.fastq" "${Raw}$String"





