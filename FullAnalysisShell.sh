# Amy Campbell
# 2/2021 repeat Iowa Wound microbiome analyses but using Qiime (standardized, equivalent database for both runs). 

###############################################
# 1. Preprocessing reads & indices using Qiime1 
###############################################
# In this section, I format dual-indexed reads in Qiime2-friendly format

# Activate virtual environment
conda activate qiime1

#-------------
# MiSeqV1V3_32
#-------------

# Set path to files for this run
MiSeqV1V3_32_Path="data/RawIntermediateFiles/MiSeqV1V3_32_raw/"
# using this : "${Raw}Undetermined_S0_L001_I1_001.fastq.gz"

# Decompress reads
gunzip "${MiSeqV1V3_32_Path}Undetermined_S0_L001_R1_001.fastq.gz"
gunzip "${MiSeqV1V3_32_Path}Undetermined_S0_L001_R2_001.fastq.gz"

# Reformat barcodes
perl -i.bak -pe 'if (m/^\@.*?\s+.*?\+.*?$/){s/\+//;}' "${MiSeqV1V3_32_Path}Undetermined_S0_L001_R1_001.fastq" "${MiSeqV1V3_32_Path}Undetermined_S0_L001_R2_001.fastq"

# Extract the barcodes from the reads for Qiime2 input
extract_barcodes.py -f "${MiSeqV1V3_32_Path}Undetermined_S0_L001_R1_001.fastq" -r "${MiSeqV1V3_32_Path}Undetermined_S0_L001_R2_001.fastq" -l 16 -o "${MiSeqV1V3_32_Path}parsed_barcodes" -c barcode_in_label

# Move all the reads and barcodes to a folder called "Input32" 
mkdir -p "${MiSeqV1V3_32_Path}Input32"
mv "${MiSeqV1V3_32_Path}parsed_barcodes/barcodes.fastq" "${MiSeqV1V3_32_Path}Input32/barcodes.fastq"
mv "${MiSeqV1V3_32_Path}Undetermined_S0_L001_R1_001.fastq" "${MiSeqV1V3_32_Path}Input32/forward.fastq"
mv "${MiSeqV1V3_32_Path}Undetermined_S0_L001_R2_001.fastq" "${MiSeqV1V3_32_Path}Input32/reverse.fastq"

# Compress the reads and barcodes for input into Qiime2 
gzip "${MiSeqV1V3_32_Path}Input32/barcodes.fastq"
gzip "${MiSeqV1V3_32_Path}Input32/forward.fastq"
gzip "${MiSeqV1V3_32_Path}Input32/reverse.fastq"


#-------------
# MiSeqV1V3_35
#-------------
MiSeqV1V3_35_Path="data/RawIntermediateFiles/MiseqV1V3_35_raw/"

# Decompress reads
gunzip "${MiSeqV1V3_35_Path}Undetermined_S0_L001_R1_001.fastq.gz"
gunzip "${MiSeqV1V3_35_Path}Undetermined_S0_L001_R2_001.fastq.gz"

# Reformat barcodes
perl -i.bak -pe 'if (m/^\@.*?\s+.*?\+.*?$/){s/\+//;}' "${MiSeqV1V3_35_Path}Undetermined_S0_L001_R1_001.fastq" "${MiSeqV1V3_35_Path}Undetermined_S0_L001_R2_001.fastq"

# Extract the barcodes from the reads for Qiime2 input
extract_barcodes.py -f "${MiSeqV1V3_35_Path}Undetermined_S0_L001_R1_001.fastq" -r "${MiSeqV1V3_35_Path}Undetermined_S0_L001_R2_001.fastq" -l 16 -o "${MiSeqV1V3_35_Path}parsed_$


# Move all the reads and barcodes to a folder called "Input35"
mkdir -p data/RawIntermediateFiles/MiSeqV1V3_35_raw/Input35
mv "${MiSeqV1V3_35_Path}Undetermined_S0_L001_R1_001.fastq" "${MiSeqV1V3_35_Path}Input35/forward.fastq"
mv "${MiSeqV1V3_35_Path}Undetermined_S0_L001_R2_001.fastq" "${MiSeqV1V3_35_Path}Input35/reverse.fastq"
mv "${MiSeqV1V3_35_Path}parsed_barcodes/barcodes.fastq" "${MiSeqV1V3_35_Path}Input35/barcodes.fastq"

gzip "${MiSeqV1V3_35_Path}Input35/forward.fastq"
gzip "${MiSeqV1V3_35_Path}Input35/reverse.fastq"
gzip "${MiSeqV1V3_35_Path}Input35/barcodes.fastq"

conda deactivate

##################################################
# 2. Demultiplexing & filtering reads using Qiime2
##################################################
conda activate qiime2

# Inspect metadata files
qiime tools inspect-metadata data/RawIntermediateFiles/Iowa35SampleMetadata.tsv 
qiime tools inspect-metadata data/RawIntermediateFiles/Iowa32SampleMetadata.tsv 

qiime metadata tabulate --m-input-file data/RawIntermediateFiles/Iowa32SampleMetadata.tsv --o-visualization data/RawIntermediateFiles/Iowa32SampleMetadata.tsv
qiime metadata tabulate --m-input-file data/RawIntermediateFiles/Iowa35SampleMetadata.tsv --o-visualization data/RawIntermediateFiles/Iowa35SampleMetadata.tsv



#--------------
# MiSeqV1V3_32
#-------------
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path "${MiSeqV1V3_32_Path}Input32/" \
  --output-path "${MiSeqV1V3_32_Path}emp-paired-end-sequences.qza"

qiime demux emp-paired \
  --i-seqs "${MiSeqV1V3_32_Path}emp-paired-end-sequences.qza" \
  --m-barcodes-file data/RawIntermediateFiles/MiSeqV1V3_32_FullMetadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences "${MiSeqV1V3_32_Path}demux.qza" \
  --p-no-golay-error-correction \
  --o-error-correction-details "{MiSeqV1V3_32_Path}demux-details.qza"

#--------------
# MiSeqV1V3_35
#-------------
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path data/RawIntermediateFiles/MiSeqV1V3_35_raw/Input32/ \
  --output-path data/RawIntermediateFiles/MiSeqV1V3_35_raw/emp-paired-end-sequences.qza

qiime demux emp-paired \
  --i-seqs "${MiSeqV1V3_35_Path}emp-paired-end-sequences.qza" \
  --m-barcodes-file data/RawIntermediateFiles/MiSeqV1V3_35_FullMetadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences "${MiSeqV1V3_35_Path}demux.qza" \
  --p-no-golay-error-correction \
  --o-error-correction-details "{MiSeqV1V3_35_Path}demux-details.qza"




