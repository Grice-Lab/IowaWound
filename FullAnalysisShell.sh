# Amy Campbell
# 09/2020 repeat Iowa Wound microbiome analyses but using Qiime (standardized, equivalent database for both runs). 

# 1. Preprocessing reads & indices using Qiime1 
###############################################
conda activate qiime1

#Format dual-indexed reads in Qiime2-friendly format

# MiSeqV1V3_32
##############
#sh ExtractBars.sh /Users/amycampbell/Iowa/Iowa2020/MiSeqV1V3_32_raw/ MiSeqV1V3_32

# MiSeqV1V3_35
sh ExtractBars.sh /Users/amycampbell/Iowa/Iowa2020/MiSeqV1V3_35_raw/ MiSeqV1V3_35

# 2. Demultiplexing & filtering reads using Qiime2
##################################################



# Linker/primer sequences for removal
# ***********************************
# 32:  TAGGTAATTGT
# 35:  GTAGAGTTTGATCCTGGCTCAG 



