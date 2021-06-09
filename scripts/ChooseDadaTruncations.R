# Amy Campbell
# 05/2021
# Choosing truncation lengths for DADA2 denoising step 
# based on FASTQC/MultiFastQC output

# Specifically, identify first position at which median quality length <25 for each of fwd and rev reads
# and for each run 

library(dplyr)
#setwd("/Users/amycampbell/Desktop/GriceLabGit/IowaWound")
# Data path
datapath = "/home/acampbe/IowaWoundData/MultiQC_PostDemuxed/"

# check sample name 
direction = function(x){
  x=as.character(x)
  substringobject = substring(x, (nchar(x) - 1), nchar(x))

    if( substringobject == "_1"){
    return("forward")
  }else{
    return("reverse")
  }
}

# Forward and reverse markings 
qualities_position = read.table(paste0(datapath, "mqc_fastqc_per_base_sequence_quality_plot_1.txt"), sep="\t", header=T)
qualities_position$Direction = sapply(qualities_position$Sample, function(x) direction(x))


qualities_fwd_32 = qualities_position %>% filter(grepl(x = Sample,pattern = "MiSeqV1V3_32")) %>% filter(Direction=="forward")
qualities_rev_32 = qualities_position %>% filter(grepl(x = Sample,pattern = "MiSeqV1V3_32")) %>% filter(Direction=="reverse")

qualities_fwd_35 = qualities_position %>% filter(grepl(x = Sample,pattern = "MiSeqV1V3_35")) %>% filter(Direction=="forward")
qualities_rev_35 = qualities_position %>% filter(grepl(x = Sample,pattern = "MiSeqV1V3_35")) %>% filter(Direction=="reverse")

#should add up to 1004 rows 
dim(qualities_fwd_32) + dim(qualities_rev_32) + dim(qualities_fwd_35) + dim(qualities_rev_35)


# Forward reads for run 32
qualities_fwd_32 = qualities_fwd_32 %>% select(-c(Sample, Direction))
apply(qualities_fwd_32, 2, median)
# forward read cutoff for run 32 should be 272 


qualities_rev_32 = qualities_rev_32 %>% select(-c(Sample, Direction))
apply(qualities_rev_32, 2, median)
# reverse read cutoff for run 32 should be 257


qualities_fwd_35 = qualities_fwd_35 %>% select(-c(Sample, Direction))
apply(qualities_fwd_35, 2, median)
# Forward read cutoff for 35 should be 262


qualities_rev_35 = qualities_rev_35 %>% select(-c(Sample, Direction))
apply(qualities_rev_35, 2, median)
# Reverse cutoff for run 35 should be 232 


