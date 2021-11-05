# Comparing old (phred score cutoff 25) and new (phred score cutoff 20) stats

# Amy Campbell
# Output of Qiime2 
# July 2021

library("dplyr")
library("stringr")
library("ggplot2")

setwd("~/Desktop/GriceLabGit/IowaWound/")

Old32Stats = read.table("/Users/amycampbell/Documents/IowaWoundData2021/DenoiseStatsOLD32.tsv", skip=1)
Old35Stats = read.table("/Users/amycampbell/Documents/IowaWoundData2021/DenoiseStatsOLD35.tsv", skip=1)

colnames(Old32Stats) = c("Sample", "input", "filtered", "pct_input_pass_filter", "denoised", "merged", "pct_input_merged", "nonchimeric", "pct_nonchimeric")
colnames(Old35Stats) = c("Sample", "input", "filtered", "pct_input_pass_filter", "denoised", "merged", "pct_input_merged", "nonchimeric", "pct_nonchimeric")

New32Stats = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoisingStats32.tsv", skip=1)
New35Stats = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoisingStats35.tsv", skip=1)
colnames(New32Stats) = c("Sample", "input_new", "filtered_new", "pct_input_pass_filter_new", "denoised_new", "merged_new", "pct_input_merged_new", "nonchimeric_new", "pct_nonchimeric_new")
colnames(New35Stats) = c("Sample", "input_new", "filtered_new", "pct_input_pass_filter_new", "denoised_new", "merged_new", "pct_input_merged_new", "nonchimeric_new", "pct_nonchimeric_new")


Stats32 = Old32Stats %>% left_join(New32Stats, by="Sample")
Stats35 = Old35Stats %>% left_join(New35Stats, by="Sample")

ggplot(Stats32,aes(x=nonchimeric, y=nonchimeric_new)) + geom_point() +geom_abline(slope=1) + ggtitle("MiSeqV1V3_32 Filtered & Merged Nonchimeric Read Counts") + xlab("Old Run") + ylab("New Run")
ggplot(Stats35,aes(x=nonchimeric, y=nonchimeric_new)) + geom_point() +geom_abline(slope=1) + ggtitle("MiSeqV1V3_35 Filtered & Merged Nonchimeric Read Counts") + xlab("Old Run") + ylab("New Run")


ggplot(Stats32,aes(x=merged, y=merged_new)) + geom_point() +geom_abline(slope=1) + ggtitle("MiSeqV1V3_32 Filtered & Merged  Read Counts") + xlab("Old Run") + ylab("New Run")
ggplot(Stats35,aes(x=merged, y=merged_new)) + geom_point() +geom_abline(slope=1) + ggtitle("MiSeqV1V3_35 Filtered & Merged Read Counts") + xlab("Old Run") + ylab("New Run")



Miseq32_27 = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoiseTest/MiSeqV1V3_32_denoising27.tsv", skip=1)
colnames(Miseq32_27) = c("Sample", "input27", "filtered27", "pct_input_pass_filter27", "denoised27", "merged27", "pct_input_merged27", "nonchimeric27", "pct_nonchimeric27")

Miseq32_25 = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoiseTest/MiSeqV1V3_32_denoising25.tsv", skip=1)
colnames(Miseq32_25) = c("Sample", "input25", "filtered25", "pct_input_pass_filter25", "denoised25", "merged25", "pct_input_merged25", "nonchimeric25", "pct_nonchimeric25")

Miseq32_23 = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoiseTest/MiSeqV1V3_32_denoising23.tsv", skip=1)
colnames(Miseq32_23) = c("Sample", "input23", "filtered23", "pct_input_pass_filter23", "denoised23", "merged23", "pct_input_merged23", "nonchimeric23", "pct_nonchimeric23")

Miseq32_20 = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoiseTest/MiSeqV1V3_32_denoising20.tsv", skip=1)
colnames(Miseq32_20) = c("Sample", "input20", "filtered20", "pct_input_pass_filter20", "denoised20", "merged20", "pct_input_merged20", "nonchimeric20", "pct_nonchimeric20")


MiSeq32Compare = Miseq32_25 %>% left_join(Miseq32_23, by="Sample")
MiSeq32Compare = MiSeq32Compare %>% left_join(Miseq32_20, by="Sample")
MiSeq32Compare = MiSeq32Compare %>% left_join(Miseq32_27, by="Sample")

m32plot = ggplot(MiSeq32Compare, aes(x=nonchimeric25,y=nonchimeric20)) + geom_point() + geom_point(data=MiSeq32Compare, aes(x=nonchimeric25,y=nonchimeric20), color="red") +
  geom_point(data=MiSeq32Compare, aes(x=nonchimeric25, y=nonchimeric23), color="blue") + 
  geom_abline(slope=1) + xlab("# Nonchimeric Reads in Sample (Phred 25)") + ylab("Nonchimeric Reads in Sample (Phred 20, 23, or 27)") +
  geom_point(data=MiSeq32Compare, aes(x=nonchimeric25, y=nonchimeric27), color="yellow")  +
  geom_label(aes(x = 10000, y = 110000, label = "Phred Cutoff:27"), fill = "darkgray", color="yellow") + 
  geom_label(aes(x = 10000, y = 120000, label = "Phred Cutoff:23"), fill = "darkgray", color="blue") + 
  geom_label(aes(x = 10000, y = 130000, label = "Phred Cutoff:20"), fill = "darkgray", color="red")+ggtitle("MiSeqV1V3_32")
  
ggsave(m32plot, file="PhredCutoffs32.png", height=10, width=10)
#+ annotate("text", label="No correction", color="red", x=1, y=4.8) + annotate("text", label="Correction", color="black", x= 2, y=4) +
 # xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") + ggtitle("QQ Plot for GWAS p-values") + annotate("text", label="Actual full set of p-values", color="blue", x=3, y=2)+
  #geom_point(data=bugwasDF_noCorrection, aes(x=expect, y=observe), color="green")




Miseq35_27 = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoiseTest/MiSeqV1V3_35_denoising27.tsv", skip=1)
colnames(Miseq35_27) = c("Sample", "input27", "filtered27", "pct_input_pass_filter27", "denoised27", "merged27", "pct_input_merged27", "nonchimeric27", "pct_nonchimeric27")


Miseq35_25 = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoiseTest/MiSeqV1V3_35_denoising25.tsv", skip=1)
colnames(Miseq35_25) = c("Sample", "input25", "filtered25", "pct_input_pass_filter25", "denoised25", "merged25", "pct_input_merged25", "nonchimeric25", "pct_nonchimeric25")

Miseq35_23 = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoiseTest/MiSeqV1V3_35_denoising23.tsv", skip=1)
colnames(Miseq35_23) = c("Sample", "input23", "filtered23", "pct_input_pass_filter23", "denoised23", "merged23", "pct_input_merged23", "nonchimeric23", "pct_nonchimeric23")

Miseq35_20 = read.table("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoiseTest/MiSeqV1V3_35_denoising20.tsv", skip=1)
colnames(Miseq35_20) = c("Sample", "input20", "filtered20", "pct_input_pass_filter20", "denoised20", "merged20", "pct_input_merged20", "nonchimeric20", "pct_nonchimeric20")


Miseq35Compare = Miseq35_25 %>% left_join(Miseq35_23, by="Sample")
Miseq35Compare = Miseq35Compare %>% left_join(Miseq35_20, by="Sample")
Miseq35Compare = Miseq35Compare %>% left_join(Miseq35_27, by="Sample")

ggplot(Miseq35Compare, aes(x=nonchimeric25,y=nonchimeric20)) + geom_point() + geom_point(data=Miseq35Compare, aes(x=nonchimeric25,y=nonchimeric20), color="red") +
  geom_point(data=Miseq35Compare, aes(x=nonchimeric25, y=nonchimeric23), color="blue") + 
  geom_abline(slope=1) + xlab("# Nonchimeric Reads in Sample (Phred 25)") + ylab("Nonchimeric Reads in Sample (Phred 20 or 23)") +  annotate("text", label="Phred Cutoff: 23", color="blue", x= 100000, y=5000) + annotate("text", label="Phred Cutoff: 20", color="red", x= 100000, y=10000) + ggtitle("MiseqV1V3_35")

m35plot = ggplot(Miseq35Compare, aes(x=nonchimeric25,y=nonchimeric20)) + geom_point() + geom_point(data=Miseq35Compare, aes(x=nonchimeric25,y=nonchimeric20), color="red") +
  geom_point(data=Miseq35Compare, aes(x=nonchimeric25, y=nonchimeric23), color="blue") + 
  geom_abline(slope=1) + xlab("# Nonchimeric Reads in Sample (Phred 25)") + ylab("Nonchimeric Reads in Sample (Phred 20, 23, or 27)") +
  geom_point(data=MiSeq32Compare, aes(x=nonchimeric25, y=nonchimeric27), color="yellow")  +
  geom_label(aes(x = 10000, y = 110000, label = "Phred Cutoff:27"), fill = "darkgray", color="yellow") + 
  geom_label(aes(x = 10000, y = 120000, label = "Phred Cutoff:23"), fill = "darkgray", color="blue") + 
  geom_label(aes(x = 10000, y = 130000, label = "Phred Cutoff:20"), fill = "darkgray", color="red") +ggtitle("MiSeqV1V3_35")
ggsave(m35plot, file="PhredCutoffs35.png", height=10, width=10)

