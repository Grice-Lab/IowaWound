# Amy Campbell
# 11/8/2021
# Analyzing # of reads in each run and how many made it through filtering 
# Also, generate metadata specifically for qiime2 diversity metrics 


setwd("~/Desktop/GriceLabGit/IowaWound/")
twocolor = c("#FFC20A", "#0C7BDC")

# Get necessary metadata
########################
patient_metadata = read.csv("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/figuring_out_metadata_5_21/GSWOUNDGRICE2015_20190221.csv")
controls_run = read.csv2("mappings/Control_Run_Info.tsv", header=T, sep=" ")
patient_mapping32 = read.csv("mappings/IA_woundpain_mapping_32_2021.csv")
patient_mapping35 = read.csv("mappings/IA_woundpain_mapping_35_2021.csv")


# Read in and format the denoising stats 
##########################################
denoise32 = read.csv2("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoisingStats32.tsv", sep="\t")
denoise32 = denoise32[-1, ]

denoise35 = read.csv2("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DenoisingStats35.tsv", sep="\t")
denoise35 = denoise35[-1, ]

sort(sapply(denoise32$merged, function(x) as.numeric(as.character(x))))
sort(sapply(denoise35$merged, function(x) as.numeric(as.character(x))))


# Merge the denoise stats with info about which samples are controls and which aren't
######################################################################################
samdata32 = denoise32 %>% select(sample.id, input, non.chimeric)
dim(samdata32)
colnames(samdata32) = c("SampleID", "TotalReads", "NonChimericFilteredReads")

samdata35 = denoise35 %>% select(sample.id, input, non.chimeric)
dim(samdata35)
colnames(samdata35) = c("SampleID", "TotalReads", "NonChimericFilteredReads")

samdata32  = samdata32 %>% left_join(controls_run, by="SampleID")
samdata35  = samdata35 %>% left_join(controls_run, by="SampleID")

# Add wound type and location
#############################
MappingSampleType32 = patient_mapping32 %>% select(X.SampleID, SubjectID)
MappingSampleType35 = patient_mapping35 %>% select(X.SampleID, SubjectID)

colnames(MappingSampleType32) = c("SampleID", "study_id")
colnames(MappingSampleType35) = c("SampleID", "study_id")

patient_metadata_subset = patient_metadata %>% select(study_id, woundloc, wound_type)
patient_metadata_subset$study_id = factor(patient_metadata_subset$study_id)
MappingSampleType32 = MappingSampleType32 %>% left_join(patient_metadata_subset, by="study_id")

# Make controls '0' for location and type 
MappingSampleType32$woundloc[is.na(MappingSampleType32$woundloc)] <- 0
MappingSampleType32$wound_type[is.na(MappingSampleType32$wound_type)] <- 0
MappingSampleType32$SampleID = factor(MappingSampleType32$SampleID)


MappingSampleType35 = MappingSampleType35 %>% left_join(patient_metadata_subset, by="study_id")

# Make controls '0' for location and type 
MappingSampleType35$woundloc[is.na(MappingSampleType35$woundloc)] <- 0
MappingSampleType35$wound_type[is.na(MappingSampleType35$wound_type)] <- 0
MappingSampleType35$SampleID = factor(MappingSampleType35$SampleID)

samdata32 = samdata32 %>% left_join(MappingSampleType32, by="SampleID")
samdata35 = samdata35 %>% left_join(MappingSampleType35, by="SampleID")

# Metadata for diversity analyses 
totalmetadata = rbind(samdata32, samdata35)
totalmetadata= totalmetadata %>% select(SampleID, ControlStatus, Run, study_id, woundloc, wound_type)
colnames(totalmetadata)  = c("sampleid", "controlstatus", "run", "studyid", "woundloc", "woundtype")
write.table(totalmetadata, file="/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv", sep = "\t", row.names = F)



melted32 = samdata32 %>% select(SampleID, TotalReads, NonChimericFilteredReads, wound_type) %>% reshape2::melt(id.vars=c("SampleID", "wound_type"))
melted32$value = sapply(melted32$value, function(x) as.numeric(as.character(x)))
readplot32 = ggplot(melted32, aes(x=SampleID, y=value, fill=variable))+ geom_bar(aes(fill=variable), stat="identity")
readplot32order = melted32 %>% filter(variable=="NonChimericFilteredReads") %>% arrange(wound_type, value)
readplot32order$wound_type

readplot32$data$SampleID = factor(readplot32$data$SampleID, levels=readplot32order$SampleID)
readplot32= readplot32 + theme_classic() + theme(axis.text.x=element_blank(), legend.position="none") + scale_fill_manual(values=twocolor) + ylab("Read Count") + ylim(c(0, 240000)) + ggtitle("Total & Filtered Non-chimeric Reads in MiSeqV1V3_32") + geom_abline(slope=0, intercept=1200, color="red") 


melted35 = samdata35 %>% select(SampleID, TotalReads, NonChimericFilteredReads, wound_type) %>% reshape2::melt(id.vars=c("SampleID", "wound_type"))
melted35$value = sapply(melted35$value, function(x) as.numeric(as.character(x)))
readplot35 = ggplot(melted35, aes(x=SampleID, y=value, fill=variable))+ geom_bar(aes(fill=variable), stat="identity")
readplot35order = melted35 %>% filter(variable=="NonChimericFilteredReads") %>% arrange(wound_type, value)
readplot35$data$SampleID = factor(readplot35$data$SampleID, levels=readplot35order$SampleID)
readplot35order$wound_type
readplot35= readplot35 + theme_classic() + theme(axis.text.x=element_blank(), legend.key.size = unit(1, 'cm'), 
                                                 legend.key.height = unit(1, 'cm'), 
                                                 legend.key.width = unit(1, 'cm'), #
                                                 legend.title = element_text(size=16), 
                                                 legend.text = element_text(size=14), axis.title.y=element_blank()) + scale_fill_manual(values=twocolor) + ylim(c(0, 240000)) + ggtitle("Total & Filtered Non-chimeric Reads in MiSeqV1V3_35") + geom_abline(slope=0, intercept=1200, color="red") #+ annotate("text", label="1200", color="red", x=100, y=-1200)

#annotate("text", label="Phred Cutoff: 23", color="blue", x= 100000, y=5000) + annotate("text", label="Phred Cutoff: 20", color="red", x= 100000, y=10000)

ggsave(gridExtra::grid.arrange(readplot32, readplot35, ncol=2), file="/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/FilteredReads.pdf", height=15, width=25)


