# Amy Campbell
# Downstream analyses comparing the cytokine and microbiome variables to the actual 
library("phyloseq")
library("dplyr")
library("stringr")
library("ggplot2")
library("ggpubr")
library("reshape2")

AddCLRDes = function(namestring){
  if(!(namestring %in% c("X", "SampleID", "StudyID", "DMMClusterAssign"))){
    return(paste0(namestring,"CLR"))
  }else{
    return(namestring)
  }
  
}

WoundMicrobiome = read.csv("~/Documents/IowaWoundData2021/PlotsForSue2022/WoundMicrobiomeDataForSEG_AEC.csv")
DominantGenera = read.csv("~/Documents/IowaWoundData2021/WoundAbundanceData_CLR_DominantGenera.csv")
CytokineData = read.csv("~/Documents/IowaWoundData2021/PlotsForSue2022/MergedData_for_SK.csv")
ClinicalData = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/GSWOUNDGRICE2015_20190221.csv")


DomGeneraCols = colnames(DominantGenera)

DomGeneraCols = sapply(DomGeneraCols, AddCLRDes)

colnames(DominantGenera) = DomGeneraCols



# Anything with technical mean CT > 35 = NA
# delta CT is techMean - 18SMean
############################################
CytokineData$Ct.Tech.Mean[CytokineData$Ct.Tech.Mean > 35] = NA
CytokineDataNonNA = CytokineData %>% filter(!is.na(Delta.Ct.Mean))
CytokineDataNonNA = CytokineDataNonNA %>% group_by(study_id, Target.Name) %>% summarise(study_id=study_id, Delta.Ct.Mean=mean(Delta.Ct.Mean)) %>% unique() %>% ungroup()
TableSamplesCytokine = table(CytokineDataNonNA$Target.Name)
CytokineDataNonNAGr20 = TableSamplesCytokine[TableSamplesCytokine >= (445*.20)]
dim(CytokineDataNonNAGr20)


CytokineDataNonNAInclude = CytokineDataNonNA %>% filter(Target.Name %in% names(CytokineDataNonNAGr20))

ClinicalDataAllIDs = data.frame(study_id = unique(ClinicalData$study_id))
ClinicalData$study_id = sapply(ClinicalData$study_id, toString)
ClinicalDataAllIDs$study_id = sapply(ClinicalDataAllIDs$study_id, toString)
for(target in names(CytokineDataNonNAGr20)){
  print(target)
  LittleSubset = CytokineDataNonNAInclude %>% filter(Target.Name==target) %>% select(study_id, Delta.Ct.Mean) 
  colnames(LittleSubset) = c("study_id",target )
  ClinicalDataAllIDs = ClinicalDataAllIDs %>% left_join(LittleSubset, by="study_id")
    
}

FullData = ClinicalDataAllIDs
DominantGenera$study_id = as.character(DominantGenera$StudyID)
FullData = FullData %>% left_join(DominantGenera, by="study_id")


FullDataDataPresent = FullData %>% select(study_id, DMMClusterAssign,names(CytokineDataNonNAGr20) )

FullDataDataPresent = FullDataDataPresent %>% left_join(ClinicalData %>% select(study_id, woundcarepain), by="study_id")


FullDataDataPresent[2:15][!is.na(FullDataDataPresent[2:15])] <- 1
FullDataDataPresent[2:15][is.na(FullDataDataPresent[2:15])] <- 0

colnames(FullDataDataPresent) = sapply(colnames(FullDataDataPresent), function(x) str_split(x, pattern="-Hs")[[1]][1] )
FullDataDataPresentMelt = FullDataDataPresent %>% reshape2::melt(id.vars=c("study_id", "woundcarepain"))
FullDataDataPresentMelt$variable = (if_else(FullDataDataPresentMelt$variable=="DMMClusterAssign", "Microbiome", as.character(FullDataDataPresentMelt$variable)))



StudyIDOrder = unique((FullDataDataPresentMelt %>% arrange(woundcarepain))$study_id)
DataAvailablePlot = ggplot(FullDataDataPresentMelt, aes(x=study_id, y=variable, fill=factor(value))) + geom_tile() + scale_fill_manual(values=c("white", "black")) + theme_classic() + theme(axis.text.x=element_text(angle=80), legend.position = "None")
DataAvailablePlot$data$study_id = factor(DataAvailablePlot$data$study_id, levels=StudyIDOrder)
DataAvailablePlot$data$variable = factor(DataAvailablePlot$data$variable , levels= c("Microbiome", setdiff(DataAvailablePlot$data$variable, "Microbiome")))
ggsave(DataAvailablePlot + theme(axis.text.x=element_blank()) + xlab("Patient") + ylab("Biological Variable"), file="~/Documents/IowaWoundData2021/PaperFigs/DataPresence.pdf", width=20, height=4)
ggsave(DataAvailablePlot + xlab("Patient") + ylab("Biological Variable"), file="~/Documents/IowaWoundData2021/PaperFigs/DataPresencePatientIDs.pdf", width=20, height=4)

View(FullDataDataPresentMelt %>% arrange(woundcarepain) %>% select(study_id, woundcarepain) %>% unique())


CytokineDataNonNA$Ct.Tech.Mean-CytokineDataNonNA$Delta.Ct.Mean

NumSamples = table(CytokineDataNonNA$Target.Name)

CytokineData100Plus = CytokineData %>% filter(Target.Name %in% names(NumSamples)[NumSamples>100]) %>% left_join(WoundMicrobiome, by="X")


ggplot(CytokineData100Plus, aes(x=Delta.Ct.Mean)) + geom_histogram() + facet_grid(.~Target.Na)

unique(CytokineData100Plus$Target.Name)

CXCL8Hs00174103_m1Data = CytokineData100Plus %>% filter(Target.Name=="CXCL8-Hs00174103_m1")


summary(aov(Delta.Ct.Mean~DMMClusterAssign, data=CXCL8Hs00174103_m1Data))
summary(aov(Delta.Ct.Mean~DMMClusterAssign, data=CXCL8Hs00174103_m1Data))

LCN2Hs01008571_m1data = CytokineData100Plus %>% filter(Target.Name=="LCN2-Hs01008571_m1") 
summary(aov(Delta.Ct.Mean~DMMClusterAssign, data=LCN2Hs01008571_m1data))
ggplot(LCN2Hs01008571_m1data, aes(x=factor(DMMClusterAssign), y=Delta.Ct.Mean))  + geom_boxplot()



# Log-fold changes for Moderate/Severe vs. Mild/None? the deltadeltaCT method? 
####################################################




