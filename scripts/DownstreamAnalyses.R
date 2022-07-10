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


dim(FullDataDataPresentMelt %>% arrange(woundcarepain) %>% select(woundcarepain, study_id) %>% unique())
DataAvailablePlot = ggplot(FullDataDataPresentMelt, aes(x=study_id, y=variable, fill=factor(value))) + geom_tile() + scale_fill_manual(values=c("white", "black")) + theme_classic() + theme(axis.text.x=element_text(angle=80, size=4), legend.position = "None")
DataAvailablePlot$data$study_id = factor(DataAvailablePlot$data$study_id, levels=StudyIDOrder)
DataAvailablePlot$data$variable = factor(DataAvailablePlot$data$variable , levels= c("Microbiome", setdiff(DataAvailablePlot$data$variable, "Microbiome")))
ggsave(DataAvailablePlot + theme(axis.text.x=element_blank()) + xlab("Patient") + ylab("Biological Variable"), file="~/Documents/IowaWoundData2021/PaperFigs/DataPresence.pdf", width=20, height=4)
ggsave(DataAvailablePlot + xlab("Patient") + ylab("Biological Variable"), file="~/Documents/IowaWoundData2021/PaperFigs/DataPresencePatientIDs.pdf", width=20, height=4)



# Log-fold changes for Moderate/Severe vs. Mild/None? the deltadeltaCT method?
# Test for differences, visualize that way if we have something to show
###############################################################################
FullData = FullData %>% left_join(ClinicalData %>% select(study_id, woundcarepain))
FullData = FullData %>% mutate(PainCatBinary = case_when( woundcarepain==0 | woundcarepain==1 | woundcarepain==2 ~"NoneMildModerate",
                                                           woundcarepain==3~ "Severe", 
                                                          ))


FullData$PainCatBinary[FullData$PainCatBinary=="NA"] = NA

listCytokines = c("ARG1-Hs00163660_m1",  "C3-Hs00163811_m1",  "C5AR1-Hs00704891_s1", "CAMP-Hs00189038_m1",  "CXCL8-Hs00174103_m1", "IL1A-Hs00174092_m1",  "IL1B-Hs01555410_m1",
"IL6-Hs00174131_m1", "LCN2-Hs01008571_m1", "MMP1-Hs00899658_m1", "MMP2-Hs01548727_m1", "MMP9-Hs00957562_m1", "TNF-Hs00174128_m1")


listPvalues = c()
listcytokinesTested = c()
listPvaluesNoneVsSevere = c()
for(cyt in listCytokines){
  testresult = (wilcox.test(FullData[,cyt] ~ FullData[,"PainCatBinary"]))
  listcytokinesTested = c(listcytokinesTested, cyt)
  listPvalues = c(listPvalues, testresult$p.value)
  
}






DataFrameCytokinePain = data.frame(pvals = listPvalues, Cytokine=listcytokinesTested)
DataFrameCytokinePain$PAdjust = p.adjust(DataFrameCytokinePain$pvals, method="BH")

# Ugh how to calculate Fold change 
####################################
geo_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

for(cyt in listCytokines){
  CtrlValues = FullData[, c(cyt, "PainCatBinary")] %>% filter(PainCatBinary=="NoneMildModerate") 
  print(CtrlValues)
  FullData[, paste0(cyt, "_CtrlMean")]= geo_mean(CtrlValues[,cyt])
  FullData[, paste0(cyt, "_TwoDeltaDeltaCT")] = 2^(-( (FullData[, cyt]) - (FullData[, paste0(cyt, "_CtrlMean")])))
  

}

View(FullData)
View(FullData %>% select(`IL6-Hs00174131_m1_CtrlMean`, PainCatBinary, `IL6-Hs00174131_m1_TwoDeltaDeltaCT`, `IL6-Hs00174131_m1`))

log2(mean((FullData %>% filter(PainCatBinary!="NoneMildModerate"))$`IL1B-Hs01555410_m1_TwoDeltaDeltaCT`, na.rm=T) /  mean((FullData %>% filter(PainCatBinary=="NoneMildModerate"))$`IL1B-Hs01555410_m1_TwoDeltaDeltaCT`, na.rm=T))
log2(mean((FullData %>% filter(PainCatBinary!="NoneMildModerate"))$`IL6-Hs00174131_m1_TwoDeltaDeltaCT`, na.rm=T) /  mean((FullData %>% filter(PainCatBinary=="NoneMildModerate"))$`IL6-Hs00174131_m1_TwoDeltaDeltaCT`, na.rm=T))
log2(mean((FullData %>% filter(PainCatBinary!="NoneMildModerate"))$`IL1B-Hs01555410_m1_TwoDeltaDeltaCT`, na.rm=T) /  mean((FullData %>% filter(PainCatBinary=="NoneMildModerate"))$`IL1B-Hs01555410_m1_TwoDeltaDeltaCT`, na.rm=T))
log2(mean((FullData %>% filter(PainCatBinary!="NoneMildModerate"))$`C5AR1-Hs00704891_s1_TwoDeltaDeltaCT`, na.rm=T) /  mean((FullData %>% filter(PainCatBinary=="NoneMildModerate"))$`C5AR1-Hs00704891_s1_TwoDeltaDeltaCT`, na.rm=T))



