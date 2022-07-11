# Amy Campbell
# Downstream analyses comparing the cytokine and microbiome variables to the actual 
library("phyloseq")
library("dplyr")
library("stringr")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("viridis")
library("gplots")

# Just distinguish "old" from "new" genera amounts (which should be the same for dominant genera included before)
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

# 2^-(DELTA CT) FOR EACH CYTOKINE
for(target in names(CytokineDataNonNAGr20)){
  print(target)
  LittleSubset = CytokineDataNonNAInclude %>% filter(Target.Name==target) %>% select(study_id, Delta.Ct.Mean) 
  colnames(LittleSubset) = c("study_id",target )
  LittleSubset[, target] = 2^-(LittleSubset[, target])
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


listkruskalCytokinesDMM = c()
listkruskalPs = c()
for(cyt in listCytokines){
  testresult = kruskal.test(FullData[,cyt] ~ FullData[,"DMMClusterAssign"])  
  listkruskalPs = c(listkruskalPs, testresult$p.value)
  listkruskalCytokinesDMM = c(listkruskalCytokinesDMM, cyt)
  
}

# Cytokines with significant KW p values adjusted
##################################################
for(cyt in listCytokines){
  if(cyt %in% c("C5AR1-Hs00704891_s1","CXCL8-Hs00174103_m1", "IL1B-Hs01555410_m1", "MMP2-Hs01548727_m1")){
    
    print(cyt)
    print(pairwise.wilcox.test(FullData[,cyt],factor(FullData[,"DMMClusterAssign"]), p.adjust.method = "BH") )
  }
}



DataFrameCytokinePain = data.frame(pvals = listPvalues, Cytokine=listcytokinesTested)
DataFrameCytokinePain$PAdjust = p.adjust(DataFrameCytokinePain$pvals, method="BH")



DataFrameCytokineDMMs = data.frame(pvalsKruskalWallis =listkruskalPs, Cytokine=listkruskalCytokinesDMM )
DataFrameCytokineDMMs$PAdjustKW = p.adjust(DataFrameCytokineDMMs$pvalsKruskalWallis, method="BH")

CytokineDataBySample = FullData %>% select(listCytokines, study_id, PainCatBinary, DMMClusterAssign, woundcarepain)

numCytokines = FullDataDataPresent %>% select(-woundcarepain,-DMMClusterAssign ) 
numCytokines$Ncyt = rowSums(numCytokines[2:ncol(numCytokines)])

MeltCytokineDataBySample

CytokineDataBySample = CytokineDataBySample %>% left_join(numCytokines %>% select(study_id, Ncyt),  by ="study_id")
CytokineDataBySample_RemoveZeroCytokine = CytokineDataBySample %>% filter(Ncyt > 0 )
MeltCytokineDataBySample = CytokineDataBySample_RemoveZeroCytokine %>% reshape2::melt(id.vars=(c("woundcarepain", "DMMClusterAssign", "PainCatBinary", "study_id", "Ncyt")))

orderHeatMap = unique((CytokineDataBySample %>% arrange(DMMClusterAssign, Ncyt))$study_id)
orderHeatMap = unique((CytokineDataBySample %>% arrange(woundcarepain, Ncyt))$study_id)

CytokineBySubjectPlot = ggplot(MeltCytokineDataBySample, aes(x=study_id, y=variable, fill=log2(value))) + geom_tile() + scale_fill_viridis(option="plasma", na.value="white")
CytokineBySubjectPlot$data$study_id = factor(CytokineBySubjectPlot$data$study_id, levels=orderHeatMap)
CytokineBySubjectPlot$data$variable = factor(CytokineBySubjectPlot$data$variable, levels=names(rev(sort(CytokineDataNonNAGr20))))

CytokineBySubjectPlot + theme_classic() + theme(axis.text.x=element_text(angle=70, hjust=.5, vjust=.5))



# DMM vs. markers
# "C5AR1-Hs00704891_s1","CXCL8-Hs00174103_m1", "IL1B-Hs01555410_m1", "MMP2-Hs01548727_m1"
FullDataContainsMicrobiome = FullData %>% filter(!is.na(DMMClusterAssign))
c5ar1plot = ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=log2(`C5AR1-Hs00704891_s1` ))) + geom_boxplot(fill="gray87") + theme_classic() + xlab("Microbial Community Type (DMM Assignment)") + ylab("-∆CT") + ggtitle("C5AR1 Expression") + theme(plot.title=element_text(hjust=.5, size=20, face="bold"))
cxcl8plot  =  ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=log2(`CXCL8-Hs00174103_m1` ))) + geom_boxplot(fill="gray87") + theme_classic() + xlab("Microbial Community Type (DMM Assignment)") + ylab("-∆CT")+ ggtitle("CXCL8 Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold"))
ilbplot  =  ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=log2(`IL1B-Hs01555410_m1` ))) + geom_boxplot(fill="gray87") + theme_classic() + xlab("Microbial Community Type (DMM Assignment)") + ylab("-∆CT")+ ggtitle("IlB Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold"))

ggplot(FullData, aes(x=log2(`CXCL8-Hs00174103_m1` ), y=log2(`IL1B-Hs01555410_m1`))) + geom_point()

cor.test(log2(FullData$`CXCL8-Hs00174103_m1` ), log2(FullData$`IL1B-Hs01555410_m1`))

FullData %>% select(listCytokines) %>% mutate_all(.funs=log2(.))


# Correlations 
JustCyt = log2(FullData %>% select(listCytokines) )
colnames(JustCyt) = sapply(colnames(JustCyt) , function(x) str_split(x, pattern="-")[[1]][1])
JustCytCor = cor(JustCyt, use="pairwise.complete.obs")

CorrelationMat$p = "NA"
for(r in 1:nrow(CorrelationMat)){
  x=CorrelationMat[r, "Var1"]
  y=CorrelationMat[r, "Var2"]
  ctest = cor.test(JustCyt[,x], JustCyt[,y])
  CorrelationMat[r, "p"] = as.numeric(ctest$p.value)
}

CorrelationMatBonfNum = nrow(CorrelationMat %>% filter(Var1!=Var2))/2
CorrelationMat$p = sapply(CorrelationMat$p, function(x) as.numeric(as.character(x)))
CorrelationMat$pAdj = (CorrelationMat$p)*CorrelationMatBonfNum
CorrelationMat = JustCytCor %>% reshape2::melt()
View(CorrelationMat)
CorrHeatMap = ggplot(CorrelationMat, aes(x=Var1, y=Var2, fill=value))+ geom_tile() + scale_fill_gradient2(high="red", low="blue", mid="white", midpoint=0) + xlab("") + ylab("") + ggtitle("Correlations Between Inflammatory\nMarker Expression (∆CT)") +theme_classic()+ theme(plot.title=element_text(size=20, hjust=.5), axis.text.x=element_text(size=15, angle=270), axis.text.y=element_text(size=15)) + labs(fill="Pearson's Correlation")

Hclustering = hclust(as.dist(1-JustCytCor))

CorrHeatMap$data$Var1 = factor(CorrHeatMap$data$Var1, levels=(Hclustering$labels)[Hclustering$order])
CorrHeatMap$data$Var2 = factor(CorrHeatMap$data$Var2, levels=(Hclustering$labels)[Hclustering$order])
CorrHeatMap
