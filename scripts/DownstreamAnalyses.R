# Amy Campbell
# Downstream analyses comparing the cytokine and microbiome variables to the actual 
#library("phyloseq")
library("dplyr")
library("stringr")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("viridis")
library("gplots")
library("rstatix")
library("correlation")
library("psych")
# # Imputed PCA functions
# library("missMDA")
# library("FactoMineR")

# Just distinguish "old" from "new" genera amounts (which should be the same for dominant genera included before)
AddCLRDes = function(namestring){
  if(!(namestring %in% c("X", "SampleID", "StudyID", "DMMClusterAssign"))){
    return(paste0(namestring,"CLR"))
  }else{
    return(namestring)
  }
  
}

WoundMicrobiome = read.csv("~/Documents/IowaWoundData2021/PlotsForSue2022/WoundMicrobiomeDataForSEG_AEC.csv")
CytokineData = read.csv("~/Documents/IowaWoundData2021/PlotsForSue2022/MergedData_for_SK.csv")
ClinicalData = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/GSWOUNDGRICE2015_20190221.csv")

ClinicalData$resting_pain_cat = case_when( ClinicalData$resting_pain == 0.0 ~ "None", 
                                             ClinicalData$resting_pain <= 3.0 & ClinicalData$resting_pain  > 0.0 ~ "Mild",
                                             ClinicalData$resting_pain > 3.0 & ClinicalData$resting_pain <=7.0 ~ "Moderate", 
                                           ClinicalData$resting_pain > 7.0  ~ "Severe"
                                             )
  
  
WoundMicrobiome$study_id = sapply(WoundMicrobiome$StudyID, as.character)

WoundDepthData = read.csv("~/Documents/IowaWoundData2021/wound_depth_covariate.csv")
##################################################################################
# 0) PREPARE CYTOKINE DATA, MERGE WITH MICROBIOME DATA, AND SUMMARIZE COMPLETENESS
##################################################################################
# correlation between pain types (resting vs. wound care)
PainRatings = ClinicalData %>% select(resting_pain_cat, woundcarepain)
PainRatings = PainRatings %>% mutate(OrdinalResting = case_when(resting_pain_cat=="None" ~ 0, 
                                                                resting_pain_cat=="Mild" ~1, 
                                                                resting_pain_cat == "Moderate" ~ 2,
                                                                resting_pain_cat == "Severe" ~3))

PainRatings = PainRatings %>% select(woundcarepain, OrdinalResting)
PainRatings$woundcarepain = factor(PainRatings$woundcarepain)
PainRatings$OrdinalResting = factor(PainRatings$OrdinalResting)

WCP_rest_cor = correlation::cor_test(PainRatings,"woundcarepain", "OrdinalResting",method="polychoric")

TablePainCareResting = table(PainRatings)
chisq.test(TablePainCareResting)
rownames(TablePainCareResting) = c("None", "Mild", "Moderate","Severe")
colnames(TablePainCareResting) = c("None", "Mild", "Moderate","Severe")
write.csv(TablePainCareResting, file="~/Documents/IowaWoundData2021/PainTypesRelationshipWhole.csv")


TablePainCareResting[,1] = TablePainCareResting[,1]/sum(TablePainCareResting[,1])
TablePainCareResting[,2] = TablePainCareResting[,2]/sum(TablePainCareResting[,2])
TablePainCareResting[,3] = TablePainCareResting[,3]/sum(TablePainCareResting[,3])
TablePainCareResting[,4] = TablePainCareResting[,4]/sum(TablePainCareResting[,4])

write.csv(TablePainCareResting, file="~/Documents/IowaWoundData2021/PainTypesRelationshipFreq.csv")



# Anything with technical mean CT > 35 = NA
# delta CT is techMean - 18SMean
############################################
CytokineData$Ct.Tech.Mean[CytokineData$Ct.Tech.Mean > 35] = NA
CytokineDataNonNA = CytokineData %>% filter(!is.na(Delta.Ct.Mean))
CytokineDataNonNA = CytokineDataNonNA %>% group_by(study_id, Target.Name) %>% summarise(study_id=study_id, Delta.Ct.Mean=mean(Delta.Ct.Mean)) %>% unique() %>% ungroup()
TableSamplesCytokine = table(CytokineDataNonNA$Target.Name)

# filter to those detected in >= 20% of samples 
CytokineDataNonNAGr20 = TableSamplesCytokine[TableSamplesCytokine >= (445*.20)]


CytokineDataNonNAInclude = CytokineDataNonNA %>% filter(Target.Name %in% names(CytokineDataNonNAGr20))

ClinicalDataAllIDs = data.frame(study_id = unique(ClinicalData$study_id))
ClinicalData$study_id = sapply(ClinicalData$study_id, toString)
ClinicalDataAllIDs$study_id = sapply(ClinicalDataAllIDs$study_id, toString)

# -(DELTA CT) FOR EACH CYTOKINE
for(target in names(CytokineDataNonNAGr20)){
  print(target)
  LittleSubset = CytokineDataNonNAInclude %>% filter(Target.Name==target) %>% select(study_id, Delta.Ct.Mean) 
  colnames(LittleSubset) = c("study_id",target )
  LittleSubset[, target] =-(LittleSubset[, target])
  ClinicalDataAllIDs = ClinicalDataAllIDs %>% left_join(LittleSubset, by="study_id")
    
}


FullData = ClinicalDataAllIDs
FullData = FullData %>% left_join(WoundMicrobiome, by="study_id")

# Summarize presence/absence of the subset of cytokines 


FullDataDataPresent = FullData %>% select(study_id, DMMClusterAssign,names(CytokineDataNonNAGr20) )

FullDataDataPresent = FullDataDataPresent %>% left_join(ClinicalData %>% select(study_id, woundcarepain, resting_pain_cat), by="study_id")


FullDataDataPresent[2:15][!is.na(FullDataDataPresent[2:15])] <- 1
FullDataDataPresent[2:15][is.na(FullDataDataPresent[2:15])] <- 0

colnames(FullDataDataPresent) = sapply(colnames(FullDataDataPresent), function(x) str_split(x, pattern="-Hs")[[1]][1] )
FullDataDataPresentMelt = FullDataDataPresent %>% reshape2::melt(id.vars=c("study_id", "woundcarepain", "resting_pain_cat"))
FullDataDataPresentMelt$variable = (if_else(FullDataDataPresentMelt$variable=="DMMClusterAssign", "Microbiome", as.character(FullDataDataPresentMelt$variable)))


StudyIDOrder = unique((FullDataDataPresentMelt %>% arrange(woundcarepain))$study_id)

dim(FullDataDataPresentMelt %>% arrange(woundcarepain) %>% select(woundcarepain, study_id) %>% unique())


DataAvailablePlot = ggplot(FullDataDataPresentMelt, aes(x=study_id, y=variable, fill=factor(value))) + geom_tile() + scale_fill_manual(values=c("white", "black")) + theme_classic() + theme(axis.text.x=element_text(angle=80, size=4), legend.position = "None")
DataAvailablePlot$data$study_id = factor(DataAvailablePlot$data$study_id, levels=StudyIDOrder)
DataAvailablePlot$data$variable = factor(DataAvailablePlot$data$variable , levels= c("Microbiome", setdiff(DataAvailablePlot$data$variable, "Microbiome")))

colorDF = FullDataDataPresentMelt %>% arrange(woundcarepain) %>% select(woundcarepain, study_id, resting_pain_cat) %>% unique()
colorDF = colorDF %>% mutate(ColorAssign = case_when(woundcarepain==0 ~ "khaki1",
                             woundcarepain==1 ~ "gold", 
                             woundcarepain==2 ~ "darkorange1",
                             woundcarepain==3~  "red3"))
colorDF = colorDF %>% mutate(ColorAssignResting = case_when(resting_pain_cat=="None" ~ "khaki1",
                                                            resting_pain_cat=="Mild" ~ "gold", 
                                                            resting_pain_cat=="Moderate" ~ "darkorange1",
                                                            resting_pain_cat=="Severe" ~  "red3"))



DataAvailablePlotResting = DataAvailablePlot + theme(axis.ticks.x = element_line(size=2, color=colorDF$ColorAssignResting))

ggsave(DataAvailablePlotResting + theme(axis.text.x=element_blank()) + xlab("Patient") + ylab("Biological Variable"), file="~/Documents/IowaWoundData2021/PaperFigs/DataPresenceRestingColorsForUse.pdf", width=20, height=4)

ggsave(DataAvailablePlot + theme(axis.text.x=element_blank()) + xlab("Patient") + ylab("Biological Variable"), file="~/Documents/IowaWoundData2021/PaperFigs/DataPresence.pdf", width=20, height=4)
ggsave(DataAvailablePlot + xlab("Patient") + ylab("Biological Variable"), file="~/Documents/IowaWoundData2021/PaperFigs/DataPresencePatientIDs.pdf", width=20, height=4)

###############################################################################
# (1) WOUND  PAIN RATINGS VS. BIOLOGICAL  VARIABLES (microbiome then cytokines)
###############################################################################

# Microbiome vs. pain
##########################

FullData_WoundCarePain  = FullData %>% left_join(ClinicalData %>% select(study_id, woundcarepain))
FullData_WoundCarePain = FullData_WoundCarePain %>% mutate(PainCatBinary = case_when( woundcarepain==0 | woundcarepain==1  ~"NoneMild",
                                                          woundcarepain==3~ "Severe", 
                                                          woundcarepain==2 ~ "Moderate" 
))




FullData_WoundCarePain$PainCatBinary[FullData_WoundCarePain$PainCatBinary=="NA"] = NA
FullData_TestSevereVsMildNone = FullData_WoundCarePain %>% filter(PainCatBinary != "Moderate")

# Common genus abundance vs. wound dressing change pain (WOUND CARE PAIN)
##########################################################################
AbundanceCLR = colnames(FullData_WoundCarePain)[grepl(colnames(FullData_WoundCarePain), pattern="_CLR")]
FullData_TestSevereVsMildNoneCLRs = FullData_WoundCarePain %>% filter(PainCatBinary != "Moderate")

pvallistGenus = c()
for(genus in AbundanceCLR){
  wiltest = (wilcox.test(FullData_TestSevereVsMildNoneCLRs[,genus] ~ FullData_TestSevereVsMildNoneCLRs[, "PainCatBinary"], exact=F) )
  pvallistGenus=append(pvallistGenus,  wiltest$p.value)
}


GenusAbundances= data.frame(Genus=AbundanceCLR, wilcox_p = pvallistGenus)
GenusAbundances$PAdj = p.adjust(GenusAbundances$wilcox_p, method="BH")

FullData_TestSevereVsMildNoneCLRs$PainCatBinary = factor(FullData_TestSevereVsMildNoneCLRs$PainCatBinary)

DataMeltedGenusAbundance = FullData_TestSevereVsMildNoneCLRs %>% select(AbundanceCLR, PainCatBinary) %>% melt(id.vars=c("PainCatBinary"))
DataMeltedGenusAbundance$Genus = sapply(DataMeltedGenusAbundance$variable, function(x) str_split(x,pattern="Abundance_CLR")[[1]][1])

GenusStats = DataMeltedGenusAbundance %>% group_by(Genus) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary")
GenusStats$y.position = GenusStats$y.position + 1
GenusPainPlot = ggplot(DataMeltedGenusAbundance, aes(x=PainCatBinary, y=value, fill=Genus)) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray34", alpha=.1)+ geom_boxplot(alpha=.4, size=.25) + facet_grid(~ Genus) +
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2) + theme_classic() + ggpubr::stat_pvalue_manual(GenusStats, label="p") +
  ylim(0, 14) + xlab("Pain Rating Category") + ggtitle("Common Genus Abundance in Wounds with \nSevere vs. None/Mild Pain Ratings") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") + ylab("CLR-transformed relative abundance") 
#colour="cyan2"

ggsave(GenusPainPlot, file="~/Documents/IowaWoundData2021/PaperFigs/GenusPain.pdf", width=10, height=7)

# Common genus abundance vs. resting pain
######################################### 

FullDataWithResting = FullData %>% left_join(ClinicalData %>% select(study_id, resting_pain_cat), by="study_id")
FullDataWithRestingNoMods = FullDataWithResting %>% filter(resting_pain_cat != "Moderate")
FullDataWithRestingNoMods$RestingBinary = if_else(FullDataWithRestingNoMods$resting_pain_cat=="Severe", "Severe", "NoneMild")
pvallistGenusResting = c()
for(genus in AbundanceCLR){
  wiltest = (wilcox.test(FullDataWithRestingNoMods[,genus] ~ FullDataWithRestingNoMods[, "RestingBinary"]))
  
  pvallistGenusResting=append(pvallistGenusResting,  wiltest$p.value)
}

GenusAbundancesRestingPain= data.frame(Genus=AbundanceCLR, wilcox_p = pvallistGenusResting)
GenusAbundancesRestingPain$PAdj = p.adjust(GenusAbundancesRestingPain$wilcox_p, method="BH")

FullDataWithRestingMelt = FullDataWithRestingNoMods %>% select(AbundanceCLR, RestingBinary) %>% melt(id.vars=c("RestingBinary"))
FullDataWithRestingMelt$Genus = sapply(FullDataWithRestingMelt$variable, function(x) str_split(x,pattern="Abundance_CLR")[[1]][1])
GenusRestingPainPlot = ggplot(FullDataWithRestingMelt, aes(x=RestingBinary, y=value, fill=Genus))+ geom_boxplot(alpha=.4) + geom_jitter(width=.2)  + facet_grid(~ Genus)  + theme_classic() +  scale_fill_brewer(palette ="Dark2" ) 

# No observed relationship between common wound genera's abundance (corynebacterium, anaerobic genera aggregated, pseudomonas, staph, strep )

# Chi-Sq for proportion of each DMM in each dressing change pain category (binary)
###################################################################################
TwoWay = FullData_TestSevereVsMildNoneCLRs %>% select(DMMClusterAssign, PainCatBinary) %>% na.omit()
TwoWay$PainCatBinary = factor(TwoWay$PainCatBinary)
TwoWay$DMMClusterAssign = factor(TwoWay$DMMClusterAssign)
chisq.test(table(TwoWay))

TwoWayAllCats = FullData_WoundCarePain %>% select(DMMClusterAssign, woundcarepain) %>% na.omit()
chisq.test(table(TwoWayAllCats))

TwoWayAllCatsTable = table(TwoWayAllCats)

write.csv(TwoWayAllCatsTable, file="~/Documents/IowaWoundData2021/PainDMMWhole.csv")

TwoWayAllCatsTable[,1] = round(TwoWayAllCatsTable[,1]/sum(TwoWayAllCatsTable[,1]), 3)
TwoWayAllCatsTable[,2] = round(TwoWayAllCatsTable[,2]/sum(TwoWayAllCatsTable[,2]),3)
TwoWayAllCatsTable[,3] = round(TwoWayAllCatsTable[,3]/sum(TwoWayAllCatsTable[,3]),3)
TwoWayAllCatsTable[,4] = round(TwoWayAllCatsTable[,4]/sum(TwoWayAllCatsTable[,4]),3)
write.csv(TwoWayAllCatsTable, file="~/Documents/IowaWoundData2021/PainDMMFreq.csv")


# Chi-Sq for proportion of each DMM in each resting pain category (binary)
##########################################################################
TwoWayResting =  FullDataWithRestingNoMods %>% select(DMMClusterAssign, RestingBinary) %>% na.omit()
chisq.test(table(TwoWayResting))
TwoWayRestingCat=  FullDataWithRestingNoMods %>% select(DMMClusterAssign, resting_pain_cat) %>% na.omit()
chisq.test(table(TwoWayRestingCat) )

TwoWayRestingAll = FullDataWithResting %>% select(DMMClusterAssign, resting_pain_cat) %>% na.omit()
TwoWayRestingAll = table(TwoWayRestingAll)

write.csv(TwoWayRestingAll, file="~/Documents/IowaWoundData2021/RestingPainDMMWhole.csv")

TwoWayRestingAll[,1] = round(TwoWayRestingAll[,1]/sum(TwoWayRestingAll[,1]), 3)
TwoWayRestingAll[,2] = round(TwoWayRestingAll[,2]/sum(TwoWayRestingAll[,2]),3)
TwoWayRestingAll[,3] = round(TwoWayRestingAll[,3]/sum(TwoWayRestingAll[,3]),3)
TwoWayRestingAll[,4] = round(TwoWayRestingAll[,4]/sum(TwoWayRestingAll[,4]),3)
write.csv(TwoWayRestingAll, file="~/Documents/IowaWoundData2021/RestingPainDMMFreq.csv")


# Wilcox for shannon and richness-based diversity metrics for severe vs. none/mild wound care pain
##################################################################################################
DiversityPain= FullData_TestSevereVsMildNoneCLRs %>% select(Genus_Shannon,Genus_Richness, PainCatBinary)
DiversityPainMelt = DiversityPain %>% select(Genus_Shannon, Genus_Richness, PainCatBinary) %>% melt(id.vars=c("PainCatBinary"))
DiversityPainMelt$Metric = DiversityPainMelt$variable
DiversityPainMeltStats = DiversityPainMelt %>% group_by(Metric) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="Pain")

# Wilcox for shannon and richness-based diversity metrics for severe vs. none/mild wound care pain
##################################################################################################
DiversityPainResting= FullDataWithRestingNoMods %>% select(Genus_Shannon,Genus_Richness, RestingBinary)
DiversityPainMeltResting = DiversityPainResting %>% select(Genus_Shannon, Genus_Richness, RestingBinary) %>% melt(id.vars=c("RestingBinary"))
DiversityPainMeltResting$Metric = DiversityPainMeltResting$variable
DiversityPainMeltRestingStats = DiversityPainMeltResting %>% group_by(Metric) %>% wilcox_test(value ~ RestingBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="Pain")


# Now, cytokines vs. pain
##########################
listCytokines = c("ARG1-Hs00163660_m1",  "C3-Hs00163811_m1",  "C5AR1-Hs00704891_s1", "CAMP-Hs00189038_m1",  "CXCL8-Hs00174103_m1", "IL1A-Hs00174092_m1",  "IL1B-Hs01555410_m1",
                  "IL6-Hs00174131_m1", "LCN2-Hs01008571_m1", "MMP1-Hs00899658_m1", "MMP2-Hs01548727_m1", "MMP9-Hs00957562_m1", "TNF-Hs00174128_m1")


# Wound care pain 
#################
FullData_TestSevereVsMildNone$Nmissing = rowSums(is.na(FullData_TestSevereVsMildNone %>% select(all_of(listCytokines))))
dim(FullData_TestSevereVsMildNone %>% filter(Nmissing < 13))
listPvalues = c()
listcytokinesTested = c()
listPvaluesNoneVsSevere = c()
kruskal_ps = c()
listpvalues_Age = c()
for(cyt in listCytokines){
  
  testresult = (wilcox.test((FullData_TestSevereVsMildNone[,cyt]) ~ FullData_TestSevereVsMildNone[,"PainCatBinary"],exact=F))
  listcytokinesTested = c(listcytokinesTested, cyt)
  listPvalues = c(listPvalues, testresult$p.value)
  kruskal_all = kruskal.test(FullData_WoundCarePain[, cyt] ~ FullData_WoundCarePain[, "woundcarepain"])
  kruskal_ps =append(kruskal_ps, kruskal_all$p.value)
  
}

DataFrameCytokinePain = data.frame(pvals = listPvalues, Cytokine=listcytokinesTested, kruskalp = kruskal_ps)
DataFrameCytokinePain$PAdjust = p.adjust(DataFrameCytokinePain$pvals, method="BH")
DataFrameCytokinePain$PAdjustKruskal = p.adjust(DataFrameCytokinePain$kruskalp, method="BH")

DataFrameCytokinePainJustWilcox = DataFrameCytokinePain %>% select(pvals, Cytokine, PAdjust)
# Resting pain 
##############

FullDataWithRestingNoMods = FullDataWithResting %>% filter(resting_pain_cat != "Moderate")
FullDataWithRestingNoMods$RestingBinary = if_else(FullDataWithRestingNoMods$resting_pain_cat=="Severe", "Severe", "NoneMild")


listPvaluesRest= c()
listcytokinesTestedRest = c()
listPvaluesNoneVsSevereRest = c()
listpvalues_Age_Rest = c()
for(cyt in listCytokines){
  
  testresult = (wilcox.test((FullDataWithRestingNoMods[,cyt]) ~ FullDataWithRestingNoMods[,"RestingBinary"],exact=F))
  listcytokinesTestedRest = c(listcytokinesTestedRest, cyt)
  listPvaluesRest = c(listPvaluesRest, testresult$p.value)

}


DataFrameCytokinePain = data.frame(pvals = listPvalues, Cytokine=listcytokinesTested, kruskalp = kruskal_ps)
DataFrameCytokinePain$PAdjust = p.adjust(DataFrameCytokinePain$pvals, method="BH")
DataFrameCytokinePain$PAdjustKruskal = p.adjust(DataFrameCytokinePain$kruskalp, method="BH")

DataFrameCytokinePain$pval_Resting = listPvaluesRest


# no wilcoxon tests yielded relationships between wound care pain, cytokines or resting pain, cytokines


# Visualizing Wound care pain vs. cytokines
#############################################
FullDataPain = FullData %>% left_join(ClinicalData %>% select(woundcarepain, study_id))
FullDataPain  = FullDataPain %>% mutate(PainCatBinary = case_when(woundcarepain %in% c(0,1) ~ "NoneMild",
                                       woundcarepain==2 ~ "Moderate", 
                                       woundcarepain==3 ~ "Severe"))

CytokineDataBySample = FullDataPain %>% select(listCytokines, study_id, PainCatBinary, DMMClusterAssign, woundcarepain)


numCytokines = FullDataDataPresent %>% select(-woundcarepain,-DMMClusterAssign ) 
numCytokines$Ncyt = rowSums(numCytokines[2: (ncol(numCytokines) - 1) ])
CytokineDataBySample = CytokineDataBySample %>% left_join(numCytokines %>% select(study_id, Ncyt),  by ="study_id")

CytokineDataBySample_RemoveZeroCytokine = CytokineDataBySample %>% filter(Ncyt > 0 )

MeltCytokineDataBySample = CytokineDataBySample_RemoveZeroCytokine %>% reshape2::melt(id.vars=(c("woundcarepain", "DMMClusterAssign", "PainCatBinary", "study_id", "Ncyt")))



orderHeatMapDMM = ((CytokineDataBySample %>% arrange(DMMClusterAssign, Ncyt))$study_id)
orderHeatMap = ((CytokineDataBySample %>% arrange(woundcarepain, Ncyt))$study_id)


names(CytokineDataNonNAGr20) = sapply(names(CytokineDataNonNAGr20), function(x) str_split(x, pattern="-")[[1]][1])
MeltCytokineDataBySample$cytokine = sapply(MeltCytokineDataBySample$variable, function(x) str_split(x, pattern="-")[[1]][1])
CytokineBySubjectPlot = ggplot(MeltCytokineDataBySample, aes(x=study_id, y=cytokine, fill=(value))) + geom_tile() + scale_fill_viridis(option="plasma", na.value="white")
CytokineBySubjectPlot$data$study_id = factor(CytokineBySubjectPlot$data$study_id, levels=orderHeatMap)

CytokineDataNonNAGr20 = sapply(CytokineDataNonNAGr20, function(x) as.numeric(as.character(x)))
CytokineBySubjectPlot$data$cytokine = factor(CytokineBySubjectPlot$data$cytokine, levels=names(rev(sort(CytokineDataNonNAGr20))))
CytokineBySubjectPlot = CytokineBySubjectPlot + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=.5, vjust=.5), axis.text.y=element_text(size=10), plot.title=element_text(face="bold", size=18, hjust=.5)) + ggtitle("Inflammatory Gene Expression by Patient (-∆CT)") + labs(fill="-∆CT")

ggsave(CytokineBySubjectPlot, file="~/Documents/IowaWoundData2021/PaperFigs/InflamVsPatient.png",width=20, height=10)
ggsave(CytokineBySubjectPlot, file="~/Documents/IowaWoundData2021/PaperFigs/InflamVsPatient.pdf",width=20, height=10)

DataFrameCytokinePainJustWilcox = DataFrameCytokinePainJustWilcox %>% select(Cytokine, pvals, PAdjust)
row.names(DataFrameCytokinePainJustWilcox) = sapply(DataFrameCytokinePainJustWilcox$Cytokine, function(x) str_split(x, pattern="-")[[1]][1])

DataFrameCytokinePainJustWilcox = DataFrameCytokinePainJustWilcox[names(rev(sort(CytokineDataNonNAGr20))),]
DataFrameCytokinePainJustWilcox$pvals = round(DataFrameCytokinePainJustWilcox$pvals, 3)
DataFrameCytokinePainJustWilcox$PAdjust = round(DataFrameCytokinePainJustWilcox$PAdjust, 3)
write.csv(DataFrameCytokinePainJustWilcox, file="~/Documents/IowaWoundData2021/CytokinesVsHealing.csv")
############################################
# (2) BIOLOGICAL VARIABLES  VS. INFLAMMATION 
#############################################

# Testing relationship between common wound genera, inflammation
###############################################################

FullDataWithOtherClinicalsInflame = FullData %>% left_join(ClinicalData %>% select(inflame,study_id), by="study_id")
DataMeltedGenusAbundanceInflammation = FullDataWithOtherClinicalsInflame %>% select(AbundanceCLR, inflame) %>% melt(id.vars=c("inflame"))
DataMeltedGenusAbundanceInflammation$Genus = sapply(DataMeltedGenusAbundanceInflammation$variable, function(x) str_split(x,pattern="Abundance_CLR")[[1]][1])
DataMeltedGenusAbundanceInflammation$inflame = factor(DataMeltedGenusAbundanceInflammation$inflame)
GenusStatsInflame = DataMeltedGenusAbundanceInflammation %>% group_by(Genus) %>% wilcox_test(value ~ inflame)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="Inflammation")

DataMeltedGenusAbundanceInflammation$inflame=factor(DataMeltedGenusAbundanceInflammation$inflame)
GenusStatsInflame$y.position = GenusStatsInflame$y.position + 1


GenusStatsInflame$xmax=2
GenusStatsInflame$xmin=1

GenusInflammationPlot = ggplot(DataMeltedGenusAbundanceInflammation, aes(x=inflame, y=value, fill=Genus)) + geom_boxplot(alpha=.4) + facet_grid(~ Genus)+
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2) + theme_classic() + ggpubr::stat_pvalue_manual(GenusStatsInflame, label="p.adj")+
  ylim(0, 14) + xlab("Inflammation (0/1)") + ggtitle("Common Genus Abundance in Wounds which are Inflamed vs. Not ") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") + ylab("CLR-transformed relative abundance") 
ggsave(GenusInflammationPlot, file="~/Documents/IowaWoundData2021/PaperFigs/GenusInflammation.pdf", width=10, height=7)


# Testing relationship between DMM Assignment, inflammation
############################################################
TwoWayInflameCluster = FullDataWithOtherClinicalsInflame %>% select(DMMClusterAssign, inflame) %>% na.omit()
chisq.test(table(TwoWayInflameCluster))
TwoWayInflameCluster = table(TwoWayInflameCluster)
TwoWayInflameCluster
write.csv(TwoWayInflameCluster, file="~/Documents/IowaWoundData2021/InflameDMMWhole.csv")


TwoWayInflameCluster[,1] = round(TwoWayInflameCluster[,1] /sum(TwoWayInflameCluster[,1]),3)
TwoWayInflameCluster[,2] = round(TwoWayInflameCluster[,2] /sum(TwoWayInflameCluster[,2]),3)
write.csv(TwoWayInflameCluster, file="~/Documents/IowaWoundData2021/InflameDMMFreq.csv")



# Now, cytokines vs. inflammation
##########################
listCytokines = c("ARG1-Hs00163660_m1",  "C3-Hs00163811_m1",  "C5AR1-Hs00704891_s1", "CAMP-Hs00189038_m1",  "CXCL8-Hs00174103_m1", "IL1A-Hs00174092_m1",  "IL1B-Hs01555410_m1",
                  "IL6-Hs00174131_m1", "LCN2-Hs01008571_m1", "MMP1-Hs00899658_m1", "MMP2-Hs01548727_m1", "MMP9-Hs00957562_m1", "TNF-Hs00174128_m1")


listPvaluesInflame = c()
listcytokinesTestedInflame = c()
listPvaluesInflame= c()
for(cyt in listCytokines){
  
  testresult = (wilcox.test((FullDataWithOtherClinicalsInflame[,cyt]) ~ FullDataWithOtherClinicalsInflame[,"inflame"],exact=F))
  
  listcytokinesTestedInflame = c(listcytokinesTestedInflame, cyt)
  listPvaluesInflame = c(listPvaluesInflame, testresult$p.value)
  
}

DataFrameCytokineInflame = data.frame( Cytokine=listcytokinesTested,pvals = listPvaluesInflame,padj = p.adjust(listPvaluesInflame, method="BH"))
DataFrameCytokineInflame$variable = sapply(DataFrameCytokineInflame$Cytokine,  function(x) (stringr::str_split(x, pattern="-"))[[1]][1])

MeltedInflammationMarkers = FullDataWithOtherClinicalsInflame %>% select(listCytokines, inflame) %>% melt(id.vars=c("inflame"))
MeltedInflammationMarkers$Cytokine = sapply(MeltedInflammationMarkers$variable, function(x) (stringr::str_split(x, pattern="-"))[[1]][1] )
MeltedInflammationMarkers$inflame = factor(if_else(MeltedInflammationMarkers$inflame==0, "Not\nInflamed", "Inflamed"))

GenusStatsInflameCytokine = MeltedInflammationMarkers %>% group_by(Cytokine) %>% wilcox_test(value ~ inflame)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="inflammation")

GenusStatsInflameCytokine$xmin = 1
GenusStatsInflameCytokine$xmax = 2
GenusStatsInflameCytokine$y.position = GenusStatsInflameCytokine$y.position - 2

colorWorkaround = (MeltedInflammationMarkers %>% select(inflame, variable) %>% unique())
colorWorkaround = colorWorkaround %>% mutate(colorAssign = if_else(inflame=="Inflamed","#984EA3","#FF7F00" ))
MeltedInflammationMarkers$inflammation =(MeltedInflammationMarkers$inflame)
boxplotinflame = ggplot(MeltedInflammationMarkers, aes(x=inflame, y=value))+ geom_boxplot(alpha=.8,fill=colorWorkaround$colorAssign, size=.25) +geom_jitter(width=.1,size=.5) + facet_grid(~Cytokine)   +
  theme_classic() + ggpubr::stat_pvalue_manual(GenusStatsInflameCytokine,label="p") +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") + ylab("-∆CT")
ggsave(boxplotinflame, file="~/Documents/IowaWoundData2021/PaperFigs/InflammationCytokines.pdf", width=16, height=10)


########################################################################################
# (3) BIOLOGICAL VARIABLES  VS. WOUND AGE(CHRONIC VS ACUTE AKA >30 days <= 30 days)
#########################################################################################

MeltedInflammationMarkers = FullDataWithOtherClinicalsInflame %>% select(listCytokines, inflame) %>% melt(id.vars=c("inflame"))
MeltedInflammationMarkers$Cytokine = sapply(MeltedInflammationMarkers$variable, function(x) (stringr::str_split(x, pattern="-"))[[1]][1] )
MeltedInflammationMarkers$inflame = factor(if_else(MeltedInflammationMarkers$inflame==0, "Not\nInflamed", "Inflamed"))

FullDataAge = FullData %>% left_join(ClinicalData %>% select(woundage,study_id, woundcarepain), by="study_id")

FullDataAgeMeltGenus= FullDataAge %>% select(AbundanceCLR, woundage) %>% melt(id.vars=c("woundage"))
FullDataAgeMeltGenus$woundage = if_else(FullDataAgeMeltGenus$woundage %in% c(1,2), "Acute", "Chronic")

FullDataAgeMeltGenus$Genus = sapply(FullDataAgeMeltGenus$variable, function(x) (stringr::str_split(x, pattern="Abundance_"))[[1]][1] )
GenusStatsAge = FullDataAgeMeltGenus %>% group_by(Genus) %>% wilcox_test(value ~ woundage)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="woundage")

GenusAgePlot = ggplot(FullDataAgeMeltGenus, aes(x=woundage, y=value, fill=Genus)) + geom_boxplot(alpha=.4, size=.25) + geom_jitter(width=.1,size=.5) + facet_grid(~Genus) +theme_classic() +
  ggpubr::stat_pvalue_manual(GenusStatsAge,label="p") +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") + ylab("-deltaCT") + scale_fill_brewer(palette="Dark2")

ggsave(GenusAgePlot,file="~/Documents/IowaWoundData2021/PaperFigs/GenusByAge.pdf", width=8, height=5)


# No association between wound age and genus shannon diversity or richness 
FullDataWithOtherClinicalDiversityAge = FullDataAge %>% select(woundage, Genus_Shannon, Genus_Richness)

FullDataWithOtherClinicalDiversityAge$woundage = if_else(FullDataWithOtherClinicalDiversityAge$woundage %in% c(1,2), "Acute", "Chronic")

wilcox.test(FullDataWithOtherClinicalDiversityAge$Genus_Shannon ~ FullDataWithOtherClinicalDiversityAge$woundage)
wilcox.test(FullDataWithOtherClinicalDiversityAge$Genus_Richness ~ FullDataWithOtherClinicalDiversityAge$woundage)


TwoWayAgeCluster = FullDataAge %>% select(DMMClusterAssign, woundage) %>% na.omit()
TwoWayAgeCluster$woundage = if_else(TwoWayAgeCluster$woundage %in% c(1, 2), "Acute", "Chronic")
chisq.test(table(TwoWayAgeCluster) )
write.csv(table(TwoWayAgeCluster), file="~/Documents/IowaWoundData2021/DMM_Age_whole.csv")

TwoWayAgeCluster = table(TwoWayAgeCluster)

TwoWayAgeCluster[,1] = round(TwoWayAgeCluster[,1]/sum(TwoWayAgeCluster[,1]), 3)
TwoWayAgeCluster[,2] = round(TwoWayAgeCluster[,2]/sum(TwoWayAgeCluster[,2]), 3)
write.csv((TwoWayAgeCluster), file="~/Documents/IowaWoundData2021/DMM_Age_freq.csv")


FullDataAge = FullDataAge %>% mutate(PainCatBinary = case_when( woundcarepain==0 | woundcarepain==1  ~"NoneMild",
                                                                woundcarepain==3~ "Severe", 
                                                                woundcarepain==2 ~ "Moderate" ))






AgeVsPain  = FullDataAge %>% select(PainCatBinary, woundage)
AgeVsPain$woundage = if_else(AgeVsPain$woundage %in% c(1, 2), "Acute", "Chronic")
AgeVsPain = AgeVsPain %>% filter(PainCatBinary != "Moderate")
chisq.test(table(AgeVsPain))

write.csv(table(AgeVsPain),  file="~/Documents/IowaWoundData2021/Pain_Age_whole_binary.csv")
TableAgePain = (table(AgeVsPain))

TableAgePain[,1] = round(TableAgePain[,1]/sum(TableAgePain[,1]),3)
TableAgePain[,2] = round(TableAgePain[,2]/sum(TableAgePain[,2]),3)
write.csv(TableAgePain,  file="~/Documents/IowaWoundData2021/Pain_Age_freq_binary.csv")


FullDataWithRestingAge = FullData %>% left_join(ClinicalData %>% select(study_id, resting_pain_cat, woundage), by="study_id")
FullDataWithRestingAge = FullDataWithRestingAge %>% filter(resting_pain_cat != "Moderate")
FullDataWithRestingAge$RestingBinary = if_else(FullDataWithRestingAge$resting_pain_cat=="Severe", "Severe", "NoneMild")

AgeVsRestPain  = FullDataWithRestingAge %>% select(RestingBinary, woundage)
AgeVsRestPain$woundage = if_else(AgeVsRestPain$woundage %in% c(1, 2), "Acute", "Chronic")
TableAgeRest = table(AgeVsRestPain)
chisq.test(table(AgeVsRestPain))
TableAgeRest[,1] = TableAgeRest[,1]/sum(TableAgeRest[,1])
TableAgeRest[,2] = TableAgePain[,2]/sum(TableAgeRest[,2])

# No association between wound age and any of the cytokines
###########################################################
FullDataAgeMeltCytokine = FullDataAge %>% select(listCytokines, woundage) %>% melt(id.vars=c("woundage"))
FullDataAgeMeltCytokine$woundage = if_else(FullDataAgeMeltCytokine$woundage %in% c(1,2), "Acute", "Chronic")
FullDataAgeMeltCytokine$Cytokine = sapply(FullDataAgeMeltCytokine$variable, function(x) (stringr::str_split(x, pattern="-"))[[1]][1] )
GenusStatsCytokinesAge = FullDataAgeMeltCytokine %>% group_by(Cytokine) %>% wilcox_test(value ~ woundage)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="age")


########################################################################################
# (4) BIOLOGICAL VARIABLES  VS. WOUND DRESSING TYPES
#########################################################################################


FullDataWithOtherClinicalsDressing = FullData %>% left_join(ClinicalData %>% select(dressingcat, study_id, woundcarepain,resting_pain_cat))

FullDataWithOtherClinicalsDressing = FullDataWithOtherClinicalsDressing %>% mutate(BinaryDressingType = if_else(dressingcat %in% c(0,2,3), "Other", "Adherent"))


FullDataWithOtherClinicalsDressing

FullDataWithOtherClinicalsDressingAnaerobes = FullDataWithOtherClinicalsDressing %>% select(AnaerobicGenusAbundance_CLR, dressingcat)
FullDataWithOtherClinicalsDressingAnaerobes$dressingcat = factor(FullDataWithOtherClinicalsDressingAnaerobes$dressingcat)

attach(FullDataWithOtherClinicalsDressingAnaerobes)
pairwise.wilcox.test( AnaerobicGenusAbundance_CLR,dressingcat, p.adjust.method = "none")
detach(FullDataWithOtherClinicalsDressingAnaerobes)

FullDataWithOtherClinicalsDressingMelt = FullDataWithOtherClinicalsDressing %>% select(AbundanceCLR, "BinaryDressingType") %>% melt(id.vars=c("BinaryDressingType"))

FullDataWithOtherClinicalsDressingMelt$Genus = sapply(FullDataWithOtherClinicalsDressingMelt$variable, function(x) str_split(x,pattern="Abundance_CLR")[[1]][1])
GenusStatsDressing = FullDataWithOtherClinicalsDressingMelt %>% group_by(Genus) %>% wilcox_test(value ~ BinaryDressingType)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="Wound  Age")
GenusStatsDressing$xmin=1
GenusStatsDressing$xmax=2
GenusDressingPlot = ggplot(FullDataWithOtherClinicalsDressingMelt, aes(x=BinaryDressingType, y=value, fill=Genus)) + geom_boxplot(alpha=.4,size=.25) + facet_grid(~ Genus)+
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2,size=.5) + theme_classic() + ggpubr::stat_pvalue_manual(GenusStatsDressing, label="p.adj") +
  ylim(0, 14) + xlab("Dressing (Adherent vs. Non-adherent\n and/or Woundvac)") + ggtitle("Common Genus Abundance in Wounds With Adherent \nvs. Non-adherent/Woundvac Dressings") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") +
  ylab("CLR-transformed relative abundance") +
  stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") 
ggsave(GenusDressingPlot, file="~/Documents/IowaWoundData2021/PaperFigs/GenusDressing.pdf", width=10, height=7) 



AllDressingTypesMelt = FullDataWithOtherClinicalsDressing %>% select(AbundanceCLR, dressingcat) %>% melt(id.vars=c("dressingcat"))
AllDressingTypesMelt$Genus = sapply(AllDressingTypesMelt$variable, function(x) str_split(x,pattern="Abundance_CLR")[[1]][1])

AllDressingTypesMelt$dressingcat = factor(AllDressingTypesMelt$dressingcat)
AllDressingTypesMelt$value =sapply(AllDressingTypesMelt$value, function(x) as.numeric(as.character(x)))
AllDressingTypesMelt$dressingcat = factor(AllDressingTypesMelt$dressingcat)
GenusStatsDressingAllTypesStat = AllDressingTypesMelt %>% group_by(Genus) %>% wilcox_test(value ~dressingcat) %>%  add_xy_position(x="Dressing type ")
GenusStatsDressingAllTypesStat = GenusStatsDressingAllTypesStat %>% mutate(xmin = case_when(group1 == 0 ~ 1,
                                                                                            group1 == 1 ~ 2, 
                                                                                            group1 == 2 ~ 3))
GenusStatsDressingAllTypesStat = GenusStatsDressingAllTypesStat %>% mutate(xmax = case_when(group2 == 1 ~ 2,
                                                                                            group2 == 2 ~ 3, 
                                                                                            group2 == 3 ~ 4))

GenusStatsDressingAllTypesStat$xmin=0
GenusStatsDressingAllTypesStat$xmax=1

AllDressingTypesMeltStaphAnaerobes = AllDressingTypesMelt %>% filter(Genus %in% c("AnaerobicGenus", "Staphylococcus") )
GenusStatsDressingAllTypesStat %>% filter(Genus %in% c("AnaerobicGenus", "Staphylococcus","Pseudomonas"))
GenusDressingPlotAllTypes = ggplot(AllDressingTypesMeltStaphAnaerobes, aes(x=dressingcat, y=value, fill=Genus)) + geom_boxplot(alpha=.4,size=.25) + facet_grid(~ Genus)+
  scale_fill_manual(values=c("#1B9E77", "#E7298A")) + geom_jitter(width=.2,size=.5) + theme_classic() +
  ylim(0, 14) + xlab("Dressing Category") + ggtitle("Common Genus Abundance in Wounds With Different Dressing Types") + 
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") +
  ylab("CLR-transformed relative abundance")  + stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray")


ggsave(GenusDressingPlotAllTypes, file="~/Documents/IowaWoundData2021/PaperFigs/GenusDressing_AllDressingTypes_StaphAnaerobes.pdf", height=7,width=8)

TwoWayDressingCluster = FullDataWithOtherClinicalsDressing %>% select(DMMClusterAssign,dressingcat ) %>% na.omit()
chisq.test(table(TwoWayDressingCluster))
TwoWayDressingClusterBinary = FullDataWithOtherClinicalsDressing %>% select(DMMClusterAssign,BinaryDressingType ) %>% na.omit()
chisq.test(table(TwoWayDressingClusterBinary))

TableDressingTypeBinary_DMMs = table(TwoWayDressingClusterBinary)
write.csv(TableDressingTypeBinary_DMMs,file="~/Documents/IowaWoundData2021/DressingType_DMM_Whole.csv")

TableDressingTypeBinary_DMMs[,1] = TableDressingTypeBinary_DMMs[,1]/sum(TableDressingTypeBinary_DMMs[,1])
TableDressingTypeBinary_DMMs[,2] = TableDressingTypeBinary_DMMs[,2]/sum(TableDressingTypeBinary_DMMs[,2])
write.csv(TableDressingTypeBinary_DMMs,file="~/Documents/IowaWoundData2021/DressingType_DMM_Freq.csv")

TwoWayDressingAll_DMM = table(TwoWayDressingCluster)
write.csv(TwoWayDressingAll_DMM,file="~/Documents/IowaWoundData2021/DressingType_All_DMM_whole.csv")

TwoWayDressingAll_DMM[,1] = round(TwoWayDressingAll_DMM[,1]/sum(TwoWayDressingAll_DMM[,1]),3)
TwoWayDressingAll_DMM[,2] = round(TwoWayDressingAll_DMM[,2]/sum(TwoWayDressingAll_DMM[,2]),3)
TwoWayDressingAll_DMM[,3] = round(TwoWayDressingAll_DMM[,3]/sum(TwoWayDressingAll_DMM[,3]),3)
TwoWayDressingAll_DMM[,4] = round(TwoWayDressingAll_DMM[,4]/sum(TwoWayDressingAll_DMM[,4]),3)

write.csv(TwoWayDressingAll_DMM,file="~/Documents/IowaWoundData2021/DressingType_All_DMM.csv")

WoundtypeDressing = ClinicalData %>% select(wound_type,dressingcat)
WoundtypeDressing = WoundtypeDressing %>% mutate(BinaryDressingType = if_else(dressingcat %in% c(0,2,3), "Other", "Adherent"))
WoundtypeDressing= WoundtypeDressing %>% select(BinaryDressingType,wound_type )
chisq.test(table(WoundtypeDressing))
WoundtypeDressing = table(WoundtypeDressing)
write.csv(WoundtypeDressing,file="~/Documents/IowaWoundData2021/DressingWoundTypeWhole.csv")


WoundtypeDressing[,1] = WoundtypeDressing[,1]/sum(WoundtypeDressing[,1])
WoundtypeDressing[,2] = WoundtypeDressing[,2]/sum(WoundtypeDressing[,2])
WoundtypeDressing[,3] = WoundtypeDressing[,3]/sum(WoundtypeDressing[,3])
WoundtypeDressing[,4] = WoundtypeDressing[,4]/sum(WoundtypeDressing[,4])
WoundtypeDressing[,5] = WoundtypeDressing[,5]/sum(WoundtypeDressing[,5])
WoundtypeDressing[,6] = WoundtypeDressing[,6]/sum(WoundtypeDressing[,6])

write.csv(TableWoundTypeDress,file="~/Documents/IowaWoundData2021/DressingWoundTypeFreq.csv")


WoundLocDressing = ClinicalData %>% select(dressingcat, woundloc)
WoundLocDressing = WoundLocDressing %>% mutate(BinaryDressingType = if_else(dressingcat %in% c(0,2,3), "Other", "Adherent"))
WoundLocDressingDF= WoundLocDressing %>% select( BinaryDressingType, woundloc)
chisq.test(table(WoundLocDressingDF))

WoundLocDressing = table(WoundLocDressingDF)
write.csv(WoundLocDressing,file="~/Documents/IowaWoundData2021/DressingWoundLocWhole.csv")

WoundLocDressing[,1] = round(WoundLocDressing[,1]/sum(WoundLocDressing[,1]),3)
WoundLocDressing[,2] = round(WoundLocDressing[,2]/sum(WoundLocDressing[,2]),3)
WoundLocDressing[,3] = round(WoundLocDressing[,3]/sum(WoundLocDressing[,3]),3)
WoundLocDressing[,4] = round(WoundLocDressing[,4]/sum(WoundLocDressing[,4]),3)
write.csv(WoundLocDressing,file="~/Documents/IowaWoundData2021/DressingWoundLocFreq.csv")



DressingFilter = FullDataWithOtherClinicalsDressing %>% select(DMMClusterAssign,dressingcat ) %>% melt(id.vars=c("dressingcat"))
FullDataWithOtherClinicalsDressing = FullDataWithOtherClinicalsDressing %>% filter(!is.na(DMMClusterAssign))
dressingplot = ggplot(FullDataWithOtherClinicalsDressing, aes(x=dressingcat, fill=factor(DMMClusterAssign))) + geom_bar(stat="count",position="fill")  +  scale_fill_brewer(palette ="Dark2" ) + theme_classic() + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))

dressingplot = dressingplot +    scale_y_continuous(labels = scales::percent) + labs(y="Percentage", x="Type of Dressing")
dressingplotTotal  =  ggplot(FullDataWithOtherClinicalsDressing, aes(x=dressingcat, fill=factor(DMMClusterAssign))) + geom_bar(stat="count")  +  scale_fill_brewer(palette ="Dark2" ) + theme_classic() + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))

ggsave(dressingplot, file="~/Documents/IowaWoundData2021/PaperFigs/DressingCategories.pdf",width=6, height=5)

DressingFilter = FullDataWithOtherClinicalsDressing %>% select(DMMClusterAssign,dressingcat ) %>% melt(id.vars=c("dressingcat")) 
FullDataWithOtherClinicalsDressing = FullDataWithOtherClinicalsDressing %>% na.omit()
dressingplot = ggplot(FullDataWithOtherClinicalsDressing, aes(x=dressingcat, fill=factor(DMMClusterAssign))) + geom_bar(stat="count",position="fill")  +  scale_fill_brewer(palette ="Dark2" ) + theme_classic() +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))
ggsave(dressingplot, file="~/Documents/IowaWoundData2021/PaperFigs/DressingCategories.pdf", width=6, height=6)


# Is the observed relationship between staphylococcus, anaerobic genera and dressing type due to wound depth?
WoundDepthData$study_id = as.character(WoundDepthData$study_id)
FullDataWithOtherClinicalsDressing = FullDataWithOtherClinicalsDressing %>% left_join(WoundDepthData, by="study_id")
correlation::cor_test(FullDataWithOtherClinicalsDressing,"wound_dp","AnaerobicGenusAbundance_CLR", method="spearman")
correlation::cor_test(FullDataWithOtherClinicalsDressing,"wound_dp","StaphylococcusAbundance_CLR", method="spearman")

ggplot(FullDataWithOtherClinicalsDressing, aes(x=wound_dp, y=AnaerobicGenusAbundance_CLR)) + geom_point()
ggplot(FullDataWithOtherClinicalsDressing, aes(x=wound_dp, y=StaphylococcusAbundance_CLR)) + geom_point()


FullDataWithOtherClinicalsDressing %>% wilcox_test(wound_dp ~ dressingcat)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="Wound  Age")
DepthByDressing = ggplot(FullDataWithOtherClinicalsDressing, aes(x=factor(dressingcat), y=wound_dp, fill=factor(dressingcat) ))+ geom_boxplot(size=.25, alpha=.4) + geom_jitter(width=.2) + scale_fill_brewer(palette="Dark2") +
  theme_classic() + xlab("Dressing Type") + ylab("Wound Depth(mm)")

ggsave(DepthByDressing,  file="~/Documents/IowaWoundData2021/PaperFigs/DressingByDepth.pdf", width=12,height=8)




FullDataWithOtherWoundType= FullData %>% left_join(ClinicalData %>% select(wound_type, study_id, woundcarepain,resting_pain_cat))

# GLM for wound depth and dressing (adherent mean = 4.048207 & 'other' mean = 2.620829 )
###################################

# https://sites.utexas.edu/sos/guided/inferential/numeric/glm/

FullDataWithOtherClinicalsDressingLocation = FullDataWithOtherClinicalsDressing %>% left_join(ClinicalData %>% select(study_id, woundloc))
FullDataWithOtherClinicalsDressingLocation$StaphylococcusAbundance_CLR

FullDataWithOtherClinicalsDressingLocation$wound_dp_cube = ((FullDataWithOtherClinicalsDressingLocation$wound_dp)^(1/3))
ModelWith=lm(StaphylococcusAbundance_CLR~ wound_dp_cube  + BinaryDressingType ,data=FullDataWithOtherClinicalsDressingLocation)
ModelWithout = lm(StaphylococcusAbundance_CLR~ wound_dp_cube ,data=FullDataWithOtherClinicalsDressingLocation)
anova(ModelWithout, ModelWith)

FullDataWithOtherClinicalsDressingNoNon = FullDataWithOtherClinicalsDressing %>% filter(dressingcat != 0)

# When you remove the 'non-adherent' category (no wound vacs), it takes away the relationship between depth and dressing type
wilcox.test(FullDataWithOtherClinicalsDressingNoNon$wound_dp ~ FullDataWithOtherClinicalsDressingNoNon$BinaryDressingType)

# it does not take away the relationship between anaerobic genera, staph genera CLR 
t.test(FullDataWithOtherClinicalsDressingNoNon$AnaerobicGenusAbundance_CLR ~ FullDataWithOtherClinicalsDressingNoNon$BinaryDressingType)
t.test(FullDataWithOtherClinicalsDressingNoNon$StaphylococcusAbundance_CLR ~ FullDataWithOtherClinicalsDressingNoNon$BinaryDressingType)


cor_test(FullDataWithOtherClinicalsDressingNoNon, x="AnaerobicGenusAbundance_CLR", y="wound_dp", method="spearman")
cor_test(FullDataWithOtherClinicalsDressingNoNon, x="StaphylococcusAbundance_CLR", y="wound_dp", method="spearman")






wilcox.test(FullDataWithOtherClinicalsDressing$wound_dp ~ FullDataWithOtherClinicalsDressing$BinaryDressingType)

ModelDressType = glm(StaphylococcusAbundance_CLR~BinaryDressingType ,data=FullDataWithOtherClinicalsDressing)

ModelDepthType = glm(StaphylococcusAbundance_CLR~ wound_dp + wound+  BinaryDressingType ,data=FullDataWithOtherClinicalsDressing)
summary(ModelDepthType)


ModelDepthType = glm(AnaerobicGenusAbundance_CLR~ wound_dp + BinaryDressingType ,data=FullDataWithOtherClinicalsDressing)
summary(ModelDepthType)

################################################
# (5) BIOLOGICAL VARIABLES  VS. WOUND  LOCATION
###############################################

# DMM Assignment
WoundLocDMM = FullData %>% left_join(ClinicalData %>% select(dressingcat, study_id, woundloc))
WoundLocDMM = table(WoundLocDMM %>% select(DMMClusterAssign, woundloc)) %>% na.omit()
chisq.test(WoundLocDMM,simulate.p.value = T)


write.csv(WoundLocDMM,file="~/Documents/IowaWoundData2021/WoundLocDMM_Whole.csv")

WoundLocDMM[,1] = round(WoundLocDMM[,1]/sum(WoundLocDMM[,1]),3)
WoundLocDMM[,2] = round(WoundLocDMM[,2]/sum(WoundLocDMM[,2]),3)
WoundLocDMM[,3] = round(WoundLocDMM[,3]/sum(WoundLocDMM[,3]),3)
WoundLocDMM[,4] = round(WoundLocDMM[,4]/sum(WoundLocDMM[,4]),3)

write.csv(WoundLocDMM,file="~/Documents/IowaWoundData2021/WoundLocDMM_Freq.csv")

WoundLoc_Genus = FullData %>% left_join(ClinicalData %>% select(study_id, woundloc))
WoundLoc_Genus_Melt = WoundLoc_Genus %>% select(AbundanceCLR, "woundloc") %>% melt(id.vars=c("woundloc"))
WoundLoc_Genus_Melt$Genus = sapply(WoundLoc_Genus_Melt$variable, function(x) str_split(x,pattern="Abundance_CLR")[[1]][1])

WoundLocStats = WoundLoc_Genus_Melt %>% group_by(Genus) %>%  wilcox_test(value ~ woundloc)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="Wound Location")


WoundLoc_Genus_Melt$woundloc = factor(WoundLoc_Genus_Melt$woundloc)
WoundLocGenusPlot = ggplot(WoundLoc_Genus_Melt, aes(x=woundloc, y=value, fill=Genus)) + geom_boxplot(alpha=.4,size=.25) + facet_grid(~ Genus)+
  scale_fill_brewer(palette="Dark2") + geom_jitter(width=.2,size=.5) + theme_classic() +
  ylim(0, 14) + xlab("Wound Location") + ggtitle("Common Genus Abundance in Wounds In Different Locations") + 
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") +
  ylab("CLR-transformed relative abundance")  + stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") 
ggsave(WoundLocGenusPlot, file="~/Documents/IowaWoundData2021/PaperFigs/DMM_WoundLoc.pdf", width=12,height=8)




################################################
# (6) BIOLOGICAL VARIABLES  VS. ONE ANOTHER
###############################################


cor_test(data=FullData,x= "AnaerobicGenusAbundance_CLR", y="StaphylococcusAbundance_CLR", method="spearman")

listCytokines = c("ARG1-Hs00163660_m1",  "C3-Hs00163811_m1",  "C5AR1-Hs00704891_s1", "CAMP-Hs00189038_m1",  "CXCL8-Hs00174103_m1", "IL1A-Hs00174092_m1",  "IL1B-Hs01555410_m1",
                  "IL6-Hs00174131_m1", "LCN2-Hs01008571_m1", "MMP1-Hs00899658_m1", "MMP2-Hs01548727_m1", "MMP9-Hs00957562_m1", "TNF-Hs00174128_m1")



# Cytokine data alone
#######################

listCytokines = c("ARG1-Hs00163660_m1",  "C3-Hs00163811_m1",  "C5AR1-Hs00704891_s1", "CAMP-Hs00189038_m1",  "CXCL8-Hs00174103_m1", "IL1A-Hs00174092_m1",  "IL1B-Hs01555410_m1",
"IL6-Hs00174131_m1", "LCN2-Hs01008571_m1", "MMP1-Hs00899658_m1", "MMP2-Hs01548727_m1", "MMP9-Hs00957562_m1", "TNF-Hs00174128_m1")
length(listCytokines)

# Cytokines vs. Microbiome data
###############################
listkruskalCytokinesDMM = c()
listkruskalPs = c()
for(cyt in listCytokines){
  testresult = kruskal.test(FullData[,cyt] ~ FullData[,"DMMClusterAssign"])  
  listkruskalPs = c(listkruskalPs, testresult$p.value)
  listkruskalCytokinesDMM = c(listkruskalCytokinesDMM, cyt)
  
}

DataFrameCytokineDMMs = data.frame(pvalsKruskalWallis =listkruskalPs, Cytokine=listkruskalCytokinesDMM )
DataFrameCytokineDMMs$PAdjustKW = p.adjust(DataFrameCytokineDMMs$pvalsKruskalWallis, method="BH")


FourCyts = FullData %>% select(c("C5AR1-Hs00704891_s1","CXCL8-Hs00174103_m1", "IL1B-Hs01555410_m1", "MMP2-Hs01548727_m1"))
GenusVars = FullData %>% select(AbundanceCLR)
DFCompare = cbind(FourCyts,GenusVars )


# four cytokine markers of interest vs. DMM
###########################################
# "C5AR1-Hs00704891_s1","CXCL8-Hs00174103_m1", "IL1B-Hs01555410_m1", "MMP2-Hs01548727_m1"
FullDataContainsMicrobiome = FullData %>% filter(!is.na(DMMClusterAssign))

colorpalette4 = RColorBrewer::brewer.pal(4, "Dark2")


C5ARStats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `C5AR1-Hs00704891_s1`) %>% wilcox_test(`C5AR1-Hs00704891_s1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p.adj<.05)
cxcl8Stats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `CXCL8-Hs00174103_m1`) %>% wilcox_test(`CXCL8-Hs00174103_m1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p.adj<.05)
ilbStats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `IL1B-Hs01555410_m1`) %>% wilcox_test(`IL1B-Hs01555410_m1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p.adj<.05)
MMP2Stats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `MMP2-Hs01548727_m1`) %>% wilcox_test(`MMP2-Hs01548727_m1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p.adj<.05)



c5ar1plot = ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`C5AR1-Hs00704891_s1` )),) + geom_boxplot(fill=colorpalette4[1]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab("-∆CT") + ggtitle("C5AR1 Expression") + theme(plot.title=element_text(hjust=.5, size=20, face="bold")) +
  ggpubr::stat_pvalue_manual(C5ARStats, label="p.adj",size=6) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

cxcl8plot  =  ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`CXCL8-Hs00174103_m1` ))) + geom_boxplot(fill=colorpalette4[2]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab("-∆CT")+ ggtitle("CXCL8 Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold")) + 
  ggpubr::stat_pvalue_manual(cxcl8Stats, label="p.adj",size=6)+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

ilbplot  =  ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`IL1B-Hs01555410_m1` ))) + geom_boxplot(fill=colorpalette4[3]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab("-∆CT")+ ggtitle("Il-1B Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold"))+
  ggpubr::stat_pvalue_manual(ilbStats, label="p.adj",size=6)+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

MMP2plot  =  ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`MMP2-Hs01548727_m1` ))) + geom_boxplot(fill=colorpalette4[4]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab("-∆CT")+ ggtitle("MMP2 Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold"))+
  ggpubr::stat_pvalue_manual(MMP2Stats, label="p.adj",size=6)+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

pdf(file="~/Documents/IowaWoundData2021/PaperFigs/FourCytokinesVsDMM.pdf", height=7,width=9)
gridExtra::grid.arrange(c5ar1plot, cxcl8plot, ilbplot, MMP2plot, ncol=2 )
dev.off()

# Genera vs. each cytokine
##########################
GeneraData = FullDataContainsMicrobiome %>% select(AbundanceCLR)
CytokineData = FullDataContainsMicrobiome %>% select(colnames(FourCyts))
CorrelationsGeneraCyt = matrix(0, nrow=ncol(GeneraData), ncol=ncol(CytokineData))
CorrelationsGeneraCytPvalues = matrix(0, nrow=ncol(GeneraData), ncol=ncol(CytokineData))


for(i in 1:ncol(GeneraData)){
  for(j in 1:ncol(CytokineData)){
    correlationtest = cor.test(GeneraData[, i], CytokineData[,j],method=c("spearman"),na.omit=T,exact=F )

    CorrelationsGeneraCyt[i,j] = correlationtest$estimate
    CorrelationsGeneraCytPvalues[i,j] = correlationtest$p.value
    
  }
}
PvalFrame = data.frame(CorrelationsGeneraCytPvalues)
rownames(PvalFrame) = colnames(GeneraData)
colnames(PvalFrame) = colnames(CytokineData)

PvalFrame$Genera = row.names(PvalFrame)
PvalFrameMelt = PvalFrame %>% reshape2::melt()
PvalFrameMelt$pvalue = PvalFrameMelt$value
PvalFrameMelt = PvalFrameMelt %>% select(Genera, variable, pvalue)

CorrEstframe = data.frame(CorrelationsGeneraCyt)
rownames(CorrEstframe) = colnames(GeneraData)
colnames(CorrEstframe) = colnames(CytokineData)
CorrEstframe$Genera=rownames(CorrEstframe)
CorrelationsMeltedCytGenera = CorrEstframe %>% reshape2::melt()

CorrelationsMeltedCytGenera = CorrelationsMeltedCytGenera %>% left_join(PvalFrameMelt, by=c("Genera", "variable"))
CorrelationsMeltedCytGenera$pvalueRounded = round(CorrelationsMeltedCytGenera$pvalue, 6)

CorrHeatMapGeneraCyt = ggplot(CorrelationsMeltedCytGenera, aes(x=variable, y=Genera, fill=value))+ geom_tile() + scale_fill_gradient2(high="red", low="blue", mid="white", midpoint=0) + xlab("") + ylab("") +
  ggtitle("Correlations Between Inflammatory Marker Expression (-∆CT)\n and Common Wound Genera ") +theme_classic()+ theme(plot.title=element_text(size=20, hjust=.5), axis.text.x=element_text(size=15, angle=270), axis.text.y=element_text(size=15)) +
  labs(fill="Spearman Correlation") + geom_text(aes(label=pvalueRounded)) + coord_equal()


ggsave(CorrHeatMapGeneraCyt, file="~/Documents/IowaWoundData2021/PaperFigs/CorrelationGeneraCytokines.pdf", width=14,height=14)


##########################################
# Cytokine correlations with one another
##########################################
JustCyt = (FullData %>% select(listCytokines) )

colnames(JustCyt) = sapply(colnames(JustCyt) , function(x) str_split(x, pattern="-")[[1]][1])
JustCytCor = cor(JustCyt, use="pairwise.complete.obs")

CorrelationMat = JustCytCor %>% melt()
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


CorrelationMat = CorrelationMat %>% mutate(significance = case_when(pAdj < .05 & pAdj >= .01 ~ "*", 
                                                                    pAdj < .01 & pAdj >= .001 ~ "**",
                                                                    pAdj < .001 & pAdj >= .0001 ~ "***",
                                                                    pAdj < .0001 ~ "****",
                                                                    TRUE ~  "")  )




CorrHeatMap = ggplot(CorrelationMat, aes(x=Var1, y=Var2, fill=value))+ geom_tile() + scale_fill_gradient2(high="red", low="blue", mid="white", midpoint=0) + xlab("") + ylab("") +
  ggtitle("Correlations Between Inflammatory\nMarker Expression (∆CT)") +theme_classic()+ theme(plot.title=element_text(size=20, hjust=.5), axis.text.x=element_text(size=15, angle=270), axis.text.y=element_text(size=15)) + labs(fill="Pearson's Correlation") +
  coord_equal()

Hclustering = hclust(as.dist(1-JustCytCor))

CorrHeatMap$data$Var1 = factor(CorrHeatMap$data$Var1, levels=(Hclustering$labels)[Hclustering$order])
CorrHeatMap$data$Var2 = factor(CorrHeatMap$data$Var2, levels=(Hclustering$labels)[Hclustering$order])

ggsave(CorrHeatMap, file="~/Documents/IowaWoundData2021/PaperFigs/CorrelationCytokinesVCytokines.pdf", width=14,height=14)


