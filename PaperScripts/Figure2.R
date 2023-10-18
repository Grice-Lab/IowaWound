
# Amy Campbell
# The heterogeneous wound microbiome varies with wound care pain, dressing practices, and inflammatory gene expression 
# Scripts to reproduce figure 2 contents 

library("dplyr")
library("stringr")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("viridis")
library("gplots")
library("rstatix")
library("gridExtra")

#library("correlation")
#library("psych")
#library("sanon")
outputfolder="/Users/amycampbell/Documents/IowaWoundData2021/PublicationFiguresOutput/"


Cytokines_And_Microbiome = read.csv("Documents/IowaWoundData2021/PublicationData/WoundMicrobiome_Cytokine_Data_Final.csv")
ClinicalData =  read.csv("Documents/IowaWoundData2021/PublicationData/ClinicalData_UpdatedPain.csv")


listCytokines = c("ARG1.Hs00163660_m1",  "C3.Hs00163811_m1",  "C5AR1.Hs00704891_s1", "CAMP.Hs00189038_m1",  "CXCL8.Hs00174103_m1", "IL1A.Hs00174092_m1",  "IL1B.Hs01555410_m1",
                  "IL6.Hs00174131_m1", "LCN2.Hs01008571_m1", "MMP1.Hs00899658_m1", "MMP2.Hs01548727_m1", "MMP9.Hs00957562_m1", "TNF.Hs00174128_m1")

color21 =c("#7570B3", "#E7298A", "#66C2A5","#E31A1C",
           "#BC80BD", "#A6D854", "#FDCDAC", "#A6761D", "#80B1D3",
           "#FF9933", "#F4CAE4", "#FFFF66", "#2E8B57", "#CCEBC5", 
           "#1F78B4", "#BEBADA","#B3B3B3","#FF6600", "#800000","#23297A", "#FB8072")
##########################################
# Figure 2C: Wound care pain vs. cytokines
##########################################

Cytokines_And_Microbiome = Cytokines_And_Microbiome %>% left_join(ClinicalData, by="study_id") %>% 
 mutate(PainCatBinary = case_when(woundcarepain %in% c(0,1) ~ "NoneMild",
                                  woundcarepain==2 ~ "Moderate",
                                  woundcarepain==3 ~ "Severe"))



Cytokines_And_Microbiome_NoModerates = Cytokines_And_Microbiome %>% filter(woundcarepain!=2)
listPvalues = c()
listcytokinesTested = c()
listPvaluesNoneVsSevere = c()
kruskal_ps = c()
listpvalues_Age = c()
for(cyt in listCytokines){
  print(cyt)
  testresult = (wilcox.test((Cytokines_And_Microbiome_NoModerates[,cyt]) ~ Cytokines_And_Microbiome_NoModerates[,"PainCatBinary"],exact=F))
  listcytokinesTested = c(listcytokinesTested, cyt)
  listPvalues = c(listPvalues, testresult$p.value)
 
}
WilcoxPvalues = data.frame(cytokine=listcytokinesTested, pvalues=listPvalues)

numCytokines = Cytokines_And_Microbiome %>% select(listCytokines) 
studyids = Cytokines_And_Microbiome$study_id
numCytokines[!is.na(numCytokines)]<-1
numCytokines[is.na(numCytokines)]<-0
numCytokines$study_id=studyids
numCytokines$Ncyt = rowSums(numCytokines[,1: (ncol(numCytokines) - 1) ])
Cytokines_And_Microbiome = Cytokines_And_Microbiome %>% left_join(numCytokines %>% select(study_id, Ncyt), by="study_id")

orderHeatMap = ((Cytokines_And_Microbiome %>% filter(Ncyt>0) %>% arrange(woundcarepain, Ncyt))$study_id)

Cytokine_studyID = Cytokines_And_Microbiome %>% select(listCytokines, study_id)
MeltCytokine_studyID = Cytokine_studyID %>% reshape2::melt(id.vars=(c("study_id")))
colnames(MeltCytokine_studyID) = c("study_id", "cytokine","value" )


CytokineBySubjectPlot = ggplot(MeltCytokine_studyID, aes(x=study_id, y=cytokine, fill=(value))) + geom_tile() + scale_fill_viridis(option="plasma", na.value="white")

cytokinesorder =  rev(c("CAMP.Hs00189038_m1","IL6.Hs00174131_m1",  "ARG1.Hs00163660_m1", "C3.Hs00163811_m1",
                    "MMP1.Hs00899658_m1","MMP2.Hs01548727_m1",  "LCN2.Hs01008571_m1","IL1A.Hs00174092_m1",
                    "TNF.Hs00174128_m1", "IL1B.Hs01555410_m1", "MMP9.Hs00957562_m1",  "CXCL8.Hs00174103_m1",
                    "C5AR1.Hs00704891_s1"))
CytokineBySubjectPlot$data$study_id = factor(CytokineBySubjectPlot$data$study_id, levels=orderHeatMap)
CytokineBySubjectPlot$data$cytokine = factor(CytokineBySubjectPlot$data$cytokine, levels=cytokinesorder)
ggsave(CytokineBySubjectPlot, file=paste0(outputfolder,"Figure2C.pdf"))






JustCyt = (Cytokines_And_Microbiome %>% select(listCytokines) )

colnames(JustCyt) = sapply(colnames(JustCyt) , function(x) str_split(x, pattern="\\.")[[1]][1])
JustCytCor = cor(JustCyt, method="spearman", use="pairwise.complete.obs")


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
  ggtitle("Correlations Between Inflammatory\nMarker Expression (∆CT)") +theme_classic()+ theme(plot.title=element_text(size=20, hjust=.5), axis.text.x=element_text(size=15, angle=270), axis.text.y=element_text(size=15)) + labs(fill="Spearman's Correlation") +
  coord_equal() #+ geom_text(aes(label=significance),color="white", size=10)

Hclustering = hclust(as.dist(1-JustCytCor))

CorrHeatMap$data$Var1 = factor(CorrHeatMap$data$Var1, levels=(Hclustering$labels)[Hclustering$order])
CorrHeatMap$data$Var2 = factor(CorrHeatMap$data$Var2, levels=(Hclustering$labels)[Hclustering$order])

ggsave(CorrHeatMap, file="~/Documents/IowaWoundData2021/PublicationFiguresOutput/Figure2A.pdf", width=14,height=14)





# Figure 2B
# Now, cytokines vs. inflammation
##########################
Cytokines_Inflammation = Cytokines_And_Microbiome %>% select(inflame, listCytokines)
listPvaluesInflame = c()
listcytokinesTestedInflame = c()
listPvaluesInflame= c()
for(cyt in listCytokines){
  
  testresult = (wilcox.test((Cytokines_Inflammation[,cyt]) ~ Cytokines_Inflammation[,"inflame"],exact=F))
  
  listcytokinesTestedInflame = c(listcytokinesTestedInflame, cyt)
  listPvaluesInflame = c(listPvaluesInflame, testresult$p.value)
  
}

DataFrameCytokineInflame = data.frame( Cytokine=listcytokinesTested,pvals = listPvaluesInflame,padj = p.adjust(listPvaluesInflame, method="bonferroni"))
DataFrameCytokineInflame$variable = sapply(DataFrameCytokineInflame$Cytokine,  function(x) (stringr::str_split(x, pattern="-"))[[1]][1])

MeltedInflammationMarkers = Cytokines_Inflammation %>% melt(id.vars=c("inflame"))
MeltedInflammationMarkers$Cytokine = sapply(MeltedInflammationMarkers$variable, function(x) (stringr::str_split(x, pattern="\\."))[[1]][1] )

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
ggsave(boxplotinflame, file="~/Documents/IowaWoundData2021/PublicationFiguresOutput/Figure2b.pdf", width=10, height=8)
ggsave(boxplotinflame, file="~/Documents/IowaWoundData2021/NewFigs_Paper_10_23/new_Figure2a.pdf", width=15, height=6)

# Reviewer #1 asked whether there was more of a relationship between cytokines and pain WITHIN each of <30 and >30 day wounds
#############################################################################################################################
Cytokines_And_Microbiome_Pre30 = Cytokines_And_Microbiome %>% filter(woundage %in% c(1,2))
Cytokines_And_Microbiome_Post30 = Cytokines_And_Microbiome %>% filter(!(woundage %in% c(1,2)))


Cytokines_And_Microbiome_Pre30 = Cytokines_And_Microbiome_Pre30 %>% filter(PainCatBinary!="Moderate")
Cytokines_And_Microbiome_Post30 = Cytokines_And_Microbiome_Post30 %>% filter(PainCatBinary!="Moderate")



listPvalues_pre30 = c()
listcytokinesTested_pre30 = c()
listPvaluesNoneVsSevere_pre30 = c()
kruskal_p_pre30 = c()
for(cyt in listCytokines){
  print(cyt)
  testresult = (wilcox.test((Cytokines_And_Microbiome_Pre30[,cyt]) ~ Cytokines_And_Microbiome_Pre30[,"PainCatBinary"],exact=F))
  listcytokinesTested_pre30 = c(listcytokinesTested_pre30, cyt)
  listPvalues_pre30 = c(listPvalues_pre30, testresult$p.value)
  
}
WilcoxPvalues_pre30 = data.frame(cytokine=listcytokinesTested_pre30, pvalues=listPvalues_pre30)
WilcoxPvalues_pre30$pAdj = p.adjust(WilcoxPvalues_pre30$pvalues,method="BH")


listPvalues_post30 = c()
listcytokinesTested_post30 = c()
listPvaluesNoneVsSevere_post30 = c()
kruskal_p_post30 = c()
for(cyt in listCytokines){
  print(cyt)
  testresult = (wilcox.test((Cytokines_And_Microbiome_Post30[,cyt]) ~ Cytokines_And_Microbiome_Post30[,"PainCatBinary"],exact=F))
  listcytokinesTested_post30 = c(listcytokinesTested_post30, cyt)
  listPvalues_post30 = c(listPvalues_post30, testresult$p.value)
  
}
WilcoxPvalues_post30 = data.frame(cytokine=listcytokinesTested_post30, pvalues=listPvalues_post30)
WilcoxPvalues_post30$pAdj = p.adjust(WilcoxPvalues_post30$pvalues,method="BH")

Cytokines_And_Microbiome_Pre30_Melt = Cytokines_And_Microbiome_Pre30 %>% select(listCytokines, PainCatBinary, wound_type) %>% melt(id.vars=c("wound_type","PainCatBinary"))
Cytokines_And_Microbiome_Post30_Melt = Cytokines_And_Microbiome_Post30 %>% select(listCytokines, PainCatBinary, wound_type) %>% melt(id.vars=c("wound_type","PainCatBinary"))


Cytokines_And_Microbiome_Pre30_Melt$Cytokine = sapply(Cytokines_And_Microbiome_Pre30_Melt$variable, function(x) str_split(string=x,pattern=".H")[[1]][1])
PainCytokine_Pre30 = Cytokines_And_Microbiome_Pre30_Melt %>% group_by(Cytokine) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary") 
PainCytokine_Pre30$xmin = 1
PainCytokine_Pre30$xmax = 2

Cytokines_And_Microbiome_Post30_Melt$Cytokine = sapply(Cytokines_And_Microbiome_Post30_Melt$variable, function(x) str_split(string=x,pattern=".H")[[1]][1])
PainCytokine_post30 = Cytokines_And_Microbiome_Post30_Melt %>% group_by(Cytokine) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary") 
PainCytokine_post30$xmin = 1
PainCytokine_post30$xmax = 2

Pre30_plot = ggplot(Cytokines_And_Microbiome_Pre30_Melt, aes(x=PainCatBinary,y=value))+ geom_boxplot() + geom_jitter(aes(color=factor(wound_type)),height=0,width=.2) + facet_grid(~Cytokine) + 
  theme_classic() + ggpubr::stat_pvalue_manual(PainCytokine_Pre30,label="p") + scale_color_manual(values=c( "#7570B3", "#FF9933", "#66C2A5","#A6D854", "#800000", "#1F78B4")) #+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") + ylab("-∆CT")


Post30_plot = ggplot(Cytokines_And_Microbiome_Post30_Melt, aes(x=PainCatBinary,y=value))+ geom_boxplot() + geom_jitter( aes(color=factor(wound_type)), height=0, width=.2) + facet_grid(~Cytokine) + 
  theme_classic() + ggpubr::stat_pvalue_manual(PainCytokine_post30,label="p")+  scale_color_manual(values=c( "#7570B3", "#FF9933", "#66C2A5","#A6D854", "#800000", "#1F78B4")) #+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") + ylab("-∆CT")
 #+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") + ylab("-∆CT")
pdf(file="~/Documents/IowaWoundData2021/NewFigs_Paper_10_23/new_FigureS2.pdf",width=14, height=14)
grid.arrange(Pre30_plot, Post30_plot, ncol=1)
dev.off()



# Is baseline inflammation different between ≤30 and >30 day old wounds?
########################################################################

MeltedInflammationMarkers = Cytokines_Inflammation %>% melt(id.vars=c("inflame"))
MeltedInflammationMarkers$Cytokine = sapply(MeltedInflammationMarkers$variable, function(x) (stringr::str_split(x, pattern="\\."))[[1]][1] )

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
ggsave(boxplotinflame, file="~/Documents/IowaWoundData2021/PublicationFiguresOutput/Figure2b.pdf", width=10, height=8)
ggsave(boxplotinflame, file="~/Documents/IowaWoundData2021/NewFigs_Paper_10_23/new_Figure2a.pdf", width=15, height=6)

