# Amy Campbell
# Figure 4 
# The heterogeneous wound microbiome varies with wound care pain, dressing practices, and inflammatory gene expression 
library("dplyr")
library("stringr")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("sanon")


outputfolder="/Users/amycampbell/Documents/IowaWoundData2021/PublicationFiguresOutput/"
Cytokines_And_Microbiome = read.csv("Documents/IowaWoundData2021/PublicationData/WoundMicrobiome_Cytokine_Data_Final.csv")
ClinicalData =  read.csv("Documents/IowaWoundData2021/PublicationData/ClinicalData_UpdatedPain.csv")
WoundDepthData = read.csv("~/Documents/IowaWoundData2021/wound_depth_covariate.csv")


AbundanceCLR = colnames(FullData)[grepl(colnames(FullData), pattern="_CLR")]



FullData = Cytokines_And_Microbiome %>% left_join(ClinicalData, by="study_id")
FullData = FullData %>% mutate(BinaryDressingType = if_else(dressingcat %in% c(2,3), "WoundVac", "NonWoundVac"))

FullData = FullData %>% left_join(WoundDepthData, by="study_id")
FullData = FullData %>% mutate(Wound_Type = case_when(wound_type==1 ~ "Pressure ulcer",
                                                      wound_type ==2 ~"Venous ulcer",
                                                      wound_type==3 ~"Arterial ulcer",
                                                      wound_type==4 ~"Surgical", 
                                                      wound_type==5 ~ "Traumatic",
                                                      wound_type==6 ~" Mixed(Traumatic/Surg)",
                                                      wound_type==7 ~ "Other"))


 
woundtypepal=c("#D55E00","#009E73","0072B2", "#E69F00","#CC79A7", "#56B4E9")

WoundDepthBoxplot = ggplot(FullData, aes(x=BinaryDressingType, y=wound_dp)) + geom_boxplot(fill="white") +
  geom_jitter(height=0, aes(color=Wound_Type)) + stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") +
  scale_color_manual(values=woundtypepal) + theme_classic() + ggpubr::stat_compare_means(method="wilcox.test") + labs(y="Wound Depth(mm)", x="Use of Woundvac", color="Wound Etiology")
ggsave(WoundDepthBoxplot,file=paste0(outputfolder, "Figure4A.pdf"),width=6,height=9)


# Subset to just surgical & mixed surgical/traumatic wounds & plot
##################################################################
Surgical_Mixed_Only = FullData %>% filter(wound_type %in% c(4,6)) %>% filter(!is.na(DMMClusterAssign))
Surgical_Mixed_Only_Melt = Surgical_Mixed_Only %>% select(AbundanceCLR, "BinaryDressingType") %>% melt(id.vars=c("BinaryDressingType"))
Surgical_Mixed_Only_Melt$Genus = sapply(Surgical_Mixed_Only_Melt$variable, function(x) str_remove(x,"Abundance_CLR"))

Plot_surgical = ggplot(Surgical_Mixed_Only_Melt, aes(x=BinaryDressingType, y=value, fill=Genus)) + geom_boxplot(alpha=.4,size=.25) + facet_grid(~ Genus)+
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2,size=.5) + theme_classic()  +
  ylim(0, 14) + xlab("Dressing (Woundvac vs. non-Woundvac)") + ggtitle("Common Genus Abundance in Wounds With Woundvac \nvs. Non-Woundvac Dressings(Surgical/Mixed Only)") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") +
  ylab("CLR-transformed relative abundance")

ggsave(Plot_surgical, file=paste0(outputfolder, "Figure4B.pdf"), width=9, height=7)

# Wound vac in surgical/mixed wounds vs. genus, stratified by sequencing run
############################################################################
Coryne_dressing = sanon(CorynebacteriumAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)
Strep_dressing = sanon(StreptococcusAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)
Staph_dressing = sanon(StaphylococcusAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)
Pseudomonas_dressing = sanon(PseudomonasAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)
Anaerobes_dressing = sanon(AnaerobicGenusAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)


# Controlling for batch using linear regression models
######################################################

modelDressingStaph = glm(StaphylococcusAbundance_CLR ~ Run+ BinaryDressingType, data=Surgical_Mixed_Only)
modelDressingPseud = glm(PseudomonasAbundance_CLR ~ Run+ BinaryDressingType, data=Surgical_Mixed_Only)

