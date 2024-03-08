# Amy Campbell
# Figure 4 + Figure S5
# The heterogeneous wound microbiome varies with wound care pain, dressing practices, and inflammatory gene expression 
library("dplyr")
library("stringr")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("sanon")


outputfolder="OutputFigures/"
Cytokines_And_Microbiome = read.csv("InputData/WoundMicrobiome_Cytokine_Data_Final.csv")
ClinicalData =  read.csv("InputData/ClinicalData_UpdatedPain.csv")
WoundDepthData = read.csv("InputData/wound_depth_covariate.csv")


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
ggsave(WoundDepthBoxplot,file=paste0(outputfolder, "OutputFigures/Figure4A.pdf"),width=6,height=9)


# Subset to just surgical & mixed surgical/traumatic wounds & plot
##################################################################
Surgical_Mixed_Only = FullData %>% filter(wound_type %in% c(4,6)) %>% filter(!is.na(DMMClusterAssign))
Surgical_Mixed_Only_Melt = Surgical_Mixed_Only %>% select(AbundanceCLR, "BinaryDressingType") %>% melt(id.vars=c("BinaryDressingType"))
Surgical_Mixed_Only_Melt$Genus = sapply(Surgical_Mixed_Only_Melt$variable, function(x) str_remove(x,"Abundance_CLR"))

Plot_surgical = ggplot(Surgical_Mixed_Only_Melt, aes(x=BinaryDressingType, y=value, fill=Genus)) + geom_boxplot(alpha=.4,size=.25) + facet_grid(~ Genus)+
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2,size=.5) + theme_classic()  +
  ylim(-1.2, 12.5) + xlab("Dressing (Woundvac vs. non-Woundvac)") + ggtitle("Common Genus Abundance in Wounds With Woundvac \nvs. Non-Woundvac Dressings(Surgical/Mixed Only)") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") +
  ylab("CLR-transformed relative abundance")

ggsave(Plot_surgical, file=paste0(outputfolder, "OutputFigures/Figure4C.pdf"), width=9, height=7)

# Wound vac in surgical/mixed wounds vs. genus, stratified by sequencing run
############################################################################
Coryne_dressing = sanon(CorynebacteriumAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)
Strep_dressing = sanon(StreptococcusAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)
Staph_dressing = sanon(StaphylococcusAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)
Pseudomonas_dressing = sanon(PseudomonasAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)
Anaerobes_dressing = sanon(AnaerobicGenusAbundance_CLR ~ grp(BinaryDressingType, ref="NonWoundVac") + strt(Run) , data=Surgical_Mixed_Only)


# Woundvac vs. not by DMM cluster
##################################
FullDataMicrobiome = FullData %>% filter(!is.na(DMMClusterAssign))
ggplot(FullDataMicrobiome,aes(x=BinaryDressingType, fill=factor(DMMClusterAssign))) + geom_bar(position="fill") + scale_fill_brewer(palette="Dark2") + theme_

FullDataMicrobiomeSurgical = FullData %>% filter(!is.na(DMMClusterAssign)) %>% filter(wound_type %in% c(4,6))
surgical_pct = ggplot(FullDataMicrobiomeSurgical,aes(x=BinaryDressingType, fill=factor(DMMClusterAssign))) + geom_bar(position="fill") + scale_fill_brewer(palette="Dark2") + theme_classic()
ggsave(surgical_pct, file="OutputFigures/Figure4B.pdf",width=10,height=6)


# Inflammation and pain in just Non-VAC trunk wounds (reviewer 3 asked)
########################################################################

# 171 wounds in these categories
FullDataMicrobiome_Trunk = FullDataMicrobiome %>% filter(woundloc==2 & BinaryDressingType=="NonWoundVac" )
FullDataMicrobiome_Trunk = FullDataMicrobiome_Trunk %>% mutate(PainCatBinary = case_when( woundcarepain==0 | woundcarepain==1  ~"NoneMild",
                                                          woundcarepain==3~ "Severe", 
                                                          woundcarepain==2 ~ "Moderate" 
))

# 111 after subsetting
FullDataMicrobiome_Trunk = FullDataMicrobiome_Trunk %>% filter(PainCatBinary!="Moderate")

Genera = c("CorynebacteriumAbundance_CLR",
           "StreptococcusAbundance_CLR" ,
           "StaphylococcusAbundance_CLR",
           "PseudomonasAbundance_CLR",
           "AnaerobicGenusAbundance_CLR" )
StratifiedPvaluesTrunk = c()
StratifiedEstimatesTrunk = c()

# Corynebacterium
Coryne_Pain_strat = sanon(CorynebacteriumAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryCorynePain = summary(Coryne_Pain_strat)
pvalues_genus_stratifiedTrunk = append(StratifiedPvaluesTrunk, SummaryCorynePain$coefficients[4])
StratifiedEstimatesTrunk = append(StratifiedEstimatesTrunk, Coryne_Pain_strat$xi)

# Streptococcus 
Streptococcus_Pain_Strat = sanon(StreptococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryStrepPain = summary(Streptococcus_Pain_Strat)
pvalues_genus_stratifiedTrunk= append(pvalues_genus_stratifiedTrunk, SummaryStrepPain$coefficients[4])
StratifiedEstimates = append(StratifiedEstimatesTrunk, Streptococcus_Pain_Strat$xi)


# Staphylococcus
Staphylococcus_Pain_strat = sanon(StaphylococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryStaphPain = summary(Staphylococcus_Pain_strat)
pvalues_genus_stratifiedTrunk = append(pvalues_genus_stratifiedTrunk, SummaryStaphPain$coefficients[4])
StratifiedEstimatesTrunk = append(StratifiedEstimatesTrunk, Staphylococcus_Pain_strat$xi)

# Pseudomonas
Pseudomonas_Pain_strat = sanon(PseudomonasAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryPseudomonasPain = summary(Pseudomonas_Pain_strat)
pvalues_genus_stratifiedTrunk = append(pvalues_genus_stratifiedTrunk, SummaryPseudomonasPain$coefficients[4])
StratifiedEstimatesTrunk = append(StratifiedEstimatesTrunk, Pseudomonas_Pain_strat$xi)

# anaerobes
Anaerobes_Pain_strat = sanon(AnaerobicGenusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryAnaerobesPain = summary(Anaerobes_Pain_strat)
pvalues_genus_stratifiedTrunk = append(pvalues_genus_stratifiedTrunk, SummaryAnaerobesPain$coefficients[4])
StratifiedEstimatesTrunk = append(StratifiedEstimatesTrunk, Anaerobes_Pain_strat$xi)


# Figure S5B
############

# 171 wounds in these categories
FullDataMicrobiome_Trunk = FullDataMicrobiome %>% filter(woundloc==2 & BinaryDressingType=="NonWoundVac" )
FullDataMicrobiome_Trunk = FullDataMicrobiome_Trunk %>% mutate(PainCatBinary = case_when( woundcarepain==0 | woundcarepain==1  ~"NoneMild",
                                                                                          woundcarepain==3~ "Severe", 
                                                                                          woundcarepain==2 ~ "Moderate" 
))

# 111 after subsetting
FullDataMicrobiome_Trunk = FullDataMicrobiome_Trunk %>% filter(PainCatBinary!="Moderate")

Genera = c("CorynebacteriumAbundance_CLR",
           "StreptococcusAbundance_CLR" ,
           "StaphylococcusAbundance_CLR",
           "PseudomonasAbundance_CLR",
           "AnaerobicGenusAbundance_CLR" )
StratifiedPvaluesTrunkInflammation = c()
StratifiedEstimatesTrunkInflammation = c()

FullDataMicrobiome_Trunk_Melt = FullDataMicrobiome_Trunk %>% select(PainCatBinary, AbundanceCLR) %>%  melt(id.vars=c("PainCatBinary"))
FullDataMicrobiome_Trunk_Melt$Genus = sapply(FullDataMicrobiome_Trunk_Melt$variable, function(x) str_remove_all(x, "Abundance_CLR"))
Plot_justTrunk_NonVAC= ggplot(FullDataMicrobiome_Trunk_Melt, aes(x=PainCatBinary, y=value, fill=Genus)) + geom_boxplot(alpha=.4,size=.25) + facet_grid(~ Genus)+
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2,height=0,size=.5) + theme_classic()  +
  ylim(-1.2, 12.5) + xlab("Pain") + ggtitle("Common Genus Abundance in Wounds With Severe vs. Mild/None Pain (Non-VAC Trunk Wounds, n=111)") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") +
  ylab("CLR-transformed relative abundance")+  stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") 

ggsave(Plot_justTrunk_NonVAC, file=paste0(outputfolder, "OutputFigures/FigureS5B.pdf"), width=9, height=7)

# Inflammation 
################
# Corynebacterium
Coryne_Pain_strat = sanon(CorynebacteriumAbundance_CLR ~ grp(inflame, ref="0") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryCorynePain = summary(Coryne_Pain_strat)
pvalues_genus_stratifiedTrunkInflammation = append(StratifiedPvaluesTrunkInflammation, SummaryCorynePain$coefficients[4])
StratifiedEstimatesTrunkInflammation = append(StratifiedEstimatesTrunkInflammation, Coryne_Pain_strat$xi)

# Streptococcus 
Streptococcus_Pain_Strat = sanon(StreptococcusAbundance_CLR ~ grp(inflame, ref="0") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryStrepPain = summary(Streptococcus_Pain_Strat)
pvalues_genus_stratifiedTrunkInflammation= append(pvalues_genus_stratifiedTrunkInflammation, SummaryStrepPain$coefficients[4])
StratifiedEstimatesInflammation = append(StratifiedEstimatesTrunkInflammation, Streptococcus_Pain_Strat$xi)

# Staphylococcus
Staphylococcus_Pain_strat = sanon(StaphylococcusAbundance_CLR ~ grp(inflame, ref="0") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryStaphPain = summary(Staphylococcus_Pain_strat)
pvalues_genus_stratifiedTrunkInflammation = append(pvalues_genus_stratifiedTrunkInflammation, SummaryStaphPain$coefficients[4])
StratifiedEstimatesTrunkInflammation = append(StratifiedEstimatesTrunkInflammation, Staphylococcus_Pain_strat$xi)

# Pseudomonas
Pseudomonas_Pain_strat = sanon(PseudomonasAbundance_CLR ~ grp(inflame, ref="0") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryPseudomonasPain = summary(Pseudomonas_Pain_strat)
pvalues_genus_stratifiedTrunkInflammation = append(pvalues_genus_stratifiedTrunkInflammation, SummaryPseudomonasPain$coefficients[4])
StratifiedEstimatesTrunkInflammation = append(StratifiedEstimatesTrunkInflammation, Pseudomonas_Pain_strat$xi)

# anaerobes
Anaerobes_Pain_strat = sanon(AnaerobicGenusAbundance_CLR ~ grp(inflame, ref="0") + strt(Run) , data=FullDataMicrobiome_Trunk)
SummaryAnaerobesPain = summary(Anaerobes_Pain_strat)
pvalues_genus_stratifiedTrunkInflammation = append(pvalues_genus_stratifiedTrunkInflammation, SummaryAnaerobesPain$coefficients[4])
StratifiedEstimatesTrunkInflammation = append(StratifiedEstimatesTrunkInflammation, Anaerobes_Pain_strat$xi)

# Figure S5A
######################
FullDataMicrobiome = FullDataMicrobiome %>% mutate(PainCatBinary = case_when( woundcarepain==0 | woundcarepain==1  ~"NoneMild",
                                                         woundcarepain==3~ "Severe", 
                                                         woundcarepain==2 ~ "Moderate" 
))
FullDataMicrobiome_BinaryPain = FullDataMicrobiome %>% filter(PainCatBinary!="Moderate")

FullDataMicrobiome_NonVAC_binPain= FullDataMicrobiome_BinaryPain %>% filter(BinaryDressingType=="NonWoundVac")
FullDataMicrobiome_VAC_binPain = FullDataMicrobiome_BinaryPain %>% filter(BinaryDressingType=="WoundVac")

Genera = c("CorynebacteriumAbundance_CLR",
           "StreptococcusAbundance_CLR" ,
           "StaphylococcusAbundance_CLR",
           "PseudomonasAbundance_CLR",
           "AnaerobicGenusAbundance_CLR" )
StratifiedPvalues_WoundVAC = c()

# Corynebacterium
Coryne_Pain_strat_VAC = sanon(CorynebacteriumAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(BinaryDressingType), data=FullDataMicrobiome_BinaryPain)
SummaryCorynePain = summary(Coryne_Pain_strat_VAC)
StratifiedPvalues_WoundVAC = append(StratifiedPvalues_WoundVAC, SummaryCorynePain$coefficients[4])

# Streptococcus 
Streptococcus_Pain_strat_VAC = sanon(StreptococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(BinaryDressingType), data=FullDataMicrobiome_BinaryPain)
SummaryStrepPain = summary(Streptococcus_Pain_strat_VAC)
StratifiedPvalues_WoundVAC= append(StratifiedPvalues_WoundVAC, SummaryStrepPain$coefficients[4])

# Staphylococcus
Staphylococcus_Pain_strat_VAC = sanon(StaphylococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(BinaryDressingType), data=FullDataMicrobiome_BinaryPain)
SummaryStaphPain = summary(Staphylococcus_Pain_strat_VAC)
StratifiedPvalues_WoundVAC = append(StratifiedPvalues_WoundVAC, SummaryStaphPain$coefficients[4])

# Pseudomonas
Pseudomonas_Pain_strat_VAC = sanon(PseudomonasAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(BinaryDressingType) , data=FullDataMicrobiome_BinaryPain)
SummaryPseudomonasPain = summary(Pseudomonas_Pain_strat_VAC)
StratifiedPvalues_WoundVAC = append(StratifiedPvalues_WoundVAC, SummaryPseudomonasPain$coefficients[4])

# anaerobes
Anaerobes_Pain_strat_VAC = sanon(AnaerobicGenusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(BinaryDressingType) , data=FullDataMicrobiome_BinaryPain)
SummaryAnaerobesPain = summary(Anaerobes_Pain_strat_VAC)
StratifiedPvalues_WoundVAC = append(StratifiedPvalues_WoundVAC, SummaryAnaerobesPain$coefficients[4])

# P values for Figure S5A
data.frame(genus=Genera, p=StratifiedPvalues_WoundVAC)

# Figure S5A
#############
FullDataMicrobiome_NonVAC_binPainMelted = FullDataMicrobiome_NonVAC_binPain %>% select(AbundanceCLR, PainCatBinary) %>% melt(id.vars=c("PainCatBinary"))
FullDataMicrobiome_NonVAC_binPainMelted$Genus = sapply(FullDataMicrobiome_NonVAC_binPainMelted$variable, function(x) str_remove(x, "Abundance_CLR"))
NonVAC_Genus_Pain = ggplot(FullDataMicrobiome_NonVAC_binPainMelted, aes(x=PainCatBinary, y=value, fill=Genus)) + geom_boxplot(alpha=.4,size=.25) + facet_grid(~ Genus)+
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2,height=0,size=.5) + theme_classic()  +
  ylim(-1.2, 12.5) + xlab("Pain") + ggtitle("Common Genus Abundance in Wounds With Severe vs. Mild/None Pain (Non-VAC Treated, n=191)") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") +
  ylab("CLR-transformed relative abundance")+  stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") 

FullDataMicrobiome_VAC_binPainMelted = FullDataMicrobiome_VAC_binPain %>% select(AbundanceCLR, PainCatBinary) %>% melt(id.vars=c("PainCatBinary"))
FullDataMicrobiome_VAC_binPainMelted$Genus = sapply(FullDataMicrobiome_VAC_binPainMelted$variable, function(x) str_remove(x, "Abundance_CLR"))
VAC_Genus_Pain = ggplot(FullDataMicrobiome_VAC_binPainMelted, aes(x=PainCatBinary, y=value, fill=Genus)) + geom_boxplot(alpha=.4,size=.25) + facet_grid(~ Genus)+
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2,height=0,size=.5) + theme_classic()  +
  ylim(-1.2, 12.5) + xlab("Pain") + ggtitle("Common Genus Abundance in Wounds With Severe vs. Mild/None Pain (VAC Treated, n=67)") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") +
  ylab("CLR-transformed relative abundance")+  stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .3, geom = "crossbar",color="gray") 
pdf("OutputFigures/S5A.pdf")
grid.arrange(VAC_Genus_Pain,NonVAC_Genus_Pain,ncol=1)
dev.off()