# Figure 3
# The heterogeneous wound microbiome varies with wound care pain, dressing practices, and inflammatory gene expression 
# Scripts to produce figure 3, Figure S3

# Wound microbiome and host factors 
library("dplyr")
library("stringr")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("viridis")
library("gplots")
library("rstatix")
library("sanon")
library("GGally")
#library("correlation")
#library("psych")
#library("sanon")
outputfolder="OutputFigures/"


Cytokines_And_Microbiome = read.csv("InputData/WoundMicrobiome_Cytokine_Data_Final.csv")
ClinicalData =  read.csv("InputData/ClinicalData_UpdatedPain.csv")


listCytokines = c("ARG1.Hs00163660_m1",  "C3.Hs00163811_m1",  "C5AR1.Hs00704891_s1", "CAMP.Hs00189038_m1",  "CXCL8.Hs00174103_m1", "IL1A.Hs00174092_m1",  "IL1B.Hs01555410_m1",
                  "IL6.Hs00174131_m1", "LCN2.Hs01008571_m1", "MMP1.Hs00899658_m1", "MMP2.Hs01548727_m1", "MMP9.Hs00957562_m1", "TNF.Hs00174128_m1")



###################################
# Figure 3A: wound genera vs. pain 
# using the stratified wilcoxon rank
# sum test
###################################
FullData = Cytokines_And_Microbiome %>% left_join(ClinicalData, by="study_id")
FullData = FullData %>% mutate(PainCatBinary = case_when( woundcarepain==0 | woundcarepain==1  ~"NoneMild",
                                               woundcarepain==3~ "Severe", 
                                               woundcarepain==2 ~ "Moderate" 
))


FullData_TestSevereVsMildNone = FullData  %>% filter(PainCatBinary != "Moderate")

# filter out those without microbiome data
FullData_TestSevereVsMildNone = FullData_TestSevereVsMildNone %>% filter(!is.na(DMMClusterAssign))
AbundanceCLR = colnames(FullData_TestSevereVsMildNone)[grepl(colnames(FullData_TestSevereVsMildNone), pattern="_CLR")]

testresultcoryne = (wilcox.test(CorynebacteriumAbundance_CLR ~ PainCatBinary, data=FullData_TestSevereVsMildNone))


Genera = AbundanceCLR
StratifiedPvalues = c()
StratifiedEstimates = c()

# Corynebacterium
Coryne_Pain_strat = sanon(CorynebacteriumAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullData_TestSevereVsMildNone)
SummaryCorynePain = summary(Coryne_Pain_strat)
pvalues_genus_stratified = append(StratifiedPvalues, SummaryCorynePain$coefficients[4])
StratifiedEstimates = append(StratifiedEstimates, Coryne_Pain_strat$xi)

# Streptococcus 
Streptococcus_Pain_Strat = sanon(StreptococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullData_TestSevereVsMildNone)
SummaryStrepPain = summary(Streptococcus_Pain_Strat)
pvalues_genus_stratified= append(pvalues_genus_stratified, SummaryStrepPain$coefficients[4])
StratifiedEstimates = append(StratifiedEstimates, Streptococcus_Pain_Strat$xi)


# Staphylococcus
Staphylococcus_Pain_strat = sanon(StaphylococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullData_TestSevereVsMildNone)
SummaryStaphPain = summary(Staphylococcus_Pain_strat)
pvalues_genus_stratified = append(pvalues_genus_stratified, SummaryStaphPain$coefficients[4])
StratifiedEstimates = append(StratifiedEstimates, Staphylococcus_Pain_strat$xi)

# Pseudomonas
Pseudomonas_Pain_strat = sanon(PseudomonasAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullData_TestSevereVsMildNone)
SummaryPseudomonasPain = summary(Pseudomonas_Pain_strat)
pvalues_genus_stratified = append(pvalues_genus_stratified, SummaryPseudomonasPain$coefficients[4])
StratifiedEstimates = append(StratifiedEstimates, Pseudomonas_Pain_strat$xi)

# anaerobes
Anaerobes_Pain_strat = sanon(AnaerobicGenusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) , data=FullData_TestSevereVsMildNone)
SummaryAnaerobesPain = summary(Anaerobes_Pain_strat)
pvalues_genus_stratified = append(pvalues_genus_stratified, SummaryAnaerobesPain$coefficients[4])
StratifiedEstimates = append(StratifiedEstimates, Anaerobes_Pain_strat$xi)

GenusAbundances = data.frame(Genus=AbundanceCLR,pvalue=pvalues_genus_stratified,estimates=StratifiedEstimates )

Genera = AbundanceCLR


GenusAbundances$padj = p.adjust(GenusAbundances$pvalue, method="BH")

GenusAbundances$group1 = "NoneMild"
GenusAbundances$group2 = "Severe"
GenusAbundances$.y.="value"
GenusAbundances$n1 = 146
GenusAbundances$n2 = 112
GenusAbundances$y.position = 12.3

GenusAbundances$xmin = 1 
GenusAbundances$xmax=2

FullData_TestSevereVsMildNone$PainCatBinary = factor(FullData_TestSevereVsMildNone$PainCatBinary)
DataMeltedGenusAbundance = FullData_TestSevereVsMildNone %>% select(AbundanceCLR, PainCatBinary) %>% melt(id.vars=c("PainCatBinary"))
DataMeltedGenusAbundance$Genus = sapply(DataMeltedGenusAbundance$variable, function(x) str_split(x,pattern="Abundance_CLR")[[1]][1])


 GenusPainPlot = ggplot(DataMeltedGenusAbundance, aes(x=PainCatBinary, y=value, fill=Genus)) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray34", alpha=.1)+ geom_boxplot(alpha=.4, size=.25) + facet_grid(~ Genus) +
    scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2) + theme_classic()  +
    ylim(-1.2, 12.5) + xlab("Pain Rating Category") + ggtitle("Common Genus Abundance in Wounds with \nSevere vs. None/Mild Pain Ratings") +
    theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") + ylab("CLR.transformed relative abundance") 
  
ggsave(GenusPainPlot, file=paste0(outputfolder, "Figure3a.pdf"))
ggsave(GenusPainPlot, file="OutputFigures/New_3a.pdf",width=12, height=6)

# Genera ~ wound duration + healing -- are Strep and Coryne still associated w/ healing when you stratify on wound duration(≤30 days vs. >30 days)?
####################################################################################################################################################
FullData_TestSevereVsMildNone_AgePain = FullData_TestSevereVsMildNone %>% mutate(WoundAgeBin = if_else(woundage %in% c(1,2), "Less30", "Gr30"))
Genera_AgePain = AbundanceCLR
StratifiedPvalues_Age_Pain = c()
StratifiedEstimates_AgePain = c()

# Corynebacterium
Coryne_Pain_strat_age = sanon(CorynebacteriumAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run ) + strt(WoundAgeBin) , data=FullData_TestSevereVsMildNone_AgePain)
SummaryCorynePain_age = summary(Coryne_Pain_strat_age)
pvalues_genus_stratified_age = append(StratifiedPvalues_Age_Pain, SummaryCorynePain_age$coefficients[4])
StratifiedEstimates_AgePain = append(StratifiedEstimates_AgePain, Coryne_Pain_strat_age$xi)

# Streptococcus 
Streptococcus_Pain_Strat_age = sanon(StreptococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(WoundAgeBin) , data=FullData_TestSevereVsMildNone_AgePain)
SummaryStrepPain_age = summary(Streptococcus_Pain_Strat_age)
pvalues_genus_stratified_age= append(pvalues_genus_stratified_age, SummaryStrepPain_age$coefficients[4])
StratifiedEstimates_AgePain = append(StratifiedEstimates_AgePain, Streptococcus_Pain_Strat_age$xi)

# Staphylococcus
Staphylococcus_Pain_strat_age = sanon(StaphylococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(WoundAgeBin) , data=FullData_TestSevereVsMildNone_AgePain)
SummaryStaphPain = summary(Staphylococcus_Pain_strat_age)
pvalues_genus_stratified_age = append(pvalues_genus_stratified_age, SummaryStaphPain$coefficients[4])
StratifiedEstimates_AgePain = append(StratifiedEstimates_AgePain, Staphylococcus_Pain_strat$xi)

# Pseudomonas
Pseudomonas_Pain_strat_age = sanon(PseudomonasAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(WoundAgeBin)  , data=FullData_TestSevereVsMildNone_AgePain)
SummaryPseudomonasPain = summary(Pseudomonas_Pain_strat_age)
pvalues_genus_stratified_age = append(pvalues_genus_stratified_age, SummaryPseudomonasPain$coefficients[4])
StratifiedEstimates_AgePain = append(StratifiedEstimates_AgePain, Pseudomonas_Pain_strat_age$xi)

# anaerobes
Anaerobes_Pain_strat_age = sanon(AnaerobicGenusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(WoundAgeBin) , data=FullData_TestSevereVsMildNone_AgePain)
SummaryAnaerobesPain = summary(Anaerobes_Pain_strat_age)
pvalues_genus_stratified_age = append(pvalues_genus_stratified_age, SummaryAnaerobesPain$coefficients[4])
StratifiedEstimates_AgePain = append(StratifiedEstimates_AgePain, Anaerobes_Pain_strat_age$xi)



GenusPainPlot = ggplot(DataMeltedGenusAbundance, aes(x=PainCatBinary, y=value, fill=Genus)) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray34", alpha=.1)+ geom_boxplot(alpha=.4, size=.25) + facet_grid(~ Genus) +
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2) + theme_classic()  +
  ylim(-1.2, 12.5) + xlab("Pain Rating Category") + ggtitle("Common Genus Abundance in Wounds with \nSevere vs. None/Mild Pain Ratings") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") + ylab("CLR.transformed relative abundance") 

MeltedByAge = FullData_TestSevereVsMildNone_AgePain %>% select(CorynebacteriumAbundance_CLR, StreptococcusAbundance_CLR, StaphylococcusAbundance_CLR, PseudomonasAbundance_CLR, AnaerobicGenusAbundance_CLR, WoundAgeBin,PainCatBinary ) %>% melt(id.vars=c("PainCatBinary","WoundAgeBin"))

MeltedByAge_Pre = MeltedByAge %>% filter(WoundAgeBin=="Less30")
MeltedByAge_Post = MeltedByAge %>% filter(WoundAgeBin=="Gr30")

Age_Pain_Pre =ggplot(MeltedByAge_Pre, aes(x=factor(PainCatBinary), y=value, fill=factor(variable))) + geom_boxplot(alpha=.4, size=.25) + geom_jitter(height=0, width=.2) + theme_classic() +scale_fill_brewer(palette ="Dark2" )+ facet_grid(.~variable) +ylim(-1.2, 12.5)

Age_Pain_Post = ggplot(MeltedByAge_Post, aes(x=factor(PainCatBinary), y=value, fill=factor(variable))) + geom_boxplot(alpha=.4, size=.25) + geom_jitter(height=0, width=.2) + theme_classic() +scale_fill_brewer(palette ="Dark2" )+ facet_grid(.~variable) +ylim(-1.2, 12.5)



MeltedByAge_Pre$Genus = sapply(MeltedByAge_Pre$variable, function(x) str_split(x, pattern="Abundance_CLR")[[1]][1])
GenusStatsPainPre = MeltedByAge_Pre %>% group_by(Genus) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary")
GenusStatsPainPre$xmin=1
GenusStatsPainPre$xmax=2

MeltedByAge_Post$Genus = sapply(MeltedByAge_Post$variable, function(x) str_split(x, pattern="Abundance_CLR")[[1]][1])
GenusStatsPainPost = MeltedByAge_Post %>% group_by(Genus) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary")
GenusStatsPainPost$xmin=1
GenusStatsPainPost$xmax=2


GenusPainPlot_Pre30 = ggplot(MeltedByAge_Pre, aes(x=PainCatBinary, y=value, fill=Genus)) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray34", alpha=.1)+ geom_boxplot(alpha=.4, size=.25) + facet_grid(~ Genus) +
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2) + theme_classic()  + stat_pvalue_manual(GenusStatsPainPre, label="p")+
  ylim(-1.2, 12.5)  + xlab("Pain Rating Category") + ggtitle("Common Genus Abundance in Wounds with \nSevere vs. None/Mild Pain Ratings (≤30 days old)") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") + ylab("CLR.transformed relative abundance") 


GenusPainPlot_Post30 = ggplot(MeltedByAge_Post, aes(x=PainCatBinary, y=value, fill=Genus)) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray34", alpha=.1)+ geom_boxplot(alpha=.4, size=.25) + facet_grid(~ Genus) +
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2) + theme_classic()  + stat_pvalue_manual(GenusStatsPainPost, label="p")+
  ylim(-1.2, 12.5)  + xlab("Pain Rating Category") + ggtitle("Common Genus Abundance in Wounds with \nSevere vs. None/Mild Pain Ratings (>30 days old)") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") + ylab("CLR.transformed relative abundance") 

pdf("OutputFigures/New_S3b.pdf", width=8, height=12)
grid.arrange(GenusPainPlot_Pre30, GenusPainPlot_Post30, ncol=1)
dev.off()

# Alternative to Figure S3B: Just make the shapes of the points based on the age category
MeltedByAge$Genus = sapply(MeltedByAge$variable, function(x) str_split(x, pattern="Abundance_CLR")[[1]][1])
GenusPainPlot_Age_Shape = ggplot(MeltedByAge, aes(x=PainCatBinary, y=value, fill=Genus)) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray34", alpha=.1)+ geom_boxplot(alpha=.4, size=.25) + facet_grid(~ Genus) +
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(aes(shape=WoundAgeBin),width=.2) + theme_classic()  + stat_pvalue_manual(GenusStatsPainPre, label="p")+
  ylim(-1.2, 12.5)  + xlab("Pain Rating Category") + ggtitle("Common Genus Abundance in Wounds with \nSevere vs. None/Mild Pain Ratings (≤30 days old)") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") + ylab("CLR.transformed relative abundance") 

GenusPainPlot_Age_Shape

# Stratified by wound location
################################
FullData_TestSevereVsMildNone_locPain = FullData_TestSevereVsMildNone# %>% mutate(woundloc = if_else(woundloc %in% c(1,2), "Less30", "Gr30"))
Genera_locPain = AbundanceCLR
StratifiedPvalues_loc_Pain = c()
StratifiedEstimates_locPain = c()

# Corynebacterium
Coryne_Pain_strat_loc = sanon(CorynebacteriumAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run ) + strt(woundloc) , data=FullData_TestSevereVsMildNone_locPain)
SummaryCorynePain_loc = summary(Coryne_Pain_strat_loc)
pvalues_genus_stratified_loc = append(StratifiedPvalues_loc_Pain, SummaryCorynePain_loc$coefficients[4])
StratifiedEstimates_locPain = append(StratifiedEstimates_locPain, Coryne_Pain_strat_loc$xi)

# Streptococcus 
Streptococcus_Pain_Strat_loc = sanon(StreptococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(woundloc) , data=FullData_TestSevereVsMildNone_locPain)
SummaryStrepPain_loc = summary(Streptococcus_Pain_Strat_loc)
pvalues_genus_stratified_loc= append(pvalues_genus_stratified_loc, SummaryStrepPain_loc$coefficients[4])
StratifiedEstimates_locPain = append(StratifiedEstimates_locPain, Streptococcus_Pain_Strat_loc$xi)

# Staphylococcus
Staphylococcus_Pain_strat_loc = sanon(StaphylococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(woundloc) , data=FullData_TestSevereVsMildNone_locPain)
SummaryStaphPain = summary(Staphylococcus_Pain_strat_loc)
pvalues_genus_stratified_loc = append(pvalues_genus_stratified_loc, SummaryStaphPain$coefficients[4])
StratifiedEstimates_locPain = append(StratifiedEstimates_locPain, Staphylococcus_Pain_strat$xi)

# Pseudomonas
Pseudomonas_Pain_strat_loc = sanon(PseudomonasAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(woundloc)  , data=FullData_TestSevereVsMildNone_locPain)
SummaryPseudomonasPain = summary(Pseudomonas_Pain_strat_loc)
pvalues_genus_stratified_loc = append(pvalues_genus_stratified_loc, SummaryPseudomonasPain$coefficients[4])
StratifiedEstimates_locPain = append(StratifiedEstimates_locPain, Pseudomonas_Pain_strat_loc$xi)

# anaerobes
Anaerobes_Pain_strat_loc = sanon(AnaerobicGenusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(woundloc) , data=FullData_TestSevereVsMildNone_locPain)
SummaryAnaerobesPain = summary(Anaerobes_Pain_strat_loc)
pvalues_genus_stratified_loc = append(pvalues_genus_stratified_loc, SummaryAnaerobesPain$coefficients[4])
StratifiedEstimates_locPain = append(StratifiedEstimates_locPain, Anaerobes_Pain_strat_loc$xi)



GenusPainPlot = ggplot(DataMeltedGenusAbundance, aes(x=PainCatBinary, y=value, fill=Genus)) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray34", alpha=.1)+ geom_boxplot(alpha=.4, size=.25) + facet_grid(~ Genus) +
  scale_fill_brewer(palette ="Dark2" ) + geom_jitter(width=.2) + theme_classic()  +
  ylim(-1.2, 12.5) + xlab("Pain Rating Category") + ggtitle("Common Genus Abundance in Wounds with \nSevere vs. None/Mild Pain Ratings") +
  theme(plot.title=element_text(hjust=.5, size=18, face="bold"), axis.text.x=element_text(size=11),strip.text.x=element_text(size=13), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.position="None") + ylab("CLR.transformed relative abundance") 



MeltedByloc = FullData_TestSevereVsMildNone_locPain %>% select(CorynebacteriumAbundance_CLR, StreptococcusAbundance_CLR, StaphylococcusAbundance_CLR, PseudomonasAbundance_CLR, AnaerobicGenusAbundance_CLR, woundloc,PainCatBinary ) %>% melt(id.vars=c("PainCatBinary","woundloc"))
MeltedByloc$Genus = sapply(MeltedByloc$variable, function(x) str_split(x, pattern="Abundance_CLR")[[1]][1])


MeltedByloc1 = MeltedByloc %>% filter(woundloc==1)
MeltedByloc2 = MeltedByloc %>% filter(woundloc==2)
MeltedByloc4 = MeltedByloc %>% filter(woundloc==4)

GenusStatsMeltedByloc1 = MeltedByloc1 %>% group_by(Genus) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary")
GenusStatsMeltedByloc2 = MeltedByloc2 %>% group_by(Genus) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary")
GenusStatsMeltedByloc4 = MeltedByloc4 %>% group_by(Genus) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary")


loc1=ggplot(MeltedByloc1, aes(x=factor(PainCatBinary), y=value, fill=Genus)) + geom_boxplot(alpha=.4, size=.25) + 
  geom_jitter(height=0, width=.2) + theme_classic() + facet_grid(.~Genus) +ylim(-1.2, 12.5) + scale_fill_brewer(palette="Dark2") +
  stat_pvalue_manual(GenusStatsMeltedByloc1,label="p")
loc2=ggplot(MeltedByloc2, aes(x=factor(PainCatBinary), y=value, fill=(Genus))) + geom_boxplot(alpha=.4, size=.25) + 
  geom_jitter(height=0, width=.2) + theme_classic() + facet_grid(.~Genus) +ylim(-1.2, 12.5)+ scale_fill_brewer(palette="Dark2")+
  stat_pvalue_manual(GenusStatsMeltedByloc2,label="p")
loc4=ggplot(MeltedByloc4, aes(x=factor(PainCatBinary), y=value, fill=factor(Genus))) + geom_boxplot(alpha=.4, size=.25) + 
  geom_jitter(height=0, width=.2) + theme_classic() + facet_grid(.~Genus) +ylim(-1.2, 12.5)+ scale_fill_brewer(palette="Dark2")+
  stat_pvalue_manual(GenusStatsMeltedByloc4,label="p")

pdf("OutputFigures/New_S3c.pdf", width=9, height=12)
grid.arrange(loc1, loc2, loc4,ncol=1)
dev.off()


# Steve suggests stratifying on both variables
###############################################
FullData_TestSevereVsMildNone_locAgePain = FullData_TestSevereVsMildNone %>% mutate(WoundAgeBin = if_else(woundage %in% c(1,2), "Less30", "Gr30"))


Genera_loc_age_Pain = AbundanceCLR
StratifiedPvalues_loc_age_Pain = c()
StratifiedEstimates_age_locPain = c()

# Corynebacterium
Coryne_Pain_strat_loc = sanon(CorynebacteriumAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run )+ strt(WoundAgeBin) + strt(woundloc) , data=FullData_TestSevereVsMildNone_locAgePain)
SummaryCorynePain_loc = summary(Coryne_Pain_strat_loc)
pvalues_genus_stratified_loc_age = append(StratifiedPvalues_loc_age_Pain, SummaryCorynePain_loc$coefficients[4])
StratifiedEstimates_age_locPain = append(StratifiedEstimates_age_locPain, Coryne_Pain_strat_loc$xi)

# Streptococcus 
Streptococcus_Pain_Strat_loc = sanon(StreptococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(WoundAgeBin) + strt(woundloc) , data=FullData_TestSevereVsMildNone_locAgePain)
SummaryStrepPain_loc = summary(Streptococcus_Pain_Strat_loc)
pvalues_genus_stratified_loc_age= append(pvalues_genus_stratified_loc_age, SummaryStrepPain_loc$coefficients[4])
StratifiedEstimates_age_locPain = append(StratifiedEstimates_age_locPain, Streptococcus_Pain_Strat_loc$xi)

# Staphylococcus
Staphylococcus_Pain_strat_loc = sanon(StaphylococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(WoundAgeBin) + strt(woundloc) , data=FullData_TestSevereVsMildNone_locAgePain)
SummaryStaphPain = summary(Staphylococcus_Pain_strat_loc)
pvalues_genus_stratified_loc_age = append(pvalues_genus_stratified_loc_age, SummaryStaphPain$coefficients[4])
StratifiedEstimates_age_locPain = append(StratifiedEstimates_age_locPain, Staphylococcus_Pain_strat$xi)

# Pseudomonas
Pseudomonas_Pain_strat_loc = sanon(PseudomonasAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(WoundAgeBin)+ strt(woundloc)  , data=FullData_TestSevereVsMildNone_locAgePain)
SummaryPseudomonasPain = summary(Pseudomonas_Pain_strat_loc)
pvalues_genus_stratified_loc_age = append(pvalues_genus_stratified_loc_age, SummaryPseudomonasPain$coefficients[4])
StratifiedEstimates_age_locPain = append(StratifiedEstimates_age_locPain, Pseudomonas_Pain_strat_loc$xi)

# anaerobes
Anaerobes_Pain_strat_loc = sanon(AnaerobicGenusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(WoundAgeBin) + strt(woundloc) , data=FullData_TestSevereVsMildNone_locAgePain)
SummaryAnaerobesPain = summary(Anaerobes_Pain_strat_loc)
pvalues_genus_stratified_loc_age = append(pvalues_genus_stratified_loc_age, SummaryAnaerobesPain$coefficients[4])
StratifiedEstimates_age_locPain = append(StratifiedEstimates_age_locPain, Anaerobes_Pain_strat_loc$xi)

# Coryne, Strep, Staph, Pseud, Anaerobes
pvalues_genus_stratified_loc_age






# Common genus abundance vs. wound dressing change pain (old method, before stratifying)
#########################################################################################
AbundanceCLR = colnames(FullData_TestSevereVsMildNone)[grepl(colnames(FullData_TestSevereVsMildNone), pattern="_CLR")]
GenusStats = DataMeltedGenusAbundance %>% group_by(Genus) %>% wilcox_test(value ~ PainCatBinary)  %>% adjust_pvalue(method = "BH") %>% add_significance()  %>%  add_xy_position(x="PainCatBinary")
GenusStats$y.position = GenusStats$y.position + 1
  
GenusPainPlot + ggpubr::stat_pvalue_manual(GenusStats, label="p")

###################################
# Figure 3B: wound genera vs. 
# inflammatory mediator genes
###################################
listCytokines = c("ARG1.Hs00163660_m1",  "C3.Hs00163811_m1",  "C5AR1.Hs00704891_s1", "CAMP.Hs00189038_m1",  "CXCL8.Hs00174103_m1", "IL1A.Hs00174092_m1",  "IL1B.Hs01555410_m1",
                  "IL6.Hs00174131_m1", "LCN2.Hs01008571_m1", "MMP1.Hs00899658_m1", "MMP2.Hs01548727_m1", "MMP9.Hs00957562_m1", "TNF.Hs00174128_m1")
listkruskalCytokinesDMM = c()
listkruskalPs = c()
for(cyt in listCytokines){
  testresult = kruskal.test(FullData[,cyt] ~ FullData[,"DMMClusterAssign"])  
  listkruskalPs = c(listkruskalPs, testresult$p.value)
  listkruskalCytokinesDMM = c(listkruskalCytokinesDMM, cyt)
  
}

DataFrameCytokineDMMs = data.frame(pvalsKruskalWallis =listkruskalPs, Cytokine=listkruskalCytokinesDMM )
Significant = DataFrameCytokineDMMs %>% filter(pvalsKruskalWallis < .05)
Cyts = FullData %>% select(Significant$Cytokine)

GenusVars = FullData %>% select(AbundanceCLR)
DFCompare = cbind(Cyts,GenusVars )

#  cytokine markers of interest vs. DMM
###########################################
FullDataContainsMicrobiome = FullData %>% filter(!is.na(DMMClusterAssign))

colorpalette6 = RColorBrewer::brewer.pal(6, "Dark2")


C5ARStats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `C5AR1.Hs00704891_s1`) %>% wilcox_test(`C5AR1.Hs00704891_s1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p<.05)
C5ARStats$signif = case_when(C5ARStats$p >=.05 ~ "NS",
                             C5ARStats$p <.05 & C5ARStats$p >= .01 ~ "*", 
                             C5ARStats$p <.01 & C5ARStats$p >= .001 ~ "**", 
                             C5ARStats$p <.001 & C5ARStats$p >= .0001 ~ "***", 
                             C5ARStats$p < .0001  ~ "****", 
                             )

cxcl8Stats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `CXCL8.Hs00174103_m1`) %>% wilcox_test(`CXCL8.Hs00174103_m1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p<.05)

cxcl8Stats$signif = case_when(cxcl8Stats$p >=.05 ~ "NS",
                              cxcl8Stats$p <.05 & cxcl8Stats$p >= .01 ~ "*", 
                             cxcl8Stats$p <.01 & cxcl8Stats$p >= .001 ~ "**", 
                             cxcl8Stats$p <.001 & cxcl8Stats$p >= .0001 ~ "***", 
                             cxcl8Stats$p < .0001  ~ "****", 
)

ilbStats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `IL1B.Hs01555410_m1`) %>% wilcox_test(`IL1B.Hs01555410_m1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p<.05)
ilbStats$signif = case_when(ilbStats$p >=.05 ~ "NS",
                            ilbStats$p <.05 & ilbStats$p >= .01 ~ "*", 
                            ilbStats$p <.01 & ilbStats$p >= .001 ~ "**", 
                            ilbStats$p <.001 & ilbStats$p >= .0001 ~ "***", 
                            ilbStats$p < .0001  ~ "****", 
)


MMP2Stats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `MMP2.Hs01548727_m1`) %>% wilcox_test(`MMP2.Hs01548727_m1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p<.05)
MMP2Stats$signif = case_when(MMP2Stats$p >=.05 ~ "NS",
                             MMP2Stats$p <.05 & MMP2Stats$p >= .01 ~ "*", 
                             MMP2Stats$p <.01 & MMP2Stats$p >= .001 ~ "**", 
                             MMP2Stats$p <.001 & MMP2Stats$p >= .0001 ~ "***", 
                             MMP2Stats$p < .0001  ~ "****", 
)



tnfStats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `TNF.Hs00174128_m1`) %>% wilcox_test(`TNF.Hs00174128_m1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p<.05)
tnfStats$signif = case_when(tnfStats$p >=.05 ~ "NS",
                            tnfStats$p <.05 & tnfStats$p >= .01 ~ "*", 
                            tnfStats$p <.01 & tnfStats$p >= .001 ~ "**", 
                            tnfStats$p <.001 & tnfStats$p >= .0001 ~ "***", 
                            tnfStats$p < .0001  ~ "****", 
)


campStats = FullDataContainsMicrobiome %>% select(DMMClusterAssign, `CAMP.Hs00189038_m1`) %>% wilcox_test(`CAMP.Hs00189038_m1` ~ DMMClusterAssign) %>% add_xy_position() %>% filter(p<.05)
campStats$signif = case_when(campStats$p >=.05 ~ "NS",
                             campStats$p <.05 & campStats$p >= .01 ~ "*", 
                             campStats$p <.01 & campStats$p >= .001 ~ "**", 
                             campStats$p <.001 & campStats$p >= .0001 ~ "***", 
                             campStats$p < .0001  ~ "****", 
)



c5ar1plot = ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`C5AR1.Hs00704891_s1` )),) + geom_boxplot(fill=colorpalette6[1]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab(".∆CT") + ggtitle("C5AR1 Expression") + theme(plot.title=element_text(hjust=.5, size=20, face="bold")) +
  ggpubr::stat_pvalue_manual(C5ARStats, label="signif",size=6) +stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

cxcl8plot  =  ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`CXCL8.Hs00174103_m1` ))) + geom_boxplot(fill=colorpalette6[2]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab(".∆CT")+ ggtitle("CXCL8 Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold")) + 
  ggpubr::stat_pvalue_manual(cxcl8Stats, label="signif",size=6)+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

ilbplot  =  ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`IL1B.Hs01555410_m1` ))) + geom_boxplot(fill=colorpalette6[3]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab(".∆CT")+ ggtitle("Il.1B Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold"))+
  ggpubr::stat_pvalue_manual(ilbStats, label="signif",size=6)+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

MMP2plot  =  ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`MMP2.Hs01548727_m1` ))) + geom_boxplot(fill=colorpalette6[4]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab(".∆CT")+ ggtitle("MMP2 Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold"))+
  ggpubr::stat_pvalue_manual(MMP2Stats, label="signif",size=6)+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

CAMPplot = ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`CAMP.Hs00189038_m1` ))) + geom_boxplot(fill=colorpalette6[5]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab(".∆CT")+ ggtitle("CAMP Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold"))+
  ggpubr::stat_pvalue_manual(campStats, label="signif",size=6)+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

TNFplot = ggplot(FullDataContainsMicrobiome, aes(x=factor(DMMClusterAssign), y=(`TNF.Hs00174128_m1` ))) + geom_boxplot(fill=colorpalette6[6]) + theme_classic() +
  xlab("Microbial Community Type (DMM Assignment)") + ylab(".∆CT")+ ggtitle("TNF.alpha Expression")+ theme(plot.title=element_text(hjust=.5, size=20, face="bold"))+
  ggpubr::stat_pvalue_manual(tnfStats, label="signif",size=6)+stat_summary(fun= "mean", fun.max= "mean", fun.min= "mean", size= .4, geom = "crossbar",color="gray68", alpha=.1)

pdf(file=paste0(outputfolder, "Figure3B.pdf"), height=7,width=12)
gridExtra::grid.arrange(c5ar1plot, cxcl8plot, ilbplot, MMP2plot, CAMPplot, TNFplot, ncol=3 )
dev.off()




###############
# Figure 3C
###############
GeneraData = FullDataContainsMicrobiome %>% select(AbundanceCLR)
CytokineData = FullDataContainsMicrobiome %>% select(Significant$Cytokine)
CorrelationsGeneraCyt = matrix(0, nrow=ncol(GeneraData), ncol=ncol(CytokineData))
CorrelationsGeneraCytPvalues = matrix(0, nrow=ncol(GeneraData), ncol=ncol(CytokineData))

Both=cbind(CytokineData, GeneraData)

pdf(file="OutputFigures/BigCorrFig.pdf", width=20,height=20)
ggpairs(Both) + theme_classic()
dev.off()

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
CorrelationsMeltedCytGenera$pvalueRounded = round(CorrelationsMeltedCytGenera$pvalue, 4)

CorrHeatMapGeneraCyt = ggplot(CorrelationsMeltedCytGenera, aes(x=variable, y=Genera, fill=value))+ geom_tile() + scale_fill_gradient2(high="red", low="blue", mid="white", midpoint=0) + xlab("") + ylab("") +
  ggtitle("Correlations Between Inflammatory Marker Expression (-∆CT)\n and Common Wound Genera ") +theme_classic()+ theme(plot.title=element_text(size=20, hjust=.5), axis.text.x=element_text(size=15, angle=270), axis.text.y=element_text(size=15)) +
  labs(fill="Spearman Correlation") + geom_text(aes(label=pvalueRounded), size=8) + coord_equal()

ggsave(CorrHeatMapGeneraCyt, file=paste0(outputfolder, "Figure3C.pdf"), width=14,height=14)


# As an aside, relationship between Strep & pain, Coryne & pain hold when stratified by woundlocation 
Staphylococcus_Pain_stratLocation = sanon(StaphylococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") +strt(Run) + strt(woundloc) , data=FullData_TestSevereVsMildNone)
Streptococcus_Pain_stratLocation = sanon(StreptococcusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(woundloc) , data=FullData_TestSevereVsMildNone)
Corynebacterium_Pain_stratLocation = sanon(CorynebacteriumAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run) + strt(woundloc) , data=FullData_TestSevereVsMildNone)
Anaerobes_Pain_stratLocation = sanon(AnaerobicGenusAbundance_CLR ~ grp(PainCatBinary, ref="NoneMild") + strt(Run)+ strt(woundloc) , data=FullData_TestSevereVsMildNone)


FullData$woundage = if_else(FullData$woundage %in% c(1,2), "Less30", "Gr30")

FullDataComplete = FullData %>% filter(!is.na(DMMClusterAssign))
FullDataComplete = FullDataComplete %>% filter(PainCatBinary!="Moderate")
Staphylococcus_Age_stratRun = sanon(StaphylococcusAbundance_CLR ~ grp(woundage, ref="Less30") + strt(Run) , data=FullDataComplete)
Streptococcus_Age_stratRun = sanon(StreptococcusAbundance_CLR ~ grp(woundage, ref="Less30") + strt(Run) , data=FullDataComplete)
Corynebacterium_Age_stratRun  = sanon(CorynebacteriumAbundance_CLR ~ grp(woundage, ref="Less30") + strt(Run) , data=FullDataComplete)
Anaerobes_Age_stratRun  = sanon(AnaerobicGenusAbundance_CLR ~ grp(woundage, ref="Less30") + strt(Run) , data=FullDataComplete)
Pseudomonas_Age_stratRun  = sanon(PseudomonasAbundance_CLR ~ grp(woundage, ref="Less30") + strt(Run) , data=FullDataComplete)



