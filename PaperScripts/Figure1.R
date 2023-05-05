# Amy Campbell
# To be submitted 2023
# The heterogeneous wound microbiome varies with wound care pain, dressing practices, and inflammatory gene expression 

#library("stringr")
#library("vegan")
#library("ggpubr")
#

#######################################
# Figure 1: Wound microbiome variables
#######################################
# In addition to doing the DMM Clustering and outputting the subsections of Figure 1, 
# This script also generates a file called WoundMicrobiomeDataForSEG_AEC.csv which 
# is used in downstream figures 



# Load packages
###############
library("dplyr")
library("phyloseq")
library("lattice")
library("parallel")
library("DirichletMultinomial")
library("reshape2")
library("ggplot2")
library("stringr")
library("microbiome")

# Functions
############
Genus_code = function(tax_table_row){
  if(tax_table_row["Genus"] !="NA"){
    return(tax_table_row["Genus"])
  }else if(tax_table_row["Family"] != "NA"){
    return(paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"))
  }else if(tax_table_row["Order"] != "NA"){
    return(paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"))
  }else if(tax_table_row["Class"]!="NA"){
    return(paste(tax_table_row["Class"],paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"), sep="_"))
  }else if(tax_table_row["Phylum"] !="NA"){
    return(paste(tax_table_row["Phylum"], paste(tax_table_row["Class"],paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"), sep="_"), sep="_"))
  }else if(tax_table_row["Phylum"]=="NA"){
    print("ouch")
    return("NULL")
  }else{
    return("NULL")
  }
}

# Load data
############
outputfolder="/Users/amycampbell/Documents/IowaWoundData2021/PublicationFiguresOutput/"
load("/Users/amycampbell/Documents/IowaWoundData2021/PublicationData/GenusLevelBatchCorrected.rda")
AnaerobeMappings = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/PublicationData/AnaerobeDB.csv")
CytokineData = read.csv("~/Documents/IowaWoundData2021/PublicationData/InflammatoryMediators.csv")

# this is important to use-- includes wound care pain classifications with the cutoffs used in the
# previous paper sfrom this dataset
ClinicalData=read.csv("Documents/IowaWoundData2021/PublicationData/ClinicalData_UpdatedPain.csv")


# Colors
#########
color21 =c("#7570B3", "#E7298A", "#66C2A5","#E31A1C",
           "#BC80BD", "#A6D854", "#FDCDAC", "#A6761D", "#80B1D3",
           "#FF9933", "#F4CAE4", "#FFFF66", "#2E8B57", "#CCEBC5", 
           "#1F78B4", "#BEBADA","#B3B3B3","#FF6600", "#800000","#23297A", "#FB8072")
OrdinationColors=c("#666666","#BEAED4","#E6AB02","#A6D854","#8DA0CB","#D95F02","#8B0A50")

randpalette18=c("#B300B3","#E6AB02",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00",
                "#681A1A","#7570B3",
                "#1F78B4","#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#FFFF33",
                "#33A02C","#E6F5C9")
# Set Random seed
##################
set.seed(19104)

# Final filtering steps on genus-level, batch-corrected 
#######################################################
BatchCorrectedPhyloseq = phyloseq::subset_taxa(BatchCorrectedPhyloseq, Order!="Chloroplast")
BatchCorrectedPhyloseq = phyloseq::subset_taxa(BatchCorrectedPhyloseq, Family!="Mitochondria")
BatchCorrectedPhyloseq = phyloseq::subset_samples(BatchCorrectedPhyloseq, study_id !="Mock")


BatchCorrectedPhyloseq@sam_data$TotalOTUsLeftBatchMitoChloro = colSums(BatchCorrectedPhyloseq@otu_table@.Data)
BatchCorrectedPhyloseq@sam_data$TotalOTUsLeftBatchMitoChloro = colSums(BatchCorrectedPhyloseq@otu_table@.Data)

# Filter to at least 1200 OTUs post-decontamination
####################################################

BatchCorrectedPhyloseq = phyloseq::subset_samples(BatchCorrectedPhyloseq, TotalOTUsLeftBatchMitoChloro>1200)
FinalSampleDist = data.frame(BatchCorrectedPhyloseq@sam_data)

BatchCorrectedPhyloseq@sam_data = sample_data(FinalSampleDist)
sample_names(BatchCorrectedPhyloseq) = (BatchCorrectedPhyloseq@sam_data$study_id)


# DMM Clustering 
################

# Aggregated to Genus level and above
# Test DMM on this
Counts = BatchCorrectedPhyloseq@otu_table

# Filter to species prevalent in at least 10 samples 
PrevalenceInfo = data.frame(Counts)
PrevalenceInfo$Prevalence = rowSums(PrevalenceInfo != 0)

# .02*(nrow(BatchCorrectedPhyloseq@sam_data))
PrevalenceInfoFilter = PrevalenceInfo %>% filter(Prevalence >=.02*(nrow(BatchCorrectedPhyloseq@sam_data) ))
TaxaKeepDMM = row.names(PrevalenceInfoFilter)

#
PhyloDMM = prune_taxa(TaxaKeepDMM, BatchCorrectedPhyloseq)

CountsDMM = t(PhyloDMM@otu_table@.Data)

CountsDMMinput = data.matrix(CountsDMM)

densityplot(log10(colSums(CountsDMMinput)), xlim=range(log10(colSums(CountsDMMinput))),xlab="Taxon representation (log 10 count)")
densityplot((colSums(CountsDMMinput)), xlim=range((colSums(CountsDMMinput))),xlab="Taxon representation (raw count)")

FitDirichlets = mclapply(1:12, dmn, count=CountsDMMinput, verbose=TRUE, seed=19104)
lplc <- sapply(FitDirichlets, laplace)


#############
# Figure 1Ci
#############
pdf(file=paste0(outputfolder,"Figure1Ci.pdf"), height=5,width=5)
plot(lplc, type = 'b', xlab = 'Dirichlet Components',ylab='Model Fit (Laplace Approximation)', main="Dirichlet Components by Laplace Model Fit (Genus Read Counts)") 
dev.off()

# get dirichlet fit object with maximum laplace-estimated likelihood (minimum negative)
bestFit <- FitDirichlets[[which.min(lplc)]]

GroupScores = data.frame(bestFit@group)

# Highest scoring component for each sample is that sample's assignment
GroupScores$assignment = sapply(1:nrow(GroupScores), function(x) which.max(GroupScores[x,]))
GroupScores$SampleID = row.names(GroupScores)
GroupScores$study_id=GroupScores$SampleID
AssignmentMapping = GroupScores %>% select(study_id,assignment )

SamDataForAssigning = data.frame(BatchCorrectedPhyloseq@sam_data)
SamDataForAssigning = SamDataForAssigning %>% left_join(AssignmentMapping, by="study_id")

BatchCorrectedPhyloseq@sam_data = sample_data(SamDataForAssigning)
rownames(BatchCorrectedPhyloseq@sam_data) = SamDataForAssigning$SampleID
colnames(BatchCorrectedPhyloseq@sam_data) = colnames(SamDataForAssigning)

###############
# Figure 1C(ii)
###############

# NMDS
########
BatchCorrectedPhyloseq@sam_data$assignment = factor(BatchCorrectedPhyloseq@sam_data$assignment)
sample_names(BatchCorrectedPhyloseq) <- BatchCorrectedPhyloseq@sam_data$SampleID
NMDSOrd <- ordinate(BatchCorrectedPhyloseq, "NMDS", "bray", weighted=T, trymax=200)

BatchCorrectedPhyloseq@sam_data$assignment = factor(BatchCorrectedPhyloseq@sam_data$assignment)
DMMordPlot = plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, color="assignment", title="Ordination of Genus-aggregated OTUs by DMM Assignment") + scale_color_manual(values=OrdinationColors)

DMMordPlot = DMMordPlot + theme_classic()

ggsave(DMMordPlot, file=paste0(outputfolder, "Figure1Cii.pdf"), width=8,height=8)



################
# Figure 1Ciii
################

sample_names(BatchCorrectedPhyloseq) = BatchCorrectedPhyloseq@sam_data$SampleID


BatchCorrectedPhyloseq@sam_data$assignment = factor(BatchCorrectedPhyloseq@sam_data$assignment)
BatchCorrectedPhyloseq@sam_data$woundage = factor(BatchCorrectedPhyloseq@sam_data$woundage)

DF_Anaerobes = data.frame(BatchCorrectedPhyloseq@tax_table@.Data)


DF_Anaerobes$X= row.names(DF_Anaerobes)

AnaerobeMappingsReduced = AnaerobeMappings %>% select(X, Anaerobe) 
DF_Anaerobes = DF_Anaerobes %>% left_join(AnaerobeMappingsReduced, by="X")
DF_Anaerobes$Genus = as.character(DF_Anaerobes$Genus)
savetaxanames = DF_Anaerobes$X

DF_Anaerobes$X=NULL
DF_Anaerobes$Species=NULL
DF_Anaerobes  = DF_Anaerobes %>% mutate(Genus = case_when(Anaerobe==1 ~ "Anaerobes", 
                                                          Anaerobe==0 ~ Genus))

DF_Anaerobes  = DF_Anaerobes %>% mutate(Genus = case_when(Anaerobe==1 ~ "Anaerobes", 
                                                          Anaerobe==0 ~ Genus))


DF_Anaerobes$Family = as.character(DF_Anaerobes$Family)
DF_Anaerobes  = DF_Anaerobes %>% mutate(Family = case_when(Anaerobe==1 ~ "Anaerobes", 
                                                           Anaerobe==0 ~ Family))

DF_Anaerobes$Order = as.character(DF_Anaerobes$Order)
DF_Anaerobes  = DF_Anaerobes %>% mutate(Order = case_when(Anaerobe==1 ~ "Anaerobes", 
                                                          Anaerobe==0 ~ Order))

DF_Anaerobes$Class = as.character(DF_Anaerobes$Class)
DF_Anaerobes  = DF_Anaerobes %>% mutate(Class = case_when(Anaerobe==1 ~ "Anaerobes", 
                                                          Anaerobe==0 ~ Class))

DF_Anaerobes$Phylum = as.character(DF_Anaerobes$Phylum)
DF_Anaerobes  = DF_Anaerobes %>% mutate(Phylum = case_when(Anaerobe==1 ~ "Anaerobes", 
                                                           Anaerobe==0 ~ Phylum))
row.names(DF_Anaerobes) = savetaxanames

AnaerobeOTUs = row.names(DF_Anaerobes %>% filter(Anaerobe == 1))
AllOTUs = data.frame(BatchCorrectedPhyloseq@tax_table)
AllOTUs$OTU = row.names(AllOTUs)

AnaerobeGenera = AllOTUs %>% filter(OTU %in% AnaerobeOTUs)

BatchCorrectedPhyloseqRanking = BatchCorrectedPhyloseq
rownames(DF_Anaerobes) = savetaxanames

BatchCorrectedPhyloseqRanking@tax_table = tax_table(DF_Anaerobes)
row.names(BatchCorrectedPhyloseqRanking@tax_table) =savetaxanames

colnames(BatchCorrectedPhyloseqRanking@tax_table) = colnames(DF_Anaerobes)

BatchCorrectedPhyloseqRanking = BatchCorrectedPhyloseqRanking %>% tax_glom(taxrank = "Genus")

PostGlomTax = data.frame(BatchCorrectedPhyloseqRanking@tax_table@.Data)
PostGlomTax$GenusAdjust = apply(PostGlomTax, 1, Genus_code)

PlotTopGenera = BatchCorrectedPhyloseqRanking
PlotTopGenera@tax_table = tax_table(PostGlomTax)
row.names(PlotTopGenera@tax_table) = row.names(PostGlomTax)
colnames(PlotTopGenera@tax_table) = colnames(PostGlomTax)


row.names(PlotTopGenera@sam_data) = BatchCorrectedPhyloseq@sam_data$SampleID

PlotTopGenera = PlotTopGenera %>% transform_sample_counts(function(x) {x/sum(x)})

PlotTopGenera1 = subset_samples(PlotTopGenera, assignment==1)
PlotTopGenera2 = subset_samples(PlotTopGenera, assignment==2)
PlotTopGenera3 = subset_samples(PlotTopGenera, assignment==3)
PlotTopGenera4 = subset_samples(PlotTopGenera, assignment==4)


topdf1 = data.frame(PlotTopGenera1@otu_table@.Data)
topdf1$rowmean = rowMeans(topdf1)
topGenera = row.names( topdf1 %>% arrange(-rowmean))[1:12]

topdf2 = data.frame(PlotTopGenera2@otu_table@.Data)
topdf2$rowmean = rowMeans(topdf2)
topGenera = append(topGenera, row.names( topdf2 %>% arrange(-rowmean))[1:12])

topdf3 = data.frame(PlotTopGenera3@otu_table@.Data)
topdf3$rowmean = rowMeans(topdf3)
topGenera = append(topGenera, row.names( topdf3 %>% arrange(-rowmean))[1:12])

topdf4 = data.frame(PlotTopGenera4@otu_table@.Data)
topdf4$rowmean = rowMeans(topdf4)
topGenera = append(topGenera, row.names( topdf4 %>% arrange(-rowmean))[1:12])


topGenera = unique(topGenera)

PlotTop12Genera = prune_taxa(topGenera, PlotTopGenera)


MeltedTop12 = PlotTop12Genera %>% psmelt() %>% select(GenusAdjust, Abundance, assignment, study_id )
MeltedTop12 = MeltedTop12 %>% group_by(GenusAdjust, assignment) %>% summarise(mean(Abundance))
MeltedTop12$abundance = MeltedTop12$`mean(Abundance)`
generapresent  = ggplot(MeltedTop12, aes(x=as.factor(assignment), y=abundance, fill=GenusAdjust)) + geom_bar(stat = "identity", position = "stack",  color = NA) + scale_fill_manual(values=(color21)) + theme_minimal() + ggtitle("Top 12 Genera in Each DMM Component (Average %)") + xlab("DMM Assignment") + ylab("% Abundance")

ggsave(generapresent, file=paste0(outputfolder, "Figure1Ciii.pdf"), height=8,width=8)


################
# Figure 1C(iv)
###############
TaxRef = data.frame(BatchCorrectedPhyloseq@tax_table@.Data)
TaxRef$GenusAdjusted = apply(TaxRef, 1, Genus_code)

TaxRef$OTU = row.names(TaxRef)

Genera = TaxRef %>% select(OTU, GenusAdjusted)
FittedVals = data.frame(bestFit@fit$Estimate)

Nclusters = ncol(FittedVals)
colnames(FittedVals) =c(1:ncol(FittedVals))
FittedVals$OTU = row.names(FittedVals)
MeltedFitted = melt(FittedVals)
colnames(MeltedFitted) = c("OTU", "Cluster", "Value")
MeltedFitted = MeltedFitted %>% left_join(Genera, by="OTU")

# # modified from https://microbiome.github.io/tutorials/DMM.html
# for (k in 1:Nclusters) {
#   dfmelt = MeltedFitted %>% filter(Cluster==k) %>% arrange(-abs(Value))
#   topOTUs = (dfmelt$OTU)[1:10]
#   print(paste0("Cluster:", k))
#   dfmelt = dfmelt %>% filter(OTU %in% topOTUs)
#   plot = ggplot(dfmelt, aes(x=Genus,y=Value)) + geom_bar(stat="identity") + coord_flip() + labs(title=paste("Top Contributors to Cluster ", k ))
#   plot$data$Genus = factor(plot$data$Genus, rev(dfmelt$Genus))
#   
#   ggsave(plot, file=paste0(outputfolder, "Figure1Ciii_", k, ".pdf") )
# }
# 


GeneraToColorDF=data.frame(GenusForColor=sort(unique(MeltedTop12$GenusAdjust)), ColorHex = color21[1:20])

Cluster1 = MeltedFitted %>% filter(Cluster==1) %>% arrange(-abs(Value))
top5OTUs1 = Cluster1 %>% filter(OTU %in%  (Cluster1$OTU)[1:5])
top5OTUs1$Genus = top5OTUs1$GenusAdjusted
top5OTUs1 = top5OTUs1 %>% mutate(GenusForColor= if_else(Genus %in% AnaerobeGenera$Genus, "Anaerobes", Genus))
top5OTUs1 = top5OTUs1 %>%left_join(GeneraToColorDF,by="GenusForColor")

palette1 = (top5OTUs1 %>% arrange(Genus))$ColorHex
plot_top5cluster1 = ggplot(top5OTUs1, aes(x=Genus,y=Value,fill=GenusForColor)) + geom_bar(stat="identity") + coord_flip() + labs(title="Top Contributors to Cluster 1" ) + scale_fill_manual(values=unique(palette1))
plot_top5cluster1$data$Genus = factor(plot_top5cluster1$data$Genus, rev(top5OTUs1$Genus)) 
plot_top5cluster1 + theme_classic() + theme(legend.position = "none")

Cluster2 = MeltedFitted %>% filter(Cluster==2) %>% arrange(-abs(Value))
top5OTUs2 = Cluster2 %>% filter(OTU %in%  (Cluster2$OTU)[1:5])
top5OTUs2$Genus = top5OTUs2$GenusAdjusted
top5OTUs2 = top5OTUs2 %>% mutate(GenusForColor= if_else(Genus %in% AnaerobeGenera$Genus, "Anaerobes", Genus))


top5OTUs2 = top5OTUs2 %>%left_join(GeneraToColorDF,by="GenusForColor")

palette2 = (top5OTUs2 %>% arrange(GenusForColor))$ColorHex
plot_top5cluster2 = ggplot(top5OTUs2, aes(x=Genus,y=Value,fill=GenusForColor)) + geom_bar(stat="identity") + coord_flip() + labs(title=paste("Top Contributors to Cluster 2" )) + scale_fill_manual(values=unique(palette2))
plot_top5cluster2$data$Genus = factor(plot_top5cluster2$data$Genus, rev(top5OTUs2$Genus))

plot_top5cluster2 + theme_classic() + theme(legend.position = "none")


Cluster3 = MeltedFitted %>% filter(Cluster==3) %>% arrange(-abs(Value))
top5OTUs3 = Cluster3 %>% filter(OTU %in%  (Cluster3$OTU)[1:5])
top5OTUs3$Genus = top5OTUs3$GenusAdjusted
top5OTUs3 = top5OTUs3 %>% mutate(GenusForColor= if_else(Genus %in% AnaerobeGenera$Genus, "Anaerobes", Genus))

top5OTUs3 = top5OTUs3 %>%left_join(GeneraToColorDF,by="GenusForColor")

palette3 = (top5OTUs3 %>% arrange(GenusForColor))$ColorHex
plot_top5cluster3 = ggplot(top5OTUs3, aes(x=Genus,y=Value,fill=GenusForColor)) + geom_bar(stat="identity") + coord_flip() + labs(title=paste("Top Contributors to Cluster 3" )) + scale_fill_manual(values=unique(palette3))
plot_top5cluster3$data$Genus = factor(plot_top5cluster3$data$Genus, rev(top5OTUs3$Genus))

plot_top5cluster3 + theme_classic() + theme(legend.position = "none")


Cluster4 = MeltedFitted %>% filter(Cluster==4) %>% arrange(-abs(Value))
top5OTUs4 = Cluster4 %>% filter(OTU %in%  (Cluster4$OTU)[1:5])
top5OTUs4$Genus = top5OTUs4$GenusAdjusted
top5OTUs4 = top5OTUs4 %>% mutate(GenusForColor= if_else(Genus %in% AnaerobeGenera$Genus, "Anaerobes", Genus))

top5OTUs4 = top5OTUs4 %>%left_join(GeneraToColorDF,by="GenusForColor")

palette4 = (top5OTUs4 %>% arrange(GenusForColor))$ColorHex
plot_top5cluster4 = ggplot(top5OTUs4, aes(x=Genus,y=Value,fill=GenusForColor)) + geom_bar(stat="identity") + coord_flip() + labs(title=paste("Top Contributors to Cluster 4" )) + scale_fill_manual(values=unique(palette4))
plot_top5cluster4$data$Genus = factor(plot_top5cluster4$data$Genus, rev(top5OTUs4$Genus))

plot_top5cluster4 + theme_classic() + theme(legend.position = "none")

pdf(file=paste0(outputfolder, "Figure1Civ.pdf"), width=7.2, height=4.3)
gridExtra::grid.arrange(plot_top5cluster1+ theme_classic() + theme(legend.position = "none"), plot_top5cluster2+ theme_classic() + theme(legend.position = "none"), plot_top5cluster3+theme_classic() + theme(legend.position = "none"), plot_top5cluster4+theme_classic() + theme(legend.position = "none"))
dev.off()

############
# Figure 1b
############
# Phylum Level 

# Order by type followed by firmicutes abundance
samplenames = sample_names(BatchCorrectedPhyloseq)
DFsampledata =data.frame(BatchCorrectedPhyloseq@sam_data) 

DFsampledata = DFsampledata %>% mutate(wound_type_string = case_when(wound_type==1 ~ "Pressure", 
                                                                           wound_type==2 ~"Venous",
                                                                           wound_type==4 ~"Surgical", 
                                                                           wound_type==5 ~ "Traumatic", 
                                                                           wound_type==6 ~ "Mixed Traumatic/Surgical",
                                                                           wound_type==7 ~ "Other"))

BatchCorrectedPhyloseq@sam_data = sample_data(DFsampledata)

df_phyla= BatchCorrectedPhyloseq %>% 
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance >.01) %>%
  group_by(Phylum)


dataframe_counts = data.frame(df_phyla)
sampnames = (unique(dataframe_counts$Sample))


df_phyla_Firm = df_phyla %>% filter(Phylum=="Firmicutes") %>% select(wound_type_string, SampleID, Abundance)
MissingFirm = unique(setdiff(df_phyla$SampleID, df_phyla_Firm$SampleID))
MissingFirm = df_phyla %>% filter(SampleID %in% MissingFirm) %>% ungroup() %>% select(wound_type_string, SampleID) %>% unique()
MissingFirm$Phylum ="Firmicutes"
MissingFirm$Abundance = 0.0
MissingFirm = MissingFirm %>% select(Phylum, wound_type_string, SampleID, Abundance)
df_phyla_Firm = rbind(MissingFirm, df_phyla_Firm)



Desired_order=unique((df_phyla_Firm %>% arrange(wound_type_string, -Abundance))$SampleID)
plotPhyla= ggplot(df_phyla, aes(x=SampleID, y=Abundance, fill=Phylum)) + geom_bar(stat="identity") + ggtitle("Phylum-level composition") + scale_fill_manual(values=(randpalette18))+ theme_minimal() + theme(axis.text.x=element_text(angle=90))
plotPhyla$data$SampleID = factor(plotPhyla$data$SampleID, levels =Desired_order)
ggsave(plotPhyla, file=paste0(outputfolder, "Figure1B.pdf"), width=20, height=6)




# # Final dataframe for SEG
##########################
FinalDF = SamDataForAssigning %>% select(SampleID, study_id, assignment)
# here
colnames(FinalDF) = c("SampleID", "StudyID", "DMMClusterAssign")

DF_Genera_Dominating5 = FinalDF

# Highest phylum per sample & its abundance
######################################################
BatchCorrectedPhyloseqPhylum = BatchCorrectedPhyloseq %>% tax_glom(taxrank="Phylum")
PhylumTaxDF = data.frame(BatchCorrectedPhyloseqPhylum@tax_table@.Data)
PhylumTaxDF$HighestPhylum = row.names(PhylumTaxDF) 
PhylumTaxDF = PhylumTaxDF %>% select(Phylum, HighestPhylum)

DF_forMax_Phylum = data.frame(t(BatchCorrectedPhyloseqPhylum@otu_table@.Data))

RowSumsPhylum = rowSums(DF_forMax_Phylum)
MaxCountPhylum = apply(DF_forMax_Phylum,1,max)

DF_forMax_Phylum$HighestPhylum = colnames(DF_forMax_Phylum)[apply(DF_forMax_Phylum,1,which.max)]
DF_forMax_Phylum$HighestPhylum = sapply(DF_forMax_Phylum$HighestPhylum, function(x) str_remove(x, "X"))

DF_forMax_Phylum$MaxPhylumAbundance = MaxCountPhylum/RowSumsPhylum

DF_forMax_Phylum$SampleID = row.names(DF_forMax_Phylum)
DF_forMax_Phylum = DF_forMax_Phylum %>% left_join(PhylumTaxDF, by="HighestPhylum")
DF_forMax_Phylum = DF_forMax_Phylum %>% select(SampleID, Phylum, MaxPhylumAbundance)

FinalDF = FinalDF %>% left_join(DF_forMax_Phylum, by="SampleID")
colnames(FinalDF) = c("SampleID", "StudyID", "DMMClusterAssign", "Most_Abundant_Phylum","Phylum_Abundance")


# Highest genus-level(or next highest) grouping per sample & its abundance
###########################################################################

DF_ForMax_Genus = data.frame(t(BatchCorrectedPhyloseq@otu_table@.Data))
BackupDF_ForMax = DF_ForMax_Genus
GenusTaxDF = data.frame(BatchCorrectedPhyloseq@tax_table@.Data)
GenusTaxDF$GenusAdjust = apply(GenusTaxDF, 1, Genus_code)

GenusTaxDF$HighestGenus = row.names(GenusTaxDF) 
GenusTaxDF = GenusTaxDF %>% select(GenusAdjust, HighestGenus)

RowSumsGenus = rowSums(DF_ForMax_Genus)
MaxCountGenus = apply(DF_ForMax_Genus,1,max)


DF_ForMax_Genus$HighestGenus = colnames(DF_ForMax_Genus)[apply(DF_ForMax_Genus,1,which.max)]
DF_ForMax_Genus$HighestGenus = sapply(DF_ForMax_Genus$HighestGenus, function(x) str_remove(x, "X"))
DF_ForMax_Genus$MaxGenusAbundance =  MaxCountGenus/RowSumsGenus

DF_ForMax_Genus$SampleID = row.names(DF_ForMax_Genus)
DF_ForMax_Genus = DF_ForMax_Genus %>% left_join(GenusTaxDF,by="HighestGenus")
DF_ForMax_Genus = DF_ForMax_Genus %>% select(SampleID, GenusAdjust, MaxGenusAbundance)
colnames(DF_ForMax_Genus) = c("SampleID", "Most_Abundant_Genus","Genus_Abundance")

FinalDF = FinalDF %>% left_join(DF_ForMax_Genus, by="SampleID")

###############################################################################################
# Take and report the genera which are "most abundant" genus in at least 5 (>1% of 406) samples
##############################################################################################
MostAbundantGenusTable = rev(sort(table(FinalDF$Most_Abundant_Genus)))
FiveorMore = names(MostAbundantGenusTable[MostAbundantGenusTable>=5])

FiveOrMorePhyloseq = BatchCorrectedPhyloseq %>%transform_sample_counts(function(x) {x/sum(x)})
FiveOrMorePhyloseqCLR = microbiome::transform(FiveOrMorePhyloseq,transform="clr" )
DF_For_FiveOrMoreGeneraCLR = data.frame((FiveOrMorePhyloseqCLR@otu_table@.Data))
TopGeneraOTUnames = data.frame(BatchCorrectedPhyloseq@tax_table@.Data) %>% filter(Genus %in% FiveorMore)


FiveOrMoreGenera = data.frame(Genus = TopGeneraOTUnames$Genus)
FiveOrMoreGenera$otu = rownames(TopGeneraOTUnames)

DF_For_FiveOrMoreGeneraCLR$otu = row.names(DF_For_FiveOrMoreGeneraCLR)
DF_For_FiveOrMoreGeneraCLR = DF_For_FiveOrMoreGeneraCLR %>% filter(otu %in% FiveOrMoreGenera$otu)
FiveOrMoreGenera = FiveOrMoreGenera %>% left_join(DF_For_FiveOrMoreGeneraCLR, by="otu")
FiveOrMoreGenera = FiveOrMoreGenera %>% select(-otu)
row.names(FiveOrMoreGenera) = FiveOrMoreGenera$Genus
FiveOrMoreGenera = FiveOrMoreGenera %>% select(-Genus)
FiveOrMoreGeneraDF = data.frame(t(FiveOrMoreGenera))
FiveOrMoreGeneraDF$SampleID = row.names(FiveOrMoreGeneraDF)
FiveOrMoreGeneraDF$SampleID = stringr::str_remove(FiveOrMoreGeneraDF$SampleID, "X")

DF_Genera_Dominating5 = FiveOrMoreGeneraDF %>% left_join(DF_Genera_Dominating5, by="SampleID")
write.csv(DF_Genera_Dominating5, file="/Users/amycampbell/Documents/IowaWoundData2021/WoundAbundanceData_CLR_DominantGenera.csv")


BatchCorrectedPhyloseqRankingPercent = BatchCorrectedPhyloseqRanking %>%transform_sample_counts(function(x) {x/sum(x)})
CLRTransformedPct = microbiome::transform(BatchCorrectedPhyloseqRankingPercent,transform="clr" )

MeltedPct = BatchCorrectedPhyloseqRankingPercent %>% psmelt()
MeltedCLR = CLRTransformedPct %>% psmelt()

# Corynebacterium 
CorynePcts = MeltedPct %>% filter(Genus=="Corynebacterium") %>% select(SampleID, Abundance)
colnames(CorynePcts) = c("SampleID", "CorynebacteriumAbundance")
CoryneCLRs = MeltedCLR %>% filter(Genus=="Corynebacterium") %>% select(SampleID, Abundance)
colnames(CoryneCLRs) = c("SampleID", "CorynebacteriumAbundance_CLR")

FinalDF = FinalDF %>% left_join(CorynePcts, by="SampleID")
FinalDF = FinalDF %>% left_join(CoryneCLRs, by="SampleID")

# Streptococcus 
StrepPcts = MeltedPct %>% filter(Genus=="Streptococcus") %>% select(SampleID, Abundance)
colnames(StrepPcts) = c("SampleID", "StreptococcusAbundance")
StrepCLRs = MeltedCLR %>% filter(Genus=="Streptococcus") %>% select(SampleID, Abundance)
colnames(StrepCLRs) = c("SampleID", "StreptococcusAbundance_CLR")

FinalDF = FinalDF %>% left_join(StrepPcts, by="SampleID")
FinalDF = FinalDF %>% left_join(StrepCLRs, by="SampleID")


# Staphylococcus  
StaphPcts = MeltedPct %>% filter(Genus=="Staphylococcus") %>% select(SampleID, Abundance)
colnames(StaphPcts) = c("SampleID", "StaphylococcusAbundance")
StaphCLRs = MeltedCLR %>% filter(Genus=="Staphylococcus") %>% select(SampleID, Abundance)
colnames(StaphCLRs) = c("SampleID", "StaphylococcusAbundance_CLR")

FinalDF = FinalDF %>% left_join(StaphPcts, by="SampleID")
FinalDF = FinalDF %>% left_join(StaphCLRs, by="SampleID")


# Pseudomonas   
PseudPcts = MeltedPct %>% filter(Genus=="Pseudomonas") %>% select(SampleID, Abundance)
colnames(PseudPcts) = c("SampleID", "PseudomonasAbundance")
PseudCLRs = MeltedCLR %>% filter(Genus=="Pseudomonas") %>% select(SampleID, Abundance)
colnames(PseudCLRs) = c("SampleID", "PseudomonasAbundance_CLR")

FinalDF = FinalDF %>% left_join(PseudPcts, by="SampleID")
FinalDF = FinalDF %>% left_join(PseudCLRs, by="SampleID")


# Anaerobes   
AnaerobePcts = MeltedPct %>% filter(Genus=="Anaerobes") %>% select(SampleID, Abundance)
colnames(AnaerobePcts) = c("SampleID", "AnaerobicGenusAbundance")
AnaerobeCLRs = MeltedCLR %>% filter(Genus=="Anaerobes") %>% select(SampleID, Abundance)
colnames(AnaerobeCLRs) = c("SampleID", "AnaerobicGenusAbundance_CLR")

FinalDF = FinalDF %>% left_join(AnaerobePcts, by="SampleID")
FinalDF = FinalDF %>% left_join(AnaerobeCLRs, by="SampleID")




# Add cytokine information
###########################
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

FinalDF$study_id = FinalDF$StudyID
FullBioData = ClinicalDataAllIDs %>% left_join(FinalDF, by="study_id")


write.csv(FullBioData, file="/Users/amycampbell/Documents/IowaWoundData2021/PublicationData/WoundMicrobiome_Cytokine_Data_Final.csv")


############
# Figure 1A
############

FullDataDataPresent = FullBioData %>% select(study_id, DMMClusterAssign,names(CytokineDataNonNAGr20) )

FullDataDataPresent = FullDataDataPresent %>% left_join(ClinicalData %>% select(study_id, UpdatedPain, resting_pain_cat), by="study_id")


FullDataDataPresent[2:15][!is.na(FullDataDataPresent[2:15])] <- 1
FullDataDataPresent[2:15][is.na(FullDataDataPresent[2:15])] <- 0

colnames(FullDataDataPresent) = sapply(colnames(FullDataDataPresent), function(x) str_split(x, pattern="-Hs")[[1]][1] )
FullDataDataPresentMelt = FullDataDataPresent %>% reshape2::melt(id.vars=c("study_id", "UpdatedPain", "resting_pain_cat"))
FullDataDataPresentMelt$variable = (if_else(FullDataDataPresentMelt$variable=="DMMClusterAssign", "Microbiome", as.character(FullDataDataPresentMelt$variable)))


StudyIDOrder = unique((FullDataDataPresentMelt %>% arrange(UpdatedPain))$study_id)

dim(FullDataDataPresentMelt %>% arrange(UpdatedPain) %>% select(UpdatedPain, study_id) %>% unique())


DataAvailablePlot = ggplot(FullDataDataPresentMelt, aes(x=study_id, y=variable, fill=factor(value))) + geom_tile() + scale_fill_manual(values=c("white", "black")) + theme_classic() + theme(axis.text.x=element_text(angle=80, size=4), legend.position = "None")
DataAvailablePlot$data$study_id = factor(DataAvailablePlot$data$study_id, levels=StudyIDOrder)
DataAvailablePlot$data$variable = factor(DataAvailablePlot$data$variable , levels= c("Microbiome", setdiff(DataAvailablePlot$data$variable, "Microbiome")))

colorDF = FullDataDataPresentMelt %>% arrange(UpdatedPain) %>% select(UpdatedPain, study_id, resting_pain_cat) %>% unique()
colorDF = colorDF %>% mutate(ColorAssign = case_when(UpdatedPain==0 ~ "khaki1",
                                                     UpdatedPain==1 ~ "gold", 
                                                     UpdatedPain==2 ~ "darkorange1",
                                                     UpdatedPain==3~  "red3"))
colorDF = colorDF %>% mutate(ColorAssignResting = case_when(resting_pain_cat=="None" ~ "khaki1",
                                                            resting_pain_cat=="Mild" ~ "gold", 
                                                            resting_pain_cat=="Moderate" ~ "darkorange1",
                                                            resting_pain_cat=="Severe" ~  "red3"))



DataAvailablePlotResting = DataAvailablePlot + theme(axis.ticks.x = element_line(size=2, color=colorDF$ColorAssignResting))

ggsave(DataAvailablePlotResting + theme(axis.text.x=element_blank()) + xlab("Patient") + ylab("Biological Variable"), file=paste0(outputfolder,"DataPresenceRestingColorsForUse_fig1A.pdf"), width=20, height=4)
ggsave(DataAvailablePlot + theme(axis.text.x=element_blank()) + xlab("Patient") + ylab("Biological Variable"), file=paste0(outputfolder, "Figure1A.pdf"), width=20, height=4)
ggsave(DataAvailablePlot + xlab("Patient") + ylab("Biological Variable"), file=paste0(outputfolder,"DataPresencePatientIDs.pdf"), width=20, height=4)




###########
# Figure S2
###########
# Plot % abundance of the top 12 mean-abundance genera in any of the 4 DMMs:
BatchCorrectedPhyloseq_Anaerobes = BatchCorrectedPhyloseq

df_genera= BatchCorrectedPhyloseqRanking %>% 
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() 
Top12_anyDMMCluster = df_genera %>% filter(OTU %in% topGenera)


PlotTop12Genera_AnyDMM_bySample = PlotTop12Genera %>% psmelt() %>% select(GenusAdjust, Abundance, assignment, study_id)

Top12GeneraPlotSample = ggplot(PlotTop12Genera_AnyDMM_bySample, aes(x=study_id, y=Abundance, fill=GenusAdjust)) + geom_bar(stat="identity") + ggtitle("Genus-level composition") + scale_fill_manual(values=(color21))+ theme_minimal() + theme(axis.text.x=element_text(angle=90))
OrderSamples = (FinalDF %>% arrange(DMMClusterAssign, -StaphylococcusAbundance, AnaerobicGenusAbundance))$StudyID
Top12GeneraPlotSample$data$study_id = factor(Top12GeneraPlotSample$data$study_id , levels=OrderSamples)

ClusterDF=(FinalDF %>% arrange(DMMClusterAssign, -StaphylococcusAbundance, AnaerobicGenusAbundance))
ClusterDF = ClusterDF %>% mutate(colortext=case_when(DMMClusterAssign==1 ~ "#4285F4", 
                                                     DMMClusterAssign==2 ~"#34A853",
                                                     DMMClusterAssign==3 ~ "#FBBC05", 
                                                     DMMClusterAssign==4 ~ "#EA4335"))
colors = (ClusterDF %>% select(StudyID, colortext) %>% unique())$colortext
Top12GeneraPlotSample = Top12GeneraPlotSample + theme(axis.text.x=element_text(color=colors))
ggsave(Top12GeneraPlotSample, file=paste0(outputfolder,"FigureS2.pdf"), width=20, height=6)





