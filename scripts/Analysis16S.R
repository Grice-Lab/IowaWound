# Amy Campbell
# Jan 2022
# Analysis using decontaminated data object output by Qiime2OutputNov.R 

# Amy Campbell
# Iowa Wound microbiome analysis 2022
library("phyloseq")
library("dplyr")
library("stringr")
library("ggplot2")
library("vegan")
library("ggpubr")
library("lattice")
library("DirichletMultinomial")
library("parallel")
library("reshape2")
library("microbiome")
# Set Random seed
##################
set.seed(19104)

# Set color palettes
##############

# 32 kind of distinguishable colors 
randpalette32= c("#E6F5C9", "#BF5B17", "#33A02C", "#B3E2CD", "#66C2A5",
                 "#66A61E", "#FFD92F", "#FFFF33", "#666666", "#6A3D9A", "#F2F2F2",
                 "#A6CEE3", "#FC8D62", "#1F78B4", "#7570B3", "#B3B3B3", "#E41A1C",
                 "#FDC086", "#FCCDE5", "#FFFFB3","#E6AB02", "#FF7F00","#E7298A",
                 "#B3DE69","#D95F02","#1B9E77", "#BC80BD", "#A6761D",
                 "#006400", "#0000B3","#681A1A", "#B300B3")

ampalette <- rev(c("#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                   "#4D4D4D", "#F0E442", "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                   "#9900FF","#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                   "#F0E442", "#4D4D4D",  "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                   "#9900FF"))
palette8colors=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


randocolor18 = c("#7570B3", "#E7298A", "#66C2A5", "#FB8072", "#BC80BD",
                 "#A6D854", "#E31A1C", "#FDCDAC", "#A6761D", "#80B1D3",
                 "#FF9933", "#F4CAE4", "#FFFF66", "#CCEBC5", "#1F78B4",
                 "#BEBADA", "#B3B3B3", "#FF6600")
#pie(rep(1,18), col=rando18_1)


rando18_1 = c("#BEAED4","#386CB0",
              "#FF7F00", "#33A02C",
              "#A6CEE3", "#A65628",
              "#66C2A5", "#666666",
              "#F4CAE4","#8B0A50",
              "#D95F02", "#FFFF99",
              "#A6761D", "#E6AB02",
              "#8DA0CB", "#6A3D9A",
              "#A6D854","#FBB4AE")

Paincolors=c("#FFFFCC", "#FFB266", "#CC6600", "#990000", "#E0E0E0")

color21 =c("#7570B3", "#E7298A", "#66C2A5","#E31A1C",
           "#BC80BD", "#A6D854", "#FDCDAC", "#A6761D", "#80B1D3",
           "#FF9933", "#F4CAE4", "#FFFF66", "#2E8B57", "#CCEBC5", 
           "#1F78B4", "#BEBADA","#B3B3B3","#FF6600", "#800000","#23297A", "#FB8072")
timepalette = c("#90E0EF","#0077B6", "#ECB02B","#E75E02","#840404" )



# blue and yellow
twocolor = c("#FFC20A", "#0C7BDC")

############
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
# (1) Phyloseq object with metadata, genus-level agglomerated & batch corrected 
load("/Users/amycampbell/Documents/IowaWoundData2021/GenusLevelBatchCorrected.rda")
AnaerobeMappings = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/AnaerobeDB.csv")

BatchCorrectedPhyloseq = phyloseq::subset_taxa(BatchCorrectedPhyloseq, Order!="Chloroplast")
BatchCorrectedPhyloseq = phyloseq::subset_taxa(BatchCorrectedPhyloseq, Family!="Mitochondria")
BatchCorrectedPhyloseq = phyloseq::subset_samples(BatchCorrectedPhyloseq, study_id !="Mock")
BatchCorrectedPhyloseq@sam_data$TotalOTUsLeftBatchMitoChloro = colSums(BatchCorrectedPhyloseq@otu_table@.Data)

BatchCorrectedPhyloseq = phyloseq::subset_samples(BatchCorrectedPhyloseq, TotalOTUsLeftBatchMitoChloro>1200)

# Plot of remaining 406 patients + their wound type, location
#############################################################
FinalSampleDist = data.frame(BatchCorrectedPhyloseq@sam_data)
FinalSampleDist = FinalSampleDist %>% mutate(wound_type_string = case_when(wound_type==1 ~ "Pressure", 
                                                                                   wound_type==2 ~"Venous",
                                                                                   wound_type==4 ~"Surgical", 
                                                                                   wound_type==5 ~ "Traumatic", 
                                                                                   wound_type==6 ~ "Mixed Traumatic/Surgical",
                                                                                   wound_type==7 ~ "Other")
                                                     
)

FinalSampleDist = FinalSampleDist %>% mutate(woundloc_string = case_when(woundloc==1 ~ "Extremity", 
                                                                                 woundloc==2 ~"Trunk",
                                                                                 woundloc==3 ~"Head/Neck", 
                                                                                 woundloc==4 ~ "Inguinal"
)

)


WoundLocRun = ggplot(FinalSampleDist, aes(x=factor(Run), fill=factor(woundloc_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound location by Run") +xlab("MiSeq Run") 
WoundTypeRun = ggplot(FinalSampleDist, aes(x=factor(Run), fill=factor(wound_type_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound type by Run") +xlab("MiSeq Run") 

WeightedUnifracDists =  phyloseq::distance(BatchCorrectedPhyloseq, method="wunifrac")
OrdinationWeightedUnifrac = ordinate(BatchCorrectedPhyloseq,"PCoA", distance=WeightedUnifracDists)
BatchCorrectedPhyloseq@sam_data$woundcarepain = factor(BatchCorrectedPhyloseq@sam_data$woundcarepain)
BatchCorrectedPhyloseq@sam_data$wound_type = factor(BatchCorrectedPhyloseq@sam_data$wound_type)
BatchCorrectedPhyloseq@sam_data$woundloc = factor(BatchCorrectedPhyloseq@sam_data$woundloc)

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

# Save fit plot as fitplotDirichlet.jpeg
jpeg("/Users/amycampbell/Documents/IowaWoundData2021/fitplotDirichlet.jpeg")
plot(lplc, type = 'b', xlab = 'Dirichlet Components',ylab='Model Fit', main="Dirichlet Components by Laplace Model Fit (Genus Read Counts)") 
dev.off()

# get dirichlet fit object with maximum laplace-estimated likelihood (minimum negative)
bestFit <- FitDirichlets[[which.min(lplc)]]

GroupScores = data.frame(bestFit@group)
 
# Highest scoring component for each sample is that sample's assignment
GroupScores$assignment = sapply(1:nrow(GroupScores), function(x) which.max(GroupScores[x,]))
GroupScores$SampleID = row.names(GroupScores)

AssignmentMapping = GroupScores %>% select(SampleID,assignment )

SamDataForAssigning = data.frame(BatchCorrectedPhyloseq@sam_data)
SamDataForAssigning = SamDataForAssigning %>% left_join(AssignmentMapping, by="SampleID")
SamDataForAssigning = SamDataForAssigning %>% mutate(wound_type_string = case_when(wound_type==1 ~ "Pressure", 
                                                                 wound_type==2 ~"Venous",
                                                                 wound_type==4 ~"Surgical", 
                                                                 wound_type==5 ~ "Traumatic", 
                                                                 wound_type==6 ~ "Mixed Traumatic/Surgical",
                                                                 wound_type==7 ~ "Other")
                                   
)

SamDataForAssigning = SamDataForAssigning %>% mutate(woundloc_string = case_when(woundloc==1 ~ "Extremity", 
                                                               woundloc==2 ~"Trunk",
                                                               woundloc==3 ~"Head/Neck", 
                                                               woundloc==4 ~ "Inguinal"
)

)




BatchCorrectedPhyloseq@sam_data = sample_data(SamDataForAssigning)
rownames(BatchCorrectedPhyloseq@sam_data) = SamDataForAssigning$SampleID
colnames(BatchCorrectedPhyloseq@sam_data) = colnames(SamDataForAssigning)

# Plot genera with top contributions to each DMM
#################################################
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

# modified from https://microbiome.github.io/tutorials/DMM.html
 for (k in 1:Nclusters) {
   dfmelt = MeltedFitted %>% filter(Cluster==k) %>% arrange(-abs(Value))
   topOTUs = (dfmelt$OTU)[1:10]
   print(paste0("Cluster:", k))
   dfmelt = dfmelt %>% filter(OTU %in% topOTUs)
   plot = ggplot(dfmelt, aes(x=Genus,y=Value)) + geom_bar(stat="identity") + coord_flip() + labs(title=paste("Top Contributors to Cluster ", k ))
   plot$data$Genus = factor(plot$data$Genus, rev(dfmelt$Genus))
   
   ggsave(plot, file=paste0(paste0("/Users/amycampbell/Documents/IowaWoundData2021/DirichletContributorsCluster",k), ".png"))
   }


NMDSOrd <- ordinate(BatchCorrectedPhyloseq, "NMDS", "bray", weighted=T, trymax=200)
BatchCorrectedPhyloseq@sam_data$assignment = factor(BatchCorrectedPhyloseq@sam_data$assignment)
BatchCorrectedPhyloseq@sam_data$woundage = factor(BatchCorrectedPhyloseq@sam_data$woundage)






# some Ordinations of interest
###############################

BatchCorrectedPhyloseq@sam_data$woundAgeCat = if_else(BatchCorrectedPhyloseq@sam_data$woundage %in% c(1,2,3), "Acute", "Chronic")
BatchCorrectedPhyloseqNoModerates = subset_samples(BatchCorrectedPhyloseq, woundcarepain !=2)


DMMordPlot = plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="assignment", title="Ordination of Genus-aggregated OTUs by DMM Assignment")  + scale_color_manual(values=rev(rando18_1[c(10, 11, 15, 17, 14, 1,8)]))
ggsave(DMMordPlot, file="/Users/amycampbell/Documents/IowaWoundData2021/DMMordinationPlot.png")

dfresults = (data.frame(BatchCorrectedPhyloseq@sam_data) %>% select(SubjectID, assignment))
resultssved = read.csv("Documents/IowaWoundData2021/PlotsForSue2022/WoundMicrobiomeDataForSEG_AEC.csv") %>% select(StudyID, DMMClusterAssign)

topgenera = read.csv("Documents/IowaWoundData2021/PlotsForSue2022/WoundMicrobiomeDataForSEG_AEC.csv") %>% select(StudyID, Most_Abundant_Genus)
# Collapse obligate anaerobes into one group & visualize
################################################################
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


MeltedTop12 = PlotTop12Genera %>% psmelt() %>% select(GenusAdjust, Abundance, assignment, study_id, )
MeltedTop12 = MeltedTop12 %>% group_by(GenusAdjust, assignment) %>% summarise(mean(Abundance))
MeltedTop12$abundance = MeltedTop12$`mean(Abundance)`
generapresent  = ggplot(MeltedTop12, aes(x=as.factor(assignment), y=abundance, fill=GenusAdjust)) + geom_bar(stat = "identity", position = "stack",  color = NA) + scale_fill_manual(values=(color21)) + theme_minimal() + ggtitle("Top 12 Genera in Each DMM Component (Average %)") + xlab("DMM Assignment") + ylab("% Abundance")
ggsave(generapresent, file="/Users/amycampbell/Documents/IowaWoundData2021/TopGeneraDMMsWithoutMock.png")

SampleDist = data.frame(PlotTopGenera@sam_data)

woundtypePlot = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(wound_type_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound Types by DMM Assignment") +xlab("DMM Assignment") 

woundlocPlot = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(woundloc_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound Locations by DMM Assignment") +xlab("DMM Assignment") 

runplot = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(Run))) +
  geom_bar() + scale_fill_manual(values=twocolor) + theme_minimal()  + ggtitle("Sequencing Run by DMM Assignment") +xlab("DMM Assignment") 

timeplot1 = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(woundage))) +
  geom_bar() + scale_fill_manual(values=timepalette) + theme_minimal()  + ggtitle("Wound age by DMM Assignment") +xlab("DMM Assignment") 

SampleDistPlot = gridExtra::grid.arrange(woundtypePlot, woundlocPlot, runplot, timeplot1)
ggsave(SampleDistPlot, width=10, height=8, file="/Users/amycampbell/Documents/IowaWoundData2021/SampleDistPlot22.pdf")


painplot1 = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(woundcarepain))) +
  geom_bar() + scale_fill_manual(values=c("#FFFFCC", "#FFB266", "#CC6600", "#990000", "#E0E0E0")) + theme_minimal()  + ggtitle("Wound pain by DMM Assignment") +xlab("DMM Assignment") +  theme(legend.position = "None")

painplot2 = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(woundcarepain))) +
  geom_bar(position="fill") + scale_fill_manual(values=c("#FFFFCC", "#FFB266", "#CC6600", "#990000", "#E0E0E0")) + theme_minimal()  + ggtitle("Wound pain by DMM Assignment") +xlab("DMM Assignment") + ylab("Proportion")

SampleDistNoModerates = SampleDist %>% filter(woundcarepain!=2)
painplot3 = ggplot(SampleDistNoModerates, aes(x=factor(assignment), fill=factor(woundcarepain))) +
  geom_bar() + scale_fill_manual(values=c("#FFFFCC", "#FFB266", "#990000", "#E0E0E0")) + theme_minimal()  + ggtitle("Wound pain by DMM Assignment") +xlab("DMM Assignment") + ylab("Proportion") +  theme(legend.position = "None")
painplot4 = ggplot(SampleDistNoModerates, aes(x=factor(assignment), fill=factor(woundcarepain))) +
  geom_bar(position="fill") + scale_fill_manual(values=c("#FFFFCC", "#FFB266", "#990000", "#E0E0E0")) + theme_minimal()  + ggtitle("Wound pain by DMM Assignment") +xlab("DMM Assignment") + ylab("Proportion")

SamplePainPlot = gridExtra::grid.arrange(painplot1, painplot2, painplot3, painplot4, ncol=2, widths=c(8, 10))
ggsave(SamplePainPlot, width=10, height=8, file="/Users/amycampbell/Documents/IowaWoundData2021/SamplePainPlot22.pdf")



save(BatchCorrectedPhyloseq, file="/Users/amycampbell/Documents/IowaWoundData2021/BatchCorrectedPhyloseqWithDMMs.rda")

countplot = ggplot(SampleDist, aes(y=TotalOTUsLeftBatchMitoChloro, x=assignment)) + geom_boxplot(fill="#CC6600") + xlab("DMM Component ") + ylab("Total #OTUs/ASVs") + theme_classic()
Shannonplot = plot_richness(BatchCorrectedPhyloseq, x="assignment", measures=c("Shannon")) + geom_boxplot(fill="#CC6600") + theme_classic() + xlab("DMM Component") + ylab("Shannon Diversity")
Richnessplot = plot_richness(BatchCorrectedPhyloseq, x="assignment", measures=c("Observed")) + geom_boxplot(fill="#CC6600") + theme_classic() + xlab("DMM Component") + ylab("Observed Genera")
gridExtra::grid.arrange(countplot, Shannonplot, Richnessplot, ncol=3)

#  Extract these diversity metrics for adding to the DF later 
RichnessDF = Richnessplot$data %>% select(SampleID, value)
colnames(RichnessDF) = c("SampleID", "Genus_Richness")
ShannonDF = Shannonplot$data %>% select(SampleID, value)
colnames(ShannonDF) =  c("SampleID", "Genus_Shannon")


# 1 = <7 days
# 2 = 8-30 days
# 3 = 31 days to 90 days
# 4 = 91 days to 1 year#
# 5 = >1 year
#90E0EF # lowest blue
#0077B6 #high blue
#ECB02B # low orange 
#E75E02 high orange
#840404 high red
#D3D3D3 gray

# Final dataframe for SEG
##########################
FinalDF = SamDataForAssigning %>% select(SampleID, study_id, assignment)
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


# Some genera of interest (Data sent to SEG in February 2022)
##############################################################

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

# Add richness info
FinalDF = FinalDF %>% left_join(RichnessDF, by="SampleID")
FinalDF = FinalDF %>% left_join(ShannonDF, by="SampleID")


write.csv(FinalDF, file="/Users/amycampbell/Documents/IowaWoundData2021/WoundMicrobiomeDataForSEG_AEC.csv")

