# Amy Campbell
# Jan 2022
# Analysis using decontaminated data object output by Qiime2OutputNov.R 


library("phyloseq")
library("dplyr")
library("stringr")
library("decontam")
library("ggplot2")
library("vegan")
library("ggpubr")
library("lattice")
library("DirichletMultinomial")

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

# blue and yellow
twocolor = c("#FFC20A", "#0C7BDC")

load("/Users/amycampbell/Documents/IowaWoundData2021/GenusLevelBatchCorrected.rda")

# BatchCorrectedPhyloseq

WeightedUnifracDists =  phyloseq::distance(BatchCorrectedPhyloseq, method="wunifrac")
OrdinationWeightedUnifrac = ordinate(BatchCorrectedPhyloseq,"PCoA", distance=WeightedUnifracDists)
BatchCorrectedPhyloseq@sam_data$woundcarepain = factor(BatchCorrectedPhyloseq@sam_data$woundcarepain)
BatchCorrectedPhyloseq@sam_data$wound_type = factor(BatchCorrectedPhyloseq@sam_data$wound_type)
BatchCorrectedPhyloseq@sam_data$woundloc = factor(BatchCorrectedPhyloseq@sam_data$woundloc)


PlotWeightedUnifracPain = plot_ordination(BatchCorrectedPhyloseq, OrdinationWeightedUnifrac, color="woundcarepain") + scale_colour_manual(values=Paincolors) + ggtitle("Principal Coordinates (Genus-level weighted UniFrac distance) Wound Pain") + theme_classic()
PlotWeightedUnifracType = plot_ordination(BatchCorrectedPhyloseq, OrdinationWeightedUnifrac, color="wound_type") + scale_colour_manual(values=rando18_1) + ggtitle("Principal Coordinates (Genus-level weighted UniFrac distance) Wound Type") + theme_classic()
PlotWeightedUnifracLoc= plot_ordination(BatchCorrectedPhyloseq, OrdinationWeightedUnifrac, color="woundloc") + scale_colour_manual(values=rando18_1) + ggtitle("Principal Coordinates (Genus-level weighted UniFrac distance) Wound Location") + theme_classic()


# DMM Clustering 
################

# Aggregated to Genus level and above

# Test DMM on this
Counts = BatchCorrectedPhyloseq@otu_table

# Filter to species prevalent in at least 10 samples 
PrevalenceInfo = data.frame(Counts)
PrevalenceInfo$Prevalence = rowSums(PrevalenceInfo != 0)

sort(PrevalenceInfo$Prevalence)

PrevalenceInfoFilter = PrevalenceInfo %>% filter(Prevalence >= 10)
TaxaKeepDMM = row.names(PrevalenceInfoFilter)

#
PhyloDMM = prune_taxa(TaxaKeepDMM, BatchCorrectedPhyloseq)

CountsDMM = t(PhyloDMM@otu_table@.Data)

CountsDMMinput = data.matrix(CountsDMM)

densityplot(log10(colSums(CountsDMMinput)), xlim=range(log10(colSums(CountsDMMinput))),xlab="Taxon representation (log 10 count)")
densityplot((colSums(CountsDMMinput)), xlim=range((colSums(CountsDMMinput))),xlab="Taxon representation (raw count)")

FitDirichlets = mclapply(1:12, dmn, count=CountsDMMinput, verbose=TRUE, seed=19104)
lplc <- sapply(FitDirichlets, laplace)
plot(lplc, type = 'b', xlab = 'Dirichlet Components',ylab='Model Fit', main="Dirichlet Components by Laplace Model Fit (Genus Read Counts)") 

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
rownames(BatchCorrectedPhyloseq@sam_data) = sample_names(Specieslevel)
colnames(BatchCorrectedPhyloseq@sam_data) = colnames(SamDataForAssigning)

# What are the most prevalent species?
# Staphylococcus_NA, Bacillus_NA, Cutibacterium_NA, Delftia_NA, Streptococcus_NA
MostPrev = row.names((PrevalenceInfo %>% arrange(-Prevalence))[1:20, ])
BatchCorrectedPhyloseq@tax_table@.Data[MostPrev, ]
hist(PrevalenceInfo$Prevalence)

PlotTopGenera = BatchCorrectedPhyloseq

#Genuslevel=Specieslevel %>% tax_glom(taxrank="Genus")

genusleveltax = data.frame(PlotTopGenera@tax_table@.Data)
write.csv(genusleveltax, file="/Users/amycampbell/Documents/IowaWoundData2021/GenusLevelTaxTable.csv")

row.names(PlotTopGenera@sam_data) = BatchCorrectedPhyloseq@sam_data$SampleID
PlotTopGenera = PlotTopGenera %>% transform_sample_counts(function(x) {x/sum(x)})

PlotTopGenera1 = subset_samples(PlotTopGenera, assignment==1)
PlotTopGenera2 = subset_samples(PlotTopGenera, assignment==2)
PlotTopGenera3 = subset_samples(PlotTopGenera, assignment==3)
PlotTopGenera4 = subset_samples(PlotTopGenera, assignment==4)
PlotTopGenera5 = subset_samples(PlotTopGenera, assignment==5)

topdf1 = data.frame(PlotTopGenera1@otu_table@.Data)
topdf1$rowmean = rowMeans(topdf1)
topGenera = row.names( topdf1 %>% arrange(-rowmean))[1:9]


topdf2 = data.frame(PlotTopGenera2@otu_table@.Data)
topdf2$rowmean = rowMeans(topdf2)
topGenera = append(topGenera, row.names( topdf2 %>% arrange(-rowmean))[1:9])

topdf3 = data.frame(PlotTopGenera3@otu_table@.Data)
topdf3$rowmean = rowMeans(topdf3)
topGenera = append(topGenera, row.names( topdf3 %>% arrange(-rowmean))[1:9])

topdf4 = data.frame(PlotTopGenera4@otu_table@.Data)
topdf4$rowmean = rowMeans(topdf4)
topGenera = append(topGenera, row.names( topdf4 %>% arrange(-rowmean))[1:9])

topdf5 = data.frame(PlotTopGenera5@otu_table@.Data)
topdf5$rowmean = rowMeans(topdf5)
topGenera = append(topGenera, row.names( topdf5 %>% arrange(-rowmean))[1:9])

topGenera = unique(topGenera)

PlotTop9Genera = prune_taxa(topGenera, PlotTopGenera)


tax_labels = data.frame(PlotTop9Genera@tax_table@.Data) %>% select("Genus")
tax_labels$Anaerobe = c(0,0,1, 0, 0, 1, 1, 1, 1,1,1,0, 0,1,1,0,1, 0)


MeltedTop9 = PlotTop9Genera %>% psmelt() %>% select(Genus, Abundance, assignment, study_id)
MeltedTop9 = MeltedTop9 %>% group_by(Genus, assignment) %>% summarise(mean(Abundance))
MeltedTop9$abundance = MeltedTop9$`mean(Abundance)`
ggplot(MeltedTop9, aes(x=as.factor(assignment), y=abundance, fill=Genus)) + geom_bar(stat = "identity", position = "stack",  color = NA) + scale_fill_manual(values=(rando18_1)) + theme_minimal() + ggtitle("Top 9 Genera in Each DMM Component (Average %)") + xlab("DMM Assignment") + ylab("% Abundance")



AnaerobicInfoIncluded = MeltedTop9 %>% left_join(tax_labels, by="Genus")
AnaerobicInfoIncluded = AnaerobicInfoIncluded %>% mutate(anaerobicGenus = case_when(Anaerobe==1~ "Anaerobes", 
                                                                                    Anaerobe==0 ~ toString(Genus[1])))
ggplot(AnaerobicInfoIncluded, aes(x=as.factor(assignment), y=abundance, fill=anaerobicGenus)) + geom_bar(stat = "identity", position = "stack",  color = NA) + scale_fill_manual(values=(randocolor18)) + theme_minimal() + ggtitle("Top 9 Genera in Each DMM Component (Average %, Anaerobes Grouped Together)") + xlab("DMM Assignment") + ylab("% Abundance")




SampleDist = data.frame(PlotTopGenera@sam_data)
ggplot(SampleDist, aes(x=factor(assignment), fill=factor(wound_type_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound Types by DMM Assignment") +xlab("DMM Assignment") 

ggplot(SampleDist, aes(x=factor(assignment), fill=factor(woundloc_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound Locations by DMM Assignment") +xlab("DMM Assignment") 



runplot = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(Run))) +
  geom_bar() + scale_fill_manual(values=c("#FFFFCC", "#FFB266", "#CC6600", "#990000", "#E0E0E0")) + theme_minimal()  + ggtitle("Run by DMM Assignment") +xlab("DMM Assignment") 

runplot = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(Run))) +
  geom_bar() + scale_fill_manual(values=twocolor) + theme_minimal()  + ggtitle("Run by DMM Assignment") +xlab("DMM Assignment") 



painplot2 = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(woundcarepain))) +
  geom_bar(position="fill") + scale_fill_manual(values=c("#FFFFCC", "#FFB266", "#CC6600", "#990000", "#E0E0E0")) + theme_minimal()  + ggtitle("Wound pain by DMM Assignment") +xlab("DMM Assignment") 



save(BatchCorrectedPhyloseq, file="/Users/amycampbell/Documents/IowaWoundData2021/BatchCorrectedPhyloseqWithDMMs.rda")







