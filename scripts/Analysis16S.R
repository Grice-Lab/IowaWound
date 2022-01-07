# Amy Campbell
# Jan 2022
# Analysis using decontaminated data object output by Qiime2OutputNov.R 


library("phyloseq")
library("dplyr")
library("stringr")
library("decontam")
library("ggplot2")
library("DESeq2")
library("vegan")
library("ggpubr")
library("RColorBrewer")
library(dplyr)
library(lattice)
library(tidyr)
library(ggplot2)
library(DirichletMultinomial)

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
# blue and yellow
twocolor = c("#FFC20A", "#0C7BDC")

load("/Users/amycampbell/Documents/IowaWoundData2021/PhyloSeqObjectPostDecontamMerged.rda")

load("/Users/amycampbell/Documents/IowaWoundData2021/PhyloSeqObjectPostDecontam.rda")

PhyloseqObjectUpdated

WeightedUnifracDists =  phyloseq::distance(PhyloseqObjectUpdated, method="wunifrac")
OrdinationWeightedUnifrac = ordinate(PhyloseqObjectUpdated,"PCoA", distance=WeightedUnifracDists)
PlotWeightedUnifracRun = plot_ordination(PhyloseqObjectUpdated, OrdinationWeightedUnifrac, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (OTU-level weighted UniFrac distance) in MiSeq Runs After Decontamination") + theme_classic()
ggsave(PlotWeightedUnifracRun, file="/Users/amycampbell/Documents/IowaWoundData2021/WeightedUnifrac.png")


# high miseq35 axis 2 point = IowaWound.Human.969
# high miseq32 axis 1 point = 120499, 120529, IowaWound.Human.1517, IowaWound.Human.1670

# all type 4-- surgical. 
# 969 is almost entirely bacteroidetes with like a little bit of proteobacteria
# 1517 and 1670 dominated by firmicutes and actinobacteriota (not unusual)
# 120499 and 120529 dominated by firmicutes as well 
WeirdSamples = subset_samples(PhyloseqObjectUpdated, SampleID %in% c("120499", "120529", "IowaWound.Human.1517", "IowaWound.Human.1670","IowaWound.Human.969"))
WeirdSamplesRelabundance = WeirdSamples %>% transform_sample_counts(function(x) {x/sum(x)})
plot_bar(WeirdSamplesRelabundance,"SampleID", fill="Phylum" )+ scale_fill_manual(values=randpalette32)
PhyloseqObjectUpdated@sam_data$wound_type = as.factor(PhyloseqObjectUpdated@sam_data$wound_type)
PhyloseqObjectUpdated@sam_data$woundloc = as.factor(PhyloseqObjectUpdated@sam_data$woundloc)
PhyloseqObjectNoWeirdSamples = subset_samples(PhyloseqObjectUpdated, !(SampleID %in% c("120499", "120529", "IowaWound.Human.1517", "IowaWound.Human.1670","IowaWound.Human.969")))

PlotWeightedUnifracRun = plot_ordination(PhyloseqObjectUpdated, OrdinationWeightedUnifrac, color="wound_type") + scale_colour_manual(values=randpalette32) + ggtitle("Principal Coordinates (OTU-level weighted UniFrac distance) in MiSeq Runs After Decontamination") + theme_classic()

# Just zooming in to see if theres any divergence between wound types when you exclude those weird points. Not really.
PlotWeightedUnifracRunNoWeird = plot_ordination(PhyloseqObjectNoWeirdSamples, OrdinationWeightedUnifrac, color="wound_type") + scale_colour_manual(values=rev(palette8colors)) + ggtitle("Principal Coordinates (OTU-level weighted UniFrac distance) in MiSeq Runs After Decontamination") + theme_classic() + xlim(-.02,.025) + ylim(-.02, .02) 
PlotWeightedUnifracRunNoWeird
PhyloseqObjectNoSurgery = subset_samples(PhyloseqObjectUpdated, wound_type!=4)
# what if you exclude surgical ? still not really  Type 7("Other") is all pretty close to the center
PlotWeightedUnifracRunNoSurg = plot_ordination(PhyloseqObjectNoSurgery, OrdinationWeightedUnifrac, color="wound_type") + scale_colour_manual(values=rev(palette8colors)) + ggtitle("Principal Coordinates Weighted Unifrac (Non-surgical wounds only) ") + theme_classic() + xlim(-.02,.025) + ylim(-.02, .02) 
ggsave(PlotWeightedUnifracRunNoSurg,file="/Users/amycampbell/Documents/IowaWoundData2021/WunifracNoSurgical.png" )

# DMM Clustering 
################

# Aggregate to species level and above
#What does it actually do with teh NA? do NA ones just get taken away??? seems bad since there are so many Staphylococcus_NAs for example
# Excluding the bad_empty keeps it from just dumping the OTUs that don't resolve to species (they get aggregated at the next highest level which is fine/what I want)
# NAstrings = c("unidentified_marine",
#               "uncultured_soil",
#               "unidentified_eubacterium",
#               "uncultured_organism", NA,
#               "NA", "uncultured_microorganism",
#               "uncultured_bacterium", "uncultured_candidate", "uncultured_compost",
#               "unidentified", "uncultured_prokaryote", "uncultured_rumen", "wastewater_metagenome")
Specieslevel = tax_glom(PhyloseqObjectUpdated, taxrank= "Species")

Counts = Specieslevel@otu_table


# Filter to species prevalent in at least 10 samples 
PrevalenceInfo = data.frame(Counts)
PrevalenceInfo$Prevalence = rowSums(PrevalenceInfo != 0)
PrevalenceInfoFilter = PrevalenceInfo %>% filter(Prevalence >= 10)
TaxaKeepDMM = row.names(PrevalenceInfoFilter)

#
PhyloDMM = prune_taxa(TaxaKeepDMM, Specieslevel)


CountsDMM = t(PhyloDMM@otu_table@.Data)

CountsDMMinput = data.matrix(CountsDMM)

densityplot(log10(colSums(CountsDMMinput)), xlim=range(log10(colSums(CountsDMMinput))),xlab="Taxon representation (log 10 count)")
densityplot((colSums(CountsDMMinput)), xlim=range((colSums(CountsDMMinput))),xlab="Taxon representation (raw count)")


# Fit dirichlet multinomial mixture models with 1,2,...,20 components and plot laplace approximation 
FitDirichlets = mclapply(1:20, dmn, count=CountsDMMinput, verbose=TRUE, seed=19104)
lplc <- sapply(FitDirichlets, laplace)
plot(lplc, type = 'b', xlab = 'Dirichlet Components',ylab='Model Fit', main="Dirichlet Components by Laplace Model Fit (ASV Counts)") 

# get dirichlet fit object with maximum laplace-estimated likelihood (minimum negative)
bestFit <- FitDirichlets[[which.min(lplc)]]

# 
GroupScores = data.frame(bestFit@group)
#colnames(GroupScores) = c("Component1", "Component2","Component3", "Component4", "Component5", "Component6")

# Highest scoring component for each sample is that sample's assignment
GroupScores$assignment = sapply(1:nrow(GroupScores), function(x) which.max(GroupScores[x,]))
GroupScores$SampleID = row.names(GroupScores)

AssignmentMapping = GroupScores %>% select(SampleID,assignment )


SamDataForAssigning = data.frame(Specieslevel@sam_data)
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




keepSamplenames = Specieslevel@sam_data$SampleID


Specieslevel@sam_data = sample_data(SamDataForAssigning)
rownames(Specieslevel@sam_data) = sample_names(Specieslevel)
colnames(Specieslevel@sam_data) = colnames(SamDataForAssigning)

# What are the most prevalent species?
# Staphylococcus_NA, Bacillus_NA, Cutibacterium_NA, Delftia_NA, Streptococcus_NA
MostPrev = row.names((PrevalenceInfo %>% arrange(-Prevalence))[1:20, ])
Specieslevel@tax_tanteble@.Data[MostPrev, ]
hist(PrevalenceInfo$Prevalence)

PlotTopGenera = Specieslevel

PlotTopGenera1 = subset_samples(PlotTopGenera, assignment==1)
PlotTopGenera2 = subset_samples(PlotTopGenera, assignment==2)
PlotTopGenera3 = subset_samples(PlotTopGenera, assignment==3)
PlotTopGenera4 = subset_samples(PlotTopGenera, assignment==4)
PlotTopGenera5 = subset_samples(PlotTopGenera, assignment==5)
PlotTopGenera6 = subset_samples(PlotTopGenera, assignment==6)



SampleDist = data.frame(PlotTopGenera@sam_data)
ggplot(SampleDist, aes(x=factor(assignment), fill=factor(wound_type_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound Types by DMM Assignment") +xlab("DMM Assignment") 

ggplot(SampleDist, aes(x=factor(assignment), fill=factor(woundloc_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound Locations by DMM Assignment") +xlab("DMM Assignment") 

painplot = ggplot(SampleDist, aes(x=factor(assignment), fill=factor(woundcarepain))) +
  geom_bar() + scale_fill_manual(values=c("#CC6600", "#990000", "#E0E0E0")) + theme_minimal()  + ggtitle("Wound pain by DMM Assignment") +xlab("DMM Assignment") 

PlotTopGenera = PlotTopGenera %>% transform_sample_counts(function(x) {x/sum(x)})


PlotTopGenera@otu_table



Specieslevel


plot_bar(Specieslevel, )
