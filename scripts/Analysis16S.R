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

randpalette18=c("#B300B3","#E6AB02",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00",
                "#681A1A","#7570B3",
                "#1F78B4","#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#FFFF33",
                "#33A02C","#E6F5C9")



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
load("/Users/amycampbell/Documents/IowaWoundData2021/PublicationData/GenusLevelBatchCorrected_Jan22.rda")

AnaerobeMappings = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/AnaerobeDB.csv")
WoundDepth = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/wound_depth_covariate.csv")

BatchCorrectedPhyloseq = phyloseq::subset_taxa(BatchCorrectedPhyloseq, Order!="Chloroplast")
BatchCorrectedPhyloseq = phyloseq::subset_taxa(BatchCorrectedPhyloseq, Family!="Mitochondria")
BatchCorrectedPhyloseq = phyloseq::subset_samples(BatchCorrectedPhyloseq, study_id !="Mock")
BatchCorrectedPhyloseq@sam_data$TotalOTUsLeftBatchMitoChloro = colSums(BatchCorrectedPhyloseq@otu_table@.Data)

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
WoundDepth$study_id = sapply(WoundDepth$study_id, as.character)


FinalSampleDist = FinalSampleDist %>% left_join(WoundDepth, by="study_id")


WoundLocRun = ggplot(FinalSampleDist, aes(x=factor(Run), fill=factor(woundloc_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound location by Run") +xlab("MiSeq Run") 
WoundTypeRun = ggplot(FinalSampleDist, aes(x=factor(Run), fill=factor(wound_type_string))) +
  geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound type by Run") +xlab("MiSeq Run") 


BatchCorrectedPhyloseq@sam_data = sample_data(FinalSampleDist)
sample_names(BatchCorrectedPhyloseq) = (BatchCorrectedPhyloseq@sam_data$study_id)

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
pdf("/Users/amycampbell/Documents/IowaWoundData2021/PaperFigs/fitplotDirichlet.pdf", height=5,width=5)
plot(lplc, type = 'b', xlab = 'Dirichlet Components',ylab='Model Fit', main="Dirichlet Components by Laplace Model Fit (Genus Read Counts)") 
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
   plot
   ggsave(plot, file=paste0(paste0("/Users/amycampbell/Documents/IowaWoundData2021/DirichletContributorsCluster",k), ".png"))
   }


dfresults = (data.frame(BatchCorrectedPhyloseq@sam_data) %>% select(SubjectID, assignment))
resultssved = read.csv("Documents/IowaWoundData2021/PlotsForSue2022/WoundMicrobiomeDataForSEG_AEC.csv") %>% select(StudyID, DMMClusterAssign)

topgenera = read.csv("Documents/IowaWoundData2021/PlotsForSue2022/WoundMicrobiomeDataForSEG_AEC.csv") %>% select(StudyID, Most_Abundant_Genus)

# Ordination visualizations
###########################

# Some sample variable coding
#############################
BatchCorrectedPhyloseq@sam_data$woundcarepain = factor(BatchCorrectedPhyloseq@sam_data$woundcarepain)
BatchCorrectedPhyloseq@sam_data$wound_type = factor(BatchCorrectedPhyloseq@sam_data$wound_type)
BatchCorrectedPhyloseq@sam_data$woundloc = factor(BatchCorrectedPhyloseq@sam_data$woundloc)
BatchCorrectedPhyloseq@sam_data$dressingcat = factor(BatchCorrectedPhyloseq@sam_data$dressingcat)
BatchCorrectedPhyloseq@sam_data$woundage = factor(if_else(BatchCorrectedPhyloseq@sam_data$woundage %in% c(1, 2), "Acute", "Chronic"))
BatchCorrectedPhyloseq@sam_data$WoundVac = if_else(BatchCorrectedPhyloseq@sam_data$dressingcat %in% c(2,3), "Yes", "No")
sample_names(BatchCorrectedPhyloseq) = BatchCorrectedPhyloseq@sam_data$SampleID



# NMDS
########
NMDSOrd <- ordinate(BatchCorrectedPhyloseq, "NMDS", "bray", weighted=T, trymax=200)
plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="Run", title="samples") + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3" ,"#E7298A", "#66A61E" ,"#E6AB02", "#A6761D" ,"#666666"))

metadf <- data.frame(sample_data(BatchCorrectedPhyloseq))
BrayDist = phyloseq::distance(BatchCorrectedPhyloseq,method = "bray")

permanovaRunBray <- adonis(BrayDist ~ Run, data = metadf)
# there are still batch effects!!


# Weighted unifrac PCOA
##########################
WeightedUnifracDists=phyloseq::distance(BatchCorrectedPhyloseq,method = "wunifrac")
OrdinationWeightedUnifrac = ordinate(BatchCorrectedPhyloseq,"PCoA", distance=WeightedUnifracDists)
plot_ordination(BatchCorrectedPhyloseq, OrdinationWeightedUnifrac, type="samples", color="Run") +scale_color_manual(values=c("#1B9E77", "#D95F02")) 

permanovaRunWunifrac<- adonis(WeightedUnifracDists ~ Run, data = metadf)

plot_ordination(BatchCorrectedPhyloseq, OrdinationWeightedUnifrac, type="samples", color="Run")

metadf <- data.frame(sample_data(BatchCorrectedPhyloseq))

unifrac.dist <- UniFrac(BatchCorrectedPhyloseq, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanovaWoundVac <- adonis(unifrac.dist ~ WoundVac, data = metadf)
permanovaWoundVac <- adonis(unifrac.dist ~ WoundVac, data = metadf)
mod <- betadisper(dis, groups)
BatchCorrectedPhyloseq@sam_data$assignment = factor(BatchCorrectedPhyloseq@sam_data$assignment)
BatchCorrectedPhyloseq@sam_data$woundage = factor(BatchCorrectedPhyloseq@sam_data$woundage)
#sample_names(BatchCorrectedPhyloseq) = BatchCorrectedPhyloseq@sam_data$study_id
BatchCorrectedPhyloseq@sam_data$WoundVac = if_else(BatchCorrectedPhyloseq@sam_data$dressingcat %in% c(2,3), "Yes", "No")

plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="Run", title="samples") + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3" ,"#E7298A", "#66A61E" ,"#E6AB02", "#A6761D" ,"#666666"))
plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="WoundVac", title="samples") + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3" ,"#E7298A", "#66A61E" ,"#E6AB02", "#A6761D" ,"#666666"))
plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="Run", title="samples") + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3" ,"#E7298A", "#66A61E" ,"#E6AB02", "#A6761D" ,"#666666"))
plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="WoundVac", title="samples") + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3" ,"#E7298A", "#66A61E" ,"#E6AB02", "#A6761D" ,"#666666"))

plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="dressingcat", title="samples") + theme_classic() + scale_color_brewer(palette = "Paired")
plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="woundage")+theme_classic()
plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="Run")+theme_classic() + scale_color_brewer(palette="Dark2")
plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="woundage")+theme_classic() + scale_color_brewer(palette="Dark2")



BatchCorrectedPhyloseq@sam_data$woundAgeCat = if_else(BatchCorrectedPhyloseq@sam_data$woundage %in% c(1,2,3), "Acute", "Chronic")
BatchCorrectedPhyloseq@sam_data$dressingcat = factor(BatchCorrectedPhyloseq@sam_data$dressingcat)

BatchCorrectedPhyloseqNoModerates = subset_samples(BatchCorrectedPhyloseq, woundcarepain !=2)

DMMordPlot = plot_ordination(BatchCorrectedPhyloseq, NMDSOrd, type="samples", color="assignment", title="Ordination of Genus-aggregated OTUs by DMM Assignment")  + scale_color_manual(values=rev(rando18_1[c(10, 11, 15, 17, 14, 1,8)]))

DMMordPlot = DMMordPlot + theme_classic()
ggsave(DMMordPlot, file="/Users/amycampbell/Documents/IowaWoundData2021/PaperFigs/DMMordinationPlot.pdf", width=8,height=8)

ggsave(DMMordPlot, file="/Users/amycampbell/Documents/IowaWoundData2021/DMMordinationPlot.png")


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

#
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
ggsave(generapresent, file="/Users/amycampbell/Documents/IowaWoundData2021/TopGeneraDMMsWithoutMock.png")
ggsave(generapresent, file="/Users/amycampbell/Documents/IowaWoundData2021/PaperFigs/TopGeneraDMMsWithoutMock.pdf",height=8,width=8)

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
save(BatchCorrectedPhyloseqRankingPercent, file="/Users/amycampbell/Documents/IowaWoundData2021/PublicationData/RelativeAbundanceGenus.rda")
save(CLRTransformedPct, file="/Users/amycampbell/Documents/IowaWoundData2021/PublicationData/CLRAbundanceGenus.rda")

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



# Plot phylum by wound type
############################
Phylum_abundance = plot_bar(BatchCorrectedPhyloseq, "SampleID","Abundance", "Phylum")


AnaerobeOTUs = row.names(DF_Anaerobes %>% filter(Anaerobe == 1))
AllOTUs = data.frame(BatchCorrectedPhyloseq@tax_table)
AllOTUs$OTU = row.names(AllOTUs)
AnaerobeGenera = AllOTUs %>% filter(OTU %in% AnaerobeOTUs)

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

pdf(file="~/Documents/IowaWoundData2021/PaperFigs/ContributorsDMMs.pdf", width=7.2, height=4.3)
gridExtra::grid.arrange(plot_top5cluster1+ theme_classic() + theme(legend.position = "none"), plot_top5cluster2+ theme_classic() + theme(legend.position = "none"), plot_top5cluster3+theme_classic() + theme(legend.position = "none"), plot_top5cluster4+theme_classic() + theme(legend.position = "none"))
dev.off()



# Phylum Level 

df_phyla= BatchCorrectedPhyloseq %>% 
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance >.01) %>%
  group_by(Phylum)


dataframe_counts = data.frame(df_phyla)
sampnames = (unique(dataframe_counts$Sample))

# Order by type followed by firmicutes abundance
df_phyla_Firm = df_phyla %>% filter(Phylum=="Firmicutes") %>% select(wound_type_string, SampleID, Abundance)
MissingFirm = unique(setdiff(df_phyla$SampleID, df_phyla_Firm$SampleID))
MissingFirm = df_phyla %>% filter(SampleID %in% MissingFirm) %>% ungroup() %>% select(wound_type_string, SampleID) %>% unique()
MissingFirm$Phylum ="Firmicutes"
MissingFirm$Abundance = 0.0
MissingFirm = MissingFirm %>% select(Phylum, wound_type_string, SampleID, Abundance)
df_phyla_Firm = rbind(MissingFirm, df_phyla_Firm)



#df_phyla_Firm %>% arrange(wound_type_string, -Firmicutes)

Desired_order=unique((df_phyla_Firm %>% arrange(wound_type_string, -Abundance))$SampleID)
plotPhyla= ggplot(df_phyla, aes(x=SampleID, y=Abundance, fill=Phylum)) + geom_bar(stat="identity") + ggtitle("Phylum-level composition") + scale_fill_manual(values=(randpalette18))+ theme_minimal() + theme(axis.text.x=element_text(angle=90))
plotPhyla$data$SampleID = factor(plotPhyla$data$SampleID, levels =Desired_order)
ggsave(plotPhyla, file="~/Documents/IowaWoundData2021/PhylumAbundance_All_By_Type.pdf", width=20, height=6)

Desired_order


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
ggsave(Top12GeneraPlotSample, file="~/Documents/IowaWoundData2021/TopGenusAbundance_All_By_Cluster.pdf", width=20, height=6)

FinalDF$AnaerobePrev = if_else(FinalDF$AnaerobicGenusAbundance>.001, 1, 0)
FinalDF$CorynePrev = if_else(FinalDF$CorynebacteriumAbundance>.001, 1, 0)
FinalDF$PseudPrev = if_else(FinalDF$PseudomonasAbundance>.001, 1, 0)
FinalDF$StaphPrev = if_else(FinalDF$StaphylococcusAbundance>.001, 1, 0)
FinalDF$StrepPrev = if_else(FinalDF$StreptococcusAbundance>.001, 1, 0)

sd(FinalDF$AnaerobicGenusAbundance)
sd(FinalDF$CorynebacteriumAbundance)
sd(FinalDF$PseudomonasAbundance)
sd(FinalDF$StaphylococcusAbundance)
sd(FinalDF$StreptococcusAbundance)

sort(table(FinalDF$Most_Abundant_Genus))

