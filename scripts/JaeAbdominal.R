# Amy Campbell
# Getting subset of for Jae of study IDs she sent of abdominal surgical wounds
library(dplyr)
library(phyloseq)
library(stringr)

# function that turns genus name of undetermined genus into 'FamilyName_NA' or Order_NA_NA' instead of just NA
# (useful for making relative abundance table)
Genus_code = function(tax_table_row){
  if(!(tax_table_row["Genus"] %in% c("uncultured", "NA"))){
    return(tax_table_row["Genus"])
  }else if(!(tax_table_row["Family"] %in% c("uncultured", "NA"))){
    return(paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"))
  }else if( !(tax_table_row["Order"] %in% c("uncultured", "NA"))){
    return(paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"))
  }else if(!(tax_table_row["Class"] %in% c("uncultured", "NA"))){
    return(paste(tax_table_row["Class"],paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"), sep="_"))
  }else if(!(tax_table_row["Phylum"] %in% c("uncultured", "NA"))){
    return(paste(tax_table_row["Phylum"], paste(tax_table_row["Class"],paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"), sep="_"), sep="_"))
  }else if(tax_table_row["Phylum"]=="NA"){
    print("ouch")
    return("NULL")
  }else{
    return("NULL")
  }
}

############
# Load data 
############
Abdominal_IDs = read.csv("~/Documents/IowaWoundData2021/Jae_AbdominalSurgicalIDs.csv")
WoundMicrobiomeData = read.csv("~/Documents/IowaWoundData2021/PlotsForSue2022/WoundMicrobiomeDataForSEG_AEC.csv")
AnaerobeMappings = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/AnaerobeDB.csv")


####################################################################################################
# The summary microbiome stats that we used for Aim 2 paper, subsetted to the abdominal surg wounds
####################################################################################################
WoundMicrobiomeDataAbdominal = WoundMicrobiomeData %>% filter(StudyID %in%  Abdominal_IDs$study_id)
write.csv(WoundMicrobiomeDataAbdominal, "~/Documents/IowaWoundData2021/JaeAbdominal/AbdominalWoundsMicrobiomeStats.csv")


#####################################################
# ASV/OTU Counts of genera (anaerobes not aggregated)
#####################################################
load("~/Documents/IowaWoundData2021/PublicationData/GenusLevelBatchCorrected_22.rda")

BatchCorrectedPhyloseq_Jae = phyloseq::subset_taxa(BatchCorrectedPhyloseq, Order!="Chloroplast")
BatchCorrectedPhyloseq_Jae = phyloseq::subset_taxa(BatchCorrectedPhyloseq_Jae, Family!="Mitochondria")
BatchCorrectedPhyloseq_Jae = phyloseq::subset_samples(BatchCorrectedPhyloseq_Jae, study_id !="Mock")

BatchCorrectedPhyloseq_Jae@sam_data$TotalOTUsLeftBatchMitoChloro = colSums(BatchCorrectedPhyloseq_Jae@otu_table@.Data)
BatchCorrectedPhyloseq_Jae@sam_data$TotalOTUsLeftBatchMitoChloro = colSums(BatchCorrectedPhyloseq_Jae@otu_table@.Data)

BatchCorrectedPhyloseq_Jae = phyloseq::subset_samples(BatchCorrectedPhyloseq_Jae, TotalOTUsLeftBatchMitoChloro>1200)
FinalSampleDist = data.frame(BatchCorrectedPhyloseq_Jae@sam_data)
FinalSampleDistJae= FinalSampleDist %>% select(Run, SampleID, study_id)

BatchCorrectedPhyloseq_Jae@sam_data = sample_data(FinalSampleDistJae)
sample_names(BatchCorrectedPhyloseq_Jae) = (BatchCorrectedPhyloseq_Jae@sam_data$study_id)

BatchCorrectedPhyloseq_Jae_Abdominal = subset_samples(BatchCorrectedPhyloseq_Jae, study_id %in% Abdominal_IDs$study_id)
save(BatchCorrectedPhyloseq_Jae_Abdominal, file="~/Documents/IowaWoundData2021/JaeAbdominal/JaeAbdominalWounds.rda")


################################################################
# Relative abundance of genera with strict anaerobes aggregated
################################################################

# This is the phyloseq object after aggregating strict anaerobes together and calculating % abundance in the otu_table 
# Object is called BatchCorrectedPhyloseqRankingPercent
load("~/Documents/IowaWoundData2021/PublicationData/RelativeAbundanceGenus.rda")
BatchCorrectedPhyloseqRankingPercent

# Genus-level abundance by patient
##################################
TaxaRows = data.frame(BatchCorrectedPhyloseqRankingPercent@tax_table@.Data)

# 'Anaerobe' variable was used for coding anaerobe y/n so we could aggregate; 
# in this data object all strict anaerobes have been collapsed into genus called 'Anaerobes'
TaxaRows$Anaerobe=NULL
TaxaRows$GenusAdjust = apply(TaxaRows, 1, Genus_code)

# relative abundance table where row names are the 'genus adjusted' 
OTUByID = data.frame(BatchCorrectedPhyloseqRankingPercent@otu_table@.Data)
row.names(OTUByID) = TaxaRows$GenusAdjust


Transposed = data.frame(t(OTUByID))
Transposed$SampleName = sapply(row.names(Transposed), function(x) str_remove(x, "X"))
TransposedJaeSamples = Transposed %>% filter(SampleName %in% WoundMicrobiomeDataAbdominal$SampleID)
TransposedJaeSamples$SampleName = NULL


# Whats the max. relative abundance of each Genus variable?
MaxRelAbundance = apply(TransposedJaeSamples, 2, max)

# Only genera found in at least .1% abundance in at least one of these 148 wounds
IncludedGenera = names(MaxRelAbundance[MaxRelAbundance>= 0.001])
IncludedGenera = sort(IncludedGenera)
TransposedJaeSamples = TransposedJaeSamples %>% select(IncludedGenera)
TransposedJaeSamples$SampleID = sapply(row.names(TransposedJaeSamples),function(x) str_remove(x, "X"))

TransposedJaeSamples = TransposedJaeSamples %>% left_join(WoundMicrobiomeData %>% select(SampleID, StudyID), by="SampleID")
TransposedJaeSamples = TransposedJaeSamples %>% select(StudyID, SampleID, IncludedGenera)
write.csv(TransposedJaeSamples, "~/Documents/IowaWoundData2021/JaeAbdominal/AllGeneraAbdominal.csv")

