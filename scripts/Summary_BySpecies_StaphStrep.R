# Amy Campbell
# "Heterogeneous wound microbiome varies..."
# Supplementary figure 6

# background:
# Reviewer 3 suggested blasting ASVs to try to get species-level classifications
# we've been cautious not to infer below genus-level information about these samples based on the V1V3 16S region,
# but here, are summarizing the overall species composition (where possible) of Strep. and Staph. in this project.

library(dplyr)
library(stringr)
library(ggplot2)

ampalette <- rev(c("#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                   "#4D4D4D", "#F0E442", "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                   "#9900FF","#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                   "#F0E442", "#4D4D4D",  "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                   "#9900FF"))
# 1. Read in all Strep-, Staph-annotated ASVs' top 5 BLASTN results to the 16S ribosomal RNA database from NCBI 
##################################################################################################################

# wound microbiome data (update path as needed)
woundmicrobiome = read.csv("/Users/amycampbell/Downloads/WoundMicrobiomeDataForSEG_AEC.csv")
studyIDsIncluded = (woundmicrobiome %>% filter(!is.na(DMMClusterAssign)))$StudyID
BlastResults = read.csv2("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/test_blast/StaphStrepOnly_Seqs/staph_strep_blast_withSubjName", sep="\t", header=F)
colnames(BlastResults) = c("ASV", "sseqid", "pident", "length","mismatch","gapopen", "qcovs","evalue","bitscore","sscinames")

# get rid of 'subspecies' or strain ID
BlastResults$SciNameSimplified = sapply(BlastResults$sscinames, function(x) paste(str_split(x, " ")[[1]][1], str_split(x, " ")[[1]][2], sep="_") )

BlastResults = BlastResults %>% filter(sseqid!="gi|2221405513|ref|NR_175559.1|")
# Iterate through the 846 unique ASVs
####################################
ASVs=unique(BlastResults$ASV)
AssignedAutomatically = data.frame()
NeedsFollowup = data.frame()

FollowupStaphs = data.frame()
for(ASVID in ASVs){
  
  SubsetFrame = BlastResults %>% filter(ASV==ASVID) %>% arrange(-bitscore)
  if("Staphylococcus_aureus" %in% SubsetFrame$SciNameSimplified){
    FollowupStaphs = rbind(FollowupStaphs, SubsetFrame)
  }
  if(length(unique(SubsetFrame$SciNameSimplified)) > 1){
    # Look for all the entries for which theres >97% identity and 99% query coverage 
    AtLeast97_99 = SubsetFrame %>% filter(pident > 97.0 & qcovs>99.0)
    if(length(unique(AtLeast97_99$SciNameSimplified))>1){
      NeedsFollowup = rbind(NeedsFollowup, SubsetFrame )
      
    }else{
      
      if(dim(AtLeast97_99)[1]==0){
        Assignment="NA"
        ASV=SubsetFrame[1,"ASV"]
        SubSeq=SubsetFrame[1,"sseqid"]
        Pident=SubsetFrame[1,"pident"]
        BitScore=SubsetFrame[1,"bitscore"]
        QueryCov=SubsetFrame[1,"qcovs"]
        
        AssignedAutomatically = rbind(AssignedAutomatically, c(ASV, Assignment, SubSeq,Pident,BitScore, QueryCov))
        
      }else{
        Assignment=AtLeast97_99[1,"SciNameSimplified"]
        ASV=AtLeast97_99[1,"ASV"]
        SubSeq=AtLeast97_99[1,"sseqid"]
        Pident=AtLeast97_99[1,"pident"]
        BitScore=AtLeast97_99[1,"bitscore"]
        QueryCov=AtLeast97_99[1,"qcovs"]
        AssignedAutomatically = rbind(AssignedAutomatically, c(ASV, Assignment, SubSeq,Pident,BitScore, QueryCov))
        
      }
         
    }

    
  }else{
    if(max(SubsetFrame$pident) > 97.0 & max(SubsetFrame$qcovs)>99.0){
      Assignment=(SubsetFrame$SciNameSimplified)[1]
     
    }else{
      Assignment="NA"
    }
  
    ASV=SubsetFrame[1,"ASV"]
    SubSeq=SubsetFrame[1,"sseqid"]
    Pident=SubsetFrame[1,"pident"]
    BitScore=SubsetFrame[1,"bitscore"]
    QueryCov=SubsetFrame[1,"qcovs"]
    
    AssignedAutomatically = rbind(AssignedAutomatically, c(ASV, Assignment, SubSeq,Pident,BitScore, QueryCov))
    
    
  }
}


# Second pass
# Of those which weren't immediately classified...
# iterate through
for(asv in unique(NeedsFollowup$ASV)){
  
  SubsetFrame = NeedsFollowup %>% filter(ASV==asv) #%>% arrange(-pident)
  SubsetFrame$pident = sapply(SubsetFrame$pident, function(x) as.numeric(as.character(x)))
  # Make a list of entries with >99% coverage the max pct identity observed
  SubsetFrame = SubsetFrame %>% filter(qcovs>99 & pident>97 )
  MaxPctID = max(SubsetFrame$pident)
  
  MaxPctRows = SubsetFrame %>% filter(pident==MaxPctID)
  
  # If theres more than one species with the same top % identity 
  # call this NA

  if(length(unique(MaxPctRows$SciNameSimplified)) > 1){
    Assignment="NA"
    ASV=SubsetFrame[1,"ASV"]
    SubSeq=SubsetFrame[1,"sseqid"]
    Pident=SubsetFrame[1,"pident"]
    BitScore=SubsetFrame[1,"bitscore"]
    QueryCov=SubsetFrame[1,"qcovs"]
    
    AssignedAutomatically = rbind(AssignedAutomatically, c(ASV, Assignment, SubSeq,Pident,BitScore, QueryCov))
    
  }else{
    # as long as this DF is not empty, assign to the top Pct Identity 
    if(dim(MaxPctRows)[1]!=0){
      Assignment=MaxPctRows[1,"SciNameSimplified"]
      ASV=SubsetFrame[1,"ASV"]
      SubSeq=SubsetFrame[1,"sseqid"]
      Pident=SubsetFrame[1,"pident"]
      BitScore=SubsetFrame[1,"bitscore"]
      QueryCov=SubsetFrame[1,"qcovs"]
      AssignedAutomatically = rbind(AssignedAutomatically, c(ASV, Assignment, SubSeq,Pident,BitScore, QueryCov))
      
      }else{
        Assignment="NA"
        ASV=SubsetFrame[1,"ASV"]
        SubSeq=SubsetFrame[1,"sseqid"]
        Pident=SubsetFrame[1,"pident"]
        BitScore=SubsetFrame[1,"bitscore"]
        QueryCov=SubsetFrame[1,"qcovs"]
        
        AssignedAutomatically = rbind(AssignedAutomatically, c(ASV, Assignment, SubSeq,Pident,BitScore, QueryCov))
        
      
    }
    
    }
  
}
colnames(AssignedAutomatically) = c("ASV", "SpeciesClassificationBLAST", "BlastNCBI_ID","PctIdentity", "Bitscore", "PctCoverage")
write.csv(AssignedAutomatically, file="~/Documents/IowaWoundData2021/Qiime2Data/test_blast/BlastAssignmentsStaphStrep.csv")


# OTU table (post-decontam but before genus-level batch correction)
load("/Users/amycampbell/Documents/IowaWoundData2021/PhyloSeqObjectPostDecontamMerged_recent.rda")

StaphStrep = subset_taxa(PhyloseqObjectUpdated,Genus %in% c("Staphylococcus", "Streptococcus"))
taxtable = data.frame(StaphStrep@tax_table) 
taxtable$OTU = row.names(taxtable)
OTUassignsStrepStaph = AssignedAutomatically %>% select(ASV, SpeciesClassificationBLAST)
taxtable$ASV = row.names(taxtable)

taxtable = taxtable %>% left_join(OTUassignsStrepStaph, by="ASV")
row.names(taxtable) = taxtable$ASV

AsTaxTable =  tax_table(taxtable)
taxa_names(AsTaxTable) = taxtable$ASV
colnames(AsTaxTable) = colnames(taxtable)
tax_table(StaphStrep)= AsTaxTable

JustStaph = subset_taxa(StaphStrep,Genus=="Staphylococcus")
RelativeAbundStaph = JustStaph %>%transform_sample_counts(function(x) {x/sum(x)})
plot_bar(RelativeAbundStaph, fill="SpeciesClassificationBLAST") + scale_fill_manual(values=ampalette)


ampaletteStrep= c(ampalette, c("#000080", "#800000", "black"))


bigpalette=c("#556B2F","#696969","#8B0000", "#808000", "#483D8B", "#008000", "#3CB371", "#008080", "#000080", "#CD5C5C", "#32CD32",
             "#DAA520", "#800080", "#FF4500","#00CED1", "#FF8C00", "#00FF00", "#00FA9A", "#dc143c", "#00BFFF", "#A020F0", "#ADFF2F", "#FF7f50",
             "#FF00FF", "#F0E68C","#FFFF54","#6495ED", "#DDA0DD", "#ADD8E6", "#7B68EE", "#7FFFD4","#FF69B4", "#FFE4C4", "#FFB6C1","#002D04")


JustStrep = subset_taxa(StaphStrep,Genus=="Streptococcus")
RelativeAbundStrep= JustStrep %>%transform_sample_counts(function(x) {x/sum(x)})
pdf("~/Documents/IowaWoundData2021/NewFigs_Paper_10_23/StrepPctAbundance.pdf", width=20, height=8)
plot_bar(RelativeAbundStrep, fill="SpeciesClassificationBLAST") + scale_fill_manual(values=bigpalette)
dev.off()

JustStrepMelted = RelativeAbundStrep %>% psmelt()
JustStaphpMelted = RelativeAbundStaph %>% psmelt()


pdf("~/Documents/IowaWoundData2021/NewFigs_Paper_10_23/StaphPctAbundance.pdf", width=20, height=8)
plot_bar(RelativeAbundStaph, fill="SpeciesClassificationBLAST") + scale_fill_manual(values=bigpalette)
dev.off()



load("~/Downloads/StaphStrepASVassigned.rda")


taxatable_staphstrep = data.frame(tax_table(StaphStrep))
taxatable_staphstrep$OTU=NULL
taxatable_staphstrep$ASV = NULL




View(taxatable_staphstrep %>% select(Species,SpeciesClassificationBLAST ))

taxatable_staphstrep$SpeciesClassificationBLAST

taxatable_staphstrep = taxatable_staphstrep %>% mutate(SpeciesAdjusted = case_when( ((Species=="NA" | Species=="uncultured_organism") & (SpeciesClassificationBLAST!="NA")) ~ SpeciesClassificationBLAST,
                                                             (!(Species=="NA" | Species=="uncultured_organism") & (SpeciesClassificationBLAST=="NA")) ~ Species,
                                                             TRUE ~ as.character(SpeciesClassificationBLAST)))
taxatable_staphstrep$Species=NULL
taxatable_staphstrep$SpeciesClassificationBLAST=NULL                      
                                

AsTaxTableAdjusted = tax_table(taxatable_staphstrep)
taxa_names(AsTaxTableAdjusted) = row.names(taxatable_staphstrep)

colnames(AsTaxTableAdjusted) = colnames(taxatable_staphstrep)

tax_table(StaphStrep) = AsTaxTableAdjusted

StaphStrepAdjusted = StaphStrep %>% tax_glom(taxrank="SpeciesAdjusted")


StaphStrepAdjusted = subset_samples(StaphStrepAdjusted, study_id %in% studyIDsIncluded)


StaphOnly = subset_taxa(StaphStrepAdjusted, Genus=="Staphylococcus")
StrepOnly = subset_taxa(StaphStrepAdjusted, Genus=="Streptococcus")


StaphOnly_Relabund = StaphOnly %>%transform_sample_counts(function(x) {x/sum(x)})


MaximumProportionStaph = apply(OTUs_Staph, 1, max)

# 23 species of Staph assigned with at least 5% abundance in at least one sample
AtLeast5 = MaximumProportionStaph[MaximumProportionStaph>0.05]
StaphOnly_Relabund_AtLeast5 = StaphOnly_Relabund %>% subset_taxa( rownames(tax_table(StaphOnly_Relabund)) %in% names(AtLeast5))

MeltedStaphAtLeast5Pct = StaphOnly_Relabund_AtLeast5 %>% psmelt()

OTUTabStaph = (data.frame(StaphOnly_Relabund_AtLeast5@otu_table))
t_OTUTabStaph = data.frame(t(OTUTabStaph))
t_OTUTabStaph[is.na(t_OTUTabStaph)] <- 0
SampleIDsOrderSepi_NA_aureus = (t_OTUTabStaph %>% arrange(X28791a195fb16e432f328fc38c217a16, X483181f68a1997d8d691fc99eb4d8919, X53fb12ba5d166f18c30e101ad3fefe25))$sampleID
SampleIDsOrderSepi_NA_aureus = sapply(SampleIDsOrderSepi_NA_aureus, function(x) str_remove(x, "X"))

t_OTUTabStaph$sampleID = colnames(OTUTabStaph)


OrderStaphEpi_NA = unique((MeltedStaphAtLeast5Pct %>% arrange(wound_type, study_id) )$study_id)

Plotstaph5pct = ggplot(MeltedStaphAtLeast5Pct, aes(x=SampleID,y=Abundance, fill=SpeciesAdjusted)) + geom_bar(stat="identity") + scale_fill_manual(values=bigpalette)+
  theme_classic() + theme(axis.text.x=element_blank())

Plotstaph5pct$data$SampleID = factor(Plotstaph5pct$data$SampleID, levels=SampleIDsOrderSepi_NA_aureus)


OTUs_Staph = data.frame(StaphOnly_Relabund@otu_table)
OTUs_Staph[is.na(OTUs_Staph)]<-0

MeanAbund = rowMeans(OTUs_Staph)

MostAbundantAvg = names(rev(sort(MeanAbund))[1:10])

OTUs_Staph_top10Abundant = subset_taxa(StaphOnly_Relabund,  rownames(tax_table(StaphOnly_Relabund)) %in% MostAbundantAvg)





StrepOnly_Relabund = StrepOnly %>%transform_sample_counts(function(x) {x/sum(x)})
OTUs_Strep = data.frame(StrepOnly_Relabund@otu_table)
OTUs_Strep[is.na(OTUs_Strep)]<-0

# strep salivarius, strep thermophilus, NA most abundant 

# 7136faeb92b8894cfb389206ac8828f9, beb20ef494d6232af88af5f79dea8650, dc4c06e76c80792a6dd6a5170e96a392   
# StrepOnly_Relabund@tax_table[rev(c("dc4c06e76c80792a6dd6a5170e96a392","beb20ef494d6232af88af5f79dea8650","7136faeb92b8894cfb389206ac8828f9")),]


MeanAbundStrep = sort(rowMeans(OTUs_Strep))
MaximumProportionStrep = apply(OTUs_Strep, 1, max)

# 36 taxa assigned with at least 5% abundance in at least one sample
AtLeast5 = MaximumProportionStrep[MaximumProportionStrep>0.05]
StrepOnly_Relabund_AtLeast5 = StrepOnly_Relabund %>% subset_taxa( rownames(tax_table(StrepOnly_Relabund)) %in% names(AtLeast5))

MeltedStrepAtLeast5Pct = StrepOnly_Relabund_AtLeast5 %>% psmelt()

OTUTabStrep = (data.frame(StrepOnly_Relabund_AtLeast5@otu_table))
t_OTUTabStrep = data.frame(t(OTUTabStrep))
t_OTUTabStrep[is.na(t_OTUTabStrep)] <- 0
SampleIDsOrderMostAbundant3 = (t_OTUTabStrep %>% arrange(X7136faeb92b8894cfb389206ac8828f9, beb20ef494d6232af88af5f79dea8650, dc4c06e76c80792a6dd6a5170e96a392))$sampleID
SampleIDsOrderMostAbundant3 = sapply(SampleIDsOrderMostAbundant3, function(x) str_remove(x, "X"))

t_OTUTabStrep$sampleID = colnames(OTUTabStrep)

bigpalette2 = bigpalette=c("#696969", "#556B2F","#8B0000", "#808000", "#483D8B", "#008000", "#3CB371", "#008080", "#000080", "#CD5C5C", "#32CD32",
                           "#DAA520", "#800080", "#FF4500","#00CED1", "#FF8C00", "#00FF00", "#00FA9A", "#dc143c", "#00BFFF", "#A020F0", "#ADFF2F", "#FF7f50",
                           "#FF00FF", "#F0E68C","#FFFF54","#6495ED", "#DDA0DD", "#ADD8E6", "#7B68EE", "#7FFFD4","#FF69B4", "#FFE4C4", "#FFB6C1","#004600", "black")

OrderStrepEpi_NA = unique((MeltedStrepAtLeast5Pct %>% arrange(wound_type, study_id) )$study_id)

PlotStrep5pct = ggplot(MeltedStrepAtLeast5Pct, aes(x=SampleID,y=Abundance, fill=SpeciesAdjusted)) + geom_bar(stat="identity") + scale_fill_manual(values=bigpalette2)+
  theme_classic() + theme(axis.text.x=element_blank())

PlotStrep5pct$data$SampleID = factor(PlotStrep5pct$data$SampleID, levels=SampleIDsOrderMostAbundant3)

pdf("~/Documents/Staph_Strep.pdf", width=20, height=20)
gridExtra::grid.arrange(PlotStrep5pct, Plotstaph5pct, ncol=1)
dev.off()
