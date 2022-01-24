# Amy Campbell
# Qiime2 output 11/2021

library("phyloseq")
library("dplyr")
library("stringr")
library("decontam")
library("ggplot2")
library("DESeq2")
library("vegan")
library("ggpubr")
library("RColorBrewer")
library("cowplot")
ampalette()
randpalette32= c("#E6F5C9", "#BF5B17", "#33A02C", "#B3E2CD", "#66C2A5",
                 "#66A61E", "#FFD92F", "#FFFF33", "#666666", "#6A3D9A", "#F2F2F2",
                 "#A6CEE3", "#FC8D62", "#1F78B4", "#7570B3", "#B3B3B3", "#E41A1C",
                 "#FDC086", "#FCCDE5", "#FFFFB3","#E6AB02", "#FF7F00","#E7298A",
                 "#B3DE69","#D95F02","#1B9E77", "#BC80BD", "#A6761D",
                 "#006400", "#0000B3","#681A1A", "#B300B3")


randpalette18=c("#B300B3","#E6AB02",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00",
                "#681A1A","#7570B3",
                "#1F78B4","#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#FFFF33",
                "#33A02C","#E6F5C9")

randpalette23 = c("#B300B3","#E6AB02", "#BF5B17",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00","#66C2A5","#66A61E",
                "#681A1A","#7570B3",
                "#1F78B4", "#FFD92F", "#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#FFFF33",
                "#33A02C","#D95F02", "#E6F5C9")


casPalette3 <- rev(c("#B85C00","#999999","#339966","#6B24B2","#56B4E9","#4D4D4D","#006600", "#CC0000", 
                     "#D119A3", "#F0E442", "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                     "#9900FF","#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                     "#F0E442", "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                     "#9900FF"))

# Set random seed to BRB zipcode :-]
set.seed(19104)

ampalette <- rev(c("#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                   "#4D4D4D", "#F0E442", "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                   "#9900FF","#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                   "#F0E442", "#4D4D4D",  "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                   "#9900FF"))
###########
# FUNCTIONS
###########
row_object <- function(taxa_id){
  # otu_taxon : of the subset of otu table columns containing just
  # otuID and taxonomy columns, a single row of that subset 
  
  defaultlist = c("NA", "NA","NA","NA", "NA", "NA", "NA")
  namelist = strsplit(toString(taxa_id), ';')
  for(i in 1:length(namelist[[1]])){
    taxaname = namelist[[1]][i]
    small_name_list = strsplit(taxaname, '__')
    mystring = small_name_list[[1]][2]
    mystring = str_remove(mystring, "\\[")
    mystring = str_remove(mystring, "\\]")
    mystring = str_remove(mystring, " ")
    defaultlist[i] = mystring
    
  }
  return(defaultlist)
}


Genus_code = function(tax_table_row){
  print(tax_table_row)
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

Species_code = function(tax_table_row){
  if(tax_table_row["Species"] !="NA"){
    return(tax_table_row["Species"])}
  else if(tax_table_row["Genus"] !="NA"){
    return(paste(tax_table_row["Genus"],tax_table_row["Species"], sep="_"))
  }else if(tax_table_row["Family"] != "NA"){
    return(paste(paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"),tax_table_row["Species"], sep="_"))
  }else if(tax_table_row["Order"] != "NA"){
    return(paste(paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"),tax_table_row["Species"], sep="_"))
  }else if(tax_table_row["Class"]!="NA"){
    return(paste(paste(tax_table_row["Class"],paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"), sep="_"),tax_table_row["Species"], sep="_"))
  }else if(tax_table_row["Phylum"] !="NA"){
    return(paste(paste(tax_table_row["Phylum"], paste(tax_table_row["Class"],paste(tax_table_row["Order"], paste(tax_table_row["Family"] , tax_table_row["Genus"], sep="_"), sep="_"), sep="_"), sep="_"),tax_table_row["Species"], sep="_"))
  }else if(tax_table_row["Phylum"]=="NA"){
    print("ouch")
    return("NULL")
  }else{
    return("NULL")
  }
}


# Returns P value of z test for differences in proportions 
Calculate_Proptest = function(row_32, row_35){
  successes32 = sum(row_32)
  successes35 = sum(row_35)
  
  N32 = length(row_32)
  N35 = length(row_35)
  
  testresult = prop.test(x=c(successes32, successes35), n=c(N32, N35))
  return(testresult$p.value)
}


# Read in data & metadata 
OTU_Table = read.csv2("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/table-with-taxonomy.tsv", header=T, sep="\t",skip=1)
runinfo = read.csv2("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/Control_Run_Info.tsv", sep=" ")
taxonomytable = read.csv2("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/taxonomy.tsv", sep="\t")
treephyseq = system.file("extdata","/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/tree.nwk",package="phyloseq")
clubgriceMetadata35 = read.csv2("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/mapping_files/run_maps/MiSeqV1V3_35.tsv", sep="\t")
twocolor = c("#FFC20A", "#0C7BDC")
patient_mapping32 = read.csv("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/IA_woundpain_mapping_32_2021.csv")
patient_mapping35 = read.csv("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/IA_woundpain_mapping_35_2021.csv")
PatientMetadata = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/GSWOUNDGRICE2015_20190221.csv")
# Process OTU table
#####################
runinfo$SampleID = sapply(runinfo$SampleID, function(x) toString(x))

OTU_Table[,2:ncol(OTU_Table)]= apply(OTU_Table[,2:ncol(OTU_Table)],2,function(x) as.numeric(as.character(x)))

totalReadsAssigned = colSums(OTU_Table[,2:ncol(OTU_Table)])

length(totalReadsAssigned[totalReadsAssigned<1200])

taxaIDs = OTU_Table$X.OTU.ID
taxadf = (data.frame(taxaIDs))
colnames(taxadf) = c("Feature.ID")
taxadf = taxadf %>% left_join(taxonomytable, by="Feature.ID")
taxadf$Taxon
taxfull = taxadf$Taxon


taxadf$Feature.ID
taxadf$Taxon

otu_IDs = data.frame(taxadf[c("Feature.ID", "Taxon")])


otu_IDs <- otu_IDs %>% group_by(Feature.ID) %>% transmute(Kingdom = (row_object(Taxon))[1],
                                                               Phylum = (row_object(Taxon))[2],
                                                               Class = (row_object(Taxon))[3],
                                                               Order = (row_object(Taxon))[4],
                                                               Family = (row_object(Taxon))[5],
                                                               Genus = (row_object(Taxon))[6],
                                                               Species = (row_object(Taxon))[7])

otu_IDs <- data.frame(otu_IDs)

# Now remove the column thats screwing up the OTU tax_glom steps
saveIDs = otu_IDs$Feature.ID
otu_IDs$Feature.ID = NULL

otu_just_ids_mat <- as.matrix(otu_IDs)
rownames(otu_just_ids_mat) = saveIDs

otu_for_phyloseq= OTU_Table
rownames(otu_for_phyloseq) = OTU_Table$X.OTU.ID
otu_for_phyloseq$X.OTU.ID=NULL
  

# Import to phyloseq
####################
OTU_tab = otu_table(otu_for_phyloseq, taxa_are_rows=TRUE)
sample_names(OTU_tab) = lapply(list(sample_names(OTU_tab)), function(x) str_replace(x, "X", ""))[[1]]
taxa = tax_table(otu_just_ids_mat)
treeobject = phyloseq::read_tree("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/tree.nwk")

row.names(runinfo) = runinfo$SampleID
sampledata= sample_data(runinfo)

# See if unifrac distances cluster by Run (before tax_glom)
###########################################################
PhyloseqObject = phyloseq(OTU_tab, taxa, sampledata, treeobject)
PhyloseqObject@sam_data$TotalOTUs = colSums(OTU_tab@.Data)

# Subset to just Bacteria first of all
PhyloseqObject = subset_taxa(PhyloseqObject, Kingdom=="Bacteria" )

PhyloseqObjectBackupInitial = PhyloseqObject
PhyloseqObjectForUniFrac = subset_samples(PhyloseqObject,TotalOTUs > 0 )


unifracdists <- phyloseq::distance(PhyloseqObjectForUniFrac, method="wunifrac")
unifracWeighted = ordinate(PhyloseqObjectForUniFrac,"PCoA", distance=unifracdists)
wUF = plot_ordination(PhyloseqObjectForUniFrac, unifracWeighted, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (ASV-level weighted UniFrac distance) in MiSeq Runs") + theme_classic()

Unweighted_unifracdists <- phyloseq::distance(PhyloseqObjectForUniFrac, method="unifrac")
unifracUnWeighted = ordinate(PhyloseqObjectForUniFrac,"PCoA", distance=Unweighted_unifracdists)
uwUF = plot_ordination(PhyloseqObjectForUniFrac, unifracUnWeighted, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (ASV-level unweighted UniFrac distance) in MiSeq Runs") + theme_classic()
gridExtra::grid.arrange(wUF, uwUF, ncol=2)

# Seems that beta divergence is significant between runs by unweighted, but not weighted, UniFrac
UniFracWeightedPreContamAdonis = adonis(unifracdists ~sample_data(PhyloseqObjectForUniFrac)$Run, permutations=9999)
UniFracUnWeightedPreContamAdonis = adonis(Unweighted_unifracdists ~sample_data(PhyloseqObjectForUniFrac)$Run, permutations=9999)

# See if unifrac distances cluster by Run (after tax_glom)
PhyloseqObjectForUniFrac_GenusGlom =  tax_glom(PhyloseqObjectForUniFrac, taxrank = "Genus")
unifracdistsGenus <- phyloseq::distance(PhyloseqObjectForUniFrac_GenusGlom, method="wunifrac")
unifracWeightedGenus = ordinate(PhyloseqObjectForUniFrac_GenusGlom,"PCoA", distance=unifracdistsGenus)
wUFgenus = plot_ordination(PhyloseqObjectForUniFrac_GenusGlom, unifracWeightedGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (Genus-level weighted UniFrac distance) in MiSeq Runs") + theme_classic()

Unweighted_unifracdistsGenus <- phyloseq::distance(PhyloseqObjectForUniFrac_GenusGlom, method="unifrac")
unifracUnWeightedGenus = ordinate(PhyloseqObjectForUniFrac_GenusGlom,"PCoA", distance=Unweighted_unifracdistsGenus)
uwUFGenus = plot_ordination(PhyloseqObjectForUniFrac_GenusGlom, unifracUnWeightedGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (Genus-level unweighted UniFrac distance) in MiSeq Runs") + theme_classic()
gridExtra::grid.arrange(wUFgenus, uwUFGenus, ncol=2)

# Once again, diff seems greater in unweighted than weighted unifrac -- maybe due to the rare species picked up by the additional reads in 32? 
UniFracWeightedPreContamGenusAdonis=adonis(unifracdistsGenus ~sample_data(PhyloseqObjectForUniFrac_GenusGlom)$Run, permutations=9999)
UniFracUnWeightedPreContamGenusAdonis=adonis(Unweighted_unifracdistsGenus ~sample_data(PhyloseqObjectForUniFrac_GenusGlom)$Run, permutations=9999)


PhyloseqObjectGenusGlom = tax_glom(PhyloseqObject, taxrank = "Genus")

RichnessPreContamASV  = plot_richness(PhyloseqObjectForUniFrac,x="Run", measures=c("Observed", "Shannon", "Simpson")) + ggtitle("Alpha Diversity (ASVs)") + geom_boxplot()
RichnessPreContamGenus =plot_richness(PhyloseqObjectGenusGlom,x="Run", measures=c("Observed", "Shannon", "Simpson")) + ggtitle("Alpha Diversity (Genera)")+ geom_boxplot()


##################
# Decontamination 
################# 
sort(totalReadsAssigned)

# Split up the two runs for decontamination purposes 
Phylo32 = subset_samples(PhyloseqObject, Run=="MiSeqV1V3_32")
Phylo35 = subset_samples(PhyloseqObject, Run=="MiSeqV1V3_35")

# (1) Identify contaminants to remove from MiSeqV1V3_32 based on decontam
#########################################################################

Controls32 = subset_samples(Phylo32, ControlStatus!="NonControl")
Controls32@sam_data$Sample_Ctrl = paste(Controls32@sam_data$SampleID, Controls32@sam_data$ControlStatus, sep="_")

# Look at the phyla present in each control (negatives + Mock community)
plot_bar(Controls32, "Sample_Ctrl", fill="Phylum") + scale_fill_manual(values=randpalette32) + ggtitle("OTUs of phyla in MiSeqV1V3_32 controls")

# Remove mock for just negative ctrls
Negatives32 = subset_samples(Controls32, SampleID!="120602")
NegativeIDs = (Negatives32@sam_data$SampleID)

# Decontam Prevalence based method on all negative controls (blanks + DNA-free water controls)
# (ID-ing OTUs/ASVs that are as prevalent or more prevalent in negative controls than in real samples)
Phylo32@sam_data$is.negative = if_else(Phylo32@sam_data$SampleID %in% NegativeIDs, TRUE, FALSE)
contaminants_prevalence32 = isContaminant(Phylo32, method="prevalence", neg="is.negative")


PrevContaminants32 = (contaminants_prevalence32 %>% filter(contaminant))
PrevContaminants32$ID = row.names(PrevContaminants32)
taxanames = row.names(taxa)
taxonomy_map_contams = data.frame(taxa)
taxonomy_map_contams$ID = row.names(taxonomy_map_contams)

contaminants_prevalence32 = PrevContaminants32 %>% left_join(taxonomy_map_contams,by="ID")

# (2) Identify ASVs in MiSeqV1V3_35 to remove based on decontam using Prevalence criteria
############################################################################################
# Decided 12/7 to remove DNAFreeWater5 because it's got an anomolously huge read count across the whole run (128059 assigned reads as compared to 50905, the next highest read count sample)
sort(totalReadsAssigned)

# Plot what's present in the negative controls (PCR controls) in MiSeqV1V3_35
Controls35 = subset_samples(Phylo35, ControlStatus!="NonControl")
Controls35@sam_data$Sample_Ctrl = paste(Controls35@sam_data$SampleID, Controls35@sam_data$ControlStatus, sep="_")
plot_bar(Controls35, "Sample_Ctrl", fill="Phylum") + scale_fill_manual(values=randpalette32) + ggtitle("OTUs of phyla in MiSeqV1V3_35 controls")



# ignore empty wells
Controls35 = subset_samples(Controls35, ControlStatus!="Empty")
waterphyla35 = plot_bar(Controls35, "SampleID", fill="Phylum") + scale_fill_manual(values=randpalette32) + ggtitle("Run 35 Water Control Phyla")

# Remove empty wells from Phylo35
Phylo35 = subset_samples(Phylo35,ControlStatus!="Empty" )

# Remove DNAFreewater5 from Phylo35 (see note above)
Phylo35 = subset_samples(Phylo35, SampleID!="DNAfreewater5")

Phylo35@sam_data$is.negative = if_else(Phylo35@sam_data$ControlStatus =="PCRControl", TRUE, FALSE)
contaminants_prevalence35 = isContaminant(Phylo35, method="prevalence", neg="is.negative")
contaminants_prevalence35$ID = row.names(contaminants_prevalence35)

contaminants_prevalence35taxa = contaminants_prevalence35 %>% left_join(taxonomy_map_contams,by="ID")
contaminantsOnly_prevalence35 = (contaminants_prevalence35taxa %>% filter(contaminant))

ToRemove35Prevalence=contaminantsOnly_prevalence35$ID

# (3) Identify ASVs in MiSeqV1V3_35  to remove based on Frequency (DNA concentration)
######################################################################################
 
# Incorporate concentration info from club grice mapping file 
sampledata35 = data.frame(Phylo35@sam_data)
clubgriceMetadata35$SampleID = clubgriceMetadata35$X.SampleID
clubgriceMetadata35 = clubgriceMetadata35 %>% select(SampleID, dna_concentration_ng_ul)
sampledata35 = sampledata35 %>% left_join(clubgriceMetadata35, by="SampleID")

sampledata35$dna_concentration_ng_ul = sapply(sampledata35$dna_concentration_ng_ul, function(x) as.numeric(as.character(x)))
Phylo35@sam_data = sample_data(sampledata35)

# Remove controls 
sample_names(Phylo35@sam_data) <- Phylo35@sam_data$SampleID
Phylo35 = subset_samples(Phylo35, ControlStatus=="NonControl")

contaminants_frequency35 = isContaminant(Phylo35, method="frequency", conc="dna_concentration_ng_ul")
contaminants_frequency35 = contaminants_frequency35 %>% filter(contaminant)
contaminants_frequency35$ID = row.names(contaminants_frequency35)

contaminants_frequency35taxa = contaminants_frequency35 %>% left_join(taxonomy_map_contams,by="ID")

ToRemove35Frequency=contaminants_frequency35taxa$ID

# (4) Remove contaminant ASVs from each run in Phyloseq Object
###############################################################

Remove35 = c(ToRemove35Frequency, ToRemove35Prevalence)
Remove32 = contaminants_prevalence32$ID

allremoves = c(Remove35, Remove32)

Run32Samples= Phylo32@sam_data$SampleID
Run35Samples=Phylo35@sam_data$SampleID

KeepTaxa = row.names(PhyloseqObject@tax_table)
KeepTaxa = setdiff(KeepTaxa, Remove32)
KeepTaxa = setdiff(KeepTaxa, Remove35)
PhyloseqObjectBackup = PhyloseqObject

removed = prune_taxa(allremoves, PhyloseqObjectBackup)
 
PhyloseqObjectBackup1 = PhyloseqObject
PhyloseqObject = prune_taxa(KeepTaxa, PhyloseqObject)

# Remove negative controls 
PhyloseqObject = subset_samples(PhyloseqObject, ControlStatus %in% c("NonControl", "Mock"))

# Filter to samples with at least 1200 assigned OTUs removes 22 samples 
# This filter removes 1 Venous ulcer, 17 surgical wounds, 2 traumatic wounds, 2 mixed traumatic + surgical 
PhyloseqObject = subset_samples(PhyloseqObject, TotalOTUs >= 1200)

# # Unifrac distances after decontamination, filtering to samples with 1200 reads and above (ASV-level)
# ######################################################################################################
# Unweighted_unifracdistsPostDecontam =  phyloseq::distance(PhyloseqObject, method="unifrac")
# unifracUnWeightedFilteredPostDecontam= ordinate(PhyloseqObject,"PCoA", distance=Unweighted_unifracdistsPostDecontam)
# uwUFFilteredPlotDecontam =  plot_ordination(PhyloseqObject, unifracUnWeightedFilteredPostDecontam, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates ( OTU-level unweighted UniFrac distance) in MiSeq Runs After Decontamination") + theme_classic()
# 
# unifracdistsFilteredOTUPostDecontam <- phyloseq::distance(PhyloseqObject, method="wunifrac")
# unifracWeightedFilteredOTUPostDecontam = ordinate(PhyloseqObject,"PCoA", distance=unifracdistsFilteredOTUPostDecontam)
# PlotunifracWeightedFilteredOTUPostDecontam = plot_ordination(PhyloseqObject, unifracWeightedFilteredOTUPostDecontam, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (OTU-level weighted UniFrac distance) in MiSeq Runs After Decontamination") + theme_classic()
# 
# # Unifrac distances after decontamination, filtering to samples with 1200 reads and above (Genus-level)
# #######################################################################################################
# PhyloseqObjectGenus = tax_glom(PhyloseqObject, taxrank = "Genus", bad_empty=c(NA, "NA"))
# 
# Unweighted_unifracdistsPostDecontamGenus =  phyloseq::distance(PhyloseqObjectGenus, method="unifrac")
# unifracUnWeightedFilteredPostDecontamGenus= ordinate(PhyloseqObjectGenus,"PCoA", distance=Unweighted_unifracdistsPostDecontamGenus)
# uwUFFilteredPlotDecontamGenus =  plot_ordination(PhyloseqObjectGenus, unifracUnWeightedFilteredPostDecontamGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates ( Genus-level unweighted UniFrac distance) in MiSeq Runs After Decontamination") + theme_classic()
# 
# unifracdistsFilteredOTUPostDecontamGenus <- phyloseq::distance(PhyloseqObjectGenus, method="wunifrac")
# unifracWeightedFilteredOTUPostDecontamGenus = ordinate(PhyloseqObjectGenus,"PCoA", distance=unifracdistsFilteredOTUPostDecontamGenus)
# PlotunifracWeightedFilteredPostDecontamGenus = plot_ordination(PhyloseqObjectGenus, unifracWeightedFilteredOTUPostDecontamGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (Genus-level weighted UniFrac distance) in MiSeq Runs After Decontamination") + theme_classic()
# 
# # PERMANOVAS 
# # OTU level
# OTULevelDF = data.frame(PhyloseqObject@sam_data)
# OTU_UnweightedAdonisPostDecontam = adonis(Unweighted_unifracdistsPostDecontam ~ Run, data=OTULevelDF, permutations=9999)
# OTU_WeightedAdonisPostDecontam = adonis(unifracdistsFilteredOTUPostDecontam ~ Run, data=OTULevelDF, permutations=9999)
# 
# #Genus level
# GenusLevelDF = data.frame(PhyloseqObjectGenus@sam_data)
# GenusUnweightedAdonisPostDecontam = adonis(Unweighted_unifracdistsPostDecontamGenus ~ Run,data=GenusLevelDF, permutations=9999)
# GenusWeightedAdonisPostDecontam = adonis(unifracdistsFilteredOTUPostDecontamGenus ~ Run, data=GenusLevelDF, permutations=9999)
# 
# 
# # Adonis permanova p-values for plots
# ######################
# PvalueWeightedPreContamASV = toString(UniFracWeightedPreContamAdonis$aov.tab$`Pr(>F)`[1])
# PvalueUnweightedPreContamASV = toString(UniFracUnWeightedPreContamAdonis$aov.tab$`Pr(>F)`[1])
# PvalueWeightedPreContamGenus = toString(UniFracWeightedPreContamGenusAdonis$aov.tab$`Pr(>F)`[1])
# PvalueUnweightedPreContamGenus = toString(UniFracUnWeightedPreContamGenusAdonis$aov.tab$`Pr(>F)`[1])
# 
# PvalueWeightedPostASV =  toString(OTU_WeightedAdonisPostDecontam$aov.tab$`Pr(>F)`[1])
# PvalueUnweightedPostASV =  toString(OTU_UnweightedAdonisPostDecontam$aov.tab$`Pr(>F)`[1])
# PvalueWeightedPostGenus = toString(GenusWeightedAdonisPostDecontam$aov.tab$`Pr(>F)`[1])
# PvalueUnweightedPostGenus = toString(GenusUnweightedAdonisPostDecontam$aov.tab$`Pr(>F)`[1])
# 
# # Pre-decontam plots for saving
# #############################
# wUF = wUF + annotate(geom="text", x=-.09, y=.01, label=paste0("permanova_P=", PvalueWeightedPreContamASV))
# uwUF = uwUF + annotate(geom="text", x=-.4, y=.2, label=paste0("permanova_P=", PvalueUnweightedPreContamASV))
# wUFgenus = wUFgenus + annotate(geom="text", x=-.6, y=.2, label=paste0("permanova_P=", PvalueWeightedPreContamGenus))
# uwUFGenus  = uwUFGenus + annotate(geom="text", x=-.22, y=.27, label=paste0("permanova_P=", PvalueUnweightedPreContamGenus))
# 
# # Post-decontam plots for saving
# #############################
# PlotunifracWeightedFilteredOTUPostDecontam = PlotunifracWeightedFilteredOTUPostDecontam + annotate(geom="text", x=-.08, y=.03, label=paste0("permanova_P=", PvalueWeightedPostASV))
# uwUFFilteredPlotDecontam = uwUFFilteredPlotDecontam + annotate(geom="text", x=-.4, y=.25, label=paste0("permanova_P=", PvalueUnweightedPostASV))
# PlotunifracWeightedFilteredPostDecontamGenus = PlotunifracWeightedFilteredPostDecontamGenus + annotate(geom="text", x=-.6, y=.22, label=paste0("permanova_P=", PvalueWeightedPostGenus))
# uwUFFilteredPlotDecontamGenus = uwUFFilteredPlotDecontamGenus + annotate(geom="text", x=-.3, y=.25, label=paste0("permanova_P=", PvalueUnweightedPostGenus))

# Save plots
##############
arrangedPreContam = gridExtra::grid.arrange(wUF, uwUF, wUFgenus, uwUFGenus)
ggsave(arrangedPreContam, file="/Users/amycampbell/Documents/IowaWoundData2021/PreDecontamUnifracDistances.pdf", height=10, width=18)

arrangedPostDecontam = gridExtra::grid.arrange(PlotunifracWeightedFilteredOTUPostDecontam, uwUFFilteredPlotDecontam, PlotunifracWeightedFilteredPostDecontamGenus, uwUFFilteredPlotDecontamGenus)
ggsave(arrangedPostDecontam, file="/Users/amycampbell/Documents/IowaWoundData2021/PostDecontamUnifracDistances.pdf", height=10, width=18)



richnessPostDecontamGenus = plot_richness(PhyloseqObjectGenus,x="Run", measures=c("Observed", "Shannon", "Simpson")) + ggtitle("Alpha Diversity (Genera) Post-Decontamination" ) + geom_boxplot()
richnessPostDecontamASV = plot_richness(PhyloseqObject,x="Run", measures=c("Observed", "Shannon", "Simpson")) + ggtitle("Alpha Diversity (ASVs) Post-Decontamination")+ geom_boxplot()
gridExtra::grid.arrange(richnessPostDecontamASV, richnessPostDecontamGenus)

# Plot mock community
#######################
# Good time to grab the mock community and look at its genus-level composition
MockComm = subset_samples(PhyloseqObject, ControlStatus=="Mock")
MockCommGenus = tax_glom(MockComm, taxrank="Genus")
MockCommGenusFilter = MockCommGenus
MockCommGenusFilter= filter_taxa(MockCommGenus, function(x) mean(x) > 0, TRUE)
GenusPlotMock = plot_bar(MockCommGenusFilter, x="SampleID", y="Abundance", fill="Genus") + scale_fill_manual(values=rev(ampalette)) + ggtitle("Mock Community Genera")
ggsave(GenusPlotMock, file="/Users/amycampbell/Documents/IowaWoundData2021/MockCommunityGenera.png")

# Add patient metadata
######################################

# Sample data minus mock community
# OTULevelDF is just data.frame(PhyloseqObject@sam_data) above
OTULevelDF = OTULevelDF %>% filter(ControlStatus!="Mock")

# Mappings from Sample to patient
patient_mapping32_select = patient_mapping32 %>% select(X.SampleID, SubjectID)
patient_mapping32_select$SampleID = patient_mapping32_select$X.SampleID

patient_mapping35_select = patient_mapping35 %>% select(X.SampleID, SubjectID)
patient_mapping35_select$SampleID = patient_mapping35_select$X.SampleID

patient_mapping = rbind(patient_mapping32_select, patient_mapping35_select)
OTULevelDF = OTULevelDF %>% left_join(patient_mapping, by="SampleID")

colnames(OTULevelDF) = c("SampleID", "ControlStatus", "Run", "TotalOTUs", "X.SampleID","study_id")

# Grab just a few 
PatientMetadata_WoundType = PatientMetadata %>% select(wound_type, woundloc, study_id)
PatientMetadata_WoundType$study_id = factor(PatientMetadata_WoundType$study_id)
OTULevelDF = OTULevelDF %>% left_join(PatientMetadata_WoundType, by="study_id")
row.names(OTULevelDF) = OTULevelDF$SampleID

# note about metadata:
# Wound types
# # 1 = Pressure ulcer
# 2 = venous ulcer
# 4 = surgical
# 5 = traumatic 
# 6 = Mixed (Traumatic + surgical)
# 7 = Other

# Map meaningful strings onto numerical codes for wound type and location 
OTULevelDF = OTULevelDF %>% mutate(wound_type_string = case_when(wound_type==1 ~ "Pressure", 
                                                                 wound_type==2 ~"Venous",
                                                                 wound_type==4 ~"Surgical", 
                                                                 wound_type==5 ~ "Traumatic", 
                                                                 wound_type==6 ~ "Mixed Traumatic/Surgical",
                                                                 wound_type==7 ~ "Other")
  
)

OTULevelDF = OTULevelDF %>% mutate(woundloc_string = case_when(woundloc==1 ~ "Extremity", 
                                                                 woundloc==2 ~"Trunk",
                                                                 woundloc==3 ~"Head/Neck", 
                                                                 woundloc==4 ~ "Inguinal"
                                                              )
                                   
)
OTULevelDF$wound_type_string = factor(OTULevelDF$wound_type_string)

OTULevelDF$woundloc_string = factor(OTULevelDF$woundloc_string)

woundtypeDistRun = ggplot(OTULevelDF, aes(x=Run, fill=wound_type_string)) + geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal()  + ggtitle("Wound Types ")


woundlocDistRun = ggplot(OTULevelDF, aes(x=Run, fill=woundloc_string)) + geom_bar() + scale_fill_manual(values=rev(randpalette32)) + theme_minimal() + ggtitle("Wound Locations")
Distributions = gridExtra::grid.arrange(woundtypeDistRun, woundlocDistRun, ncol=2)
ggsave(Distributions, file="/Users/amycampbell/Documents/IowaWoundData2021/WoundLocTypeByRun.png")
PatientMetadata$study_id = factor(PatientMetadata$study_id)
patient_mapping$study_id = patient_mapping$SubjectID
patient_mapping = patient_mapping %>% left_join(PatientMetadata, by="study_id")


sampledataphyloseqobj = data.frame(PhyloseqObject@sam_data)
sampledataphyloseqobj = sampledataphyloseqobj %>% left_join(patient_mapping, by="SampleID")
PhyloseqObject@sam_data = sample_data(sampledataphyloseqobj)
sample_names(PhyloseqObject) = sampledataphyloseqobj$SampleID


# Average together OTU counts for 2 duplicated patients 
########################################################
PhyloseqObject@sam_data$study_id = factor(PhyloseqObject@sam_data$study_id)

PhyloseqObject@sam_data$study_id = sapply(PhyloseqObject@sam_data$study_id, function(x) toString(x))

# Get just those patients' samples 
PhyloseqObject_2728 = subset_samples(PhyloseqObject, study_id %in% c("27", "28"))
PhyloseqObjectMerged = merge_samples(PhyloseqObject_2728, "study_id")

PhyloseqObjectMerged@sam_data$SampleID = c("120361_120383","120362_120384")
PhyloseqObjectMerged@sam_data$X.SampleID = c("120361_120383","120362_120384")
sample_names(PhyloseqObjectMerged) = PhyloseqObjectMerged@sam_data$SampleID


PhyloseqObjectMerged@sam_data$SubjectID = factor(PhyloseqObjectMerged@sam_data$study_id)
PhyloseqObjectMerged@sam_data$Run = "MiSeqV1V3_32"
PhyloseqObjectMerged@sam_data$ControlStatus="NonControl" 
PhyloseqObjectMerged@sam_data$othwtype="" 

# everything except these two samples 
PhyloseqObjectUpdated = subset_samples(PhyloseqObject, !(study_id %in% c("27", "28")))

# now add them back in 
PhyloseqObjectUpdated = merge_phyloseq(PhyloseqObjectUpdated, PhyloseqObjectMerged)

save(PhyloseqObjectUpdated, file="PhyloSeqObjectPostDecontamMerged.rda")

# Phylum level abundance 
########################
PhyloseqObjectUpdated =  subset_taxa(PhyloseqObjectUpdated, Phylum!="NA")

df_phyla_prenorm = PhyloseqObjectUpdated %>% 
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance >.01) %>%
  group_by(Phylum)


dataframe_counts = data.frame(df_phyla_prenorm)
sampnames = (unique(dataframe_counts$Sample))

dataframe_counts = dataframe_counts %>% filter(Phylum=="Firmicutes")

# samples with 0 Firmicutes (or less than 1%; weird)
zero = setdiff(sampnames,dataframe_counts$Sample)
zero_32 =zero[!grepl("IowaWound", zero)]
zero_35 =zero[grepl("IowaWound", zero)]

dataframe_counts_32 = dataframe_counts %>% filter(Run=="MiSeqV1V3_32") %>% arrange(Abundance)
samples32 = append(zero_32, dataframe_counts_32$Sample)

dataframe_counts_35 = dataframe_counts %>% filter(Run=="MiSeqV1V3_35") %>% arrange(Abundance)
samples35 = append(zero_35, dataframe_counts_35$Sample)
Desired_order= append(samples32, samples35)
View(df_phyla)

plotPhylaPreNorm= ggplot(df_phyla_prenorm, aes(x=SampleID, y=Abundance, fill=Phylum)) + geom_bar(stat="identity") + ggtitle("Phylum-level composition by run (Pre-Normalization)") + scale_fill_manual(values=(randpalette18))+ theme_minimal() + theme(axis.text.x=element_text(angle=90))

plotPhylaPreNorm$data$SampleID = factor(plotPhylaPreNorm$data$Sample, levels =Desired_order)

plotPhylaPreNorm + theme(axis.text.x = element_blank())

# 
# # NA strings:
# NAstrings = c("unidentified_marine",
# "uncultured_soil",
# "unidentified_eubacterium",
# "uncultured_organism", NA,
# "NA", "uncultured_microorganism",
# "uncultured_bacterium", "uncultured_candidate", "uncultured_compost",
# "unidentified", "uncultured_prokaryote", "uncultured_rumen", "wastewater_metagenome")

# taxes = data.frame(PhyloseqObjectUpdated@tax_table@.Data)
# sort(unique(taxes$Species))
# Specieslevel = tax_glom(PhyloseqObjectUpdated, taxrank= "Species", bad_empty = NAstrings)
# taxesspecies = data.frame(Specieslevel@tax_table@.Data)
# write.csv(taxesspecies, file="SpeciesLevelTaxTable.csv")
# 

# ###############################################
# # Testing out effects of variance stabilization 
# ###############################################
# 
# # for purposes of VST the factor doesn't do much  (confirmed this empirically by )
# deseqobj_fromphylo = phyloseq_to_deseq2(PhyloseqObjectUpdated, ~Run)
# 
# # Calculate geometric means without running into divide by zero errors
# # Via Phyloseq vignette (https://rdrr.io/bioc/phyloseq/src/inst/doc/phyloseq-mixture-models.R)
# get_geometric_means <- function(x){
#   mean = exp(sum(log(x[x>0]), na.rm=TRUE) / length(x))
#   return(mean)
# }
# gmeans = apply(counts(deseqobj_fromphylo), 1, get_geometric_means)
# estimates_sizes = estimateSizeFactors(deseqobj_fromphylo, geoMeans = gmeans)
# disp_estimates = estimateDispersions(estimates_sizes)
# 
# v_stabilized_counts= getVarianceStabilizedData(disp_estimates)
# v_stabilized_countsForZeroing = v_stabilized_counts
# v_stabilized_countsForZeroing[which(v_stabilized_countsForZeroing < 0)] <- 0
# 
# PhyloseqObjectUpdatedNormalized = PhyloseqObjectUpdated
# PhyloseqObjectUpdatedNormalized@otu_table = otu_table(v_stabilized_countsForZeroing, taxa_are_rows = T)
# df_phyla_postnorm = PhyloseqObjectUpdatedNormalized %>% 
#   tax_glom(taxrank = "Phylum") %>%
#   transform_sample_counts(function(x) {x/sum(x)}) %>%
#   psmelt() %>% 
#   filter(Abundance >.01) %>%
#   group_by(Phylum)
# 
# plotPhylaPostNorm= ggplot(df_phyla_postnorm, aes(x=SampleID, y=Abundance, fill=Phylum)) + geom_bar(stat="identity") + ggtitle("Phylum-level composition by run (Post-Normalization)") + scale_fill_manual(values=(randpalette23))+ theme_minimal() + theme(axis.text.x=element_text(angle=90))
# 
# plotPhylaPostNorm$data$SampleID = factor(plotPhylaPostNorm$data$Sample, levels =Desired_order)
# 
# extractlegend=cowplot::get_legend(plotPhylaPostNorm)
# plotPhylaPostNorm + theme(axis.text.x = element_blank())
# 
# Pre_Post = gridExtra::grid.arrange(plotPhylaPreNorm +theme(axis.text.x = element_blank(), legend.position = "none") ,extractlegend, plotPhylaPostNorm + theme(axis.text.x = element_blank(), legend.position = "none"), nrow=2, widths=c(.8,.1))
# ggsave(Pre_Post, file="/Users/amycampbell/Documents/IowaWoundData2021/PrePostNorm.pdf", height=12, width=25)

#################################################
# Testing batch effect correction
#################################################

SpeciesLevelGlom = PhyloseqObjectUpdated %>% tax_glom(taxrank="Species")# %>% transform_sample_counts(function(x) {x/sum(x)})
GenusLevelGlom = PhyloseqObjectUpdated %>% tax_glom(taxrank="Genus") #
# 2a41defcfcab7155dcb50a0c646ab68b
SpeciesLevelGlomPct = SpeciesLevelGlom %>% transform_sample_counts(function(x) {x/sum(x)})
GenusLevelGlomPct = GenusLevelGlom %>% transform_sample_counts(function(x) {x/sum(x)})

# Genus-level batch correction 
GenusLevelBatchCorrect = data.frame(GenusLevelGlomPct@otu_table@.Data)

# Genera with <.1% abundance are considered missing
GenusLevelBatchCorrect[GenusLevelBatchCorrect < .001] <- 0
# Genera with >.1% abundance are considered present 
GenusLevelBatchCorrect[GenusLevelBatchCorrect >= .001] <- 1

#Divide into the two runs 
GenusLevelDF32 = GenusLevelBatchCorrect %>% select(!contains("IowaWound.Human."))
GenusLevelDF35 = GenusLevelBatchCorrect %>% select(contains("IowaWound.Human."))

dataframe_for_recordingGenus= data.frame(otuID = row.names(GenusLevelDF32))
PvectorGenus= sapply(1:nrow(dataframe_for_recordingGenus), function(x) Calculate_Proptest(GenusLevelDF32[x,], GenusLevelDF35[x, ]))
dataframe_for_recordingGenus$Pvector = PvectorGenus
dataframe_for_recordingGenus$adjPvector = PvectorGenus*nrow(dataframe_for_recordingGenus)

rm_Genera = (dataframe_for_recordingGenus %>% filter(adjPvector< .05))$otuID

taxagenusGlom = data.frame(GenusLevelGlom@tax_table@.Data)
taxagenusGlom$otuID = row.names(taxagenusGlom)

# Sneathia, Mycoplasma, Neisseria, and Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium to be removed
taxagenusGlom %>% filter(otuID %in% rm_Genera)

KeepTaxaTestGenusBatch = setdiff(row.names(GenusLevelGlom@tax_table@.Data), rm_Genera)
TestUnifracPostBatchEffectGenus = prune_taxa(KeepTaxaTestGenusBatch, GenusLevelGlom)

# Since 120560 only had Mycoplasma after decontam, it's now an empty sample and has to be removed 
TestUnifracPostBatchEffectGenus@sam_data$TotalOTUsLeftBatch = colSums(TestUnifracPostBatchEffectGenus@otu_table@.Data)
TestUnifracPostBatchEffectGenus = subset_samples(TestUnifracPostBatchEffectGenus, TotalOTUsLeftBatch>1200)

UW_testUnifrac_PostBatchGenus =  phyloseq::distance(TestUnifracPostBatchEffectGenus, method="unifrac")
W_testUnifrac_PostBatchGenus =   phyloseq::distance(TestUnifracPostBatchEffectGenus, method="wunifrac")

UW_ordination_PostbatchGenus = ordinate(TestUnifracPostBatchEffectGenus,"PCoA", distance=UW_testUnifrac_PostBatchGenus)
W_ordination_PostbatchGenus = ordinate(TestUnifracPostBatchEffectGenus,"PCoA", distance=W_testUnifrac_PostBatchGenus)

plot_UW_postbatchGenus = plot_ordination(TestUnifracPostBatchEffectGenus, UW_ordination_PostbatchGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates ( Genus-level unweighted UniFrac distance) in MiSeq Runs After Batch Correction") + theme_classic()
plot_W_postbatchGenus = plot_ordination(TestUnifracPostBatchEffectGenus, W_ordination_PostbatchGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates ( Genus-level weighted UniFrac distance) in MiSeq Runs After Batch Correction") + theme_classic()

PostBatchDFGenusTest = data.frame(TestUnifracPostBatchEffectGenus@sam_data)

PostBatch_UnweightedAdonisGenusTest = adonis(UW_testUnifrac_PostBatchGenus ~ Run, data=PostBatchDFGenusTest, permutations=9999)
PostBatch_WeightedAdonisGenusTest = adonis(W_testUnifrac_PostBatchGenus ~ Run, data=PostBatchDFGenusTest, permutations=9999)











# Test DMM on this
TestUnifracPostBatchEffectGenus
Counts = TestUnifracPostBatchEffectGenus@otu_table


# Filter to species prevalent in at least 10 samples 
PrevalenceInfo = data.frame(Counts)
PrevalenceInfo$Prevalence = rowSums(PrevalenceInfo != 0)

sort(PrevalenceInfo$Prevalence)

PrevalenceInfoFilter = PrevalenceInfo %>% filter(Prevalence >= 10)
TaxaKeepDMM = row.names(PrevalenceInfoFilter)

#
PhyloDMM = prune_taxa(TaxaKeepDMM, TestUnifracPostBatchEffectGenus)

CountsDMM = t(PhyloDMM@otu_table@.Data)

CountsDMMinput = data.matrix(CountsDMM)

densityplot(log10(colSums(CountsDMMinput)), xlim=range(log10(colSums(CountsDMMinput))),xlab="Taxon representation (log 10 count)")
densityplot((colSums(CountsDMMinput)), xlim=range((colSums(CountsDMMinput))),xlab="Taxon representation (raw count)")
library(DirichletMultinomial)


# Fit dirichlet multinomial mixture models with 1,2,...,20 components and plot laplace approximation 
FitDirichlets = mclapply(1:12, dmn, count=CountsDMMinput, verbose=TRUE, seed=19104)
lplc <- sapply(FitDirichlets, laplace)
plot(lplc, type = 'b', xlab = 'Dirichlet Components',ylab='Model Fit', main="Dirichlet Components by Laplace Model Fit (Genus Read Counts)") 

# get dirichlet fit object with maximum laplace-estimated likelihood (minimum negative)
bestFit <- FitDirichlets[[which.min(lplc)]]
















OTULevelGlom = PhyloseqObjectUpdated %>% transform_sample_counts(function(x) {x/sum(x)})

OTULevelBatchCorrect = data.frame(OTULevelGlom@otu_table@.Data)
OTULevelBatchCorrect[OTULevelBatchCorrect < .0001] <- 0
OTULevelBatchCorrect[OTULevelBatchCorrect >= .0001] <- 1

OTULevelDF32 = OTULevelBatchCorrect %>% select(!contains("IowaWound.Human."))
OTULevelDF35 = OTULevelBatchCorrect %>% select(contains("IowaWound.Human."))

dataframe_for_recordingOTUs = data.frame(otuID = row.names(OTULevelDF32))
PvectorOTU = sapply(1:nrow(dataframe_for_recordingOTUs), function(x) Calculate_Proptest(OTULevelDF32[x,], OTULevelDF35[x, ]))
dataframe_for_recordingOTUs$Pvector = PvectorOTU
dataframe_for_recordingOTUs$PvectorAdj = dataframe_for_recordingOTUs$Pvector * nrow(dataframe_for_recordingOTUs)
rm_OTUs = (dataframe_for_recordingOTUs %>% filter(PvectorAdj< .01))$otuID
OTULevelGlom@tax_table@.Data[rm_OTUs,]


KeepTaxaTestOTUBatch = setdiff(row.names(PhyloseqObjectUpdated@otu_table@.Data), rm_OTUs)
PhyloseqObjectUpdatedTestRemove = prune_taxa(KeepTaxaTestOTUBatch, PhyloseqObjectUpdated)


# Test Unifrac with Batch Effect Correction (OTU level )
TestUnifracPostBatchEffect = PhyloseqObjectUpdatedTestRemove
UW_testUnifrac_PostBatch =  phyloseq::distance(TestUnifracPostBatchEffect, method="unifrac")
W_testUnifrac_PostBatch <- phyloseq::distance(TestUnifracPostBatchEffect, method="wunifrac")


UW_ordination_Postbatch = ordinate(TestUnifracPostBatchEffect,"PCoA", distance=UW_testUnifrac_PostBatch)
W_ordination_Postbatch = ordinate(TestUnifracPostBatchEffect,"PCoA", distance=W_testUnifrac_PostBatch)

plot_UW_postbatch = plot_ordination(TestUnifracPostBatchEffect, UW_ordination_Postbatch, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates ( OTU-level unweighted UniFrac distance) in MiSeq Runs After Batch Correction") + theme_classic()
plot_W_postbatch = plot_ordination(TestUnifracPostBatchEffect, W_ordination_Postbatch, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates ( OTU-level weighted UniFrac distance) in MiSeq Runs After Batch Correction") + theme_classic()


p_val_UW = PostBatch_UnweightedAdonis$aov.tab$`Pr(>F)`
p_val_W = PostBatch_WeightedAdonis$aov.tab$`Pr(>F)`

plot_UW_postbatch = plot_UW_postbatch + annotate("text", x=-.1, y=.2, label=paste0("PermanovaP=", toString(p_val_UW[1])))
plot_W_postbatch = plot_W_postbatch + annotate("text", x=-.06, y=.04, label=paste0("PermanovaP=", toString(p_val_W[1])))

PostBatchDF = data.frame(TestUnifracPostBatchEffect@sam_data)

PostBatch_UnweightedAdonis = adonis(UW_testUnifrac_PostBatch ~ Run, data=PostBatchDF, permutations=9999)
PostBatch_WeightedAdonis = adonis(W_testUnifrac_PostBatch ~ Run, data=PostBatchDF, permutations=9999)

# Test Unifrac with OTU-level Batch effect (Genus level Phylo)

GenusPostBatchEffect = PhyloseqObjectUpdatedTestRemove %>% tax_glom(taxrank = "Genus")
UW_testUnifrac_PostBatchGenus =  phyloseq::distance(GenusPostBatchEffect, method="unifrac")
W_testUnifrac_PostBatchGenus <- phyloseq::distance(GenusPostBatchEffect, method="wunifrac")


UW_ordination_PostbatchGenus = ordinate(GenusPostBatchEffect,"PCoA", distance=UW_testUnifrac_PostBatchGenus)
W_ordination_PostbatchGenus = ordinate(GenusPostBatchEffect,"PCoA", distance=W_testUnifrac_PostBatchGenus)

plot_UW_postbatchGenus = plot_ordination(GenusPostBatchEffect, UW_ordination_PostbatchGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates ( Genus-level unweighted UniFrac distance) in MiSeq Runs After Batch Correction") + theme_classic()
plot_W_postbatchGenus = plot_ordination(GenusPostBatchEffect, W_ordination_PostbatchGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates ( Genus-level weighted UniFrac distance) in MiSeq Runs After Batch Correction") + theme_classic()

PostBatchDFGenus = data.frame(GenusPostBatchEffect@sam_data)

PostBatch_UnweightedAdonisGenus = adonis(UW_testUnifrac_PostBatchGenus ~ Run, data=PostBatchDFGenus, permutations=9999)
PostBatch_WeightedAdonisGenus = adonis(W_testUnifrac_PostBatchGenus ~ Run, data=PostBatchDFGenus, permutations=9999)
PostBatch_UnweightedAdonisGenus
PostBatch_WeightedAdonisGenus


p_val_UW = PostBatch_UnweightedAdonis$aov.tab$`Pr(>F)`
p_val_W = PostBatch_WeightedAdonis$aov.tab$`Pr(>F)`

plot_UW_postbatch = plot_UW_postbatch + annotate("text", x=-.1, y=.2, label=paste0("PermanovaP=", toString(p_val_UW[1])))
plot_W_postbatch = plot_W_postbatch + annotate("text", x=-.06, y=.04, label=paste0("PermanovaP=", toString(p_val_W[1])))





