# Amy Campbell
# Output of Qiime2 
# July 2021

library("phyloseq")
library("dplyr")
library("stringr")
library("decontam")
library("ggplot2")
library("DESeq2")
library("vegan")
library("ggpubr")
setwd("~/Desktop/GriceLabGit/IowaWound/")
ampalette <- rev(c("#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                   "#4D4D4D", "#F0E442", "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                   "#9900FF","#B85C00","#999999","#339966","#6B24B2","#56B4E9","#D119A3","#006600", "#CC0000", 
                   "#F0E442", "#4D4D4D",  "#CC99FF","#663300","#33CC33", "#0072B2", "#FF9900", 
                   "#9900FF"))
twocolor = c("#FFC20A", "#0C7BDC")
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
  if(tax_table_row[1,c("Genus")] !="NA"){
    return(tax_table_row[1,c("Genus")])
  }else if(tax_table_row[1,c("Family")] != "NA"){
    return(paste(tax_table_row[1,c("Family")] , tax_table_row[1, c("Genus")], sep="_"))
  }else if(tax_table_row[1, c("Order")] != "NA"){
    return(paste(tax_table_row[1, c("Order")], paste(tax_table_row[1,c("Family")] , tax_table_row[1, c("Genus")], sep="_"), sep="_"))
  }else if(tax_table_row[1, c("Class")]!="NA"){
    return(paste(tax_table_row[1, c("Class")],paste(tax_table_row[1, c("Order")], paste(tax_table_row[1,c("Family")] , tax_table_row[1, c("Genus")], sep="_"), sep="_"), sep="_"))
  }else if(tax_table_row[1, c("Phylum")]!="NA"){
    return(paste(tax_table_row[1, c("Phylum")], paste(tax_table_row[1, c("Class")],paste(tax_table_row[1, c("Order")], paste(tax_table_row[1,c("Family")] , tax_table_row[1, c("Genus")], sep="_"), sep="_"), sep="_"), sep="_"))
  }else if(tax_table_row[1, c("Phylum")]=="NA"){
    print("ouch")
    return("NULL")
  }else{
    return("NULL")
  }
}

# Read in data & metadata 
OTU_Table = read.csv2("data/table.from_biom_w_taxonomy.txt", header=T, sep="\t",skip=1)
controls_run = read.csv2("mappings/Control_Run_Info.tsv", header=T, sep=" ")
patient_mapping32 = read.csv("mappings/IA_woundpain_mapping_32_2021.csv")
patient_mapping35 = read.csv("mappings/IA_woundpain_mapping_35_2021.csv")
patient_metadata = read.csv("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/figuring_out_metadata_5_21/GSWOUNDGRICE2015_20190221.csv")

patient_mapping32$SampleID = as.factor(patient_mapping32$X.SampleID)
patient_mapping35$SampleID = as.factor(patient_mapping35$X.SampleID)

denoise32 = read.csv2("data/DenoiseStats32.tsv", sep="\t")
denoise32 = denoise32[-1, ]

denoise35 = read.csv2("data/DenoiseStats35.tsv", sep="\t")
denoise35 = denoise35[-1, ]
 
taxonomy_column = OTU_Table$taxonomy 
taxonomy_map = OTU_Table %>% select(X.OTU.ID, taxonomy)

taxonomy_map <- taxonomy_map %>% group_by(X.OTU.ID) %>% transmute(Kingdom = (row_object(taxonomy))[1],
                                                               Phylum = (row_object(taxonomy))[2],
                                                               Class = (row_object(taxonomy))[3],
                                                               Order = (row_object(taxonomy))[4],
                                                               Family = (row_object(taxonomy))[5],
                                                               Genus = (row_object(taxonomy))[6],
                                                               Species = (row_object(taxonomy))[7])

taxonomy_map <- data.frame(taxonomy_map)

# Make taxonomy matrix for phyloseq
####################################

save_OTUs = taxonomy_map$X.OTU.ID
taxonomy_map$X.OTU.ID = NULL 
taxonomy_matrix = as.matrix(taxonomy_map)
row.names(taxonomy_matrix) = save_OTUs
taxonomy_matrix = tax_table(taxonomy_matrix)

OTU_Table$taxonomy =NULL 
save_rows =  OTU_Table$X.OTU.ID
OTU_Table$X.OTU.ID = NULL
OTU_Table = mapply(function(x) as.numeric(as.character(x)),OTU_Table)
rownames(OTU_Table) = save_rows
OTU_tab = otu_table(OTU_Table, taxa_are_rows=TRUE)

sample_names(OTU_tab) = lapply(list(sample_names(OTU_tab)), function(x) str_replace(x, "X", ""))[[1]]

# Make sample data objects for input into phyloseq
###################################################
Samples32 = controls_run %>% filter(Run == "MiSeqV1V3_32")
Samples35 = controls_run %>% filter(Run == "MiSeqV1V3_35")

sampledata32 = sample_data(Samples32)
sample_names(sampledata32) <- Samples32$SampleID

sampledata35 = sample_data(Samples35)
sample_names(sampledata35) <- Samples35$SampleID

sampledataboth = sample_data(controls_run)
sample_names(sampledataboth) = controls_run$SampleID

# Phyloseq objects from runs 32 and 35
######################################
phylo32 = phyloseq(OTU_tab, taxonomy_matrix, sampledata32)
phylo35 = phyloseq(OTU_tab, taxonomy_matrix, sampledata35)
phylofull = phyloseq(OTU_tab, taxonomy_matrix, sampledataboth)

# Filter all cyanobacteria
##########################
phylo32 = subset_taxa(phylo32,Phylum !="Cyanobacteria")
phylo35 = subset_taxa(phylo35,Phylum !="Cyanobacteria")
phylofull = subset_taxa(phylofull,Phylum !="Cyanobacteria")

# Filter all non-bacteria
phylo32 = subset_taxa(phylo32, Kingdom == "Bacteria")
phylo35 = subset_taxa(phylo35, Kingdom == "Bacteria")
phylofull = subset_taxa(phylofull, Kingdom == "Bacteria")

# Look at non-chimeric reads in each sample, characterized & uncharacterized
#############################################################################
denoise35 = denoise35 %>% select(sample.id, filtered, non.chimeric, input)
denoise32 = denoise32 %>% select(sample.id, filtered, non.chimeric, input)

# Look at proportions of 
phylo32@sam_data$Reads = colSums(phylo32@otu_table@.Data)
phylo35@sam_data$Reads = colSums(phylo35@otu_table@.Data)

samdata32 = data.frame(phylo32@sam_data) %>% filter(ControlStatus=="NonControl")
samdata35 = data.frame(phylo35@sam_data) %>% filter(ControlStatus=="NonControl")

colnames(denoise32) = c("SampleID", "Filtered", "NonChimeric","input")
colnames(denoise35) = c("SampleID", "Filtered", "NonChimeric","input")

samdata32 = samdata32 %>% left_join(denoise32, by="SampleID")
samdata35 = samdata35 %>% left_join(denoise35, by="SampleID")

samdata32$Filtered = sapply(samdata32$Filtered, function(x) as.numeric(as.character(x)))
samdata35$Filtered = sapply(samdata35$Filtered, function(x) as.numeric(as.character(x)))

samdata32$NonChimeric = sapply(samdata32$NonChimeric, function(x) as.numeric(as.character(x)))
samdata35$NonChimeric = sapply(samdata35$NonChimeric, function(x) as.numeric(as.character(x)))

samdata32$Chimeric = samdata32$Filtered - samdata32$NonChimeric
samdata35$Chimeric = samdata35$Filtered - samdata35$NonChimeric

patient_metadata$SubjectID = as.factor(patient_metadata$study_id)
samdata32 = samdata32 %>% left_join(patient_mapping32, by="SampleID")
samdata32 = samdata32 %>% left_join(patient_metadata, by = "SubjectID")

samdata35 = samdata35 %>% left_join(patient_mapping35, by="SampleID")
samdata35 = samdata35 %>% left_join(patient_metadata, by = "SubjectID")

# Plots of proportion chimera
chimeras32melt = samdata32 %>% select(SampleID, Chimeric, NonChimeric, wound_type) %>% reshape2::melt(id.vars=c("SampleID", "wound_type"))
chimeras32 = ggplot(chimeras32melt, aes(x=SampleID, y=value, fill=variable))+ geom_bar(aes(fill=variable), stat="identity")
chimeras32_order = chimeras32melt %>% filter(variable=="NonChimeric") %>% arrange(wound_type, value)
chimeras32$data$SampleID = factor(chimeras32$data$SampleID, levels=chimeras32_order$SampleID)
chimeras32 = chimeras32 + theme_classic() + scale_fill_manual(values=twocolor) +theme(axis.text.x=element_blank()) + xlab("") + ggtitle("Run 32 Chimeric Reads") + ylim(c(0, 170000))


chimeras35melt = samdata35 %>% select(SampleID, Chimeric, NonChimeric, wound_type) %>% reshape2::melt(id.vars=c("SampleID", "wound_type"))
chimeras35 = ggplot(chimeras35melt, aes(x=SampleID, y=value, fill=variable))+ geom_bar(aes(fill=variable), stat="identity")
chimeras35_order = chimeras35melt %>% filter(variable=="NonChimeric") %>% arrange(wound_type, value)
chimeras35$data$SampleID = factor(chimeras35$data$SampleID, levels=chimeras35_order$SampleID)
chimeras35 = chimeras35 + theme_classic() + scale_fill_manual(values=twocolor) +theme(axis.text.x=element_blank()) + xlab("") + ylim(c(0, 170000))+ ggtitle("Run 35 Chimeric Reads")
ggsave(gridExtra::grid.arrange(chimeras32, chimeras35), file="ChimericReads.pdf")

# Plot of proportion reads mapped 
samdata32$Uncharacterized = samdata32$NonChimeric - samdata32$Reads
samdata35$Uncharacterized = samdata35$NonChimeric - samdata35$Reads
samdata32$Characterized =  samdata32$Reads
samdata35$Characterized =  samdata35$Reads

characterized32 = samdata32 %>% select(SampleID, Uncharacterized, Characterized, wound_type) %>% reshape2::melt(id.vars=c("SampleID", "wound_type"))
characterized35 = samdata35 %>% select(SampleID, Uncharacterized, Characterized, wound_type) %>% reshape2::melt(id.vars=c("SampleID", "wound_type"))

characterizedplot32 = ggplot(characterized32, aes(fill=variable,y=value,x=SampleID)) + geom_bar(stat="identity") + theme_classic() + scale_fill_manual(values=twocolor) +theme(axis.text.x=element_blank()) + xlab("") + ylim(c(0, 170000))+ ggtitle("Run 32 Characterized Non-chimeric Reads")
characterizedplot35 = ggplot(characterized35, aes(fill=variable, y=value, x=SampleID)) + geom_bar(stat="identity") + theme_classic() + scale_fill_manual(values=twocolor) +theme(axis.text.x=element_blank()) + xlab("") + ylim(c(0, 170000))+ ggtitle("Run 35 Characterized Non-chimeric Reads")

characterizedplot32Order = (characterized32 %>% filter(variable=="Characterized") %>% arrange(wound_type, value))$SampleID
characterizedplot32$data$SampleID = factor(characterizedplot32$data$SampleID, levels=characterizedplot32Order)

characterizedplot35Order = (characterized35 %>% filter(variable=="Characterized") %>% arrange(wound_type, value))$SampleID
characterizedplot35$data$SampleID = factor(characterizedplot35$data$SampleID, levels=characterizedplot35Order)

ggsave(gridExtra::grid.arrange(characterizedplot32, characterizedplot35), file="UncharacterizedReads.pdf")



# Actual filtration step 
# Filtered out
samdata32$input = sapply(samdata32$input, function(x) as.numeric(as.character(x)))
samdata32$Removed = samdata32$input - samdata32$Filtered
samdata35$input = sapply(samdata35$input, function(x) as.numeric(as.character(x)))
samdata35$Removed = samdata35$input - samdata35$Filtered
samdata32$Kept = samdata32$Filtered
samdata35$Kept = samdata35$Filtered

Filtered32 = samdata32 %>% select(SampleID, Kept, Removed, wound_type) %>% reshape2::melt(id.vars=c("SampleID", "wound_type"))
Filtered35 = samdata35 %>% select(SampleID, Kept, Removed, wound_type) %>% reshape2::melt(id.vars=c("SampleID", "wound_type"))


filterplot32 =  ggplot(Filtered32, aes(fill=variable,y=value,x=SampleID)) + geom_bar(stat="identity") + theme_classic() + scale_fill_manual(values=twocolor) +theme(axis.text.x=element_blank()) + xlab("") + ggtitle("Run 32 Raw Read Filtering/Denoising") + ylim(c(0,250000))
filterplot32$data$variable = factor(filterplot32$data$variable,levels=c("Removed", "Kept"))

filterplot35 =  ggplot(Filtered35, aes(fill=variable,y=value,x=SampleID)) + geom_bar(stat="identity") + theme_classic() + scale_fill_manual(values=twocolor) +theme(axis.text.x=element_blank()) + xlab("") + ggtitle("Run 35 Raw Read Filtering/Denoising")+ ylim(c(0,250000))
filterplot32$data$variable = factor(filterplot32$data$variable,levels=c("Removed", "Kept"))
filterplot35$data$variable = factor(filterplot35$data$variable,levels=c("Removed", "Kept"))

filterplot32_Order = (filterplot32$data %>% filter(variable=="Kept") %>% arrange(wound_type, value))$SampleID
filterplot35_Order = (filterplot35$data %>% filter(variable=="Kept") %>% arrange(wound_type, value))$SampleID

filterplot32$data$SampleID = factor(filterplot32$data$SampleID, levels=filterplot32_Order)
filterplot35$data$SampleID = factor(filterplot35$data$SampleID, levels=filterplot35_Order)

ggsave(gridExtra::grid.arrange(filterplot32, filterplot35, ncol=2), file="ReadFiltering.pdf", height=7, width=15)

# Compare the final pre-decontam reads together by run and wound type 
bigframe = rbind(samdata32, samdata35)
bigframe = bigframe %>% mutate(WoundDescription = case_when(wound_type == 1 ~ "PressureUlcer",
                                                            wound_type==2 ~ "VenousUlcer",
                                                            wound_type==3 ~"ArterialUlcer",
                                                            wound_type==4 ~ "Surgical", 
                                                            wound_type==5 ~ "Traumatic",
                                                            wound_type==6 ~ "MixedTraumaticSurgical",
                                                            TRUE ~ "Other"))
ggplot(bigframe, aes(x=WoundDescription, y=Reads, group=WoundDescription)) + geom_boxplot() + facet_grid(.~Run)



# Filter non-control samples to those with minimum 1000 reads prior to decontam
###############################################################################

phylo32_original = sample_names(phylo32)
phylo35_original = sample_names(phylo35)

phylo32@sam_data$Reads <- colSums(phylo32@otu_table@.Data)
phylo35@sam_data$Reads <- colSums(phylo35@otu_table@.Data)


# Before any decontamination
#############################
phylofull@sam_data$Reads = colSums(phylofull@otu_table@.Data)
phylofull = subset_samples(phylofull, ( Reads >1000))

phylofullOTUdf = data.frame(rowSums(phylofull@otu_table@.Data))
colnames(phylofullOTUdf) = c("Sum")
phylofullOTUdfKeep = phylofullOTUdf %>% filter(Sum > 10)
phylofull = prune_taxa(row.names(phylofullOTUdfKeep), phylofull)
phylofull_ordination = ordinate(phylofull, "PCoA", "bray")
colSums(phylofull@otu_table@.Data)
phylofullPlot= plot_ordination(phylofull, phylofull_ordination, color="Run") + 
  ggtitle("Bray Curtis-based PCoA Coordinates 1 and 2 (Before Decontamination)") + scale_color_manual(values=twocolor) + geom_point(size=3) 

# Adding genus code 
pre_Decontam_tax = phylofull@tax_table
pre_Decontam_tax[is.na(pre_Decontam_tax)] <- "NA"
tax_table(phylofull) <- pre_Decontam_tax
PreDecontamGenusCode= sapply(rownames(phylofull@tax_table), function(x) Genus_code(phylofull@tax_table[x, ]))
physeq_tax_preDecontam = data.frame(phylofull@tax_table)
physeq_tax_preDecontam$GenusCode = PreDecontamGenusCode
phylofull@tax_table = tax_table(physeq_tax_preDecontam)
colnames(phylofull@tax_table) = colnames(physeq_tax_preDecontam)
rownames(phylofull@tax_table) = rownames(physeq_tax_preDecontam)
phylofullGenusCode = phylofull %>% tax_glom(taxrank="GenusCode")
phylofull_ordination_GenusCode = ordinate(phylofullGenusCode, "PCoA", "bray")
phylofullPlotGenuscode= plot_ordination(phylofullGenusCode, phylofull_ordination_GenusCode, color="Run") + 
  ggtitle("Bray Curtis-based PCoA Coordinates 1 and 2 (Before Decontamination; Aggregated to Genus)") + scale_color_manual(values=twocolor) + geom_point(size=3) 


# 217 - 208 = 9 removed from run 32
# 233 - 155 = 78 removed from run 35
# 363 positive samples left 

phylo32 <- subset_samples(phylo32, !(ControlStatus=="NonControl" & Reads < 1000))
phylo35 <- subset_samples(phylo35,  !(ControlStatus=="NonControl" & Reads < 1000))

DepthFiltered32 = setdiff(phylo32_original, sample_names(phylo32))
DepthFiltered35 = setdiff(phylo35_original, sample_names(phylo35))

DepthFilteredPatients32 = patient_mapping32 %>% filter(X.SampleID %in% DepthFiltered32)
DepthFilteredPatientMetadata32 = patient_metadata %>% filter(study_id %in% DepthFilteredPatients32$SubjectID)
# one traumatic wound (5), one mixed (6), 7 surgical (4) in 32 
# 64 surgical (4), 

DepthFilteredPatients35 = patient_mapping35 %>% filter(X.SampleID %in% DepthFiltered35)
DepthFilteredPatientMetadata35 = patient_metadata %>% filter(study_id %in% DepthFilteredPatients35$SubjectID)

patient_metadata=patient_metadata %>% mutate(Removed = if_else(study_id %in% DepthFilteredPatients35$SubjectID | study_id %in% DepthFilteredPatients32$SubjectID,"Remove", "Keep"))
patient_metadata = patient_metadata %>% mutate(WoundDescription = case_when(wound_type == 1 ~ "PressureUlcer",
                                                         wound_type==2 ~ "VenousUlcer",
                                                         wound_type==3 ~"ArterialUlcer",
                                                         wound_type==4 ~ "Surgical", 
                                                         wound_type==5 ~ "Traumatic",
                                                         wound_type==6 ~ "MixedTraumaticSurgical",
                                                         TRUE ~ "Other"))

depthremoved1 = ggplot(patient_metadata, aes(Removed)) + aes(x=Removed, fill=WoundDescription) + geom_bar() + scale_fill_manual(values=rev(ampalette)) + xlab("") + ylab("Count total") + ggtitle("Depth-removed samples by\ntotal count")

depthremoved2 = ggplot(patient_metadata, aes(Removed)) + aes(x=Removed, fill=WoundDescription) + geom_bar(position="fill") + scale_fill_manual(values=rev(ampalette)) + xlab("") + ylab("Proportion total") + ggtitle("Depth-removed samples by proportion")

ggsave(gridExtra::grid.arrange(depthremoved1, depthremoved2, ncol=2), file="depth_removed.pdf", width=11, height=7)

################################################################
# Process decontamination for MiSeqV1V3_32 --> phylo32_filtered
################################################################

# (1) Identify contaminants to remove from MiSeqV1V3_32 based on decontam
#########################################################################
Controls32 = subset_samples(phylo32, ControlStatus!="NonControl")
Controls32@sam_data$Sample_Ctrl = paste(Controls32@sam_data$SampleID, Controls32@sam_data$ControlStatus, sep="_")
plot_bar(Controls32, "Sample_Ctrl", fill="Phylum") + scale_fill_manual(values=ampalette) + ggtitle("OTUs of phyla in MiSeqV1V3_32 controls")

# Remove mock for just negative ctrls
Negatives32 = subset_samples(Controls32, SampleID!="120602")
NegativeIDs = (Negatives32@sam_data$SampleID)

# Decontam Prevalence based method on all negative controls (blanks + DNA-free water controls)
phylo32@sam_data$is.negative = if_else(phylo32@sam_data$SampleID %in% NegativeIDs, TRUE, FALSE)
contaminants_prevalence = isContaminant(phylo32, method="prevalence", neg="is.negative")

contaminants_prevalence$ID = row.names(contaminants_prevalence)
taxonomy_map_contams = taxonomy_map
taxonomy_map_contams$ID = save_OTUs
contaminants_prevalence = contaminants_prevalence %>% left_join(taxonomy_map_contams,by="ID")

# when you include all DNA-free water controls and blank controls 
contaminants_allNegativeControls = contaminants_prevalence %>% filter(contaminant)

# Firmicutes genera
Negatives32_Firmicutes = subset_taxa(Negatives32, Phylum=="Firmicutes")
Negatives32_Firmicutes_abundance  = transform_sample_counts(Negatives32_Firmicutes, function(x) x / sum(x) )
Negatives32_Firmicutes_abundance@otu_table@.Data[,"120603"] <- 0
Negatives32_Firmicutes_abundance_keep = filter_taxa(Negatives32_Firmicutes_abundance, function(x) mean(x) > 1e-5, TRUE)

Negatives32_Firmicutes = prune_taxa(row.names(Negatives32_Firmicutes_abundance_keep@tax_table),Negatives32_Firmicutes)
plot_bar(Negatives32_Firmicutes, "Sample_Ctrl", fill="Genus") + scale_fill_manual(values=ampalette) + ggtitle("OTUs of Firmicutes genera in negative controls")

# (2) Further identify OTUs present at >.1% in both blank swab extraction controls to remove from BOTH runs
###########################################################################################################
ExtractionControls = subset_samples(phylo32, ControlStatus=="ExtractionControl")
ExtractionControlsRelAbundance = transform_sample_counts(ExtractionControls, function(x) x / sum(x) )
ExtractionControlsRelAbundanceDF = data.frame(ExtractionControlsRelAbundance@otu_table@.Data)
ExtractionControlsRelAbundanceDF$OTU = row.names(ExtractionControlsRelAbundance@otu_table@.Data)
ExtractionControlsRelAbundanceDF = ExtractionControlsRelAbundanceDF %>% mutate(Extraction1Present = if_else(X120363 > .001, 1, 0), Extraction2Present = if_else(X120586 > .001, 1, 0))
ExtractionControls_Positive = ExtractionControlsRelAbundanceDF %>% filter(Extraction1Present+Extraction2Present == 2)

# What are we removing?
taxonomy_map_removes = taxonomy_map
taxonomy_map_removes$OTU = save_OTUs
ExtractionControls_Positive = ExtractionControls_Positive %>% left_join(taxonomy_map_removes, by="OTU")
intersect(contaminants_prevalence$ID, ExtractionControls_Positive$OTU)

# Filter based on findings
##########################
Remove32_Prevalence = contaminants_allNegativeControls$ID
TaxaToKeep = setdiff(taxa_names(phylo32),contaminants_allNegativeControls$ID)
phylo32_filtered = subset_samples(phylo32, ControlStatus == "NonControl")
phylo32_filtered = prune_taxa(TaxaToKeep, phylo32_filtered)


# One Staphylococcus OTU, one cutibacterium OTU, one Ralstonia OTU; All three of these were ID'd as contaminants by Decontam as well. 
Remove_from_Both = ExtractionControls_Positive$OTU

################################################################
# Process decontamination for MiSeqV1V3_35 --> phylo35_filtered
################################################################

# DNA-free water only
Controls35 = subset_samples(phylo35, ControlStatus!="NonControl")
# ignore empty wells
Controls35 = subset_samples(Controls35, ControlStatus!="Empty")
waterphyla35 = plot_bar(Controls35, "SampleID", fill="Phylum") + scale_fill_manual(values=ampalette) + ggtitle("Run 35 Water Control Phyla")

Controls35_Firmicutes = subset_taxa(Controls35, Phylum=="Firmicutes")
Controls35_FirmicutesAbundance =  transform_sample_counts(Controls35_Firmicutes, function(x) x / sum(x) )
Controls35_FirmicutesAbundance@otu_table@.Data[,c("DNAfreewater1", "DNAfreewater2", "DNAfreewater3", "DNAfreewater4")] <- 0
Controls35_Firmicutes_abundance_keep = filter_taxa(Controls35_FirmicutesAbundance, function(x) mean(x) > 1e-5, TRUE)
Controls35_Firmicutes = prune_taxa(taxa_names(Controls35_Firmicutes_abundance_keep), Controls35_Firmicutes)

watergenerafirm35 = plot_bar(Controls35_Firmicutes, "SampleID", fill="Genus") + scale_fill_manual(values=ampalette) + ggtitle("Run 35 Water Control Genera in Firmicutes")

ggsave(gridExtra::grid.arrange(waterphyla35, watergenerafirm35, ncol=2), file="Controls35_contents.pdf",width=13, height=7)

# Prevalence-based contaminant identification
#############################################
phylo35 = subset_samples(phylo35,ControlStatus!="Empty" )
phylo35@sam_data$is.negative = if_else(phylo35@sam_data$ControlStatus =="PCRControl", TRUE, FALSE)
contaminants_prevalence35 = isContaminant(phylo35, method="prevalence", neg="is.negative")
contaminants_prevalence35$ID = row.names(contaminants_prevalence35)
contaminants_prevalence35 = contaminants_prevalence35 %>% left_join(taxonomy_map_contams,by="ID")
OTUs_Remove_prevalence35 = (contaminants_prevalence35 %>% filter(contaminant))$ID

Remove35_Prevalence = OTUs_Remove_prevalence35

# Concentration-based contaminant identification
################################################
phylo35@sam_data$SampleID
ConcentrationInfo = read.csv2("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/mapping_files/run_maps/MiSeqV1V3_35.tsv", sep="\t")
ConcentrationInfo$SampleID = ConcentrationInfo$X.SampleID 
ConcentrationInfo = ConcentrationInfo %>% select(SampleID, final_dna_concentration_ng_ul)
phylo35DF = data.frame(phylo35@sam_data)
phylo35DF = phylo35DF %>% left_join(ConcentrationInfo, by="SampleID")
phylo35@sam_data$Concentration = sapply(phylo35DF$final_dna_concentration_ng_ul, function(x) as.numeric(as.character(x)))
phylo35 = subset_samples(phylo35, ControlStatus=="NonControl")
contaminants_frequency35 <- isContaminant(phylo35, method="frequency", conc="Concentration")
contaminants_frequency35 = contaminants_frequency35 %>% filter(contaminant)
contaminants_frequency35$ID = row.names(contaminants_frequency35)
contaminants_frequency35 = contaminants_frequency35 %>% left_join(taxonomy_map_contams, by="ID")

# none of the contaminants are in common with the prevalence-based ones
intersect(row.names(contaminants_frequency35), OTUs_Remove_prevalence35)

OTUs_Remove_concentration35 = contaminants_frequency35$ID

# Assemble list of OTUs to remove
# First combine prevalence-based and concentration-based contaminants
OTUs_Remove35 = append(OTUs_Remove_prevalence35, OTUs_Remove_concentration35)
# Then add the ones ID'd in the extraction controls from run 32
OTUs_Remove35 = unique(append(OTUs_Remove35, Remove_from_Both))
#make list of OTUs to keep 
OTUs_keep35 = setdiff(taxa_names(phylo35), OTUs_Remove35)

phylo35_filtered = prune_taxa(OTUs_keep35, phylo35)


###############################################
# Join together and do final prevalence filter
###############################################

# Join
phylo_joined = merge_phyloseq(phylo32_filtered, phylo35_filtered)
otunames = rownames(phylo_joined@otu_table)
sampnames = colnames(phylo_joined@otu_table)

# Strict filtering
###################
# Lists of removed OTUs
# Test stricter phylojoin version where you remove all OTUs removed from either from both 
Remove35_Prevalence
Remove32_Prevalence
OTUs_Remove_concentration35
AllContaminants = c(OTUs_Remove_concentration35, c(Remove35_Prevalence, Remove32_Prevalence))

# phylo_joined_STRICT = prune_taxa(setdiff(taxa_names(phylo_joined@otu_table), AllContaminants), phylo_joined)
# phylo_joined_relabundance_STRICT = transform_sample_counts(phylo_joined_STRICT, function(x) x / sum(x) )
# OTU_table_full_STRICT= data.frame(phylo_joined_relabundance_STRICT@otu_table)
# OTU_table_full_STRICT[is.na(OTU_table_full_STRICT)] <- 0 
# OTU_table_full_STRICT[OTU_table_full_STRICT < .001 ] <- 0
# OTU_table_full_STRICT$Prevalence <- apply(OTU_table_full_STRICT, 1, function(x) sum(x>0))
# OTU_table_full_STRICT  = subset(OTU_table_full_STRICT, Prevalence >=1 )
# OTU_table_full_STRICT$Prevalence = NULL
# OTUs_keep_prevalence_STRICT = rownames(OTU_table_full_STRICT)
# phylo_joined_STRICT = prune_taxa(OTUs_keep_prevalence_STRICT, phylo_joined_STRICT)
# 
# braycurtis_run_STRICT = phyloseq::distance(phylo_joined_STRICT, method = "bray")
# ordinationbray_STRICT = ordinate(phylo_joined_STRICT, method="PCoA", distance=braycurtis_run_STRICT)
# plot_ordination(phylo_joined_STRICT, ordinationbray_STRICT, color="Run") + theme(aspect.ratio=1) + scale_color_manual(values=twocolor)
# plot_ordination(phylo_joined_STRICT, ordinationbray_STRICT, color="Run") + theme(aspect.ratio=1) + scale_color_manual(values=wound_type)


phylo_joined_relabundance = transform_sample_counts(phylo_joined, function(x) x / sum(x) )
OTU_table_full = data.frame(phylo_joined_relabundance@otu_table)
OTU_table_full[is.na(OTU_table_full)] <- 0 
OTU_table_full[OTU_table_full < .001 ] <- 0


OTU_table_full$Prevalence <- apply(OTU_table_full, 1, function(x) sum(x>0))


# Filter to OTUs present in at least .1% in at least 1 sample
OTU_table_full  = subset(OTU_table_full, Prevalence >=1 )
OTU_table_full$Prevalence = NULL
OTUs_keep_prevalence = rownames(OTU_table_full)

phylo_joined = prune_taxa(OTUs_keep_prevalence, phylo_joined)


phylo_joined_backup = phylo_joined
phylo_joined_backup@sam_data$Reads_After_Decontam = colSums(phylo_joined@otu_table)
phylo_joined_lowcounts = subset_samples(phylo_joined_backup, Reads_After_Decontam <1000)
phylo_joined_lowcounts_Actino = subset_taxa(phylo_joined_lowcounts,Phylum=="Actinobacteriota")
phylo_joined_lowcounts_Firmicutes= subset_taxa(phylo_joined_lowcounts,Phylum=="Firmicutes")

plot_bar(phylo_joined_lowcounts_Actino, fill="Family")
plot_bar(phylo_joined_lowcounts_Firmicutes, fill="Family")



phylojoined_backup_prefilter = phylojoined

phylo_joined@sam_data$RemainingReads = colSums(phylo_joined@otu_table@.Data)
plotting_DF = data.frame(phylo_joined@sam_data)
patient_mapping32_subset = patient_mapping32 %>% select(SampleID,SubjectID)
patient_mapping35_subset = patient_mapping35 %>% select(SampleID, SubjectID)
patientmapping = rbind(patient_mapping32_subset, patient_mapping35_subset)

woundTypeInfo = patient_metadata %>% select(SubjectID, wound_type)
patientmapping = patientmapping %>% left_join(woundTypeInfo, by="SubjectID")
plotting_DF = plotting_DF %>% left_join(patientmapping, by="SampleID")


plotting_DF = plotting_DF %>% mutate(WoundDescription = case_when(wound_type == 1 ~ "PressureUlcer",
                                                                  wound_type==2 ~ "VenousUlcer",
                                                                  wound_type==3 ~"ArterialUlcer",
                                                                  wound_type==4 ~ "Surgical", 
                                                                  wound_type==5 ~ "Traumatic",
                                                                  wound_type==6 ~ "MixedTraumaticSurgical",
                                                                  TRUE ~ "Other"))

plottingDF = ggplot(plotting_DF, aes(x=WoundDescription, y=RemainingReads, group=WoundDescription, fill=WoundDescription)) + geom_boxplot() + facet_grid(.~Run) + scale_fill_manual(values=rev(ampalette))

plot_richness(phylo_joined, measures=c("Shannon"), x="Run") + theme_classic()

# Beta diversity between the runs; 
# Can't reject null that dispersions are the same between the runs

beta <- betadisper(braycurtis_run, phylo_joined@sam_data$Run)
permutest(beta)


braycurtis_run = phyloseq::distance(phylo_joined, method = "bray")
ordinationbray = ordinate(phylo_joined, method="PCoA", distance=braycurtis_run)
plot_ordination(phylo_joined, ordinationbray, color="Run") + theme(aspect.ratio=1) + scale_color_manual(values=twocolor)
plot_ordination(phylo_joined, ordinationbray, color="Run") + theme(aspect.ratio=1) + scale_color_manual(values=wound_type)

  

plot_richness(phylo_joined, x="Run", measures=c("Observed", "Simpson", "Shannon", "InvSimpson"))+ theme_classic() +  geom_boxplot() + stat_compare_means(method = "t.test", label.x=1.2)


 ggsave(plottingDF, file="ComparingFinalOTUcounts_RunType.pdf", width=14, height=7)

deseqobj_full= phyloseq_to_deseq2(phylo_joined, ~1) 
save(deseqobj_full, file="data/DESeqObject21.rda")

deseqobj = deseqobj_full
get_geometric_means <- function(x){
  mean = exp(sum(log(x[x>0]), na.rm=TRUE) / length(x))
   return(mean)
 }
 
matrixcounts = counts(deseqobj)
# 
gmeans = apply(counts(deseqobj), 1, get_geometric_means)
# 
 estimates_sizes = estimateSizeFactors(deseqobj, geoMeans = gmeans)
# 
 disp_estimates = estimateDispersions(estimates_sizes)
 
 v_stabilized = getVarianceStabilizedData(disp_estimates)
 

 save(v_stabilized, file="data/V_stabilized_IowaWound_Combined_2021.rda")
 
####################
# Post-normalization
####################
 
# Set < 0 count values to 0 after the log transformation
v_stabilized[v_stabilized < 0] <-  0 
 
# Feed normalized data into the phyloseq object
phylo_joined_postNorm = phylo_joined
phylo_joined_postNorm@otu_table = otu_table(v_stabilized, taxa_are_rows = T)

pre_norm_relative = transform_sample_counts(phylo_joined, function(x) x / sum(x) )

pre_norm = plot_bar(pre_norm_relative , "SampleID", fill="Phylum") + scale_fill_manual(values=ampalette)

other_phyla = setdiff(unique(pre_norm$data$Phylum), c("Bacteroidota", "Proteobacteria","Actinobacteriota","Firmicutes"))
phyla_order=c( other_phyla , "Bacteroidota","Proteobacteria","Actinobacteriota", "Firmicutes")

pre_norm$data$Phylum = factor(pre_norm$data$Phylum,levels=phyla_order)


pre_norm_data = phylo_joined %>% 
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% 
  psmelt() %>% 
  group_by(Phylum) 

PRENormPalette = unique(rev(c("#B85C00","#999999","#339966","#6B24B2","#4D4D4D", "#CC0000", 
                       "#CC99FF","#663300","#33CC33", "#0072B2",  
                       "#9900FF","#B85C00","#339966","#6B24B2","#56B4E9","#992900","#006600", 
                       "#4D4D4D", "#0072B2","#FF9900")))

Palette17 = c("blue4","aquamarine4","coral1", "darkgoldenrod1", "brown4","darkgreen", "darkolivegreen3","dodgerblue2", "darksalmon","mediumpurple4","lightseagreen", "tan2", "paleturquoise", "yellow1", "maroon4", "palevioletred2", "orange4")

prenorm_plot = ggplot(pre_norm_data, aes(x=Sample, y=Abundance, fill=Phylum)) + geom_bar(stat="identity")  +
  ggtitle("Phylum-Level Composition of Pre-normalized MiSeq_32 and 35 Samples Combined (Grouped by MiSeq Run)") + theme_classic()+
  theme(axis.text.x=element_blank(), legend.title=element_text(size=14), legend.text=element_text(size=12), legend.key.size = unit(.8, "cm")) + 
  scale_fill_manual(values=rev(Palette17)) 

SampleOrder = unique((prenorm_plot$data %>% filter(Phylum=="Firmicutes") %>% arrange(Run, Abundance))$Sample)

prenorm_plot$data$Sample = factor(prenorm_plot$data$Sample, levels=SampleOrder)
prenorm_plot$data$Phylum = factor(prenorm_plot$data$Phylum, levels=phyla_order)


post_norm_data = phylo_joined_postNorm %>% 
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% 
  psmelt() %>% 
  group_by(Phylum) 


postnorm_plot = ggplot(post_norm_data, aes(x=Sample, y=Abundance, fill=Phylum)) + geom_bar(stat="identity")  +
  ggtitle("Phylum-Level Composition of Normalized MiSeq_32 and 35 Samples Combined (Grouped by MiSeq Run)") + theme_classic()+
  theme(axis.text.x=element_blank(), legend.title=element_text(size=14), legend.text=element_text(size=12), legend.key.size = unit(.8, "cm")) + 
  scale_fill_manual(values=rev(Palette17)) 
postnorm_plot$data$Sample = factor(postnorm_plot$data$Sample, levels=SampleOrder)
postnorm_plot$data$Phylum = factor(postnorm_plot$data$Phylum, levels=phyla_order)

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  postnorm_plot 
)

postnorm_plot = postnorm_plot + theme(legend.position="none")
prenorm_plot = prenorm_plot + theme(legend.position="none")
gridExtra::grid.arrange(prenorm_plot, postnorm_plot, nrow=2)

TwoPlots = cowplot::plot_grid(prenorm_plot, postnorm_plot, nrow=2)
combinedplots = cowplot::plot_grid(TwoPlots, legend, rel_widths = c(20, 2))  

prenormOrdination= ordinate(phylo_joined, "PCoA", "bray")
PreNormBrayPlot = plot_ordination(phylo_joined, prenormOrdination, color="Run") + 
  ggtitle("Bray Curtis-based PCoA Coordinates 1 and 2 (Pre-normalization)") + scale_color_manual(values=twocolor) + geom_point(size=2) 


postnormOrdination= ordinate(phylo_joined_postNorm, "PCoA", "bray")
PostNormBrayPlot = plot_ordination(phylo_joined_postNorm, postnormOrdination, color="Run") + 
  ggtitle("Bray Curtis-based PCoA Coordinates 1 and 2(after VST normalization)") + scale_color_manual(values=twocolor) + geom_point(size=3) 

gridExtra::grid.arrange(PreNormBrayPlot, PostNormBrayPlot, nrow=2)


genus_level =  phylo_joined_postNorm %>% 
  tax_glom(taxrank = "Genus") 
postnormOrdinationGenus= ordinate(genus_level, "PCoA", "bray")
PostNormBrayPlotGenus = plot_ordination(genus_level, postnormOrdinationGenus, color="Run") + 
  ggtitle("Bray Curtis-based PCoA Coordinates 1 and 2(genus-aggregated)") + scale_color_manual(values=twocolor) + geom_point(size=2) 

gridExtra::grid.arrange(PostNormBrayPlot, PostNormBrayPlotGenus)



# Adding genus code 
taxtablepost = phylo_joined_postNorm@tax_table
taxtablepost[is.na(taxtablepost)] <- "NA"
tax_table(phylo_joined_postNorm) <- taxtablepost
Vstabilized_GenusCode= sapply(rownames(phylo_joined_postNorm@tax_table), function(x) Genus_code(phylo_joined_postNorm@tax_table[x, ]))
physeq_tax = data.frame(phylo_joined_postNorm@tax_table)
physeq_tax$GenusCode = Vstabilized_GenusCode
phylo_joined_postNorm@tax_table = tax_table(physeq_tax)
colnames(phylo_joined_postNorm@tax_table) = colnames(physeq_tax)
rownames(phylo_joined_postNorm@tax_table) = rownames(physeq_tax)


# Adding genus code 
taxtablepre = phylo_joined@tax_table
taxtablepre[is.na(taxtablepre)] <- "NA"
tax_table(phylo_joined) <- taxtablepre
PreNormGenusCode= sapply(rownames(phylo_joined@tax_table), function(x) Genus_code(phylo_joined@tax_table[x, ]))
physeq_tax_pre = data.frame(phylo_joined@tax_table)
physeq_tax_pre$GenusCode = PreNormGenusCode
phylo_joined@tax_table = tax_table(physeq_tax_pre)
colnames(phylo_joined@tax_table) = colnames(physeq_tax_pre)
rownames(phylo_joined@tax_table) = rownames(physeq_tax_pre)




genus_code_level =  phylo_joined_postNorm %>% 
  tax_glom(taxrank = "GenusCode") 
postnormOrdinationGenusCode= ordinate(genus_code_level, "PCoA", "bray")
PostNormBrayPlotGenusCode = plot_ordination(genus_code_level, postnormOrdinationGenus, color="Run") + 
  ggtitle("Bray Curtis-based PCoA Coordinates 1 and 2(genus-aggregated after VST Normalization)") + scale_color_manual(values=twocolor) + geom_point(size=3) 

gridExtra::grid.arrange(PostNormBrayPlot, PostNormBrayPlotGenusCode)



genus_code_level_prop =  phylo_joined_postNorm %>% 
  tax_glom(taxrank = "GenusCode")  %>% 
  transform_sample_counts(function(x) {x/sum(x)})
postnormOrdinationGenusCode= ordinate(genus_code_level_prop, "PCoA", "bray")
PostNormBrayPlotGenusCodeProp = plot_ordination(genus_code_level_prop, postnormOrdinationGenusCode, color="Run") + 
  ggtitle("Bray Curtis-based PCoA Coordinates 1 and 2(genus-aggregated after VST Normalization)") + scale_color_manual(values=twocolor) + geom_point(size=3) 

gridExtra::grid.arrange(PostNormBrayPlot, PostNormBrayPlotGenusCode)

post_norm_unaggregatedPCA = prcomp(t(phylo_joined_postNorm@otu_table))
post_norm_aggregatedPCA = prcomp(t(genus_code_level))


phylo_joined_postNorm_32Only = subset_samples(phylo_joined_postNorm, Run=="MiSeqV1V3_32")
phylo_joined_postNorm_35Only = subset_samples(phylo_joined_postNorm, Run=="MiSeqV1V3_35")




#Pre_Batch_Plot + xlim(c(-.4,.45)) + ylim(c(-.4, .4))
#  gmeans_runs = apply(counts(deseqobjRun), 1, get_geometric_means)
#  # 
#  estimates_sizes_run = estimateSizeFactors(deseqobjRun, geoMeans = gmeans_runs)
#  # 
#  disp_estimates_run = estimateDispersions(estimates_sizes_run)
#  
#  v_stabilized_run = getVarianceStabilizedData(disp_estimates_run)

# # Compare what happens with and without a design factor
# run = unlist(as.list(v_stabilized_run))
# regular = unlist(as.list(v_stabilized))
# 
# df_plot = data.frame(run)
# df_plot$regular = regular
# ggplot(df_plot, aes(x=run, y=regular)) + geom_point()

