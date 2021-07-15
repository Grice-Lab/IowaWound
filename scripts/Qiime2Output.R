# Amy Campbell
# Output of Qiime2 
# July 2021

library("phyloseq")
library("dplyr")
library("stringr")
library("decontam")
library("ggplot2")
setwd("~/Desktop/GriceLabGit/IowaWound/")
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

# Read in data & metadata 
OTU_Table = read.csv2("data/table.from_biom_w_taxonomy.txt", header=T, sep="\t",skip=1)
controls_run = read.csv2("mappings/Control_Run_Info.tsv", header=T, sep=" ")
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


################################################################
# Process decontamination for MiSeqV1V3_32 --> phylo32_filtered
################################################################

# (1) Identify contaminants to remove from MiSeqV1V3_32 based on decontam
#########################################################################
Controls32 = subset_samples(phylo32, ControlStatus!="NonControl")
Controls32@sam_data$Sample_Ctrl = paste(Controls32@sam_data$SampleID, Controls32@sam_data$ControlStatus, sep="_")
plot_bar(Controls32, "Sample_Ctrl", fill="Phylum") + scale_fill_manual(values=ampalette)

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

Negatives32_Firmicutes = prune_taxa(row.names(Negatives32_Firmicutes_abundance_REMOVE@tax_table),Negatives32_Firmicutes)
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
plot_bar(Controls35, "SampleID", fill="Phylum") + scale_fill_manual(values=ampalette) + ggtitle("Run 35 Water Control Phyla")

Controls35_Firmicutes = subset_taxa(Controls35, Phylum=="Firmicutes")
Controls35_FirmicutesAbundance =  transform_sample_counts(Controls35_Firmicutes, function(x) x / sum(x) )
Controls35_FirmicutesAbundance@otu_table@.Data[,c("DNAfreewater1", "DNAfreewater2", "DNAfreewater3", "DNAfreewater4")] <- 0
Controls35_Firmicutes_abundance_keep = filter_taxa(Controls35_FirmicutesAbundance, function(x) mean(x) > 1e-5, TRUE)
Controls35_Firmicutes = prune_taxa(taxa_names(Controls35_Firmicutes_abundance_keep), Controls35_Firmicutes)

plot_bar(Controls35_Firmicutes, "SampleID", fill="Genus") + scale_fill_manual(values=ampalette) + ggtitle("Run 35 Water Control Genera in Firmicutes")

# Prevalence-based contaminant identification
#############################################
phylo35 = subset_samples(phylo35,ControlStatus!="Empty" )
phylo35@sam_data$is.negative = if_else(phylo35@sam_data$ControlStatus =="PCRControl", TRUE, FALSE)
contaminants_prevalence35 = isContaminant(phylo35, method="prevalence", neg="is.negative")
contaminants_prevalence35$ID = row.names(contaminants_prevalence35)
contaminants_prevalence35 = contaminants_prevalence35 %>% left_join(taxonomy_map_contams,by="ID")
OTUs_Remove_prevalence35 = (contaminants_prevalence35 %>% filter(contaminant))$ID

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

deseqobj_full = phyloseq_to_deseq2(phylo_joined, ~SampleID)

save(deseqobj_full, file="data/DESeqObject21.rda")



