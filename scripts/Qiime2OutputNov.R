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

# Set random seed to BRB zipcode :-]
set.seed(19104)

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
OTU_Table = read.csv2("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/table-with-taxonomy.tsv", header=T, sep="\t",skip=1)
runinfo = read.csv2("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/Control_Run_Info.tsv", sep=" ")
taxonomytable = read.csv2("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/taxonomy.tsv", sep="\t")
treephyseq = system.file("extdata","/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/tree.nwk",package="phyloseq")
twocolor = c("#FFC20A", "#0C7BDC")

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
PhyloseqObjectForUniFrac = subset_samples(PhyloseqObject,TotalOTUs > 0 )


unifracdists <- phyloseq::distance(PhyloseqObjectForUniFrac, method="wunifrac")
unifracWeighted = ordinate(PhyloseqObjectForUniFrac,"PCoA", distance=unifracdists)
wUF = plot_ordination(PhyloseqObjectForUniFrac, unifracWeighted, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (ASV-level weighted UniFrac distance) in MiSeq Runs") + theme_classic()

Unweighted_unifracdists <- phyloseq::distance(PhyloseqObjectForUniFrac, method="unifrac")
unifracUnWeighted = ordinate(PhyloseqObjectForUniFrac,"PCoA", distance=Unweighted_unifracdists)
uwUF = plot_ordination(PhyloseqObjectForUniFrac, unifracUnWeighted, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (ASV-level unweighted UniFrac distance) in MiSeq Runs") + theme_classic()
gridExtra::grid.arrange(wUF, uwUF, ncol=2)

# Seems that beta divergence is significant between runs by unweighted, but not weighted, UniFrac
adonis(unifracdists ~sample_data(PhyloseqObjectForUniFrac)$Run)
adonis(Unweighted_unifracdists ~sample_data(PhyloseqObjectForUniFrac)$Run)

# See if unifrac distances cluster by Run (after tax_glom)
PhyloseqObjectForUniFrac_GenusGlom =  tax_glom(PhyloseqObjectForUniFrac, taxrank = "Genus")
unifracdistsGenus <- phyloseq::distance(PhyloseqObjectForUniFrac_GenusGlom, method="wunifrac")
unifracWeightedGenus = ordinate(PhyloseqObjectForUniFrac_GenusGlom,"PCoA", distance=unifracdistsGenus)
wUFgenus = plot_ordination(PhyloseqObjectForUniFrac_GenusGlom, unifracWeightedGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (Genus-level weighted UniFrac distance) in MiSeq Runs") + theme_classic()

Unweighted_unifracdistsGenus <- phyloseq::distance(PhyloseqObjectForUniFrac_GenusGlom, method="unifrac")
unifracUnWeightedGenus = ordinate(PhyloseqObjectForUniFrac_GenusGlom,"PCoA", distance=Unweighted_unifracdistsGenus)
uwUFGenus = plot_ordination(PhyloseqObjectForUniFrac_GenusGlom, unifracUnWeightedGenus, color="Run") + scale_colour_manual(values=twocolor) + ggtitle("Principal Coordinates (Genus-level unweighted UniFrac distance) in MiSeq Runs") + theme_classic()
gridExtra::grid.arrange(wUFgenus, uwUFGenus, ncol=2)

adonis(unifracdistsGenus ~sample_data(PhyloseqObjectForUniFrac_GenusGlom)$Run)
adonis(Unweighted_unifracdistsGenus ~sample_data(PhyloseqObjectForUniFrac_GenusGlom)$Run)











rownames(OTU_Table) = taxfull
OTU_Table$Taxonomy = taxfull
transposedOTU = data.frame(t(OTU_Table[,2:ncol(OTU_Table)]))
colnames(transposedOTU) = taxaIDs
totalNonZeros = rowSums(transposedOTU!=0) 

transposedOTU$sums = rowSums(transposedOTU)
transposedOTU$Sample=rownames(transposedOTU)
ggplot(transposedOTU, aes(x=sums)) + geom_histogram(binwidth=1000, fill="aquamarine", color="black") +scale_x_continuous( breaks=seq(0,160000,5000)) + theme_classic() + theme(axis.text.x=element_text(angle=270)) + xlab("Total ASV Counts") + ylab("Frequency") + ggtitle("Assignments Across Runs")
transposedOTU$nonzeros =totalNonZeros 

transposedOTU$SampleID = row.names(transposedOTU)
transposedOTU$Run
transposedOTU$SampleID = sapply(transposedOTU$SampleID, function(x) str_replace(x, "X", ""))
transposedOTU = transposedOTU %>% left_join(runinfo, by="SampleID")
ggplot(transposedOTU, aes(x=nonzeros)) + geom_histogram(binwidth=5, fill="salmon", color="black") +scale_x_continuous( breaks=seq(0,320, 20)) + theme_classic() + theme(axis.text.x=element_text(angle=270)) + xlab("ASV Richness") + ylab("Frequency") + ggtitle("Assignments Across Runs(Richness)") + facet_grid(.~Run)
ggplot(transposedOTU, aes(x=sums)) + geom_histogram(binwidth=1000, fill="aquamarine", color="black") +scale_x_continuous( breaks=seq(0,160000,5000)) + theme_classic() + theme(axis.text.x=element_text(angle=270)) + xlab("Total ASV Counts") + ylab("Frequency") + ggtitle("Assignments Across Runs")+ facet_grid(.~Run)

transposedOTU32 = transposedOTU %>% filter(Run=="MiSeqV1V3_32")
median(transposedOTU32$sums)
median(transposedOTU32$nonzeros)

transposedOTU35 = transposedOTU %>% filter(Run=="MiSeqV1V3_35")
median(transposedOTU35$sums)
median(transposedOTU35$nonzeros)


controls_run = read.csv2("mappings/Control_Run_Info.tsv", header=T, sep=" ")
patient_mapping32 = read.csv("mappings/IA_woundpain_mapping_32_2021.csv")
patient_mapping35 = read.csv("mappings/IA_woundpain_mapping_35_2021.csv")
patient_metadata = read.csv("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/figuring_out_metadata_5_21/GSWOUNDGRICE2015_20190221.csv")

patient_mapping32$SampleID = as.factor(patient_mapping32$X.SampleID)
patient_mapping35$SampleID = as.factor(patient_mapping35$X.SampleID)
