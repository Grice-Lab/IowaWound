library(dplyr)
library(tidyr)

# Include all samples from MiSeqV1V3_32 
#######################################
ClubGriceMetadata = read.csv("data/RawIntermediateFiles/MiSeqV1V3_32_raw/MiSeqV1V3_32.tsv", sep="\t")
UpdatedCuratedMetadata = read.csv("data/RawIntermediateFiles/Iowa32SampleMetadata.tsv", sep="\t")
UpdatedCuratedMetadata_Subset = UpdatedCuratedMetadata %>% select(c(X.SampleID, SubjectID))

# Test that barcodes are equivalent between sets  
testframe =  UpdatedCuratedMetadata %>% left_join(ClubGriceMetadata, by="X.SampleID")
testframe$BarcodeSequence.y = factor(testframe$BarcodeSequence.y, levels = unique(testframe$BarcodeSequence.y))
(testframe$BarcodeSequence.x == as.factor(testframe$BarcodeSequence.y))

newdataframe = ClubGriceMetadata %>% left_join(UpdatedCuratedMetadata_Subset, by="X.SampleID")
newdataframe$SubjectID = factor(newdataframe$SubjectID, levels=c(levels(newdataframe$SubjectID), "Exclude"))
newdataframe$SubjectID[is.na(newdataframe$SubjectID)] <- "Exclude"
write.table(newdataframe, quote=F, sep='\t', file="data/RawIntermediateFiles/MiSeqV1V3_32_FullMetadata.tsv", row.names = F)

 
# Include all samples from MiSeqV1V3_35
ClubGriceMetadata35 = read.csv("data/RawIntermediateFiles/MiSeqV1V3_35.tsv", sep="\t")
UpdatedCuratedMetadata35 = read.csv("data/RawIntermediateFiles/Iowa35SampleMetadata.tsv", sep="\t")
colnames(UpdatedCuratedMetadata35)
UpdatedCuratedMetadata35_Subset = UpdatedCuratedMetadata35 %>% select(c(X.SampleID, SubjectID))
newdataframe35 = ClubGriceMetadata35 %>% left_join(UpdatedCuratedMetadata35_Subset, by = "X.SampleID")
newdataframe35 = newdataframe35 %>% select(c(X.SampleID, SubjectID.y, library_concentration_ng_ul, final_dna_concentration_ng_ul, dna_concentration_ng_ul, project, Plate, WellPosition, BarcodeSequence, LinkerPrimerSequence, ReversePrimerSequence ))
colnames(newdataframe35) = c("#SampleID", "SubjectID", "library_concentration_ng_ul","final_dna_concentration_ng_ul",  "dna_concentration_ng_ul", "project", "Plate", "WellPosition", "BarcodeSequence", "LinkerPrimerSequence", "ReversePrimerSequence" )
write.table(newdataframe35, quote=F, sep='\t', file="data/RawIntermediateFiles/MiSeqV1V3_35_FullMetadata.tsv", row.names = F)


# View(newdataframe35)
