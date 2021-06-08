# Amy Campbell
# 06/2021
# Making sense of metadata and making samples-to-keep.tsv files 

library(dplyr)
library(stringr)
metadata32 = read.csv2("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/mapping_files/metadata/IowaWound_metadata.tsv", sep = '\t', header=T)
Amy_metadataOct20 = read.csv("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/mapping_files/run_maps/IA_woundpain_mappingfile_AEC_UpdateOct20.csv", header=T)



# the metadata EAG originally emailed me in april 2019
#######################################################
EAG_April19 = read.csv("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/figuring_out_metadata_5_21/IA_woundpain_mappingfile_SentFromEAG_4-19.csv", header=T)

Total_32_Samples_April19 = EAG_April19 %>% filter(RunName == "MiSeqV1V3_32" )
Total_35_Samples_April19 = EAG_April19 %>% filter(RunName == "MiSeqV1V3_35" )

# EAG missing 
# 174, 203, 225,370, 425, 494, 521, 540, 691, 704, 825, 1051

Missing_Accordingto_Julia = c(174, 20, 370, 425, 494, 521, 540, 691, 704, 825, 1051, 225)

# Metadata sent from EAG 5/19 (Gardner Data, no run info)
#########################################################
EAG_May19 = read.csv("/Users/amycampbell/Desktop/GriceLabGit/IowaWound/figuring_out_metadata_5_21/GSWOUNDGRICE2015_20190221.csv", header=T)
dim(EAG_May19) # Should be 445 subjects represented total (unlikely)


# Run 32
#########
# Run map 
runmap32 = read.csv2("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/mapping_files/run_maps/MiSeqV1V3_32.tsv", sep = '\t', header=T)
metadata32 = metadata32 %>% select(X.SampleID, Patient.ID)
runmap32reduced = runmap32 %>% select(X.SampleID, Control)
mapped_metadata32 = runmap32reduced %>% left_join(metadata32, by ="X.SampleID")

# metadata from Club_Grice 

CGMetadata = read.csv2("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/mapping_files/metadata/IowaWound_metadata.tsv", sep = "\t", header=T)
Missing_From_April19 = setdiff(CGMetadata$Description, Total_32_Samples_April19$Description)
CGMetadata_Missing_From_April19 = CGMetadata %>% filter(Description %in% Missing_From_April19)


# These are studyIDs of patients whose samples were run on the MiSeq but who didn't 
# end up in the metadata sheet EAG gave me in april 2019 (because they weren't in )
unique(CGMetadata_Missing_From_April19$Patient.ID)
# 2     4     5     6     7     12    14    15    16    19    20    21    22    Blank


# Only ones in CGMetadata_Missing_From_April19 that should actually be included are the blanks (Sample_120363, Sample_120586)
setdiff(CGMetadata_Missing_From_April19$Patient.ID, EAG_May19$study_id)


# Amy_metadataOct20
##################################
Amy_metadataOct20_Run32 = Amy_metadataOct20 %>% filter(RunName=="MiSeqV1V3_32")
dim(Amy_metadataOct20_Run32)

setdiff(Amy_metadataOct20_Run32$Description, EAG_April19$Description) # only differences between that and EAG's metadata from 04/19 are "Sample_120586", "Sample_120363" which are blank swabs from that run (bless)
setdiff(Total_32_Samples_April19$Description, Amy_metadataOct20_Run32$Description)
Amy_metadataOct20_Run32 %>% filter(Description %in% c("Sample_120586", "Sample_120363"))

# buuut I should also include the mock community (120602). 

# all the patients missing from EAG's April 2019 metadata are also missing from GSWOUNDGRICE2015_20190221.csv which indicates to me that
# they dropped from the study after their samples were submitted. They should be excluded.  
Patients_MissingFromApril19 = c(2, 4, 5, 6 , 7, 12, 14, 15, 16, 19, 20, 21, 22)

# Mark for exclusion any non-control "Patient.ID"'s 
Mark_For_Exclusion = unique(setdiff(CGMetadata$Patient.ID, EAG_May19$study_id))
Mark_For_Exclusion = Mark_For_Exclusion[!(Mark_For_Exclusion %in% c("Blank", "Water", "Mock"))]

# Study IDs (aka patient IDs) marked for exclusion based on absence from metadata
IncludedSamples32  = (CGMetadata %>% filter(!(Patient.ID %in% Mark_For_Exclusion)))

length(unique(IncludedSamples32$Patient.ID))

# 27 and 28 have two samples run on the MiSeq
sort(table(IncludedSamples32$Patient.ID))
IncludedSamples32$X.SampleID

write.table(IncludedSamples32$X.SampleID, quote=F, row.names=F, col.names=F, file="/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/samples-to-keep-32.tsv")


# Run 35
#########

runmap35 = read.csv2("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/mapping_files/run_maps/MiSeqV1V3_35.tsv", sep = '\t', header=T)

Amy_metadataOct20_Run35 = Amy_metadataOct20 %>% filter(RunName=="MiSeqV1V3_35")

setdiff(Amy_metadataOct20_Run35$Description, Total_35_Samples_April19$Description)

setdiff(Amy_metadataOct20_Run35$Description, Total_35_Samples_April19$Description)
setdiff( Total_35_Samples_April19$Description, Amy_metadataOct20_Run35$Description)

Iowa35 = runmap35 %>% filter(project=="Iowa")
missing_from_April19 = setdiff(Iowa35$Description,Total_35_Samples_April19$Description)

# patients to be excluded because they aren't in the up to date metadata 
Patients_Missing_from_UTD_metadata = setdiff(Iowa35$SubjectID, EAG_May19$study_id)
IncludeSamples35 = runmap35 %>% filter(project=="Iowa" | is.na(project)) %>% filter(X.SampleID != "Mysterysample1")
IncludeSamples35 = IncludeSamples35 %>% filter(!(SubjectID %in% Patients_Missing_from_UTD_metadata))
View(IncludeSamples35)

write.table(IncludeSamples35$X.SampleID, quote=F, row.names=F, col.names=F, file="/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/samples-to-keep-35.tsv")

# "IowaWound.Human.1252", "IowaWound.Human.947", "IowaWound.Human.1501","IowaWound.Human.1576", "IowaWound.Human.1311", "IowaWound.Human.1495"

# 1576 is in the metadata but not EAG's sample list 
# IowaWound.Human.1576
# Also MiSeqV1V3_35_barcode_IowaWound.Human.1258?
# KEEP 

# All subjects represented in EAG's sample list of run 35 from 04-19
Subjects35_April19 = Total_35_Samples_April19$SubjectID
Subjects35_April19 = Subjects35_April19[!(Subjects35_April19=="Water") & !is.na(Subjects35_April19)]
Subjects35_April19 = sapply(Subjects35_April19, function(x) as.numeric(as.character(x)))
AllSubjects35_Total = c(Subjects35_April19, 1576)
# remove those that shouldnt be here 
AllSubjects35_Total = (AllSubjects35_Total[! AllSubjects35_Total %in% c(1252, 947, 1501, 1630, 1566, 1311, 1495)])


exclude_35 = missing_from_May19_Metadata35

exclude_35_all= setdiff(Iowa35$SubjectID,  EAG_May19$study_id)

AllSubjects32_april2019 = sapply(AllSubjects32_april2019, function(x) as.numeric(as.character(x)))
AllSubjectsrepresented = append(AllSubjects35_Total, AllSubjects32_april2019)
intersect(AllSubjects35_Total, AllSubjects32_april2019)
length(AllSubjectsrepresented)

AllSubjectsrepresentedUnique = unique(AllSubjectsrepresented)
total_missing = setdiff(EAG_May19$study_id, AllSubjectsrepresentedUnique)
total_missing

Missing_Accordingto_Julia = c(174, 20, 370, 425, 494, 521, 540, 691, 704, 825, 1051, 225)
Missing_accordingto_EAG = c(174,203, 225,370, 425,494, 521, 540, 691, 704, 825, 1051)

setdiff(total_missing,Missing_Accordingto_Julia )
setdiff(Missing_Accordingto_Julia, total_missing )

# Ones that EAG never asked about 432
setdiff(total_missing,Missing_accordingto_EAG )
setdiff(Missing_accordingto_EAG, total_missing )

View(Iowa35)
