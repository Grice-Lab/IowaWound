# Amy Campbell
# 09/09/2020 
# Merging Iowa metadata with the qPCR results Simon gave me

library("reshape2")
library("dplyr")

setwd( "/Users/amycampbell/Desktop/GriceLabGit/IowaWound/")
qpcrdata = read.csv("data/Compiled_IWS_18S_Analyzed_Ct35.csv")
metadata = read.csv("data/IowaWoundPainMetaData.csv")

qpcrdata = qpcrdata[c("Study.ID", "Target.Name", "Ct.Tech.Mean", "Delta.Ct.Mean")]
qpcrdata$study_id = qpcrdata$Study.ID

metadata$study_id = sapply(metadata$study_id, function(x) as.character(x))
qpcrdata$Study.ID = NULL


qpcrdata = qpcrdata %>% left_join(metadata, by = "study_id")
write.csv(qpcrdata, file="MergedData_for_SK.csv")

