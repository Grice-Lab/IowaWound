# Exploring Clinical Data
library(dplyr)a
library(ggplot2)

ClinicalData = read.csv("/Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/GSWOUNDGRICE2015_20190221.csv")

ClinicalData$woundcarepain

ClinicalData = ClinicalData %>%  mutate(PainCatBinary = case_when( woundcarepain==0 | woundcarepain==1  ~"NoneMild",
                                                                   woundcarepain==3~ "Severe", 
                                                                   woundcarepain==2 ~ "Moderate" 
))
ClinicalData = ClinicalData %>% filter(PainCatBinary !="Moderate")


ggplot(ClinicalData, aes(y=pain_dress_change, x=dressingcat, group=dressingcat)) + geom_boxplot()
