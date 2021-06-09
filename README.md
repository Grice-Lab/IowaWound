# IowaWound
Iowa Wound Pain analyses 2020

# (1) Make manifest .tsv’s for important already-demultiplexed reads

## Input
- Raw demultiplexed data path from Qi’s pipeline in IowaWoundData/MiSeqV1V3_35/demultiplexed, IowaWoundData/MiSeqV1V3_32/demultiplexed
## Outputs
-  GH:IowaWound/mappings/Manifest35.tsv
-  GH:IowaWound/mappings/Manifest35.tsv	

## Scripts

- [make_manifest_32.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/make_manifest_35.sh)
- [make_manifest_35.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/make_manifest_32.sh)
    
# (2) Import data into Qiime2 

## Input
- Raw demultiplexed data path from Qi’s pipeline in IowaWoundData/MiSeqV1V3_35/demultiplexed, IowaWoundData/MiSeqV1V3_32/demultiplexed
-  GH:IowaWound/mappings/Manifest35.tsv, GH:IowaWound/mappings/Manifest32.tsv
- 
## Outputs
- IowaWoundData/MiSeqV1V3_35/paired-end-demux35.qza
- IowaWoundData/MiSeqV1V3_32/paired-end-demux32.qza

## Scripts

[import_data.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/import_data.sh)
