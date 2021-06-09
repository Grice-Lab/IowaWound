# IowaWound
Iowa Wound Pain analyses 2021

# (1) Make manifest .tsv’s for important already-demultiplexed reads

## Input
- Raw demultiplexed data path from Qi’s pipeline in IowaWoundData/MiSeqV1V3_35/demultiplexed, IowaWoundData/MiSeqV1V3_32/demultiplexed
## Outputs
-  GH:IowaWound/mappings/Manifest35.tsv
-  GH:IowaWound/mappings/Manifest35.tsv	

## Scripts

- [make_manifest_32.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/make_manifest_35.sh)
- [make_manifest_35.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/make_manifest_32.sh)
    
# (2) Import data into Qiime2 and summarize input

## Input
- Raw demultiplexed data path from Qi’s pipeline in IowaWoundData/MiSeqV1V3_35/demultiplexed, IowaWoundData/MiSeqV1V3_32/demultiplexed
-  GH:IowaWound/mappings/Manifest35.tsv, GH:IowaWound/mappings/Manifest32.tsv

## Outputs
- IowaWoundData/MiSeqV1V3_35/paired-end-demux35.qza
- IowaWoundData/MiSeqV1V3_32/paired-end-demux32.qza

## Scripts

- [import_data.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/import_data.sh)
- [summarize_input.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/summarize_input.sh)

# (3) FastQC to get quality metrics 

## Input
- Raw demultiplexed data path from Qi’s pipeline in IowaWoundData/MiSeqV1V3_35/demultiplexed, IowaWoundData/MiSeqV1V3_32/demultiplexed

## Outputs
- Full MultiQC output to IowaWoundData/MultiQC_PostDemuxed

## Scripts

- [fastqc.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/fastqc.sh)

# (3) FastQC to get quality metrics 

## Input
- Raw demultiplexed data path from Qi’s pipeline in IowaWoundData/MiSeqV1V3_35/demultiplexed, IowaWoundData/MiSeqV1V3_32/demultiplexed

## Outputs
- Full MultiQC output to IowaWoundData/MultiQC_PostDemuxed

## Scripts

- [fastqc.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/fastqc.sh)

# (4) Selection of truncation lengths and denoising with DADA

## Input
- IowaWoundData/MultiQC_PostDemuxed/mqc_fastqc_per_base_sequence_quality_plot_1.txt 
- IowaWoundData/MiSeqV1V3_35/paired-end-demux35.qza
- IowaWoundData/MiSeqV1V3_32/paired-end-demux32.qza
- denoising.sh parameters for MiSeqV1V3_32 chosen in [ChooseDadaTruncations.R](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/ChooseDadaTruncations.R):
- - –p-trunc-len-f 272
- - -p-trunc-len-r 257
- denoising.sh parameters for MiSeqV1V3_35 chosen in [ChooseDadaTruncations.R](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/ChooseDadaTruncations.R):
- - –p-trunc-len-f 262
- - –p-trunc-len-r 232 

## Outputs
- IowaWoundData/MiSeqV1V3_32/rep-seqs32.qza
- IowaWoundData/MiSeqV1V3_32/table32.qza
- IowaWoundData/MiSeqV1V3_32/denoising-stats32.qza
- IowaWoundData/MiSeqV1V3_35/rep-seqs35.qza
- IowaWoundData/MiSeqV1V3_35/table35.qza
- IowaWoundData/MiSeqV1V3_35/denoising-stats35.qza

## Scripts
- [ChooseDadaTruncations.R](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/ChooseDadaTruncations.R)
- [denoise.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/denoise.sh)

# (5) Visualizations of denoising output 
## Input
- IowaWoundData/MiSeqV1V3_32/rep-seqs32.qza
- IowaWoundData/MiSeqV1V3_32/table32.qza
- IowaWoundData/MiSeqV1V3_32/denoising-stats32.qza
- IowaWoundData/MiSeqV1V3_35/rep-seqs35.qza
- IowaWoundData/MiSeqV1V3_35/table35.qza
- IowaWoundData/MiSeqV1V3_35/denoising-stats35.qza
## Output
- IowaWoundData/MiSeqV1V3_32/rep-seqs32.qzv
- IowaWoundData/MiSeqV1V3_32/table32.qzv
- IowaWoundData/MiSeqV1V3_32/denoising-stats32.qzv
- IowaWoundData/MiSeqV1V3_35/rep-seqs35.qzv
- IowaWoundData/MiSeqV1V3_35/table35.qzv
- IowaWoundData/MiSeqV1V3_35/denoising-stats35.qzv

## Scripts
[summarize_denoise.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/summarize_denoise.sh)


# (6) Make phylogeny 
## Input
- IowaWoundData/MiSeqV1V3_32/rep-seqs32.qza 
- IowaWoundData/MiSeqV1V3_35/rep-seqs35.qza 

## Outputs
-  IowaWoundData/MiSeqV1V3_32/aligned-rep-seqs32.qza
-  IowaWoundData/MiSeqV1V3_35/aligned-rep-seqs35.qza
-  IowaWoundData/MiSeqV1V3_32/masked-aligned-rep-seqs32.qza
-  IowaWoundData/MiSeqV1V3_35/masked-aligned-rep-seqs35.qza
-  IowaWoundData/MiSeqV1V3_32/unrooted-tree32.qza
-  IowaWoundData/MiSeqV1V3_35/unrooted-tree35.qza
-  IowaWoundData/MiSeqV1V3_32/rooted-tree32.qza
-  IowaWoundData/MiSeqV1V3_32/rooted-tree35.qza
## Scripts

- [make_phylogeny.sh](https://github.com/Grice-Lab/IowaWound/blob/master/scripts/make_phylogeny.sh)
    
