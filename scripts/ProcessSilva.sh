#!/bin/bash
# Amy Campbell
# 07/2021
# Making a 97% identity OTU database for the V1V3 region from silva 
# vrsion 138 
# for closed reference OTU assignment of ASVs output by deblur 

source /home/acampbe/software/miniconda3/bin/activate Qiime2Env


# Based on tutorial https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494

qiime tools import \
    --type 'FeatureData[SILVATaxonomy]' \
    --input-path /home/acampbe/DownloadedDatabases/Silva/tax_slv_ssu_138.txt \
    --output-path /home/acampbe/DownloadedDatabases/Silva/taxranks-silva-138-ssu-nr99.qza \

qiime tools import \
    --type 'FeatureData[SILVATaxidMap]' \
    --input-path /home/acampbe/DownloadedDatabases/Silva/taxmap_slv_ssu_ref_nr_138.txt \
    --output-path /home/acampbe/DownloadedDatabases/Silva/taxmap-silva-138-ssu-nr99.qza \

qiime tools import \
    --type 'Phylogeny[Rooted]' \
    --input-path /home/acampbe/DownloadedDatabases/Silva/tax_slv_ssu_138.tre \
    --output-path /home/acampbe/DownloadedDatabases/Silva/taxtree-silva-138-nr99.qza \

qiime tools import \
    --type 'FeatureData[RNASequence]' \
    --input-path /home/acampbe/DownloadedDatabases/Silva/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta \
    --output-path /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs.qza

qiime rescript parse-silva-taxonomy \
    --i-taxonomy-tree /home/acampbe/DownloadedDatabases/Silva/taxtree-silva-138-nr99.qza \
    --i-taxonomy-map /home/acampbe/DownloadedDatabases/Silva/taxmap-silva-138-ssu-nr99.qza \
    --i-taxonomy-ranks /home/acampbe/DownloadedDatabases/Silva/taxranks-silva-138-ssu-nr99.qza \
    --p-include-species-labels \
    --o-taxonomy /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-tax.qza

qiime rescript cull-seqs \
    --i-sequences /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs.qza \
    --o-clean-sequences /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs-cleaned.qza

qiime rescript dereplicate \
    --i-sequences /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs-cleaned.qza  \
    --i-taxa /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-tax.qz \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-tax-derep-uniq.qza

# 27F-534R primers
qiime feature-classifier extract-reads \
    --i-sequences /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer AGAGTTTGATCCTGGCTCAG \
    --p-r-primer ATTACCGCGGCTGCTGG \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs-27f-534r.qza

qiime rescript dereplicate \
    --i-sequences /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs-27f-534r.qza \
    --i-taxa /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-tax-derep-uniq.qza \
    --p-rank-handles 'silva' \
    --p-mode 'lca' \
    --p-perc-identity 97 \
    --o-dereplicated-sequences /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs-27f-534r-lca.qza \
    --o-dereplicated-taxa  /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-tax-27f-534r-derep-lca.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-seqs-27f-534r-lca.qza \
  --i-reference-taxonomy /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-tax-27f-534r-derep-lca.qza \
  --o-classifier /home/acampbe/DownloadedDatabases/Silva/silva-138-ssu-nr99-27f-534r-classifier.qza


#/home/acampbe/DownloadedDatabases/Silva/
#/home/acampbe/DownloadedDatabases/Silva/
