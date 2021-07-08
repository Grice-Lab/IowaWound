#!/bin/bash
# Amy Campbell
# 07/2021
# Making a 97% identity OTU database for the V1V3 region from silva 
# vrsion 138 
# for closed reference OTU assignment of ASVs output by deblur 

#source /home/acampbe/software/miniconda3/bin/activate Qiime2Env


# Based on tutorial https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494

qiime tools import \
    --type 'FeatureData[SILVATaxonomy]' \
    --input-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/tax_slv_ssu_138.txt \
    --output-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/taxranks-silva-138-ssu-nr99.qza \

qiime tools import \
    --type 'FeatureData[SILVATaxidMap]' \
    --input-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/taxmap_slv_ssu_ref_nr_138.txt \
    --output-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/taxmap-silva-138-ssu-nr99.qza \

qiime tools import \
    --type 'Phylogeny[Rooted]' \
    --input-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/tax_slv_ssu_138.tre \
    --output-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/taxtree-silva-138-nr99.qza \

qiime tools import \
    --type 'FeatureData[RNASequence]' \
    --input-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta \
    --output-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs.qza

qiime rescript parse-silva-taxonomy \
    --i-taxonomy-tree /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/taxtree-silva-138-nr99.qza \
    --i-taxonomy-map /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/taxmap-silva-138-ssu-nr99.qza \
    --i-taxonomy-ranks /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/taxranks-silva-138-ssu-nr99.qza \
    --p-include-species-labels \
    --o-taxonomy /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-tax.qza

qiime rescript cull-seqs \
    --i-sequences /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs.qza \
    --o-clean-sequences /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs-cleaned.qza

qiime rescript dereplicate \
    --i-sequences /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs-cleaned.qza  \
    --i-taxa /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-tax.qza \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-tax-derep-uniq.qza

# 27F-534R primers
qiime feature-classifier extract-reads \
    --i-sequences /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer AGAGTTTGATCCTGGCTCAG \
    --p-r-primer ATTACCGCGGCTGCTGG \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs-27f-534r.qza

qiime rescript dereplicate \
    --i-sequences /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs-27f-534r.qza \
    --i-taxa /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-tax-derep-uniq.qza \
    --p-rank-handles 'silva' \
    --p-mode 'lca' \
    --p-perc-identity .97 \
    --o-dereplicated-sequences /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs-27f-534r-lca.qza \
    --o-dereplicated-taxa  /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-tax-27f-534r-derep-lca.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-seqs-27f-534r-lca.qza \
  --i-reference-taxonomy /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-tax-27f-534r-derep-lca.qza \
  --o-classifier /Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/silva-138-ssu-nr99-27f-534r-classifier.qza


#/Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/
#/Users/amycampbell/Desktop/GriceLabGit/IowaWound/silva/
