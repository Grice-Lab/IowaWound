# Amy Campbell
# Nov 2021
# diversity metrics for combined 

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/rooted-tree-merged.qza \
  --i-table /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/mergedtable.qza \
  --p-sampling-depth 1200 \
  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
  --output-dir /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results

# Alpha diversity
#################
qiime diversity alpha-group-significance \
  --i-alpha-diversity /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
  --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/faith-pd-group-significance.qzv

#qiime diversity alpha-group-significance \
  --i-alpha-diversity /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/evenness_vector.qza \
  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
  --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/evenness-group-significance.qzv

#qiime diversity alpha-group-significance \
  --i-alpha-diversity /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/shannon_vector.qza \
  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
  --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/shannon_group-significance.qzv


# Beta diversity
################

qiime diversity beta-group-significance \
  --i-distance-matrix /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
  --m-metadata-column run \
  --o-visualization  /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/unweighted-unifrac-run-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
  --m-metadata-column run \
  --o-visualization  /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/braycurt-run-significance.qzv \
  --p-pairwise


##first, use the unweighted unifrac data as input
#qiime emperor plot \
#  --i-pcoa /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
#  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
#  --p-custom-axes run \
#  --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/unweighted-unifrac-emperor-Run.qzv

#now repeat with bray curtis
#qiime emperor plot \
#  --i-pcoa /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/bray_curtis_pcoa_results.qza \
#  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
#  --p-custom-axes run \
#  --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/bray-curtis-emperor-Run.qzv



#qiime emperor plot \
#  --i-pcoa /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
#  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
#  --p-custom-axes woundtype \
#  --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/unweighted-unifrac-emperor-woundtype.qzv

#now repeat with bray curtis
#qiime emperor plot \
#  --i-pcoa /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/bray_curtis_pcoa_results.qza \
#  --m-metadata-file /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/DiversityMetadata.tsv \
#  --p-custom-axes woundtype \
#  --o-visualization /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/core-metrics-results/bray-curtis-emperor-woundtype.qzv


