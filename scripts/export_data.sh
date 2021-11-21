# Amy Campbell
# Updated 11-2021

# export from qiime
###################
# this will output /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/tree.nwk
qiime tools export \
  --input-path /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/unrooted-tree-merged.qza \
  --output-path /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput


mkdir -p /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput
qiime tools export --input-path /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/mergedtable.qza --output-path /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput
qiime tools export --input-path /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/taxonomy.qza --output-path /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput

# make backup
#############
cp /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/taxonomy.tsv /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/taxonomy_backup.tsv

# change first line
####################
tail -n +2  /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/taxonomy.tsv > /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/biom-taxonomy.tsv

echo "#OTUID\ttaxonomy\tconfidence" | cat -  /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/biom-taxonomy.tsv > temp && mv temp  /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/biom-taxonomy.tsv

# Add metadata
biom add-metadata -i /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/feature-table.biom -o /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/table-with-taxonomy.biom --observation-metadata-fp /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/biom-taxonomy.tsv --sc-separated taxonomy

biom convert -i /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/table-with-taxonomy.biom -o /Users/amycampbell/Documents/IowaWoundData2021/Qiime2Data/RInput/table-with-taxonomy.tsv --to-tsv
