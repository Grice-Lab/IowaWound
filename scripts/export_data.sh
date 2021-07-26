
# export from qiime
###################
qiime tools export --input-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/table-cr-97.qza --output-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported
qiime tools export --input-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/taxonomy_v1v3.qza --output-path /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported

# make backup
#############
cp /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported/taxonomy.tsv /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported/biom-taxonomy_v1v3.tsv

# change first line
####################
tail -n +2  /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported/taxonomy.tsv > /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported/biom-taxonomy_v1v3.tsv

echo "#OTUID\ttaxonomy\tconfidence" | cat -  /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported/biom-taxonomy_v1v3.tsv > temp && mv temp  /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported/biom-taxonomy_v1v3.tsv

# Add metadata
biom add-metadata -i /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported/feature-table.biom -o /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/table-with-taxonomy.biom --observation-metadata-fp /Users/amycampbell/Desktop/GriceLabGit/IowaWound/data/exported/biom-taxonomy_v1v3.tsv --sc-separated taxonomy

